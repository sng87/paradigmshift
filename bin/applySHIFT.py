#!/usr/bin/env python
"""
applySHIFT.py

Author: Sam Ng
Last Updated: 2014-08-04
"""
import math, os, random, re, string, sys, types
from copy import deepcopy

import pandas
import networkx

from optparse import OptionParser
from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack

## default variables
base_directory = os.path.dirname(os.path.abspath(__file__))
paradigm_executable = 'paradigm'

## as classes
class ParadigmSetup:
    """
    Stores relevant information for preparing a Paradigm run [Paradigm-Shift specific]
    Dependencies: logger, returnColumns, returnRows, readList
    """
    def __init__(self, directory, include_samples = None, nulls = 30, batch_size = 50, public = False):
        self.directory = directory.rstrip('/')
        self.nulls = nulls
        self.batch_size = batch_size
        self.public = public
        (self.config, self.params) = (None, None)
        assert(os.path.exists('%s/config.txt' % (directory)))
        assert(os.path.exists('%s/params.txt' % (directory)))
        self.config = '%s/config.txt' % (directory)
        self.params = '%s/params.txt' % (directory)
        (self.imap, self.dogma) = (None, None)
        for file in os.listdir(directory):
            if file.endswith('.imap'):
                self.imap = '%s/%s' % (directory, file)
            elif file.endswith('.dogma'):
                self.dogma = '%s/%s' % (directory, file)
        (self.pathway) = (None)
        assert(os.path.exists('%s/clusterFiles' % (directory)))
        for file in os.listdir('%s/clusterFiles' % (directory)):
            if file.endswith('pathway.tab'):
                assert(self.pathway == None)
                self.pathway = '%s/clusterFiles/%s' % (directory, file)
        assert(self.pathway != None)
        (self.genome, self.mrna, self.protein, self.active, self.ipl) = ([], [], [], [], [])
        f = open(self.config, 'r')
        for line in f:
            if line.isspace():
                continue
            if line.startswith('evidence'):
                tokens = line.lstrip('evidence [').rstrip(']').split(',')
                attach_node = None
                attach_file = None
                for token in tokens:
                    parts = token.split('=')
                    if parts[0] == 'node':
                        attach_node = parts[1]
                    elif parts[0] == 'suffix':
                        attach_file = parts[1]
                if attach_node == 'genome':
                    self.genome.append('%s/clusterFiles/%s' % (directory, attach_file))
                elif attach_node == 'mRNA':
                    self.mrna.append('%s/clusterFiles/%s' % (directory, attach_file))
                elif attach_node == 'protein':
                    self.protein.append('%s/clusterFiles/%s' % (directory, attach_file))
                elif attach_node == 'active':
                    self.active.append('%s/clusterFiles/%s' % (directory, attach_file))
        f.close()
        assert(len(self.mrna) > 0)
        if os.path.exists('%s/merge_merged_unfiltered.tab' % (directory)):
            self.ipl.append('%s/merge_merged_unfiltered.tab' % (directory))
        elif os.path.exists('%s/merge_merged.tab' % (directory)):
            self.ipl.append('%s/merge_merged.tab' % (directory))
        (self.features) = (None)
        for file in self.genome + self.mrna:
            if self.features == None:
                self.features = returnColumns(file)
            else:
                self.features = list(set(self.features) & set(returnColumns(file)))
        self.features.sort()
        (self.samples) = (None)
        for file in self.genome + self.mrna + self.protein + self.active:
            if self.samples == None:
                self.samples = returnRows(file)
            else:
                self.samples = list(set(self.samples) & set(returnRows(file)))
        if include_samples:
            self.samples = list(set(self.samples) & set(readList(include_samples)))
        self.samples.sort()
        ## we really should do a check to make sure the sample and feature spaces match
        ## for all files, since this is a sticking point to many downstream methods

class Model:
    """
    Stores models
    """
    def __init__(self, model_directory):
        self.random_seed = 0
        self.focus_node = model_directory.split('/')[-1]
        self.focus_genes = self.focus_node.split('_')
        self.model_directory = model_directory
        self.directory = '%s_model' % (self.focus_node)
    def setSeed(self, seed_input = None):
        if seed_input == None:
            self.random_seed = random.randint(0,999999999)
        else:
            if os.path.exists(seed_input):
                f = open(seed_input, 'r')
                self.random_seed = int(f.readline().rstrip())
                f.close()
            else:
                self.random_seed = int(seed_input)
    def printSeed(self):
        o = open('seed.log', 'w')
        o.write('%s\n' % (self.random_seed))
        o.close()

class Pathway:
    """
    Paradigm compatible pathway class [2014-6-9]
    Dependencies: logger
    """
    def __init__(self, input, pid = None):
        self.nodes = {}
        self.interactions = {}
        self.pid = pid
        if type(input) == types.TupleType:
            self.nodes = deepcopy(input[0])
            self.interactions = deepcopy(input[1])
        elif type(input) == types.StringType:
            (self.nodes, self.interactions) = self.readSPF(input)
        else:
            logger('ERROR: invalid input for pathway import (%s)\n' % (input), die = True)
    def readSPF(self, input_file):
        nodes = {}
        interactions = {}
        f = open(input_file, 'r')
        for line in f:
            if line.isspace():
                continue
            pline = line.rstrip().split('\t')
            if len(pline) == 2:
                nodes[pline[1]] = pline[0]
            elif len(pline) == 3:
                if pline[0] not in interactions:
                    interactions[pline[0]] = {}
                if pline[1] not in interactions[pline[0]]:
                    interactions[pline[0]][pline[1]] = pline[2]
                else:
                    interactions[pline[0]][pline[1]] += ';' + pline[2]
            else:
                logger('ERROR: line length not 2 or 3 (%s)\n' % (line), die = True)
        f.close()
        return(nodes, interactions)
    def writeSPF(self, output_file, reverse = False):
        o = open(output_file, 'w')
        for node in self.nodes:
            o.write('%s\t%s\n' % (self.nodes[node], node))
        for source in self.interactions:
            for target in self.interactions[source]:
                for interaction in self.interactions[source][target].split(';'):
                    o.write('%s\t%s\t%s\n' % (source, target, interaction))
        o.close()
    def writeSIF(self, output_file):
        o = open(output_file, 'w')
        for source in self.interactions:
            for target in self.interactions[source]:
                for interaction in self.interactions[source][target].split(';'):
                    o.write('%s\t%s\t%s\n' % (source, interaction, target))
        o.close()
    def networkx(self):
        networkx_graph = networkx.MultiDiGraph()
        for node in self.nodes:
            networkx_graph.add_node(node, type = self.nodes[node])
        for source in self.interactions:
            for target in self.interactions[source]:
                for interaction in self.interactions[source][target].split(';'):
                    networkx_graph.add_edge(source, target, interaction = interaction)
        return(networkx_graph)
    def reverse(self):
        reversed_interactions = {}
        for source in self.interactions:
            for target in self.interactions[source]:
                if target not in reversed_interactions:
                    reversed_interactions[target] = {}
                reversed_interactions[target][source] = self.interactions[source][target]
        return(reversed_interactions)
    def getShortestPaths(self, source, target, max_distance = None):
        shortest_paths = []
        all_walks = [[source]]
        while len(shortest_paths) == 0:
            next_walks = []
            while len(all_walks) > 0:
                current_walk = all_walks.pop()
                if current_walk[-1] not in self.interactions:
                    continue
                for intermediate in self.interactions[current_walk[-1]]:
                    if intermediate in current_walk:
                        continue
                    if intermediate == target:
                        complete_walk = current_walk + [intermediate]
                        shortest_paths.append([(complete_walk[index],
                                              self.interactions[complete_walk[index]][complete_walk[index + 1]],
                                              complete_walk[index + 1])
                                              for index in range(len(complete_walk) - 1)])
                    next_walks.append(current_walk + [intermediate])
            if len(next_walks) == 0:
                break
            if max_distance is not None:
                if len(next_walks[0]) >= max_distance + 1:
                    break
            all_walks = deepcopy(next_walks)
        return(shortest_paths)
    def getAllPaths(self, source, target, max_distance):
        all_paths = []
        all_walks = [[source]]
        for distance in range(1, max_distance + 1):
            next_walks = []
            while len(all_walks) > 0:
                current_walk = all_walks.pop()
                if current_walk[-1] not in self.interactions:
                    continue
                for intermediate in self.interactions[current_walk[-1]]:
                    if intermediate in current_walk:
                        continue
                    if intermediate == target:
                        complete_walk = current_walk + [intermediate]
                        all_paths.append([(complete_walk[index], 
                                         self.interactions[complete_walk[index]][complete_walk[index + 1]],
                                         complete_walk[index + 1])
                                         for index in range(len(complete_walk) - 1)])
                    next_walks.append(current_walk + [intermediate])
            if len(next_walks) == 0:
                break
            all_walks = deepcopy(next_walks)
        return(all_paths)
    def getNeighbors(self, node, max_distance):
        reversed_interactions = self.reverse()
        seen_nodes = set([node])
        border_nodes = [node]
        frontier_nodes = []
        for distance in range(max_distance):
            while len(border_nodes) > 0:
                current_node = border_nodes.pop()
                if current_node in self.interactions:
                    for target in self.interactions[current_node]:
                        if target not in seen_nodes:
                            seen_nodes.update([target])
                            frontier_nodes.append(target)
                if current_node in reversed_interactions:
                    for source in reversed_interactions[current_node]:
                        if source not in seen_nodes:
                            seen_nodes.update([source])
                            frontier_nodes.append(source)
            border_nodes = deepcopy(frontier_nodes)
            frontier_nodes = []
        return(list(seen_nodes))
    def subsetPathway(self, subsetted_nodes):
        subsetted_pathway = Pathway( ({}, {}) )
        for source in subsetted_nodes:
            if source not in subsetted_pathway.nodes:
                subsetted_pathway.nodes[source] = self.nodes[source]
            if source in self.interactions:
                for target in self.interactions[source]:
                    if target not in subsetted_nodes:
                        continue
                    if target not in subsetted_pathway.nodes:
                        subsetted_pathway.nodes[target] = self.nodes[target]
                    if source not in subsetted_pathway.interactions:
                        subsetted_pathway.interactions[source] = {}
                    subsetted_pathway.interactions[source][target] = self.interactions[source][target]
        return(subsetted_pathway)
    def appendPathway(self, append_pathway):
        for node in append_pathway.nodes:
            if node not in self.nodes:
                self.nodes[node] = append_pathway.nodes[node]
        for source in append_pathway.interactions:
            assert(source in self.nodes)
            for target in append_pathway.interactions[source]:
                assert(target in self.nodes)
                if source not in self.interactions:
                    self.interactions[source] = {}
                if target not in self.interactions[source]:
                    self.interactions[source][target] = append_pathway.interactions[source][target]
                else:
                    self.interactions[source][target] = ';'.join(list(set(self.interactions[source][target].split(';')) |
                                                                      set(append_pathway.interactions[source][target].split(';'))))

## as functions
def logger(message, file = None, die = False):
    """
    Writes messages to standard error [2014-3-1]
    """
    if file is None:
        sys.stderr.write(message)
    else:
        o = open(file, 'a')
        o.write(message)
        o.close()
    if die:
        sys.exit(1)

def returnRows(input_file, sep = '\t', index = 0, header = True):
    """
    Returns the rows of a file without loading it into memory [2014-3-1]
    Dependencies: logger
    """
    rows = []
    f = open(input_file, 'r')
    if header:
        line = f.readline()
        if line.isspace():
            logger('ERROR: encountered a blank header\n', die = True)
    for line in f:
        if line.isspace():
            continue
        rows.append(line.rstrip().split(sep)[index])
    f.close()
    return(rows)

def returnColumns(input_file, sep = '\t', index = 0):
    """
    Returns the columns of a file without loading it into memory [2014-3-1]
    Dependencies: logger
    """
    f = open(input_file, 'r')
    if index > 0:
        for n in range(index):
            f.readline()
    line = f.readline().rstrip()
    if line.isspace():
        logger('ERROR: encountered a blank header\n', die = True)
    f.close()
    return(line.split(sep)[1:])

def readList(input_file, header = False):
    """
    Reads in a simple one column list [2014-3-1]
    """
    input_list = []
    f = open(input_file, 'r')
    if header:
        f.readline()
    for line in f:
        if line.isspace():
            continue
        input_list.append(line.rstrip())
    f.close()
    return(input_list)

def getBatchCount(cohort_size, batch_size = 15):
    """
    Batching function for public Paradigm executable [Paradigm-Shift specific]
    """
    return(int((cohort_size + batch_size - 1)/batch_size))

def getChunks(sample_list, batch_size = 15):
    """
    Chunking function for public Paradigm executable [Paradigm-Shift specific]
    """
    for index in xrange(0, len(sample_list), batch_size):
        yield sample_list[index:index + batch_size]

def generateData(focus_genes, upstream_features, downstream_features, allow_features, include_samples, data_files, nulls = 0, random_seed = 1):
    """
    Generates data files for Paradigm-Shift runs [Paradigm-Shift specific]
    """
    random.seed(random_seed)
    for file in data_files:
        data_frame = pandas.read_csv(file, sep = '\t', index_col = 0)
        output_features = list(set(data_frame.columns) & (set(upstream_features) | set(downstream_features)))
        output_samples = include_samples
        data_frame[output_features].loc[output_samples].transpose().to_csv('data/transposed_%s' % (file.split('/')[-1]), sep = '\t', na_rep = 'NA', index_label = 'genes')
        output_features = list(set(data_frame.columns) & set(upstream_features))
        output_samples = include_samples
        data_frame[output_features].loc[output_samples].to_csv('data/up_%s' % (file.split('/')[-1]), sep = '\t', na_rep = 'NA', index_label = 'samples')
        output_features = list((set(data_frame.columns) & set(downstream_features)) - set(focus_genes))
        output_samples = include_samples
        data_frame[output_features].loc[output_samples].to_csv('data/down_%s' % (file.split('/')[-1]), sep = '\t', na_rep = 'NA', index_label = 'samples')
    for null in range(1, nulls + 1):
        feature_pool = list((set(upstream_features) | set(downstream_features)) - set(focus_genes))
        permute_pool = list(set(allow_features) - set(feature_pool) - set(focus_genes))
        permute_map = {focus_gene : focus_gene for focus_gene in focus_genes}
        for permute in zip(feature_pool, random.sample(permute_pool, len(feature_pool))):
            permute_map[permute[0]] = permute[1]
        for file in data_files:
            data_frame = pandas.read_csv(file, sep = '\t', index_col = 0)
            output_features = list(set(data_frame.columns) & set(upstream_features))
            output_samples = include_samples
            permute_features = [permute_map[feature] for feature in output_features]
            permute_data_frame = data_frame[permute_features].loc[output_samples].copy()
            permute_data_frame.columns = output_features
            permute_data_frame.to_csv('data/up_N%s_%s' % (null, file.split('/')[-1]), sep = '\t', na_rep = 'NA', index_label = 'samples')
            output_features = list((set(data_frame.columns) & set(downstream_features)) - set(focus_genes))
            output_samples = include_samples
            permute_features = [permute_map[feature] for feature in output_features]
            permute_data_frame = data_frame[permute_features].loc[output_samples].copy()
            permute_data_frame.columns = output_features
            permute_data_frame.to_csv('data/down_N%s_%s' % (null, file.split('/')[-1]), sep = '\t', na_rep = 'NA', index_label = 'samples')

def generateBatchedData(focus_genes, upstream_features, downstream_features, allow_features, include_samples, data_files, nulls = 0, batch_size = 50, random_seed = 1):
    """
    Generates data files for Paradigm-Shift runs [Paradigm-Shift specific]
    Dependencies: getBatchCount, getChunks
    """
    random.seed(random_seed)
    for file in data_files:
        data_frame = pandas.read_csv(file, sep = '\t', index_col = 0)
        output_features = list(set(data_frame.columns) & (set(upstream_features) | set(downstream_features)))
        output_samples = include_samples
        data_frame[output_features].loc[output_samples].transpose().to_csv('data/transposed_%s' % (file.split('/')[-1]), sep = '\t', na_rep = 'NA', index_label = 'genes')
        batches = getBatchCount(len(include_samples), batch_size = batch_size)
        chunker = getChunks(include_samples, batch_size = batch_size)
        for b in range(batches):
            current_chunk = chunker.next()
            output_features = list(set(data_frame.columns) & set(upstream_features))
            output_samples = current_chunk
            data_frame[output_features].loc[output_samples].to_csv('data/up_b%s_%s_%s' % (b, batches, file.split('/')[-1]), sep = '\t', na_rep = 'NA', index_label = 'samples')
            output_features = list((set(data_frame.columns) & set(downstream_features)) - set(focus_genes))
            output_samples = current_chunk
            data_frame[output_features].loc[output_samples].to_csv('data/down_b%s_%s_%s' % (b, batches, file.split('/')[-1]), sep = '\t', na_rep = 'NA', index_label = 'samples')
    for null in range(1, nulls + 1):
        feature_pool = list((set(upstream_features) | set(downstream_features)) - set(focus_genes))
        permute_pool = list(set(allow_features) - set(feature_pool) - set(focus_genes))
        permute_map = {focus_gene : focus_gene for focus_gene in focus_genes}
        for permute in zip(feature_pool, random.sample(permute_pool, len(feature_pool))):
            permute_map[permute[0]] = permute[1]
        for file in data_files:
            data_frame = pandas.read_csv(file, sep = '\t', index_col = 0)
            batches = getBatchCount(len(include_samples), batch_size = batch_size)
            chunker = getChunks(include_samples, batch_size = batch_size)
            for b in range(batches):
                current_chunk = chunker.next()
                output_features = list(set(data_frame.columns) & set(upstream_features))
                output_samples = current_chunk
                permute_features = [permute_map[feature] for feature in output_features]
                permute_data_frame = data_frame[permute_features].loc[output_samples].copy()
                permute_data_frame.columns = output_features
                permute_data_frame .to_csv('data/up_N%s_b%s_%s_%s' % (null, b, batches, file.split('/')[-1]), sep = '\t', na_rep = 'NA', index_label = 'samples')
                output_features = list((set(data_frame.columns) & set(downstream_features)) - set(focus_genes))
                output_samples = current_chunk
                permute_features = [permute_map[feature] for feature in output_features]
                permute_data_frame = data_frame[permute_features].loc[output_samples].copy()
                permute_data_frame.columns = output_features
                permute_data_frame .to_csv('data/down_N%s_b%s_%s_%s' % (null, b, batches, file.split('/')[-1]), sep = '\t', na_rep = 'NA', index_label = 'samples')

## jt classes
class jtCmd(Target):
    def __init__(self, command, directory, file = None):
        Target.__init__(self, time=1000)
        self.command = command
        self.directory = directory
        self.file = file
    def run(self):
        os.chdir(self.directory)
        if self.file:
            o = open(self.file, 'a')
            o.write('%s\n' % (self.command))
            o.close()
        os.system(self.command)

class queueModels(Target):
    def __init__(self, model_list, paradigm_setup, directory):
        Target.__init__(self, time=10000)
        self.model_list = model_list
        self.paradigm_setup = paradigm_setup
        self.directory = directory
    def run(self):
        os.chdir(self.directory)
        
        ## queue models
        if not os.path.exists('analysis'):
            os.mkdir('analysis')
        if len(self.model_list) > 0:
            model = self.model_list[0]
            if not os.path.exists('analysis/%s' % (model.directory)):
                os.mkdir('analysis/%s' % (model.directory))
                #prepareModel(model, self.paradigm_setup, self.directory).run()
                self.addChildTarget(prepareModel(model,
                                                 self.paradigm_setup,
                                                 self.directory))
            else:
                pass
            #queueModels(self.model_list[1:], self.paradigm_setup, self.directory).run()
            self.setFollowOnTarget(queueModels(self.model_list[1:],
                                               self.paradigm_setup,
                                               self.directory))

class prepareModel(Target):
    def __init__(self, model, paradigm_setup, directory):
        Target.__init__(self, time=10000)
        self.model = model
        self.paradigm_setup = paradigm_setup
        self.directory = directory
    def run(self):
        ps_directory = '%s/analysis/%s' % (self.directory, self.model.directory)
        os.chdir(ps_directory)
        
        ## load model
        assert((os.path.exists('%s/upstream_pathway.tab' % (self.model.model_directory))) and (os.path.exists('%s/downstream_pathway.tab' % (self.model.model_directory))))
        selected_upstream = Pathway('%s/upstream_pathway.tab' % (self.model.model_directory))
        selected_downstream = Pathway('%s/downstream_pathway.tab' % (self.model.model_directory))
        
        ## generate data and run paradigm
        os.mkdir('data')
        data_files = self.paradigm_setup.genome + self.paradigm_setup.mrna + self.paradigm_setup.protein + self.paradigm_setup.active
        upstream_features = list(set(selected_upstream.nodes.keys()) & set(self.paradigm_setup.features))
        downstream_features = list(set(selected_downstream.nodes.keys()) & set(self.paradigm_setup.features))
        if self.paradigm_setup.public:
            generateBatchedData(self.model.focus_genes,
                                upstream_features,
                                downstream_features,
                                self.paradigm_setup.features,
                                self.paradigm_setup.samples,
                                data_files,
                                nulls = self.paradigm_setup.nulls,
                                batch_size = self.paradigm_setup.batch_size,
                                random_seed = self.model.random_seed)
        else:
            generateData(self.model.focus_genes,
                         upstream_features,
                         downstream_features,
                         self.paradigm_setup.features,
                         self.paradigm_setup.samples,
                         data_files,
                         nulls = self.paradigm_setup.nulls,
                         random_seed = self.model.random_seed)
        selected_upstream.writeSPF('upstream_pathway.tab')
        selected_downstream.writeSPF('downstream_pathway.tab')
        self.addChildTarget(runParadigm(self.model,
                                        self.paradigm_setup,
                                        ps_directory))
        self.setFollowOnTarget(computeMShifts(self.model,
                                              self.paradigm_setup,
                                              self.directory))

class runParadigm(Target):
    def __init__(self, analysis, paradigm_setup, directory, nulls = None):
        Target.__init__(self, time=10000)
        self.analysis = analysis
        self.paradigm_setup = paradigm_setup
        self.directory = directory
        if nulls is None:
            self.nulls = self.paradigm_setup.nulls
        else:
            self.nulls = nulls
    def run(self):
        os.chdir(self.directory)
        os.mkdir('paradigm')
        logger('Running Paradigm inference ...\n', file = 'progress.log')
        
        ## copy paradigm files
        os.system('cp %s .' % (self.paradigm_setup.config))
        os.system('cp %s .' % (self.paradigm_setup.params))
        if self.paradigm_setup.dogma:
            os.system('cp %s .' % (self.paradigm_setup.dogma))
        if self.paradigm_setup.imap:
            os.system('cp %s .' % (self.paradigm_setup.imap))
        ## run Paradigm (observed and nulls)
        batches = getBatchCount(len(self.paradigm_setup.samples), batch_size = self.paradigm_setup.batch_size)
        if batches == 0:
            self.addChildTarget(jtCmd('%s -p upstream_pathway.tab -c config.txt -b %s -o %s' % (paradigm_executable, 'data/up_', 'paradigm/%s_upstream.fa' % (self.analysis.focus_node)), self.directory, file = 'jobs.list'))
            self.addChildTarget(jtCmd('%s -p downstream_pathway.tab -c config.txt -b %s -o %s' % (paradigm_executable, 'data/down_', 'paradigm/%s_downstream.fa' % (self.analysis.focus_node)), self.directory, file = 'jobs.list'))
            for null in range(1, self.nulls + 1):
                self.addChildTarget(jtCmd('%s -p upstream_pathway.tab -c config.txt -b %s -o %s' % (paradigm_executable, 'data/up_N%s_' % (null), 'paradigm/N%s_%s_upstream.fa' % (null, self.analysis.focus_node)), self.directory, file = 'jobs.list'))
                self.addChildTarget(jtCmd('%s -p downstream_pathway.tab -c config.txt -b %s -o %s' % (paradigm_executable, 'data/down_N%s_' % (null), 'paradigm/N%s_%s_downstream.fa' % (null, self.analysis.focus_node)), self.directory, file = 'jobs.list'))
        elif self.paradigm_setup.public:
            os.mkdir('outputFiles')
            for b in range(batches):
                self.addChildTarget(jtCmd('%s -p upstream_pathway.tab -c config.txt -b %s -o %s' % (paradigm_executable, 'data/up_b%s_%s_' % (b, batches), 'outputFiles/%s_upstream_b%s_%s.fa' % (self.analysis.focus_node, b, batches)), self.directory, file = 'jobs.list'))
                self.addChildTarget(jtCmd('%s -p downstream_pathway.tab -c config.txt -b %s -o %s' % (paradigm_executable, 'data/down_b%s_%s_' % (b, batches), 'outputFiles/%s_downstream_b%s_%s.fa' % (self.analysis.focus_node, b, batches)), self.directory, file = 'jobs.list'))
                for null in range(1, self.nulls + 1):
                    self.addChildTarget(jtCmd('%s -p upstream_pathway.tab -c config.txt -b %s -o %s' % (paradigm_executable, 'data/up_N%s_b%s_%s_' % (null, b, batches), 'outputFiles/N%s_%s_upstream_b%s_%s.fa' % (null, self.analysis.focus_node, b, batches)), self.directory, file = 'jobs.list'))
                    self.addChildTarget(jtCmd('%s -p downstream_pathway.tab -c config.txt -b %s -o %s' % (paradigm_executable, 'data/down_N%s_b%s_%s_' % (null, b, batches), 'outputFiles/N%s_%s_downstream_b%s_%s.fa' % (null, self.analysis.focus_node, b, batches)), self.directory, file = 'jobs.list'))
        else:
            os.mkdir('outputFiles')
            for b in range(batches):
                self.addChildTarget(jtCmd('%s -p upstream_pathway.tab -c config.txt -b %s -o %s -s %s,%s' % (paradigm_executable, 'data/up_', 'outputFiles/%s_upstream_b%s_%s.fa' % (self.analysis.focus_node, b, batches), b, batches), self.directory, file = 'jobs.list'))
                self.addChildTarget(jtCmd('%s -p downstream_pathway.tab -c config.txt -b %s -o %s -s %s,%s' % (paradigm_executable, 'data/down_', 'outputFiles/%s_downstream_b%s_%s.fa' % (self.analysis.focus_node, b, batches), b, batches), self.directory, file = 'jobs.list'))
                for null in range(1, self.nulls + 1):
                    self.addChildTarget(jtCmd('%s -p upstream_pathway.tab -c config.txt -b %s -o %s -s %s,%s' % (paradigm_executable, 'data/up_N%s_' % (null), 'outputFiles/N%s_%s_upstream_b%s_%s.fa' % (null, self.analysis.focus_node, b, batches), b, batches), self.directory, file = 'jobs.list'))
                    self.addChildTarget(jtCmd('%s -p downstream_pathway.tab -c config.txt -b %s -o %s -s %s,%s' % (paradigm_executable, 'data/down_N%s_' % (null), 'outputFiles/N%s_%s_downstream_b%s_%s.fa' % (null, self.analysis.focus_node, b, batches), b, batches), self.directory, file = 'jobs.list'))
        self.setFollowOnTarget(collectParadigm(self.analysis, self.paradigm_setup, self.directory, nulls = self.nulls))

class collectParadigm(Target):
    def __init__(self, analysis, paradigm_setup, directory, nulls = None):
        Target.__init__(self, time=10000)
        self.analysis = analysis
        self.paradigm_setup = paradigm_setup
        self.directory = directory
        if nulls is None:
            self.nulls = self.paradigm_setup.nulls
        else:
            self.nulls = nulls
    def run(self):
        os.chdir(self.directory)
        
        batches = getBatchCount(len(self.paradigm_setup.samples), batch_size = self.paradigm_setup.batch_size)
        for b in range(batches):
            os.system('cat outputFiles/%s_upstream_b%s_%s.fa >> paradigm/%s_upstream.fa' % (self.analysis.focus_node, b, batches, self.analysis.focus_node))
            os.system('cat outputFiles/%s_downstream_b%s_%s.fa >> paradigm/%s_downstream.fa' % (self.analysis.focus_node, b, batches, self.analysis.focus_node))
            for null in range(1, self.nulls + 1):
                os.system('cat outputFiles/N%s_%s_upstream_b%s_%s.fa >> paradigm/N%s_%s_upstream.fa' % (null, self.analysis.focus_node, b, batches, null, self.analysis.focus_node))
                os.system('cat outputFiles/N%s_%s_downstream_b%s_%s.fa >> paradigm/N%s_%s_downstream.fa' % (null, self.analysis.focus_node, b, batches, null, self.analysis.focus_node))
        if os.path.exists('outputFiles'):
            os.system('rm -rf outputFiles')


class computeMShifts(Target):
    def __init__(self, model, paradigm_setup, directory):
        Target.__init__(self, time=10000)
        self.model = model
        self.paradigm_setup = paradigm_setup
        self.directory = directory
    def run(self):
        ps_directory = '%s/analysis/%s' % (self.directory, self.model.directory)
        os.chdir(ps_directory)
        
        ## read in Paradigm inferences
        assert(os.path.exists('paradigm/%s_upstream.fa' % (self.model.focus_node)))
        assert(os.path.exists('paradigm/%s_downstream.fa' % (self.model.focus_node)))
        upstream_ipls = readParadigm('paradigm/%s_upstream.fa' % (self.model.focus_node))[1]
        downstream_ipls = readParadigm('paradigm/%s_downstream.fa' % (self.model.focus_node))[1]
        null_upstream_ipls = {}
        null_downstream_ipls = {}
        for null in range(1, self.paradigm_setup.nulls + 1):
            assert(os.path.exists('paradigm/N%s_%s_upstream.fa' % (null, self.model.focus_node)))
            assert(os.path.exists('paradigm/N%s_%s_downstream.fa' % (null, self.model.focus_node)))
            null_upstream_ipls[null] = readParadigm('paradigm/N%s_%s_upstream.fa' % (null, self.model.focus_node))[1]
            null_downstream_ipls[null] = readParadigm('paradigm/N%s_%s_downstream.fa' % (null, self.model.focus_node))[1]
        
        ## compute raw p-shifts
        raw_shifts = {'real' : {}}
        for null in range(1, self.paradigm_setup.nulls + 1):
            raw_shifts['null%s' % (null)] = {}
        for sample in self.paradigm_setup.samples:
            assert((sample in downstream_ipls.columns) and (sample in upstream_ipls.columns))
            raw_shifts['real'][sample] = (downstream_ipls[sample][self.model.focus_node] - upstream_ipls[sample][self.model.focus_node])
            for null in range(1, self.paradigm_setup.nulls + 1):
                raw_shifts['null%s' % (null)][sample] = (null_downstream_ipls[null][sample][self.model.focus_node] - null_upstream_ipls[null][sample][self.model.focus_node])
        pandas.DataFrame(raw_shifts).to_csv('all_shifts.tab', sep = '\t', index_label = 'id')

def as_main():
    ## check for fresh run
    if os.path.exists('.jobTree'):
        logger('WARNING: .jobTree directory found, remove it first to start a fresh run\n')
    
    ## parse arguments
    parser = OptionParser(usage = '%prog [options] paradigm_directory model_directory')
    Stack.addJobTreeOptions(parser)
    parser.add_option('--jobFile', help='Add as child of jobFile rather than new jobTree')
    parser.add_option('-n', '--nulls', dest='nulls', default=0)
    parser.add_option('-b', '--batchsize', dest='batch_size', default=50)
    parser.add_option('-y', '--public', action='store_true', dest='paradigm_public', default=False)
    parser.add_option('-g', '--galaxy', dest='galaxy_run', action='store_true', default=False)
    options, args = parser.parse_args()
    logger('Using Batch System : %s\n' % (options.batchSystem))
    
    if len(args) == 1:
        if args[0] == 'clean':
            command = 'rm -rf .jobTree analysis'
            logger(command)
            os.system(command)
            sys.exit(0)
    
    assert(len(args) == 2)
    paradigm_directory = os.path.abspath(args[0])
    model_directory = os.path.abspath(args[1])
    
    ## set Galaxy Executables
    global paradigm_executable, circleplot_executable
    if options.galaxy_run:
        paradigm_executable = os.path.join(base_directory, 'paradigm')
    
    ## set Paradigm files
    paradigm_setup = ParadigmSetup(paradigm_directory,
                                   nulls = int(options.nulls),
                                   batch_size = int(options.batch_size),
                                   public = options.paradigm_public)
    
    ## gather models
    model_list = []
    for file in os.listdir(model_directory):
        if os.path.exists('%s/%s/upstream_pathway.tab' % (model_directory, file)) and os.path.exists('%s/%s/downstream_pathway.tab' % (model_directory, file)):
            model = Model('%s/%s' % (model_directory, file))
            model.setSeed()
            model_list.append(model)
    
    ## run
    s = Stack(queueModels(model_list,
                          paradigm_setup,
                          os.getcwd().rstrip('/')))
    if options.jobFile:
        s.addToJobFile(options.jobFile)
    else:
        if options.jobTree == None:
            options.jobTree = './.jobTree'
        
        failed = s.startJobTree(options)
        #failed = ''
        #queueModels(model_list, paradigm_setup, os.getcwd().rstrip('/')).run()
        if failed:
            logger('%d jobs failed\n' % failed)
        else:
            os.system('rm -rf .lastjobTree')
            os.system('mv .jobTree .lastjobTree')

if __name__ == '__main__':
    from applySHIFT import *
    as_main()