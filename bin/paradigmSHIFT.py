#!/usr/bin/env python
"""
paradigmSHIFT.py
    by Sam Ng
"""
import math, os, random, re, shutil, string, sys, types
from copy import deepcopy

import pandas
import networkx

from optparse import OptionParser
from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack

## executables
bin_directory = os.path.dirname(os.path.abspath(__file__))
# paradigm_executable = os.path.join(bin_directory, 'paradigm')
paradigm_executable = 'paradigm'
circleplot_executable = os.path.join(bin_directory, 'circlePlot.py')

## classes
class ParadigmSetup:
    """
    Stores relevant information for preparing Paradigm runs [paradigmSHIFT.py specific]
    Dependencies: returnColumns, returnRows, readList
    """
    def __init__(self, directory, include_samples = None, null_size = 30, batch_size = 50, pathway_file = None, public = False):
        assert(os.path.exists('%s/clusterFiles' % (directory)))
        self.directory = directory
        self.null_size = null_size
        self.batch_size = batch_size
        self.public = public
        ## locate config and params files
        (self.config, self.params) = (None, None)
        assert(os.path.exists('%s/config.txt' % (directory)))
        assert(os.path.exists('%s/params.txt' % (directory)))
        self.config = '%s/config.txt' % (directory)
        self.params = '%s/params.txt' % (directory)
        ## locate imap and dogma files (optional)
        (self.imap, self.dogma) = (None, None)
        for file in os.listdir(directory):
            if file.endswith('.imap'):
                self.imap = '%s/%s' % (directory, file)
            elif file.endswith('.dogma'):
                self.dogma = '%s/%s' % (directory, file)
        ## locate pathway file or set user defined pathway file (optional)
        (self.pathway) = (None)
        if pathway_file is None:
            for file in os.listdir('%s/clusterFiles' % (directory)):
                if file.endswith('pathway.tab'):
                    assert(self.pathway == None)
                    self.pathway = '%s/clusterFiles/%s' % (directory, file)
        else:
            assert(os.path.exists(pathway_file))
            self.pathway = os.path.abspath(pathway_file) 
        assert(self.pathway != None)
        ## locate data files from config file and clusterFiles directory
        self.data = []
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
                    self.data.append(attach_file)
                    self.genome.append('%s/clusterFiles/%s' % (directory, attach_file))
                elif attach_node == 'mRNA':
                    self.data.append(attach_file)
                    self.mrna.append('%s/clusterFiles/%s' % (directory, attach_file))
                elif attach_node == 'protein':
                    self.data.append(attach_file)
                    self.protein.append('%s/clusterFiles/%s' % (directory, attach_file))
                elif attach_node == 'active':
                    self.data.append(attach_file)
                    self.active.append('%s/clusterFiles/%s' % (directory, attach_file))
        f.close()
        assert(len(self.mrna) > 0)
        if os.path.exists('%s/merge_merged_unfiltered.tab' % (directory)):
            self.ipl.append('%s/merge_merged_unfiltered.tab' % (directory))
        elif os.path.exists('%s/merge_merged.tab' % (directory)):
            self.ipl.append('%s/merge_merged.tab' % (directory))
        ## determine rows/columns from data files or user defined (optional)
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
        ## convert params file to version 1 if public
        if self.public:
            self.convertParams(self.params, self.data)
        #### should add a check that the feature and sample spaces match for all data
        #### files, which should be true for any valid Paradigm run
        def convertParams(params_file, data_files):
            f = open(params_file, "r")
            file_header = f.readline().rstrip()
            attribute_map = {}
            for part in file_header.lstrip("> parameters ").lstrip("> ").split(" "):
                parts = part.split("=")
                attribute_map[parts[0]] = parts[1]
            if "version" not in attribute_map:
                attribute_map["version"] = "2"
            if "version" == "1":
                f.close()
            elif "version" == "2":
                data_indices = []
                table_index = 0
                header_map = {}
                params_map = {}
                for line in f:
                    if line.startswith(">"):
                        table_index += 1
                        params_map[table_index] = []
                        (child, parents) = line.rstrip().split(" ")[-1].split("=")
                        if child in data_files:
                            header_map[table_index] = "> child='%s' edge1='-obs>'\n" % (child)
                            data_indices.append(table_index)
                        else:
                            header_map[table_index] = "> child='%s' edge1='-obs>'\n" % (child)
                    else:
                        params_map[table_index].append(line.rstrip())
                f.close()
                shutil.move(params_file, params_file + ".v%s" % attribute_map["version"])
                attribute_map["version"] = "1"
                o = open(params_file, "w")
                o.write("> %s\n" % (" ".join(["=".join(list(item)) for item in attribute_map.items()])))
                for data_index in data_indices:
                    o.write(header_map[data_index])
                    for index, param in enumerate(params_map[data_index]):
                        first_index = index % 3
                        second_index = (index - first_index)/3
                        o.write("%s\t%s\t%s\n" % (first_index, second_index, param))
                o.close()

class Parameters:
    """
    Stores parameters used for the Paradigm-Shift algorithm [paradigmSHIFT.py specific]
    """
    def __init__(self):
        self.random_seed = 0
        self.min_alterations = 20
        self.in_parallel = False
        self.n_rounds = 1
        self.m_folds = 5
        self.search_distance = 2
        self.threshold = [0.84]
        self.cost = [0.0]
        self.selection_method = ['tt']
        self.model_directory = None
        self.cross_validation = True
        self.cross_validation_two_sided = False
        self.cross_validation_threshold = 0.55
        self.separation_method = 'tt'
        self.report_directory = 'report'
    def importConfig(self, config_file):
        f = open(config_file, 'r')
        for line in f:
            if line.isspace():
                continue
            pline = line.rstrip().split('\t')
            if line.startswith('search_distance'):
                self.search_distance = int(pline[1])
            elif line.startswith('threshold'):
                self.threshold = [float(item) for item in pline[1].split(',')]
            elif line.startswith('cost'):
                self.cost = [float(item) for item in pline[1].split(',')]
            elif line.startswith('selection') or line.startswith('selection_method'):
                self.selection_method = [item for item in pline[1].split(',')]
            elif line.startswith('model'):
                self.model_directory = pline[1].rstrip('/')
            elif line.startswith('cv') or line.startswith('cross_validation'):
                self.cross_validation = bool(int(pline[1]))
            elif line.startswith('msep') or line.startswith('separation_method'):
                self.separation_method = pline[1]
            elif line.startswith('report'):
                self.report_directory = pline[1].rstrip('/')
        f.close()
    def setSeed(self, seed_input = None):
        if seed_input is None:
            self.random_seed = random.randint(0, 999999999)
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

class Alterations:
    """
    Stores the alterations for a particular analysis run [paradigmSHIFT.py specific]
    """
    def __init__(self, analysis_name, focus_genes, all_samples, positive_samples, negative_samples = None):
        self.analysis_name = analysis_name
        self.focus_genes = focus_genes
        self.focus_node = '_'.join(self.focus_genes)
        self.positive_samples = list(set(positive_samples) & set(all_samples))
        if negative_samples:
            self.negative_samples = list(set(negative_samples) & set(all_samples))
        else:
            self.negative_samples = list(set(all_samples) - set(positive_samples))
        self.directory = re.compile('[\W_]+').sub('_', analysis_name)
    def getProportion(self):
        return(float(len(self.positive_samples))/float(len(self.positive_samples) + len(self.negative_samples)))

class Pathway:
    """
    Paradigm compatible pathway class [2014-12-14]
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
            raise Exception("invalid input for pathway import (%s)\n" % (input))
    def readSPF(self, input_file):
        nodes = {}
        interactions = {}
        f = open(input_file, "r")
        for line in f:
            if line.isspace():
                continue
            pline = line.rstrip().split("\t")
            if len(pline) == 2:
                nodes[pline[1]] = pline[0]
            elif len(pline) == 3:
                if pline[0] not in interactions:
                    interactions[pline[0]] = {}
                if pline[1] not in interactions[pline[0]]:
                    interactions[pline[0]][pline[1]] = pline[2]
                else:
                    interactions[pline[0]][pline[1]] += ";" + pline[2]
            else:
                raise Exception("expected length 2 or 3 on line %s\n" % (line))
        f.close()
        return(nodes, interactions)
    def writeSPF(self, output_file, reverse = False):
        o = open(output_file, "w")
        for node in self.nodes:
            o.write("%s\t%s\n" % (self.nodes[node], node))
        for source in self.interactions:
            for target in self.interactions[source]:
                for interaction in self.interactions[source][target].split(";"):
                    o.write("%s\t%s\t%s\n" % (source, target, interaction))
        o.close()
    def writeSIF(self, output_file):
        o = open(output_file, "w")
        for source in self.interactions:
            for target in self.interactions[source]:
                for interaction in self.interactions[source][target].split(";"):
                    o.write("%s\t%s\t%s\n" % (source, interaction, target))
        o.close()
    def networkx(self):
        networkx_graph = networkx.MultiDiGraph()
        for node in self.nodes:
            networkx_graph.add_node(node, type = self.nodes[node])
        for source in self.interactions:
            for target in self.interactions[source]:
                for interaction in self.interactions[source][target].split(";"):
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
    def getShortestPaths(self, source, target, search_distance = None):
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
            if search_distance is not None:
                if len(next_walks[0]) >= search_distance + 1:
                    break
            all_walks = deepcopy(next_walks)
        return(shortest_paths)
    def getAllPaths(self, source, target, search_distance):
        all_paths = []
        all_walks = [[source]]
        for distance in range(1, search_distance + 1):
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
    def getNeighbors(self, node, search_distance):
        reversed_interactions = self.reverse()
        seen_nodes = set([node])
        border_nodes = [node]
        frontier_nodes = []
        for distance in range(search_distance):
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
                    self.interactions[source][target] = ";".join(list(set(self.interactions[source][target].split(";")) |
                                                                      set(append_pathway.interactions[source][target].split(";"))))

## functions
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

def returnRows(input_file, sep = "\t", index = 0, header = True):
    """
    Returns the rows of a file without loading it into memory [2014-12-14]
    """
    rows = []
    f = open(input_file, "r")
    if header:
        line = f.readline()
        if line.isspace():
            raise Exception("encountered a blank header\n")
    for line in f:
        if line.isspace():
            continue
        rows.append(line.rstrip().split(sep)[index])
    f.close()
    return(rows)

def returnColumns(input_file, sep = "\t", index = 0):
    """
    Returns the columns of a file without loading it into memory [2014-12-14]
    """
    f = open(input_file, "r")
    if index > 0:
        for n in range(index):
            f.readline()
    line = f.readline().rstrip()
    if line.isspace():
        raise Exception("encountered a blank header\n")
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

def getFullNeighborhood(focus_genes, global_pathway, search_distance = 2):
    """
    Simplifies complex interactions to match the Paradigm-Shift inference model [Paradigm-Shift specific]
    """
    def globComplexes(focus_gene, global_pathway):
        group_index = {0 : [focus_gene]}
        if focus_gene not in global_pathway.interactions:
            return(group_index)
        focus_complexes = []
        for target in global_pathway.interactions[focus_gene]:
            if global_pathway.interactions[focus_gene][target] == 'component>':
                focus_complexes.append(target)
        seen_complexes = set()
        grouped_complexes = []
        while len(focus_complexes) > 0:
            current_complex = focus_complexes.pop(0)
            if current_complex in seen_complexes:
                continue
            seen_complexes.update([current_complex])
            current_group = [current_complex]
            border_nodes = [current_complex]
            while len(border_nodes) > 0:
                current_node = border_nodes.pop(0)
                if current_node in global_pathway.interactions:
                    for target in global_pathway.interactions[current_node]:
                        if global_pathway.interactions[current_node][target] == 'component>' and global_pathway.nodes[target] == 'complex' and target not in seen_complexes:
                            current_group.append(target)
                            border_nodes.append(target)
                            seen_complexes.update([target])
            grouped_complexes.append(deepcopy(current_group))
        for index, group in enumerate(grouped_complexes):
            group_index[index + 1] = deepcopy(group)
        return(group_index)
    def searchUpstream(focus_gene, focus_group, global_pathway, reversed_interactions, search_distance = 2):
        upstream_pathway = Pathway( ({focus_gene : global_pathway.nodes[focus_gene]}, {}) )
        seen_nodes = set(focus_group)
        seen_nodes.update([focus_gene])
        border_nodes = deepcopy(focus_group)
        frontier_nodes = []
        for distance in range(1, search_distance + 1):
            while(len(border_nodes) > 0):
                current_node = border_nodes.pop(0)
                if current_node in reversed_interactions:
                    for source in reversed_interactions[current_node]:
                        if source == focus_gene:
                            continue
                        elif global_pathway.nodes[source] == 'abstract':
                            continue
                        elif distance == 1:
                            if reversed_interactions[current_node][source] == 'component>':
                                if global_pathway.nodes[source] == 'complex' and source not in seen_nodes:
                                    border_nodes.append(source)
                                    seen_nodes.update([source])
                                elif global_pathway.nodes[source] == 'protein' or global_pathway.nodes[source] == 'family':
                                    upstream_pathway.nodes[source] = global_pathway.nodes[source]
                                    if source not in upstream_pathway.interactions:
                                        upstream_pathway.interactions[source] = {}
                                    upstream_pathway.interactions[source][focus_gene] = '-a>'
                                    frontier_nodes.append(source)
                                    seen_nodes.update([source])
                            elif reversed_interactions[current_node][source].startswith('-a') or reversed_interactions[current_node][source].startswith('-t'):
                                upstream_pathway.nodes[source] = global_pathway.nodes[source]
                                if source not in upstream_pathway.interactions:
                                    upstream_pathway.interactions[source] = {}
                                upstream_pathway.interactions[source][focus_gene] = reversed_interactions[current_node][source]
                                frontier_nodes.append(source)
                                seen_nodes.update([source])
                        elif source in upstream_pathway.nodes:
                            upstream_pathway.interactions[source][current_node] = reversed_interactions[current_node][source]
                        else:
                            upstream_pathway.nodes[source] = global_pathway.nodes[source]
                            if source not in upstream_pathway.interactions:
                                upstream_pathway.interactions[source] = {}
                            upstream_pathway.interactions[source][current_node] = reversed_interactions[current_node][source]
                            frontier_nodes.append(source)
                            seen_nodes.update([source])
            border_nodes = deepcopy(frontier_nodes)
        return(upstream_pathway)
    def searchDownstream(focus_gene, focus_group, global_pathway, reversed_interactions, search_distance = 2):
        downstream_pathway = Pathway( ({focus_gene : global_pathway.nodes[focus_gene]}, {}) )
        seen_nodes = set(focus_group)
        seen_nodes.update([focus_gene])
        border_nodes = deepcopy(focus_group)
        frontier_nodes = []
        for distance in range(1, search_distance + 1):
            while(len(border_nodes) > 0):
                current_node = border_nodes.pop(0)
                if current_node in global_pathway.interactions:
                    for target in global_pathway.interactions[current_node]:
                        if target == focus_gene:
                            continue
                        elif global_pathway.nodes[target] == 'abstract':
                            continue
                        elif distance == 1:
                            if global_pathway.interactions[current_node][target] == 'component>':
                                continue
                            elif global_pathway.interactions[current_node][target].startswith('-a') or global_pathway.interactions[current_node][target].startswith('-t'):
                                downstream_pathway.nodes[target] = global_pathway.nodes[target]
                                if focus_gene not in downstream_pathway.interactions:
                                    downstream_pathway.interactions[focus_gene] = {}
                                downstream_pathway.interactions[focus_gene][target] = global_pathway.interactions[current_node][target]
                                frontier_nodes.append(target)
                                seen_nodes.update([target])
                        elif target in downstream_pathway.nodes:
                            if current_node not in downstream_pathway.interactions:
                                downstream_pathway.interactions[current_node] = {}
                            downstream_pathway.interactions[current_node][target] = global_pathway.interactions[current_node][target]
                        else:
                            downstream_pathway.nodes[target] = global_pathway.nodes[target]
                            if current_node not in downstream_pathway.interactions:
                                downstream_pathway.interactions[current_node] = {}
                            downstream_pathway.interactions[current_node][target] = global_pathway.interactions[current_node][target]
                            frontier_nodes.append(target)
                            seen_nodes.update([target])
            border_nodes = deepcopy(frontier_nodes)
        return(downstream_pathway)
    
    reversed_interactions = global_pathway.reverse()
    group_index = {}
    group_pathways = {}
    for focus_gene in focus_genes:
        group_index[focus_gene] = globComplexes(focus_gene, global_pathway)
        group_pathways[focus_gene] = {}
        for index in group_index[focus_gene]:
            group_pathways[focus_gene][index] = [deepcopy(searchUpstream(focus_gene,
                                                                         group_index[focus_gene][index],
                                                                         global_pathway,
                                                                         reversed_interactions,
                                                                         search_distance = search_distance)),
                                                 deepcopy(searchDownstream(focus_gene,
                                                                           group_index[focus_gene][index],
                                                                           global_pathway,
                                                                           reversed_interactions,
                                                                           search_distance = search_distance))]
    return(group_index, group_pathways)

def getRelevantPaths(focus_genes, upstream_pathway, downstream_pathway, search_distance = 2):
    """
    Identifies relevant paths between features and the focus gene [Paradigm-Shift specific]
    """
    upstream_path_map = {}
    downstream_path_map = {}
    for focus_gene in focus_genes:
        for node in upstream_pathway.nodes:
            if node == focus_gene:
                continue
            all_paths = upstream_pathway.getAllPaths(node, focus_gene, search_distance)
            legal_paths = []        ## upstream paths must start with a protein
            for path in all_paths:
                if upstream_pathway.nodes[path[0][0]] == 'protein':
                    legal_paths.append(path)
            if len(legal_paths) > 0:
                if node not in upstream_path_map:
                    upstream_path_map[node] = deepcopy(legal_paths)
                else:
                    upstream_path_map[node] += deepcopy(legal_paths)
        for node in downstream_pathway.nodes:
            if node == focus_gene:
                continue
            all_paths = downstream_pathway.getAllPaths(focus_gene, node, search_distance)
            legal_paths = []        ## downstream paths must end with a protein and a transcriptional edge
            for path in all_paths:
                if downstream_pathway.nodes[path[-1][2]] == 'protein' and path[-1][1].startswith('-t'):
                    legal_paths.append(path)
            if len(legal_paths) > 0:
                if node not in downstream_path_map:
                    downstream_path_map[node] = deepcopy(legal_paths)
                else:
                    downstream_path_map[node] += deepcopy(legal_paths)
    return(upstream_path_map, downstream_path_map)

def computeFeatureRanks(score_map):
    rank_map = {}
    for attachment in score_map:
        rank_map[attachment] = {}
        ranked_features = []
        for feature in score_map[attachment]:
            try:
                fval = float(score_map[attachment][feature])
                if fval != fval:
                    raise ValueError
                ranked_features.append(feature)
            except ValueError:
                rank_map[attachment][feature] = 0
        ranked_features.sort(lambda x, y: cmp(abs(score_map[attachment][x]), abs(score_map[attachment][y])))
        for index, feature in enumerate(ranked_features):
            rank_map[attachment][feature] = float(index + 1)/len(ranked_features)
    return(rank_map)

def computeFeatureScores(data_map, positive_samples, negative_samples, upstream_pathway_map, downstream_pathway_map, method = 'variance'):
    def getScoreByVariance(data_frame, all_samples):
        score_map = {}
        for feature in data_frame.columns:
            score_map[feature] = computeMean(list(data_frame[feature].loc[all_samples]), return_sd = True, sample_sd = True)[1]
        return(score_map)
    def getScoreByWelchsTT(data_frame, positive_samples, negative_samples):
        score_map = {}
        for feature in data_frame.columns:
            score_map[feature] = computeWelchsT(list(data_frame[feature].loc[positive_samples]),
                                                list(data_frame[feature].loc[negative_samples]),
                                                alpha = 0.0,
                                                return_df = False)
        return(score_map)
	# def scoreSVM(matData, posSamples, negSamples):
		# import liblinear
		# import liblinearutil
		# svmLabels = []
		# svmData = []
		# featureMap = {}
		# for sample in (posSamples + negSamples):
			# if sample in posSamples:
				# svmLabels.append(1)
			# else:
				# svmLabels.append(-1)
			# svmData.append({})
			# for i, feature in enumerate(matData.keys()):
				# if i+1 not in featureMap:
					# featureMap[i+1] = feature
				# try:
					# svmData[-1][i+1] = float(matData[feature][sample])
				# except ValueError:
					# svmData[-1][i+1] = 0.0
		# prob = liblinear.problem(svmLabels, svmData)
		# param = liblinear.parameter('-s 3 -c 5 -q')
		# liblinearutil.save_model('model_file', liblinearutil.train(prob, param))
		# # m = liblinearutil.train(prob, param)
		# # testLabels = [] ## like svmLabels
		# # testData = [] ## like svmData
		# # for i in range(len(testLabels)):
		# #     x0, max_idx = gen_feature_nodearray(testData[i])
		# #     label = liblinear.predict(m, x0)
		# #     if label == testLabels[i]:
		# #         print "correct!"
		# weights = rList("model_file")[6:]
		# scoreMap = {}
		# for feature in featureMap.keys():
			# scoreMap[featureMap[feature]] = float(weights[feature-1])
		# return (scoreMap)
    def getScoreByGelNet(data_frame, positive_samples, negative_samples, upstream_pathway_map, downstream_pathway_map, l1 = 0.95, l2 = 0.05):
        """
        Refactoring of https://github.com/ucscCancer/pathway_tools/blob/master/scripts/diffuse.py [613be7d5baad339e8ddc852be6e10baff0cf8f9c]
        """
        from array import array
        from numpy import dot, genfromtxt, identity
        from scipy.sparse import coo_matrix
        from scipy.sparse.linalg import expm
        
        class SciPYKernel:
            def __init__(self):
                self.laplacian = None
                self.index2node = None
                self.kernel = None
                self.labels = None
            def readKernel(self, kernel_file):
                self.kernel = coo_matrix(genfromtxt(kernel_file, delimiter = '\t')[1:, 1:])
                f = open(kernel_file, 'r')
                self.labels = f.readline().rstrip().split('\t')[1:]
                f.close()
            def makeKernel(self, networkx_graph, diffusion_time = 0.1):
                ## parse the network, build the graph laplacian
                edges, nodes, node_out_degrees = self.parseGraph(networkx_graph)
                num_nodes = len(nodes)
                node_order = list(nodes)
                index2node = {}
                node2index = {}
                for i in range(0, num_nodes):
                    index2node[i] = node_order[i]
                    node2index[node_order[i]] = i
                ## construct the diagonals
                row = array('i')
                col = array('i')
                data = array('f')
                for i in range(0, num_nodes):
                    ## diag entries: out degree
                    degree = 0
                    if index2node[i] in node_out_degrees:
                        degree = node_out_degrees[index2node[i]]
                    ## append to the end
                    data.insert(len(data), degree)
                    row.insert(len(row), i)
                    col.insert(len(col), i)
                ## add edges
                for i in range(0, num_nodes):
                    for j in range(0, num_nodes):
                        if i == j:
                            continue
                        if (index2node[i], index2node[j]) not in edges:
                            continue
                        ## append index to i-th row, j-th column
                        row.insert(len(row), i)
                        col.insert(len(col), j)
                        ## -1 for laplacian edges
                        data.insert(len(data), -1)
                ## graph laplacian
                L = coo_matrix((data,(row, col)), shape=(num_nodes,num_nodes)).tocsc()
                self.laplacian = L
                self.index2node = index2node
                self.kernel = expm(-1*diffusion_time*L)
                self.labels = node_order
            def printLaplacian(self):
                cx = self.laplacian.tocoo()
                for i,j,v in zip(cx.row, cx.col, cx.data):
                    a = self.index2node[i]
                    b = self.index2node[j]
                    print '\t'.join([a,b,str(v)])
            def parseGraph(self, networkx_graph):
                edges = set()
                nodes = set()
                degrees = {}
                for source in networkx_graph.edge:
                    for target in networkx_graph.edge[source]:
                        if (source, target) in edges:
                            continue
                        edges.add((source, target))
                        edges.add((target, source))
                        nodes.add(source)
                        nodes.add(target)
                        if source not in degrees:
                            degrees[source] = 0
                        if target not in degrees:
                            degrees[target] = 0
                        degrees[source] += 1
                        degrees[target] += 1
                return (edges, nodes, degrees)
            def getRow(self, source):
                cx = heat_kernel.kernel.tocoo()
                edges = {}
                for i, j, v in zip(cx.row, cx.col, cx.data):
                    a = self.index2node[i]
                    b = self.index2node[j]
                    edges[(a, b)] = str(v)
                heats = {}
                for target in sorted(self.labels):
                    if (source, target) in edges:
                        heats[target] = float(edges[(source, target)])
                    else:
                        heats[target] = 0.0
                return(heats)
            def writeKernel(self, output_file):
                o = open(output_file, 'w')
                cx = self.kernel.tocoo()
                edges = {}
                for i, j, v in zip(cx.row, cx.col, cx.data):
                    a = self.index2node[i]
                    b = self.index2node[j]
                    edges[(a, b)] = str(v)
                ## iterate through rows
                ## sort labels in alphabetical order
                o.write('Key\t' + '\t'.join(sorted(self.labels)) + '\n')
                for nodeA in sorted(self.labels):
                    printstr = nodeA
                    # through columns
                    for nodeB in sorted(self.labels):
                        if (nodeA, nodeB) in edges:
                            printstr += '\t' + edges[(nodeA, nodeB)]
                        else:
                            printstr += '\t0'
                    o.write(printstr + '\n')
                o.close()
            def kernelMultiplyOne(self, vector):
                array = []
                ## loop over gene names in the network kernel: add the starting value if 
                ## it's present in the supplied input vector
                for label in self.labels:
                    if label in vector:
                        array.append(vector[label])
                    else:
                        array.append(0)
                ## take the dot product
                value = self.kernel*array
                return_vec = {}
                idx = 0
                for label in self.labels:
                    return_vec[label] = float(value[idx])
                    idx += 1
                return return_vec
            @staticmethod
            def getAngle(v1, v2):
                arry1 = []
                arry2 = []
                for key in v1:
                    arry1.append(float(v1[key]))
                    arry2.append(float(v2[key]))
                mag_1 = math.sqrt(dot(arry1,arry1))
                mag_2 = math.sqrt(dot(arry2,arry2))
                cos_theta = dot(arry1,arry2)/(mag_1*mag_2)
                return math.acos(cos_theta)
            def diffuse(self, vector, reverse=False):
                ## reverse is not used: heat diffusion is undirected
                diffused_vector = self.kernelMultiplyOne(vector)
                return diffused_vector
        
        class_map = {}
        for sample in positive_samples:
            class_map[sample] = 1
        for sample in negative_samples:
            class_map[sample] = 0
        class_vector = pandas.DataFrame(pandas.Series(class_map), columns = ['class'])
        sample_intersection = list(set(class_vector.index) & set(data_frame.index))
        sample_intersection.sort()
        upstream_score_map = {}
        for focus_gene in upstream_pathway_map:
            heat_kernel = SciPYKernel()
            heat_kernel.makeKernel(upstream_pathway_map[focus_gene].networkx())
            heat_vector = pandas.DataFrame(pandas.Series(heat_kernel.getRow(focus_gene)), columns = ['heat'])
            weight_vector = ((heat_vector/max(heat_vector.icol(0))) - 1.0).abs()**8
            weight_vector.columns = ['weight']
            feature_intersection = list(set(weight_vector.index) & set(data_frame.columns) - set([focus_gene]))
            feature_intersection.sort()
            penalty_matrix = pandas.DataFrame(identity(len(feature_intersection)))
            penalty_matrix.columns = feature_intersection
            penalty_matrix.index = feature_intersection
            data_frame[feature_intersection].loc[sample_intersection].to_csv('data.X.matrix', sep = '\t', index_label = 'id', na_rep = 0)
            class_vector.loc[sample_intersection].to_csv('class.y.vector', sep = '\t', index_label = 'id')
            weight_vector.loc[feature_intersection].to_csv('weight.d.vector', sep = '\t', index_label = 'id')
            penalty_matrix.to_csv('penalty.P.matrix', sep = '\t', index_label = 'id')
            os.system('Rscript ~/bin/paradigmshift/bin/gelnet.R')
            f = open('gelnet.w.vector', 'r')
            for line in f:
                if line.isspace():
                    continue
                pline = line.rstrip().split('\t')
                if pline[0] not in upstream_score_map:
                    upstream_score_map[pline[0]] = abs(float(pline[1]))
                else:
                    upstream_score_map[pline[0]] = max(upstream_score_map[pline[0]], abs(float(pline[1])))
            f.close()
        downstream_score_map = {}
        for focus_gene in downstream_pathway_map:
            heat_kernel = SciPYKernel()
            heat_kernel.makeKernel(downstream_pathway_map[focus_gene].networkx())
            heat_vector = pandas.DataFrame(pandas.Series(heat_kernel.getRow(focus_gene)), columns = ['heat'])
            weight_vector = ((heat_vector/max(heat_vector.icol(0))) - 1.0).abs()**8
            weight_vector.columns = ['weight']
            feature_intersection = list(set(weight_vector.index) & set(data_frame.columns) - set([focus_gene]))
            feature_intersection.sort()
            penalty_matrix = pandas.DataFrame(identity(len(feature_intersection)))
            penalty_matrix.columns = feature_intersection
            penalty_matrix.index = feature_intersection
            data_frame[feature_intersection].loc[sample_intersection].groupby(level=0, axis=0).mean().to_csv('data.X.matrix', sep = '\t', index_label = 'id', na_rep = 0)
            class_vector.loc[sample_intersection].to_csv('class.y.vector', sep = '\t', index_label = 'id')
            weight_vector.loc[feature_intersection].to_csv('weight.d.vector', sep = '\t', index_label = 'id')
            penalty_matrix.to_csv('penalty.P.matrix', sep = '\t', index_label = 'id')
            os.system('Rscript ~/bin/paradigmshift/bin/gelnet.R')
            f = open('gelnet.w.vector', 'r')
            for line in f:
                if line.isspace():
                    continue
                pline = line.rstrip().split('\t')
                if pline[0] not in downstream_score_map:
                    downstream_score_map[pline[0]] = abs(float(pline[1]))
                else:
                    downstream_score_map[pline[0]] = max(downstream_score_map[pline[0]], abs(float(pline[1])))
            f.close()
        combined_score_map = {}
        for feature in upstream_score_map:
            if feature in downstream_score_map:
                combined_score_map[feature] = 0.0
            else:
                combined_score_map[feature] = upstream_score_map[feature]
        for feature in downstream_score_map:
            if feature in upstream_score_map:
                combined_score_map[feature] = 0.0
            else:
                combined_score_map[feature] = downstream_score_map[feature]
        return(combined_score_map)
    
    all_samples = list(set(positive_samples) | set(negative_samples))
    score_map = {}
    for attachment in data_map:
        if method == 'variance':
            score_map[attachment] = getScoreByVariance(data_map[attachment], all_samples)
        elif method == 'tt':
            score_map[attachment] = getScoreByWelchsTT(data_map[attachment], positive_samples, negative_samples)
        elif method == 'gelnet':
            score_map[attachment] = getScoreByGelNet(data_map[attachment], positive_samples, negative_samples, upstream_pathway_map, downstream_pathway_map, l1 = 0.95, l2 = 0.05)
    score_frame = pandas.DataFrame(score_map)
    score_frame.to_csv('feature.scores', sep = '\t', na_rep = 'NA')
    
    if method in ['gelnet']:
        return(score_map)
    else:
        return(computeFeatureRanks(score_map))

def collapsePath(path, max_distance = 2):
    if len(path) <= max_distance:
        return(path)
    collapsed_node = "//".join([edge[2] for edge in path[:len(path) - max_distance + 1]])
    if [edge[1][-1] for edge in path[:len(path) - max_distance + 1]].count("|") % 2 == 0:
        if len(path[len(path) - max_distance][1]) == 3 and path[len(path) - max_distance][1].startswith("-t"):
            collapsed_interaction = "-t>"
        else:
            collapsed_interaction = "-a>"
    else:
        if len(path[len(path) - max_distance][1]) == 3 and path[len(path) - max_distance][1].startswith("-t"):
            collapsed_interaction = "-t|"
        else:
            collapsed_interaction = "-a|"
    collapsed_path = [(path[0][0], collapsed_interaction, collapsed_node), (collapsed_node, path[len(path) - max_distance + 1][1], path[len(path) - max_distance + 1][2])] + path[len(path) - max_distance + 2:]
    return(collapsed_path)

def getSelectedNeighborhood(focus_node, focus_genes, data_map, positive_samples, negative_samples, upstream_path_map, downstream_path_map, upstream_pathway_map, downstream_pathway_map, base_upstream_pathway, base_downstream_pathway, threshold = 0.84, cost = 0.0, method = 'variance' , max_distance = 2):
    def getUpstreamValue(feature, value_map):
        value_list = []
        if 'mrna' in value_map:
            if feature in value_map['mrna']:
                try:
                    fval = abs(float(value_map['mrna'][feature]))
                    if fval != fval:
                        raise ValueError
                except ValueError:
                    fval = 0
            else:
                fval = 0
            value_list.append(fval)
        if 'active' in value_map:
            if feature in value_map['active']:
                try:
                    fval = abs(float(value_map['active'][feature]))
                    if fval != fval:
                        raise ValueError
                except ValueError:
                    fval = 0
            else:
                fval = 0
            value_list.append(fval)
        return(max(value_list))
    def getDownstreamValue(feature, value_map):
        value_list = []
        if 'mrna' in value_map:
            if feature in value_map['mrna']:
                try:
                    fval = abs(float(value_map['mrna'][feature]))
                    if fval != fval:
                        raise ValueError
                except ValueError:
                    fval = 0
            else:
                fval = 0
            value_list.append(fval)
        return(max(value_list))
    
    score_map = computeFeatureScores(data_map, positive_samples, negative_samples, upstream_pathway_map, downstream_pathway_map, method = method)
    logger('## upstream neighborhood\n', file = 'selection.log')
    selected_upstream_pathway = Pathway( ({}, {}) )
    for focus_gene in focus_genes:
        selected_upstream_pathway.nodes[focus_gene] = base_upstream_pathway.nodes[focus_gene]
    if len(focus_genes) > 1:
        selected_upstream_pathway.nodes[focus_node] = 'abstract'
        for focus_gene in focus_genes:
            selected_upstream_pathway.interactions[focus_gene] = {}
            selected_upstream_pathway.interactions[focus_gene][focus_node] = '-ap>'
    selected_upstream_features = []
    ranked_upstream_features = []
    for feature in upstream_path_map:
        if getUpstreamValue(feature, score_map) > 0:
            ranked_upstream_features.append(feature)
    ranked_upstream_features.sort(lambda x, y: cmp(getUpstreamValue(y, score_map), getUpstreamValue(x, score_map)))
    if len(ranked_upstream_features) > 0:
        while (len(selected_upstream_features) < 4) or (getUpstreamValue(ranked_upstream_features[0], score_map) > threshold + cost*len(selected_upstream_features)):
            current_feature = ranked_upstream_features.pop(0)
            logger('%s\t%s\t%s\t%s\n' % (current_feature, getUpstreamValue(current_feature, score_map), len(selected_upstream_features), upstream_path_map[current_feature]), file = 'selection.log')
            selected_upstream_features.append(current_feature)
            for path in upstream_path_map[current_feature]:
                collapsed_path = collapsePath(path, max_distance = max_distance)
                for edge in collapsed_path:
                    if edge[0] not in selected_upstream_pathway.nodes:
                        if edge[0] not in base_upstream_pathway.nodes:
                            selected_upstream_pathway.nodes[edge[0]] = "abstract"
                        else:
                            selected_upstream_pathway.nodes[edge[0]] = base_upstream_pathway.nodes[edge[0]]
                    if edge[2] not in selected_upstream_pathway.nodes:
                        if edge[2] not in base_upstream_pathway.nodes:
                            selected_upstream_pathway.nodes[edge[2]] = "abstract"
                        else:
                            selected_upstream_pathway.nodes[edge[2]] = base_upstream_pathway.nodes[edge[2]]
                    if edge[0] not in selected_upstream_pathway.interactions:
                        selected_upstream_pathway.interactions[edge[0]] = {}
                    selected_upstream_pathway.interactions[edge[0]][edge[2]] = edge[1]
            if len(ranked_upstream_features) == 0:
                break
    logger('## downstream neighborhood\n', file = 'selection.log')
    selected_downstream_pathway = Pathway( ({}, {}) )
    for focus_gene in focus_genes:
        selected_downstream_pathway.nodes[focus_gene] = base_downstream_pathway.nodes[focus_gene]
    if len(focus_genes) > 1:
        selected_downstream_pathway.nodes[focus_node] = 'abstract'
        selected_downstream_pathway.interactions[focus_node] = {}
        for focus_gene in focus_genes:
            selected_downstream_pathway.interactions[focus_node][focus_gene] = '-ap>'
    selected_downstream_features = []
    ranked_downstream_features = []
    for feature in downstream_path_map:
        if getDownstreamValue(feature, score_map) > 0:
            ranked_downstream_features.append(feature)
    ranked_downstream_features.sort(lambda x, y: cmp(getDownstreamValue(y, score_map), getDownstreamValue(x, score_map)))
    if len(ranked_downstream_features) > 0:
        while (len(selected_downstream_features) < 4) or (getDownstreamValue(ranked_downstream_features[0], score_map) > threshold + cost*len(selected_downstream_features)):
            current_feature = ranked_downstream_features.pop(0)
            logger('%s\t%s\t%s\t%s\n' % (current_feature, getDownstreamValue(current_feature, score_map), len(selected_downstream_features), downstream_path_map[current_feature]), file = 'selection.log')
            selected_downstream_features.append(current_feature)
            for path in downstream_path_map[current_feature]:
                collapsed_path = collapsePath(path, max_distance = max_distance)
                for edge in collapsed_path:
                    if edge[0] not in selected_downstream_pathway.nodes:
                        if edge[0] not in base_downstream_pathway.nodes:
                            selected_downstream_pathway.nodes[edge[0]] = "abstract"
                        else:
                            selected_downstream_pathway.nodes[edge[0]] = base_downstream_pathway.nodes[edge[0]]
                    if edge[2] not in selected_downstream_pathway.nodes:
                        if edge[2] not in base_downstream_pathway.nodes:
                            selected_downstream_pathway.nodes[edge[2]] = "abstract"
                        else:
                            selected_downstream_pathway.nodes[edge[2]] = base_downstream_pathway.nodes[edge[2]]
                    if edge[0] not in selected_downstream_pathway.interactions:
                        selected_downstream_pathway.interactions[edge[0]] = {}
                    selected_downstream_pathway.interactions[edge[0]][edge[2]] = edge[1]
            if len(ranked_downstream_features) == 0:
                break
    return(selected_upstream_pathway, selected_downstream_pathway, ((len(selected_upstream_features) >= 4) and (len(selected_downstream_features) >= 4)))

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

def generateData(focus_genes, upstream_features, downstream_features, allow_features, include_samples, data_files, null_size = 0, random_seed = 1):
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
    for null in range(1, null_size + 1):
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

def generateBatchedData(focus_genes, upstream_features, downstream_features, allow_features, include_samples, data_files, null_size = 0, batch_size = 50, random_seed = 1):
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
    for null in range(1, null_size + 1):
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

def readParadigm(input_file, include_features = None):
    """
    Read in a Paradigm format output file [2014-3-10]
    """
    likelihood_map = {}
    ipl_map = {}
    f = open(input_file, 'r')
    for line in f:
        if line.isspace():
            continue
        if line.startswith(">"):
            pline = re.split("[= ]", line.rstrip())
            sample = pline[1]
            likelihood_map[sample] = float(pline[3])
            ipl_map[sample] = {}
        else:
            pline = line.rstrip().split('\t')
            feature = pline[0]
            if include_features != None:
                if feature not in include_features:
                    continue
            ipl_map[sample][feature] = float(pline[1])
    f.close()
    ipl_data = pandas.DataFrame(ipl_map)
    return(likelihood_map, ipl_data)

def convertListToFloat(input_list):
    """
    Converts a list of strings to floats and removes non-float values [2014-3-1]
    """
    float_list = []
    for element in input_list:
        try:
            float_val = float(element)
            if float_val != float_val:
                raise ValueError
            float_list.append(float_val)
        except ValueError:
            continue
    return(float_list)

def computeMean(input_list, null = 'NA', return_sd = False, sample_sd = True):
    """
    Computes mean and optionally the sd [2014-3-1]
    Dependencies: convertListToFloat
    """
    float_list = convertListToFloat(input_list)
    if len(float_list) == 0:
        mean = null
        sd = null
    else:
        mean = sum(float_list)/float(len(float_list))
        sd = 0.0
        for element in float_list:
            sd += (element - mean)**2
        if len(float_list) > 1:
            if sample_sd:
                sd = math.sqrt(sd/(len(float_list) - 1))
            else:
                sd = math.sqrt(sd/len(float_list))
        else:
            sd = 0.0
    if return_sd:
        return(mean, sd)
    else:
        return(mean)

def normalizeDataFrame(data_frame, include_samples = None):
    """
    Normalizes a pandas dataframe per feature based on the stats of a set of samples [2014-3-11]
    Dependencies: convertListToFloat, computeMean
    """
    normalize_data_frame = data_frame.copy()
    if include_samples:
        normalize_samples = include_samples
    else:
        normalize_samples = data_frame.columns
    for feature in normalize_data_frame.index:
        (mean, sd) = computeMean(list(data_frame[normalize_samples].loc[feature]), return_sd = True)
        for sample in normalize_data_frame.columns:
            normalize_data_frame[sample][feature] = (data_frame[sample][feature] - mean)/sd
    return(normalize_data_frame)

def getConfusion(index, ranked_features, classification_map):
    """
    Computes the confusion matrix given a split on a ranked list [2014-3-11]
    """
    (tp, fn, fp, tn) = (0, 0, 0, 0)
    for feature in ranked_features[:index]:
        if classification_map[feature] == 1:
            tp += 1
        elif classification_map[feature] == 0:
            fp += 1
    for feature in ranked_features[index:]:
        if classification_map[feature] == 1:
            fn += 1
        elif classification_map[feature] == 0:
            tn += 1
    return(tp, fn, fp, tn)

def computeAUC(ranked_features, score_map, classification_map):
    """
    Computes the AUC and points for an plotting [2014-3-11]
    Dependencies: getConfusion
    """
    auc = 0.0
    points = []
    x = 0.0
    index = 0
    while (index <= len(ranked_features)):
        (tp, fn, fp, tn) = getConfusion(index, ranked_features, classification_map)
        try:
            tpr = tp/float(tp+fn)
            fpr = fp/float(fp+tn)
        except ZeroDivisionError:
            return('NA', [])
        points.append( (fpr, tpr) )
        if fpr > x:
            dx = fpr - x
            x = fpr
            auc += dx*tpr
        index += 1
        if index < len(ranked_features):
            while (score_map[ranked_features[index - 1]] == score_map[ranked_features[index]]):
                index += 1
                if index == len(ranked_features):
                    break
        elif index == len(ranked_features):
            pass
        else:
            break
    return(auc, points)

def computeWelchsT(values_1, values_2, alpha = 0.0, return_df = False):
    """
    Computes the t_statistic for a the Welch's T-test for unpaired distributions of unequal variance [2014-3-11]
    Would be better to use scipy.stats.ttest_ind
    Dependencies: computeMean
    """
    (mean_1, sd_1) = computeMean(values_1, return_sd = True)
    (mean_2, sd_2) = computeMean(values_2, return_sd = True)
    try:
        t_statistic = (mean_1 - mean_2)/(math.sqrt((sd_1**2)/len(values_1) + (sd_2**2)/len(values_2)) + alpha)
    except ZeroDivisionError:
        t_statistic = 'NA'
    except TypeError:
        t_statistic = 'NA'
    if return_df:
        df = (((sd_1**2)/len(values_1) + (sd_2**2)/len(values_2))**2)/((((sd_1**2)/len(values_1))**2)/(len(values_1) - 1) + (((sd_2**2)/len(values_2))**2)/(len(values_2) - 1))
        return(t_statistic, df)
    else:
        return(t_statistic)

def computeSeparation(positive_samples, negative_samples, pshift_map, method = 'tt'):
    """
    Computes the separation statistic for a Paradigm-Shift run [Paradigm-Shift specific]
    Dependencies: computeWelchsT, computeMean
    """
    positive_values = []
    negative_values = []
    for sample in positive_samples + negative_samples:
        if pshift_map[sample] == 'NA':
            continue
        if sample in positive_samples:
            positive_values.append(pshift_map[sample])
        elif sample in negative_samples:
            negative_values.append(pshift_map[sample])
    if method == 'tt':
        separation_statistic = computeWelchsT(positive_values, negative_values)
    else:
        separation_statistic = 'NA'
    return(separation_statistic)

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

class queueAnalyses(Target):
    def __init__(self, analysis_list, paradigm_setup, global_pathway, parameters, directory, last_report_list = [], run_analyses = None, total_analyses = None):
        Target.__init__(self, time=10000)
        self.analysis_list = analysis_list
        self.paradigm_setup = paradigm_setup
        self.global_pathway = global_pathway
        self.parameters = parameters
        self.directory = directory
        self.last_report_list = last_report_list
        self.run_analyses = run_analyses
        self.total_analyses = total_analyses
    def run(self):
        os.chdir(self.directory)
        
        ## queue analyses
        report_list = deepcopy(self.last_report_list)
        if self.total_analyses is None:
            self.run_analyses = 0
            self.total_analyses = len(self.analysis_list)
        if not os.path.exists('analysis'):
            os.mkdir('analysis')
        if len(self.analysis_list) > 0:
            analysis = self.analysis_list[0]
            if not os.path.exists('analysis/%s' % (analysis.directory)):
                logger('Running analysis on %s (%s/%s)\n' % (analysis.analysis_name, self.run_analyses + 1, self.total_analyses), file = 'analysis/progress.log')
                os.mkdir('analysis/%s' % (analysis.directory))
                report_list.append(analysis.directory)
                self.addChildTarget(branchFolds(analysis,
                                                self.paradigm_setup,
                                                self.global_pathway,
                                                self.parameters,
                                                self.directory))
            else:
                logger('Already performed analysis on %s (%s/%s)\n' % (analysis.analysis_name, self.run_analyses + 1, self.total_analyses), file = 'analysis/progress.log')
            self.setFollowOnTarget(queueAnalyses(self.analysis_list[1:],
                                                 self.paradigm_setup,
                                                 self.global_pathway,
                                                 self.parameters,
                                                 self.directory,
                                                 last_report_list = report_list,
                                                 run_analyses = self.run_analyses + 1,
                                                 total_analyses = self.total_analyses))
        else:
            if self.parameters.report_directory != None:
                self.setFollowOnTarget(makeReport(report_list,
                                                  self.parameters.report_directory,
                                                  self.directory))

class branchAnalyses(Target):
    def __init__(self, analysis_list, paradigm_setup, global_pathway, parameters, directory, run_analyses = None, total_analyses = None):
        Target.__init__(self, time=10000)
        self.analysis_list = analysis_list
        self.paradigm_setup = paradigm_setup
        self.global_pathway = global_pathway
        self.parameters = parameters
        self.directory = directory
        self.run_analyses = run_analyses
        self.total_analyses = total_analyses
    def run(self):
        os.chdir(self.directory)
        
        ## branch analyses
        report_list = []
        if self.total_analyses is None:
            self.run_analyses = 0
            self.total_analyses = len(self.analysis_list)
        if not os.path.exists('analysis'):
            os.mkdir('analysis')
        for analysis in self.analysis_list:
            if not os.path.exists('analysis/%s' % (analysis.directory)):
                logger('Running analysis on %s (%s/%s)\n' % (analysis.analysis_name, self.run_analyses + 1, self.total_analyses), file = 'analysis/progress.log')
                os.mkdir('analysis/%s' % (analysis.directory))
                report_list.append(analysis.directory)
                self.addChildTarget(branchFolds(analysis,
                                                self.paradigm_setup,
                                                self.global_pathway,
                                                self.parameters,
                                                self.directory))
            else:
                logger('Already performed analysis on %s (%s/%s)\n' % (analysis.analysis_name, self.run_analyses + 1, self.total_analyses), file = 'analysis/progress.log')
            self.run_analyses += 1
        if self.parameters.report_directory != None:
            self.setFollowOnTarget(makeReport(report_list,
                                              self.parameters.report_directory,
                                              self.directory))

class branchFolds(Target):
    def __init__(self, analysis, paradigm_setup, global_pathway, parameters, directory):
        Target.__init__(self, time=10000)
        self.analysis = analysis
        self.paradigm_setup = paradigm_setup
        self.global_pathway = global_pathway
        self.parameters = parameters
        self.directory = directory
    def run(self):
        ps_directory = '%s/analysis/%s' % (self.directory, self.analysis.directory)
        os.chdir(ps_directory)
        random.seed(self.parameters.random_seed+123454321)
        
        ## cross validation
        if self.parameters.cross_validation:
            logger('Performing cross-validation ...\n', file = 'progress.log')
            fold_map = {}
            for round in range(1, self.parameters.n_rounds + 1):
                fold_map[round] = {}
                for fold in range(1, self.parameters.m_folds + 1):
                    fold_map[round][fold] = []
                select_positive = deepcopy(self.analysis.positive_samples)
                select_negative = deepcopy(self.analysis.negative_samples)
                while len(select_positive) + len(select_negative) > 0:
                    for fold in range(1, self.parameters.m_folds + 1):
                        if len(select_positive) > 0:
                            fold_map[round][fold].append(select_positive.pop(random.randint(0, len(select_positive) - 1)))
                        elif len(select_negative) > 0:
                            fold_map[round][fold].append(select_negative.pop(random.randint(0, len(select_negative) - 1)))
            for round in range(1, self.parameters.n_rounds + 1):
                for fold in range(1, self.parameters.m_folds + 1):
                    fold_index = (round - 1)*self.parameters.m_folds + fold
                    os.mkdir('fold%s' % (fold_index))
                    self.addChildTarget(branchParameters(fold_index,
                                                         self.analysis,
                                                         fold_map[round][fold],
                                                         self.paradigm_setup,
                                                         self.global_pathway,
                                                         self.parameters,
                                                         self.directory))
        
        ## run final
        os.mkdir('final')
        self.setFollowOnTarget(branchParameters(0,
                                                self.analysis,
                                                self.paradigm_setup.samples,
                                                self.paradigm_setup,
                                                self.global_pathway,
                                                self.parameters,
                                                self.directory))

class branchParameters(Target):
    def __init__(self, fold, analysis, training_samples, paradigm_setup, global_pathway, parameters, directory):
        Target.__init__(self, time=10000)
        self.fold = fold
        self.analysis = analysis
        self.training_samples = training_samples
        self.paradigm_setup = paradigm_setup
        self.global_pathway = global_pathway
        self.parameters = parameters
        self.directory = directory
    def run(self):
        if self.fold == 0:
            ps_directory = '%s/analysis/%s/final' % (self.directory, self.analysis.directory)
            os.chdir(ps_directory)
            logger('Constructing final model ...\n', file = 'progress.log')
        else:
            ps_directory = '%s/analysis/%s/fold%s' % (self.directory, self.analysis.directory, self.fold)
            os.chdir(ps_directory)
        
        ## average cross validation aucs for final run
        if self.fold == 0:
            auc_list = []
            auc_lines = []
            auc_average = None
            for round in range(1, self.parameters.n_rounds + 1):
                for fold in range(1, self.parameters.m_folds + 1):
                    fold_index = (round - 1)*self.parameters.m_folds + fold
                    if os.path.exists('../fold%s/auc.tab' % (fold_index)):
                        f = open('../fold%s/auc.tab' % (fold_index), 'r')
                        line = f.readline()
                        f.close()
                        (auc_train, auc_test, auc_params) = line.rstrip().split('\t')
                    else:
                        (auc_train, auc_test, auc_params) = ('---', '---', '---')
                    auc_list.append(auc_test)
                    auc_lines.append('%s\t%s\t%s\t%s' % (fold_index, auc_train, auc_test, auc_params))
            auc_average = computeMean(auc_list)
            o = open('../average_auc.tab', 'w')
            o.write('> %s\tMean(AUC) = %s\n' % (self.analysis.analysis_name, auc_average))
            o.write('# fold\ttrain\ttest\tparameters\n')
            o.write('%s\n' % ('\n'.join(auc_lines)))
            o.close()
            if self.parameters.cross_validation:
                if auc_average == 'NA':
                    os.chdir('..')
                    shutil.rmtree('final')
                    for round in range(1, self.parameters.n_rounds + 1):
                        for fold in range(1, self.parameters.m_folds + 1):
                            fold_index = (round - 1)*self.parameters.m_folds + fold
                            shutil.rmtree('fold%s' % (fold_index))
                    logger('Terminating analysis (did not pass cross-validation) ...\n', file = 'progress.log')
                    return
                elif auc_average < self.parameters.cross_validation_threshold:
                    os.chdir('..')
                    shutil.rmtree('final')
                    for round in range(1, self.parameters.n_rounds + 1):
                        for fold in range(1, self.parameters.m_folds + 1):
                            fold_index = (round - 1)*self.parameters.m_folds + fold
                            shutil.rmtree('fold%s' % (fold_index))
                    logger('Terminating analysis (did not pass cross-validation) ...\n', file = 'progress.log')
                    return
        
        ## branch parameters
        for threshold in self.parameters.threshold:
            for cost in self.parameters.cost:
                for method in self.parameters.selection_method:
                    current_parameters = [threshold, cost, method]
                    os.mkdir('param_%s' % ('_'.join([str(parameter) for parameter in current_parameters])))
                    self.addChildTarget(selectNeighborhood(self.fold,
                                                           self.analysis,
                                                           self.training_samples,
                                                           self.paradigm_setup,
                                                           self.global_pathway,
                                                           current_parameters,
                                                           self.parameters,
                                                           self.directory))
        
        ## evaluate models
        self.setFollowOnTarget(compareParameters(self.fold,
                                                 self.analysis,
                                                 self.parameters,
                                                 self.directory))

class selectNeighborhood(Target):
    def __init__(self, fold, analysis, training_samples, paradigm_setup, global_pathway, current_parameters, parameters, directory):
        Target.__init__(self, time=10000)
        self.fold = fold
        self.analysis = analysis
        self.training_samples = training_samples
        self.paradigm_setup = paradigm_setup
        self.global_pathway = global_pathway
        self.current_parameters = current_parameters
        self.parameters = parameters
        self.directory = directory
    def run(self):
        if self.fold == 0:
            ps_directory = '%s/analysis/%s/final/param_%s' % (self.directory, self.analysis.directory, '_'.join([str(parameter) for parameter in self.current_parameters]))
            os.chdir(ps_directory)
        else:
            ps_directory = '%s/analysis/%s/fold%s/param_%s' % (self.directory, self.analysis.directory, self.fold, '_'.join([str(parameter) for parameter in self.current_parameters]))
            os.chdir(ps_directory)
        logger('Selecting neighborhood ...\n', file = 'progress.log')
        
        ## use trained model if it exists
        selected_upstream = None
        selected_downstream = None
        if self.parameters.model_directory is not None:
            model_path = ''
            if self.parameters.model_directory.startswith('/'):
                model_path = '%s/%s' % (self.parameters.model_directory, self.analysis.focus_node)
            else:
                model_path = '%s/%s/%s' % (self.directory, self.parameters.model_directory, self.analysis.focus_node)
            if (os.path.exists('%s/upstream_pathway.tab' % (model_path))) and (os.path.exists('%s/downstream_pathway.tab' % (model_path))):
                logger('Using trained model ...\n', file = 'progress.log')
                selected_upstream = Pathway('%s/upstream_pathway.tab' % (model_path))
                selected_downstream = Pathway('%s/downstream_pathway.tab' % (model_path))
                selection_pass = True
        
        ## get upstream and downstream base neighborhoods per complex, then combine
        ## this allows for a more sophisticated breakdown of downstream targets by complex
        ## but for now it remains agnostic
        if (selected_upstream is None) and (selected_downstream is None):
            (group_index, group_pathways) = getFullNeighborhood(self.analysis.focus_genes,
                                                                self.global_pathway,
                                                                search_distance = self.parameters.search_distance)
            upstream_pathway_map = {}
            downstream_pathway_map = {}
            for focus_gene in self.analysis.focus_genes:
                upstream_pathway_map[focus_gene] = Pathway( ({focus_gene : self.global_pathway.nodes[focus_gene]}, {}) )
                downstream_pathway_map[focus_gene] = Pathway( ({focus_gene : self.global_pathway.nodes[focus_gene]}, {}) )
                for index in group_pathways[focus_gene]:
                    upstream_pathway_map[focus_gene].appendPathway(group_pathways[focus_gene][index][0])
                    downstream_pathway_map[focus_gene].appendPathway(group_pathways[focus_gene][index][1])
                upstream_pathway_map[focus_gene].writeSPF('upstream_base.%s.tab' % (focus_gene))
                downstream_pathway_map[focus_gene].writeSPF('downstream_base.%s.tab' % (focus_gene))
                # upstream_pathway_map[focus_gene].writeSIF('upstream_base.%s.sif' % (focus_gene))
                # downstream_pathway_map[focus_gene].writeSIF('downstream_base.%s.sif' % (focus_gene))
            
            ## identify all interaction paths relevant to the Paradigm-Shift task
            base_upstream_pathway = Pathway( ({}, {}) )
            base_downstream_pathway = Pathway( ({}, {}) )
            for focus_gene in self.analysis.focus_genes:
                base_upstream_pathway.appendPathway(upstream_pathway_map[focus_gene])
                base_downstream_pathway.appendPathway(downstream_pathway_map[focus_gene])
            (upstream_path_map, downstream_path_map) = getRelevantPaths(self.analysis.focus_genes,
                                                                        base_upstream_pathway,
                                                                        base_downstream_pathway,
                                                                        search_distance = self.parameters.search_distance)
            
            ## score and select features
            data_map = {}
            assert(len(self.paradigm_setup.mrna) > 0)
            data_map['mrna'] = pandas.read_csv(self.paradigm_setup.mrna[0], sep = '\t', index_col = 0)
            if len(self.paradigm_setup.active) > 0:
                data_map['active'] = pandas.read_csv(self.paradigm_setup.active[0], sep = '\t', index_col = 0)
            positive_samples = list(set(self.training_samples) & set(self.analysis.positive_samples))
            negative_samples = list(set(self.training_samples) & set(self.analysis.negative_samples))
            (selected_upstream, selected_downstream, selection_pass) = getSelectedNeighborhood(self.analysis.focus_node,
                                                   self.analysis.focus_genes,
                                                   data_map,
                                                   positive_samples,
                                                   negative_samples,
                                                   upstream_path_map,
                                                   downstream_path_map,
                                                   upstream_pathway_map,
                                                   downstream_pathway_map,
                                                   base_upstream_pathway,
                                                   base_downstream_pathway,
                                                   threshold = self.current_parameters[0],
                                                   cost = self.current_parameters[1],
                                                   method = self.current_parameters[2])
        
        if not selection_pass:
            o = open('auc.tab', 'w')
            o.write('---\t---\n')
            o.close()
        else:
            os.mkdir('data')
            data_files = self.paradigm_setup.genome + self.paradigm_setup.mrna + self.paradigm_setup.protein + self.paradigm_setup.active
            upstream_features = list(set(selected_upstream.nodes.keys()) & set(self.paradigm_setup.features))
            downstream_features = list(set(selected_downstream.nodes.keys()) & set(self.paradigm_setup.features))
            if self.paradigm_setup.public:
                if self.fold == 0:
                    generateBatchedData(self.analysis.focus_genes,
                                        upstream_features,
                                        downstream_features,
                                        self.paradigm_setup.features,
                                        self.paradigm_setup.samples,
                                        data_files,
                                        null_size = self.paradigm_setup.null_size,
                                        batch_size = self.paradigm_setup.batch_size,
                                        random_seed = self.parameters.random_seed + self.fold)
                else:
                    generateBatchedData(self.analysis.focus_genes,
                                        upstream_features,
                                        downstream_features,
                                        self.paradigm_setup.features,
                                        self.paradigm_setup.samples,
                                        data_files,
                                        null_size = 0,
                                        batch_size = self.paradigm_setup.batch_size,
                                        random_seed = self.parameters.random_seed + self.fold)
            else:
                if self.fold == 0:
                    generateData(self.analysis.focus_genes,
                                 upstream_features,
                                 downstream_features,
                                 self.paradigm_setup.features,
                                 self.paradigm_setup.samples,
                                 data_files,
                                 null_size = self.paradigm_setup.null_size,
                                 random_seed = self.parameters.random_seed + self.fold)
                else:
                    generateData(self.analysis.focus_genes,
                                 upstream_features,
                                 downstream_features,
                                 self.paradigm_setup.features,
                                 self.paradigm_setup.samples,
                                 data_files,
                                 null_size = 0,
                                 random_seed = self.parameters.random_seed + self.fold)
            selected_upstream.writeSPF('upstream_pathway.tab')
            selected_downstream.writeSPF('downstream_pathway.tab')
            if self.fold == 0:
                self.addChildTarget(runParadigm(self.analysis,
                                                self.paradigm_setup,
                                                ps_directory))
            else:
                self.addChildTarget(runParadigm(self.analysis,
                                                self.paradigm_setup,
                                                ps_directory,
                                                null_size = 0))
            self.setFollowOnTarget(computeShifts(self.fold,
                                                 self.analysis,
                                                 self.training_samples,
                                                 self.paradigm_setup,
                                                 self.current_parameters,
                                                 self.parameters,
                                                 self.directory))

class runParadigm(Target):
    def __init__(self, analysis, paradigm_setup, directory, null_size = None):
        Target.__init__(self, time=10000)
        self.analysis = analysis
        self.paradigm_setup = paradigm_setup
        self.directory = directory
        if null_size is None:
            self.null_size = self.paradigm_setup.null_size
        else:
            self.null_size = null_size
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
        ## run Paradigm (real and null)
        batches = getBatchCount(len(self.paradigm_setup.samples), batch_size = self.paradigm_setup.batch_size)
        if batches == 0:
            self.addChildTarget(jtCmd('%s -p upstream_pathway.tab -c config.txt -b %s -o %s' % (paradigm_executable, 'data/up_', 'paradigm/%s_upstream.fa' % (self.analysis.focus_node)), self.directory, file = 'jobs.list'))
            self.addChildTarget(jtCmd('%s -p downstream_pathway.tab -c config.txt -b %s -o %s' % (paradigm_executable, 'data/down_', 'paradigm/%s_downstream.fa' % (self.analysis.focus_node)), self.directory, file = 'jobs.list'))
            for null in range(1, self.null_size + 1):
                self.addChildTarget(jtCmd('%s -p upstream_pathway.tab -c config.txt -b %s -o %s' % (paradigm_executable, 'data/up_N%s_' % (null), 'paradigm/N%s_%s_upstream.fa' % (null, self.analysis.focus_node)), self.directory, file = 'jobs.list'))
                self.addChildTarget(jtCmd('%s -p downstream_pathway.tab -c config.txt -b %s -o %s' % (paradigm_executable, 'data/down_N%s_' % (null), 'paradigm/N%s_%s_downstream.fa' % (null, self.analysis.focus_node)), self.directory, file = 'jobs.list'))
        elif self.paradigm_setup.public:
            os.mkdir('outputFiles')
            for b in range(batches):
                self.addChildTarget(jtCmd('%s -p upstream_pathway.tab -c config.txt -b %s -o %s' % (paradigm_executable, 'data/up_b%s_%s_' % (b, batches), 'outputFiles/%s_upstream_b%s_%s.fa' % (self.analysis.focus_node, b, batches)), self.directory, file = 'jobs.list'))
                self.addChildTarget(jtCmd('%s -p downstream_pathway.tab -c config.txt -b %s -o %s' % (paradigm_executable, 'data/down_b%s_%s_' % (b, batches), 'outputFiles/%s_downstream_b%s_%s.fa' % (self.analysis.focus_node, b, batches)), self.directory, file = 'jobs.list'))
                for null in range(1, self.null_size + 1):
                    self.addChildTarget(jtCmd('%s -p upstream_pathway.tab -c config.txt -b %s -o %s' % (paradigm_executable, 'data/up_N%s_b%s_%s_' % (null, b, batches), 'outputFiles/N%s_%s_upstream_b%s_%s.fa' % (null, self.analysis.focus_node, b, batches)), self.directory, file = 'jobs.list'))
                    self.addChildTarget(jtCmd('%s -p downstream_pathway.tab -c config.txt -b %s -o %s' % (paradigm_executable, 'data/down_N%s_b%s_%s_' % (null, b, batches), 'outputFiles/N%s_%s_downstream_b%s_%s.fa' % (null, self.analysis.focus_node, b, batches)), self.directory, file = 'jobs.list'))
        else:
            os.mkdir('outputFiles')
            for b in range(batches):
                self.addChildTarget(jtCmd('%s -p upstream_pathway.tab -c config.txt -b %s -o %s -s %s,%s' % (paradigm_executable, 'data/up_', 'outputFiles/%s_upstream_b%s_%s.fa' % (self.analysis.focus_node, b, batches), b, batches), self.directory, file = 'jobs.list'))
                self.addChildTarget(jtCmd('%s -p downstream_pathway.tab -c config.txt -b %s -o %s -s %s,%s' % (paradigm_executable, 'data/down_', 'outputFiles/%s_downstream_b%s_%s.fa' % (self.analysis.focus_node, b, batches), b, batches), self.directory, file = 'jobs.list'))
                for null in range(1, self.null_size + 1):
                    self.addChildTarget(jtCmd('%s -p upstream_pathway.tab -c config.txt -b %s -o %s -s %s,%s' % (paradigm_executable, 'data/up_N%s_' % (null), 'outputFiles/N%s_%s_upstream_b%s_%s.fa' % (null, self.analysis.focus_node, b, batches), b, batches), self.directory, file = 'jobs.list'))
                    self.addChildTarget(jtCmd('%s -p downstream_pathway.tab -c config.txt -b %s -o %s -s %s,%s' % (paradigm_executable, 'data/down_N%s_' % (null), 'outputFiles/N%s_%s_downstream_b%s_%s.fa' % (null, self.analysis.focus_node, b, batches), b, batches), self.directory, file = 'jobs.list'))
        self.setFollowOnTarget(collectParadigm(self.analysis, self.paradigm_setup, self.directory, null_size = self.null_size))

class collectParadigm(Target):
    def __init__(self, analysis, paradigm_setup, directory, null_size = None):
        Target.__init__(self, time=10000)
        self.analysis = analysis
        self.paradigm_setup = paradigm_setup
        self.directory = directory
        if null_size is None:
            self.null_size = self.paradigm_setup.null_size
        else:
            self.null_size = null_size
    def run(self):
        os.chdir(self.directory)
        
        batches = getBatchCount(len(self.paradigm_setup.samples), batch_size = self.paradigm_setup.batch_size)
        for b in range(batches):
            os.system('cat outputFiles/%s_upstream_b%s_%s.fa >> paradigm/%s_upstream.fa' % (self.analysis.focus_node, b, batches, self.analysis.focus_node))
            os.system('cat outputFiles/%s_downstream_b%s_%s.fa >> paradigm/%s_downstream.fa' % (self.analysis.focus_node, b, batches, self.analysis.focus_node))
            for null in range(1, self.null_size + 1):
                os.system('cat outputFiles/N%s_%s_upstream_b%s_%s.fa >> paradigm/N%s_%s_upstream.fa' % (null, self.analysis.focus_node, b, batches, null, self.analysis.focus_node))
                os.system('cat outputFiles/N%s_%s_downstream_b%s_%s.fa >> paradigm/N%s_%s_downstream.fa' % (null, self.analysis.focus_node, b, batches, null, self.analysis.focus_node))
        if os.path.exists('outputFiles'):
            os.system('rm -rf outputFiles')

class computeShifts(Target):
    def __init__(self, fold, analysis, training_samples, paradigm_setup, current_parameters, parameters, directory, alpha = 0.05):
        Target.__init__(self, time=10000)
        self.fold = fold
        self.analysis = analysis
        self.training_samples = training_samples
        self.paradigm_setup = paradigm_setup
        if self.fold == 0:
            self.null_size = self.paradigm_setup.null_size
        else:
            self.null_size = 0
        self.current_parameters = current_parameters
        self.parameters = parameters
        self.directory = directory
        self.alpha = alpha
    def run(self):
        if self.fold == 0:
            ps_directory = '%s/analysis/%s/final/param_%s' % (self.directory, self.analysis.directory, '_'.join([str(parameter) for parameter in self.current_parameters]))
            os.chdir(ps_directory)
        else:
            ps_directory = "%s/analysis/%s/fold%s/param_%s" % (self.directory, self.analysis.directory, self.fold, '_'.join([str(parameter) for parameter in self.current_parameters]))
            os.chdir(ps_directory)
        logger('Computing P-Shifts ...\n', file = 'progress.log')
        
        ## define sample groups
        training_all = self.training_samples
        training_positive = list(set(self.analysis.positive_samples) & set(training_all))
        training_negative = list(set(self.analysis.negative_samples) & set(training_all))
        testing_all = list((set(self.paradigm_setup.samples) - set(self.training_samples)) & (set(self.analysis.positive_samples) | set(self.analysis.negative_samples)))
        testing_positive = list(set(self.analysis.positive_samples) & set(testing_all))
        testing_negative = list(set(self.analysis.negative_samples) & set(testing_all))
        
        ## read in Paradigm inferences
        assert(os.path.exists('paradigm/%s_upstream.fa' % (self.analysis.focus_node)))
        assert(os.path.exists('paradigm/%s_downstream.fa' % (self.analysis.focus_node)))
        upstream_ipls = readParadigm('paradigm/%s_upstream.fa' % (self.analysis.focus_node))[1]
        downstream_ipls = readParadigm('paradigm/%s_downstream.fa' % (self.analysis.focus_node))[1]
        upstream_ipls.to_csv('upstream_paradigm.tab', sep = '\t')
        downstream_ipls.to_csv('downstream_paradigm.tab', sep = '\t')
        # normalized_upstream_ipls = normalizeDataFrame(upstream_ipls.loc[[self.analysis.focus_node]], include_samples = training_negative)
        # normalized_downstream_ipls = normalizeDataFrame(downstream_ipls.loc[[self.analysis.focus_node]], include_samples = training_negative)
        null_upstream_ipls = {}
        null_downstream_ipls = {}
        # normalized_null_upstream_ipls = {}
        # normalized_null_downstream_ipls = {}
        for null in range(1, self.null_size + 1):
            assert(os.path.exists('paradigm/N%s_%s_upstream.fa' % (null, self.analysis.focus_node)))
            assert(os.path.exists('paradigm/N%s_%s_downstream.fa' % (null, self.analysis.focus_node)))
            null_upstream_ipls[null] = readParadigm('paradigm/N%s_%s_upstream.fa' % (null, self.analysis.focus_node))[1]
            null_downstream_ipls[null] = readParadigm('paradigm/N%s_%s_downstream.fa' % (null, self.analysis.focus_node))[1]
            # normalized_null_upstream_ipls[null] = normalizeDataFrame(null_upstream_ipls[null].loc[[self.analysis.focus_node]], include_samples = training_negative)
            # normalized_null_downstream_ipls[null] = normalizeDataFrame(null_downstream_ipls[null].loc[[self.analysis.focus_node]], include_samples = training_negative)
        
        ## compute raw and normalized p-shifts
        raw_shifts = {'real' : {}}
        for null in range(1, self.null_size + 1):
            raw_shifts['null%s' % (null)] = {}
        for sample in self.paradigm_setup.samples:
            assert((sample in downstream_ipls.columns) and (sample in upstream_ipls.columns))
            raw_shifts['real'][sample] = (downstream_ipls[sample][self.analysis.focus_node] - upstream_ipls[sample][self.analysis.focus_node])
            for null in range(1, self.null_size + 1):
                raw_shifts['null%s' % (null)][sample] = (null_downstream_ipls[null][sample][self.analysis.focus_node] - null_upstream_ipls[null][sample][self.analysis.focus_node])
        pandas.DataFrame(raw_shifts).to_csv('all_shifts.tab', sep = '\t', index_label = 'id')
        o = open('pshift.tab', 'w')
        o.write("> %s\tP-Shifts:Table\n" % (self.analysis.focus_node))
        o.write("# sample\tclass\tP-Shift\n")
        for sample in self.paradigm_setup.samples:
            if sample in training_positive + testing_positive:
                o.write("%s\t+\t%s\n" % (sample, raw_shifts['real'][sample]))
            elif sample in training_negative + testing_negative:
                o.write("%s\t-\t%s\n" % (sample, raw_shifts['real'][sample]))
            else:
                o.write("%s\t?\t%s\n" % (sample, raw_shifts['real'][sample]))
        o.close()
        o = open('wildtype_shifts.tab', 'w')
        o.write("sample\tP-Shift\n")
        for sample in self.paradigm_setup.samples:
            if sample in training_negative + testing_negative:
                o.write("%s\t%s\n" % (sample, raw_shifts['real'][sample]))
        o.close()
        #### output sample by real + null matrix
        # o = open('normalized_pshift.tab', 'w')
        # o.write("> %s\tNormalized_P-Shifts:Table\n" % (self.analysis.focus_node))
        # o.write("# sample\tclass\tP-Shift\n")
        # for sample in self.paradigm_setup.samples:
        #     if sample in training_positive + testing_positive:
        #         o.write("%s\t+\t%s\n" % (sample, normalized_shifts[sample]))
        #     else:
        #         o.write("%s\t-\t%s\n" % (sample, normalized_shifts[sample]))
        # o.close()
        raw_centered_shifts = {}
        # normalized_centered_shifts = {}
        (raw_negative_mean, raw_negative_sd) = computeMean([raw_shifts['real'][sample] for sample in training_negative], return_sd = True)
        (raw_positive_mean, raw_positive_sd) = computeMean([raw_shifts['real'][sample] for sample in training_positive], return_sd = True)
        # (normalized_negative_mean, normalized_negative_sd) = computeMean([normalized_shifts[sample] for sample in training_negative], return_sd = True)
        # (normalized_positive_mean, normalized_positive_sd) = computeMean([normalized_shifts[sample] for sample in training_positive], return_sd = True)
        for sample in self.paradigm_setup.samples:
            raw_centered_shifts[sample] = (raw_shifts['real'][sample] - raw_negative_mean)/raw_negative_sd
            # normalized_centered_shifts[sample] = (normalized_shifts[sample] - normalized_negative_mean)/normalized_negative_sd
        o = open('pshift.centered.tab', 'w')
        o.write("> %s\tP-Shifts:Table\n" % (self.analysis.focus_node))
        o.write("# sample\tclass\tP-Shift\n")
        for sample in self.paradigm_setup.samples:
            if sample in training_positive + testing_positive:
                o.write("%s\t+\t%s\n" % (sample, raw_centered_shifts[sample]))
            else:
                o.write("%s\t-\t%s\n" % (sample, raw_centered_shifts[sample]))
        o.close()
        # o = open('normalized_pshift.centered.tab', 'w')
        # o.write("> %s\tNormalized_P-Shifts:Table\n" % (self.analysis.focus_node))
        # o.write("# sample\tclass\tP-Shift\n")
        # for sample in self.paradigm_setup.samples:
        #     if sample in training_positive + testing_positive:
        #         o.write("%s\t+\t%s\n" % (sample, normalized_centered_shifts[sample]))
        #     else:
        #         o.write("%s\t-\t%s\n" % (sample, normalized_centered_shifts[sample]))
        # o.close()
        
        ## compute auc from stats
        absolute_shifts = {}
        classification_map = {}
        for sample in self.paradigm_setup.samples:
            absolute_shifts[sample] = abs(raw_centered_shifts[sample])
            if sample in training_positive + testing_positive:
                classification_map[sample] = 1
            else:
                classification_map[sample] = 0
        training_samples_sorted = deepcopy(training_all)
        if self.parameters.cross_validation_two_sided:
            training_samples_sorted.sort(lambda x, y: cmp(absolute_shifts[y], absolute_shifts[x]))
            training_auc = computeAUC(training_samples_sorted, absolute_shifts, classification_map)[0]
        else:
            if raw_positive_mean >= raw_negative_mean:
                training_samples_sorted.sort(lambda x, y: cmp(raw_shifts['real'][y], raw_shifts['real'][x]))
            elif raw_positive_mean < raw_negative_mean:
                training_samples_sorted.sort(lambda x, y: cmp(raw_shifts['real'][x], raw_shifts['real'][y]))
            training_auc = computeAUC(training_samples_sorted, raw_shifts['real'], classification_map)[0]
        if self.fold != 0:
            testing_samples_sorted = deepcopy(testing_all)
            if self.parameters.cross_validation_two_sided:
                testing_samples_sorted.sort(lambda x, y: cmp(absolute_shifts[y], absolute_shifts[x]))
                testing_auc = computeAUC(testing_samples_sorted, absolute_shifts, classification_map)[0]
            else:
                if raw_positive_mean >= raw_negative_mean:
                    testing_samples_sorted.sort(lambda x, y: cmp(raw_shifts['real'][y], raw_shifts['real'][x]))
                elif raw_positive_mean < raw_negative_mean:
                    testing_samples_sorted.sort(lambda x, y: cmp(raw_shifts['real'][x], raw_shifts['real'][y])) 
                testing_auc = computeAUC(testing_samples_sorted, raw_shifts['real'], classification_map)[0]
        else:
            testing_auc = '---'
        o = open('auc.tab', 'w')
        o.write('%s\t%s\n' % (training_auc, testing_auc))
        o.close()
        
        ## compute m-separation and significance
        if self.fold == 0:
            o = open('positive.scores', 'w')
            o.write('%s\n' % ('\n'.join([str(raw_shifts['real'][sample]) for sample in training_positive])))
            o.close()
            o = open('negative.scores', 'w')
            o.write('%s\n' % ('\n'.join([str(raw_shifts['real'][sample]) for sample in training_negative])))
            o.close()
            mseparation_map = {}
            mseparation_map['real'] = computeSeparation(training_positive,
                                                        training_negative,
                                                        raw_shifts['real'],
                                                        method = self.parameters.separation_method)
            for null in range(1, self.null_size + 1):
                mseparation_map['null%s' % (null)] = computeSeparation(training_positive,
                                                                       training_negative,
                                                                       raw_shifts['null%s' % (null)],
                                                                       method = self.parameters.separation_method)
            o = open('real.scores', 'w')
            o.write('%s\n' % (mseparation_map['real']))
            o.close()
            o = open('null.scores', 'w')
            o.write('%s\n' % ('\n'.join([str(mseparation_map['null%s' % (null)]) for null in range(1, self.null_size + 1)])))
            o.close()
            real_scores = [mseparation_map['real']]
            null_scores = [mseparation_map['null%s' % (null)] for null in range(1, self.null_size + 1)]
            (null_mean, null_sd) = computeMean(null_scores, return_sd = True)
            significance_score = (real_scores[0] - null_mean)/(null_sd + self.alpha)
            o = open('significance.tab', 'w')
            o.write('# Focus_Genes\tNegative_Samples\tPositive_Samples\tM-Separation\tZ-Score\n')
            o.write('%s\t%s\t%s\t%s\t%s\n' % (self.analysis.focus_node, len(training_negative), len(training_positive), real_scores[0], significance_score))
            o.close()
        else:
            o = open('significance.tab', 'w')
            o.write('# Focus_Genes\tNegative_Samples\tPositive_Samples\tM-Separation\tZ-Score\n')
            o.write('%s\t%s\t%s\t%s\t%s\n' % (self.analysis.focus_node, len(training_negative), len(training_positive), '---', '---'))
            o.close()

        ## output files for circleplots
        if self.fold == 0:
            o = open('up.features', 'w')
            for feature in set(upstream_ipls.index) - set(self.analysis.focus_genes):
                o.write('%s\n' % (feature))
            o.close()
            o = open('alteration.features', 'w')
            o.write('%s\n' % ('\n'.join(self.analysis.focus_genes)))
            if len(self.analysis.focus_genes) > 0:
                o.write('%s\n' % (self.analysis.focus_node))
            o.close()
            o = open('down.features', 'w')
            for feature in set(downstream_ipls.index) - set(self.analysis.focus_genes):
                o.write('%s\n' % (feature))
            o.close()
            o = open('include.samples', 'w')
            for sample in self.paradigm_setup.samples:
                o.write('%s\n' % (sample))
            o.close()
            o = open('alteration.circle', 'w')
            o.write('id')
            for sample in self.paradigm_setup.samples:
                o.write('\t%s' % (sample))
            o.write('\n*')
            for sample in self.paradigm_setup.samples:
                if sample in self.analysis.positive_samples:
                    o.write('\t1')
                elif sample in self.analysis.negative_samples:
                    o.write('\t0')
                else:
                    o.write('\t0.5')
            o.write('\n')
            o.close()
            os.system('median-center.py data/transposed_%s expression.circle' % (self.paradigm_setup.mrna[0].split('/')[-1]))
            o = open('shift.circle', 'w')
            o.write('id')
            for sample in self.paradigm_setup.samples:
                o.write('\t%s' % (sample))
            o.write('\n*')
            for sample in self.paradigm_setup.samples:
                o.write('\t%s' % (raw_shifts['real'][sample]))
            o.write('\n')
            o.close()
            o = open('color.map', 'w')
            o.write('> 1\n0\t255.255.255\n1\t0.0.0\n0.5\t150.150.150\n')
            o.close()

class compareParameters(Target):
    def __init__(self, fold, analysis, parameters, directory):
        Target.__init__(self, time=10000)
        self.fold = fold
        self.analysis = analysis
        self.parameters = parameters
        self.directory = directory
    def run(self):
        if self.fold == 0:
            ps_directory = '%s/analysis/%s/final' % (self.directory, self.analysis.directory)
            os.chdir(ps_directory)
            logger('Identifying the best parameters by training validation ...\n', file = 'progress.log')
        else:
            ps_directory = '%s/analysis/%s/fold%s' % (self.directory, self.analysis.directory, self.fold)
            os.chdir(ps_directory)
        
        ## report AUCs for model with best training validation
        best_training_auc = 0.0
        best_testing_auc = 0.0
        best_parameters = None
        for threshold in self.parameters.threshold:
            for cost in self.parameters.cost:
                for method in self.parameters.selection_method:
                    current_parameters = [threshold, cost, method]
                    f = open('param_%s/auc.tab' % ('_'.join([str(parameter) for parameter in current_parameters])), 'r')
                    (current_training_auc, current_testing_auc) = f.readline().rstrip().split('\t')
                    f.close()
                    try:
                        current_training_auc = float(current_training_auc)
                    except ValueError:
                        continue
                    try:
                        current_testing_auc = float(current_testing_auc)
                    except ValueError:
                        pass
                    if current_training_auc > best_training_auc:
                        best_training_auc = current_training_auc
                        best_testing_auc = current_testing_auc
                        best_parameters = [str(threshold), str(cost), str(method)]
        if self.fold != 0:
            o = open('auc.tab', 'w')
            if best_training_auc == 0.0:
                o.write('---\t---\t---\n')
            else:
                o.write('%s\t%s\t%s\n' % (best_training_auc, best_testing_auc, ','.join(best_parameters)))
            o.close()
        else:
            if best_parameters is not None:
                self.setFollowOnTarget(generateOutput(self.analysis, best_parameters, self.parameters, self.directory))

class generateOutput(Target):
    def __init__(self, analysis, best_parameters, parameters, directory):
        Target.__init__(self, time=10000)
        self.analysis = analysis
        self.best_parameters = best_parameters
        self.parameters = parameters
        self.directory = directory
    def run(self):
        ps_directory = '%s/analysis/%s' % (self.directory, self.analysis.directory)
        os.chdir(ps_directory)
        logger('Generating output files ...\n', file = 'progress.log')
        
        ## copy tables
        os.system('cp final/param_%s/significance.tab significance.tab' % ('_'.join(self.best_parameters)))
        os.system('cp final/param_%s/pshift.tab pshift.tab' % ('_'.join(self.best_parameters)))
        # os.system('cp final/param_%s/normalized_pshift.tab normalized_pshift.tab' % ('_'.join(self.best_parameters)))
        
        ## output m-separation and significance plots
        os.system('mseparation.R %s final/param_%s/positive.scores final/param_%s/negative.scores' % (self.analysis.focus_node,
                                                                                          '_'.join(self.best_parameters),
                                                                                          '_'.join(self.best_parameters)))
        os.system('significance.R %s final/param_%s/real.scores final/param_%s/null.scores' % (self.analysis.focus_node,
                                                                                          '_'.join(self.best_parameters),
                                                                                          '_'.join(self.best_parameters)))
        
        ## output circleplots
        os.mkdir('img')
        os.chdir('final/param_%s' % ('_'.join(self.best_parameters)))
        os.system('%s -m color.map -o \"%s;alteration.circle,shift.circle\" -s include.samples -f alteration.features ../../img/ alteration.circle expression.circle shift.circle' % (circleplot_executable, self.analysis.focus_node))
        os.system('%s -m color.map -o \"%s;alteration.circle,shift.circle\" -s include.samples -f up.features ../../img/ alteration.circle expression.circle' % (circleplot_executable, self.analysis.focus_node))
        os.system('%s -m color.map -o \"%s;alteration.circle,shift.circle\" -s include.samples -f down.features ../../img/ alteration.circle expression.circle' % (circleplot_executable, self.analysis.focus_node))
        os.chdir('../..')
        
        ## output sif
        combined_pathway = Pathway( ({}, {}) )
        combined_pathway.appendPathway(Pathway('final/param_%s/upstream_pathway.tab' % ('_'.join(self.best_parameters))))
        combined_pathway.appendPathway(Pathway('final/param_%s/downstream_pathway.tab' % ('_'.join(self.best_parameters))))
        combined_pathway.writeSIF('pshift_%s.sif' % (self.analysis.focus_node))
        
        ## cleanup unneeded files
        for round in range(1, self.parameters.n_rounds + 1):
            for fold in range(1, self.parameters.m_folds + 1):
                fold_index = (round - 1)*self.parameters.m_folds + fold
                if os.path.exists('fold%s' % (fold_index)):
                    shutil.rmtree('fold%s' % (fold_index))
        logger('Completed analysis ...\n', file = 'progress.log')

class makeReport(Target):
    def __init__(self, report_list, report_directory, directory):
        Target.__init__(self, time=10000)
        self.report_list = report_list
        self.report_directory = report_directory
        self.directory = directory
    def run(self):
        os.chdir(self.directory)
        logger('... Done\n', file = 'analysis/progress.log')
        
        if not os.path.exists(self.report_directory):
            os.mkdir('%s' % (self.report_directory))
        
        for analysis in self.report_list:
            if os.path.exists('%s/%s' % (self.report_directory, analysis)):
                continue
            os.mkdir('%s/%s' % (self.report_directory, analysis))
            os.system('cp analysis/%s/pshift.tab %s/%s' % (analysis, self.report_directory, analysis))
            os.system('cp analysis/%s/significance.tab %s/%s' % (analysis, self.report_directory, analysis))
            os.system('cp analysis/%s/*.sif %s/%s' % (analysis, self.report_directory, analysis))
            os.system('cp -r analysis/%s/img %s/%s' % (analysis, self.report_directory, analysis))
            os.system('cp analysis/%s/*.pdf %s/%s' % (analysis, self.report_directory, analysis))
        
        ## cytoscape-web
        # for gene in self.includeFeatures:
        #     if os.path.exists('analysis/%s/sig.tab' % (gene)):
        #         tableFiles = []
        #         tableFiles.append('analysis/%s/sig.tab' % (gene))
        #         tableFiles.append('msepPlot:analysis/%s/%s.msep.pdf' % (gene, gene))
        #         tableFiles.append('backgroundPlot:analysis/%s/%s.background.pdf' %
        #                                                                      (gene, gene))
        #         tableFiles.append('analysis/%s/avgAUC.tab' % (gene))
        #         tableFiles.append('analysis/%s/pshift.tab' % (gene))
        #         system('pathmark-report.py -t %s analysis/%s %s' %
        #                                      (','.join(tableFiles), gene, self.reportDir))
        #         system('cp analysis/%s/pshift* %s' % (gene, self.reportDir))

def ps_main():
    ## check for fresh run
    if os.path.exists('.jobTree'):
        logger('ERROR: .jobTree directory found, remove it first to start a fresh run\n', die = True)
    
    ## parse arguments
    parser = OptionParser(usage = '%prog [options] paradigm_directory analysis_file')
    Stack.addJobTreeOptions(parser)
    parser.add_option('--jobFile', help = 'Add as child of jobFile rather than new jobTree')
    parser.add_option('-c', '--config', dest = 'config_file', default = None, help = '')
    parser.add_option('-s', '--samples', dest = 'include_samples', default = None, help = '')
    parser.add_option('-f', '--features', dest = 'include_features', default = None, help = '')
    parser.add_option('-n', '--null_size', dest = 'null_size', default = 30, help = '')
    parser.add_option('-b', '--batch_size', dest = 'batch_size', default = 50, help = '')
    parser.add_option('-p', '--pathway', dest = 'pathway_file', default = None, help = '')
    parser.add_option('-y', '--public', action = 'store_true', dest = 'paradigm_public', default = False, help = '')
    parser.add_option('-z', '--seed', dest = 'seed', default = None, help = '')
    options, args = parser.parse_args()
    logger('Using Batch System : %s\n' % (options.batchSystem))
    
    if len(args) == 1:
        if args[0] == 'clean':
            command = 'rm -rf .jobTree analysis report'
            logger(command)
            os.system(command)
            sys.exit(0)
    
    assert(len(args) == 2)
    paradigm_directory = os.path.abspath(args[0])
    analysis_file = args[1]
    
    ## set Paradigm files
    paradigm_setup = ParadigmSetup(paradigm_directory.rstrip('/'),
                                   include_samples = options.include_samples,
                                   null_size = int(options.null_size),
                                   batch_size = int(options.batch_size),
                                   pathway_file = options.pathway_file,
                                   public = options.paradigm_public)
    
    ## set pathway
    global_pathway = Pathway(paradigm_setup.pathway)
    
    ## set parameters
    parameters = Parameters()
    if options.config_file:
        parameters.importConfig(options.config_file)
    parameters.setSeed(seed_input = options.seed)
    parameters.printSeed()
    
    ## read alterations
    if options.include_features:
        include_features = readList(options.include_features)
    else:
        include_features = None
    altered_list = []
    f = open(analysis_file, 'r')
    for line in f:
        if line.isspace():
            continue
        pline = line.rstrip().split('\t')
        if len(pline) == 3:
            altered = Alterations(pline[0],
                                  pline[1].split(','),
                                  paradigm_setup.samples,
                                  pline[2].split(','),
                                  negative_samples = None)
        elif len(pline) == 4:
            altered = Alterations(pline[0],
                                  pline[1].split(','),
                                  paradigm_setup.samples,
                                  pline[2].split(','),
                                  negative_samples = pline[3].split(','))
        if len(set(altered.focus_genes) & set(global_pathway.nodes)) == len(altered.focus_genes):
            if include_features != None:
                if len(set(altered.focus_genes) & set(include_features)) == 0:
                    continue
            if len(altered.positive_samples) >= parameters.min_alterations:
                altered_list.append(altered)
    f.close()
    
    ## run
    if parameters.in_parallel:
        s = Stack(branchAnalyses(altered_list,
                                 paradigm_setup,
                                 global_pathway,
                                 parameters,
                                 os.getcwd().rstrip('/')))
    else:
        s = Stack(queueAnalyses(altered_list,
                                paradigm_setup,
                                global_pathway,
                                parameters,
                                os.getcwd().rstrip('/')))
    if options.jobFile:
        s.addToJobFile(options.jobFile)
    else:
        if options.jobTree == None:
            options.jobTree = './.jobTree'
        
        failed = s.startJobTree(options)
        if failed:
            logger('%d jobs failed\n' % failed)
        else:
            os.system('rm -rf .lastjobTree')
            os.system('mv .jobTree .lastjobTree')

if __name__ == '__main__':
    from paradigmSHIFT import *
    ps_main()
