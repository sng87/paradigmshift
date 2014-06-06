#!/usr/bin/env python
"""
summarizeMAF.py

Author: Sam Ng
Last Updated: 2014-03-22
"""
import math, os, random, re, string, sys, types
from copy import deepcopy
from optparse import OptionParser

truncation_mutations = ["Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation", "Splice_Site"]
missense_mutations = ["Missense_Mutation"]
any_mutations = truncation_mutations + missense_mutations

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

def readMAF(input_file):
    """
    Reads a MAF format file [2014-3-1]
    Dependencies: logger
    """
    mutation_classification = {}
    f = open(input_file, 'r')
    line = f.readline().rstrip()
    if line.isspace():
        logger('ERROR: encountered a blank header\n', die = True)
    pline = line.split('\t')
    hugo_column = -1
    tumor_column = -1
    class_column = -1
    for index, column in enumerate(pline):
        if column == 'Hugo_Symbol':
            hugo_column = index
        elif column == 'Tumor_Sample_Barcode':
            tumor_column = index
        elif column == 'Variant_Classification':
            class_column = index
    assert(hugo_column != -1)
    assert(tumor_column != -1)
    assert(class_column != -1)
    samples = []
    for line in f:
        if line.isspace():
            continue
        pline = line.rstrip().split('\t')
        if pline[hugo_column] not in mutation_classification:
            mutation_classification[pline[hugo_column]] = {}
        if pline[tumor_column] not in mutation_classification[pline[hugo_column]]:
            mutation_classification[pline[hugo_column]][pline[tumor_column]] = []
        if pline[class_column] not in mutation_classification[pline[hugo_column]][pline[tumor_column]]:
            mutation_classification[pline[hugo_column]][pline[tumor_column]].append(pline[class_column])
        if pline[tumor_column] not in samples:
            samples.append(pline[tumor_column])
    f.close()
    return(samples, mutation_classification)

def formatSample(sample):
    """this function may need to be modified depending on the project"""
    return(sample[0:12])

def main():
    ## parse arguments
    parser = OptionParser(usage = '%prog [options] paradigm_directory analysis_file')
    options, args = parser.parse_args()
    
    assert(len(args) == 1)
    maf_file = args[0]
    
    ## read maf
    (maf_samples, maf_data) = readMAF(maf_file)
    
    ## write mutation events
    mutation_set_map = {}
    maf_features = maf_data.keys()
    maf_features.sort(lambda x, y: cmp(len(maf_data[y].keys()), len(maf_data[x].keys())))
    f = open('mutation_events.tab', 'w')
    for gene in maf_features:
        mutated_samples = []
        for sample in maf_data[gene]:
            if len(set(maf_data[gene][sample]) & set(any_mutations)) > 0:
                mutated_samples.append(formatSample(sample))
        if len(mutated_samples) == 0:
            continue
        f.write('%s_mutation\t%s\t%s\n' % (gene, gene, ','.join(mutated_samples)))
        mutation_set_map[gene] = deepcopy(mutated_samples)
    f.close()
    f = open('include.samples', 'w')
    f.write('%s\n' % ('\n'.join([formatSample(sample) for sample in maf_samples])))
    f.close()
    
    ## write mutation.list_t
    f = open('mutation.list_t', 'w')
    for gene in maf_features:
        if gene not in mutation_set_map:
            continue
        f.write('%s\t%s\n' % (gene, '\t'.join(mutation_set_map[gene])))
    f.close()

if __name__ == "__main__":
    main()
