#!/usr/bin/env python
"""
summarizeMAF.py
    by Sam Ng
"""
import math, os, random, re, string, sys, types
from copy import deepcopy

from optparse import OptionParser

#### To Do:
####     - Rewrite this to handle multiple MAFs

## default variables
truncation_calls = ["Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation", "Splice_Site"]
missense_calls = ["Missense_Mutation"]
noncoding_calls = ["RNA"]
include_calls = truncation_calls + missense_calls
exclude_calls = ["Silent"]
assert(len(set(include_calls) & set(exclude_calls)) == 0)

def readMAF(input_file, gene_column_name = "Hugo_Symbol", sample_column_name = "Tumor_Sample_Barcode", call_column_name = "Variant_Classification"):
    """
    Reads a MAF format file [2014-6-24]
    """
    mutation_classification = {}
    f = open(input_file, "r")
    line = f.readline().rstrip()
    while line.startswith("#"):
        line = f.readline().rstrip()
    assert(not line.isspace())
    pline = line.split("\t")
    hugo_column = -1
    tumor_column = -1
    class_column = -1
    for index, column in enumerate(pline):
        if column == gene_column_name:
            hugo_column = index
        elif column == sample_column_name:
            tumor_column = index
        elif column == call_column_name:
            class_column = index
    assert(hugo_column != -1)
    assert(tumor_column != -1)
    assert(class_column != -1)
    samples = []
    for line in f:
        if line.isspace():
            continue
        pline = line.rstrip().split("\t")
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

def formatSample(sample, length = 12):
    """this function may need to be modified depending on the project [SummarizeMAF Specific]"""
    return(sample[0:length])

def main():
    ## parse arguments
    parser = OptionParser(usage = "%prog [options] maf_file")
    parser.add_option("-i", "--include", dest="include_calls", default=None)
    parser.add_option("-e", "--exclude", dest="exclude_calls", default=None)
    parser.add_option("-l", "--length", dest="sample_length", default="12")
    parser.add_option("-o", "--output", dest="output_name", default="mutation")
    options, args = parser.parse_args()
    
    assert(len(args) == 1)
    maf_file = args[0]
    
    ## set options
    global include_calls, exclude_calls
    if options.include_calls:
        include_calls = options.include_calls.split(",")
    if options.exclude_calls:
        exclude_calls = options.exclude_calls.split(",")
    
    ## read maf
    (maf_samples, maf_data) = readMAF(maf_file)
    
    ## write samples list
    f = open("include.samples", "w")
    f.write("%s\n" % ("\n".join([formatSample(sample, length = int(options.sample_length)) for sample in maf_samples])))
    f.close()
    
    ## write tab
    mutation_set_map = {}
    maf_features = maf_data.keys()
    maf_features.sort(lambda x, y: cmp(len(maf_data[y].keys()), len(maf_data[x].keys())))
    f = open("%s_events.tab" % (options.output_name), "w")
    for gene in maf_features:
        mutated_samples = []
        for sample in maf_data[gene]:
            if options.include_calls:
                if len(set(maf_data[gene][sample]) & set(include_calls)) > 0:
                     mutated_samples.append(formatSample(sample, length = int(options.sample_length)))
            elif options.exclude_calls:
                if len(set(maf_data[gene][sample]) - set(exclude_calls)) > 0:
                     mutated_samples.append(formatSample(sample, length = int(options.sample_length)))
            else:
                if len(set(maf_data[gene][sample]) & set(include_calls)) > 0:
                     mutated_samples.append(formatSample(sample, length = int(options.sample_length)))
        if len(mutated_samples) > 0:
            f.write("%s_mutation\t%s\t%s\n" % (gene, gene, ",".join(mutated_samples)))
            mutation_set_map[gene] = deepcopy(mutated_samples)
    f.close()
    
    ## write list_t
    f = open("%s_events.list_t" % (options.output_name), "w")
    for gene in maf_features:
        if gene not in mutation_set_map:
            continue
        f.write("%s\t%s\n" % (gene, "\t".join(mutation_set_map[gene])))
    f.close()

if __name__ == "__main__":
    main()
