"""
metagene_plot.py

Generates metagene plots for protein-coding and noncoding genes
from ChIP-seq coverage data stored in pickle files.

Inputs:
- GTF file of genes
- Pickle files for IP and Input
- Flanking size around TSS/TES
- Output plot file path
- Optional gene set file (True/False and path)
- Optional protein-coding and noncoding gene lists

Outputs:
- PDF plot showing normalized coverage depth upstream, across gene, and downstream

Dependencies:
- Python 3
- numpy, matplotlib, csv
"""

import csv
import numpy as np
import math
from collections import defaultdict
import matplotlib.pyplot as plt
import sys
import pickle
import os

# ----------------------------
# Input arguments
# ----------------------------
gtfFilePath = sys.argv[1]
pickleFilePathIP = sys.argv[2]
pickleFilePathInput = sys.argv[3]
flankingSize = int(sys.argv[4])
plotFilePath = sys.argv[5]
inputGeneSet = sys.argv[6].lower() == 'true'

# Optional gene set file
if inputGeneSet:
    geneSetFilePath = sys.argv[7]
    if not os.path.exists(geneSetFilePath):
        raise FileNotFoundError(f"Gene set file not found: {geneSetFilePath}")
    geneSet = [row[0] for row in csv.reader(open(geneSetFilePath), delimiter='\t')]
else:
    geneSet = None

# Protein-coding and noncoding gene lists
pcoding_gene_file = "Pcoding-Q4.txt"
noncoding_gene_file = "nonCode-Q4.txt"

for f in [pcoding_gene_file, noncoding_gene_file]:
    if not os.path.exists(f):
        raise FileNotFoundError(f"Gene list file not found: {f}")

pcoding_genes = [line.strip() for line in open(pcoding_gene_file)]
noncoding_genes = [line.strip() for line in open(noncoding_gene_file)]

# ----------------------------
# Read GTF file
# ----------------------------
if not os.path.exists(gtfFilePath):
    raise FileNotFoundError(f"GTF file not found: {gtfFilePath}")

gtfContents = list(csv.reader(open(gtfFilePath), delimiter='\t'))
chromosomes = list(set([line[0] for line in gtfContents]))
genesByChromosome = {chrom: defaultdict(tuple) for chrom in chromosomes}

for line in gtfContents:
    if line[2] == 'transcript':
        chrom = line[0]
        geneID = [entry for entry in line[8].split(';') if 'gene_id' in entry][0].split('"')[1]
        strand = line[6]
        start = int(line[3])
        end = int(line[4])
        genesByChromosome[chrom][geneID] = (strand, start, end)

# ----------------------------
# Functions
# ----------------------------
def get_metagene_values(pickleFilePath):
    geneCoverages = pickle.load(open(pickleFilePath, 'rb'))
    geneOrder, geneVals, upVals, downVals = [], [], [], []

    for chrom in chromosomes:
        for geneID in genesByChromosome[chrom]:
            geneLength = len(geneCoverages[chrom][geneID]) - 2*flankingSize
            interval = geneLength / 100
            intervalBreaks = [interval * i for i in range(101)]
            intervalStartEndVals = [(math.ceil(intervalBreaks[i]) + flankingSize + 1,
                                     math.floor(intervalBreaks[i+1]) + 1 + flankingSize + 1)
                                    for i in range(100)]
            intervalCoverages = [np.mean([geneCoverages[chrom][geneID][pos] for pos in range(start, end)])
                                 for start, end in intervalStartEndVals]
            upstreamCoverage = geneCoverages[chrom][geneID][:flankingSize+1]
            downstreamCoverage = geneCoverages[chrom][geneID][geneLength + flankingSize -1 :]
            if not any(math.isnan(v) for v in intervalCoverages):
                geneVals.append(intervalCoverages)
                upVals.append(upstreamCoverage)
                downVals.append(downstreamCoverage)
                geneOrder.append(geneID)

    return geneOrder, {'geneVals': geneVals, 'upVals': upVals, 'downVals': downVals}


def average_metagene_values(geneOrder, metageneValues, geneSet):
    geneIndices = [i for i, g in enumerate(geneOrder) if g in geneSet]
    geneValsArray = np.array([metageneValues['geneVals'][i] for i in geneIndices])
    upValsArray = np.array([metageneValues['upVals'][i] for i in geneIndices])
    downValsArray = np.array([metageneValues['downVals'][i] for i in geneIndices])
    return np.mean(geneValsArray, axis=0), np.mean(upValsArray, axis=0), np.mean(downValsArray, axis=0)


# ----------------------------
# Run metagene analysis
# ----------------------------
geneOrderIP, metageneIP = get_metagene_values(pickleFilePathIP)
geneOrderInput, metageneInput = get_metagene_values(pickleFilePathInput)

pcoding_ip_vals, pcoding_ip_up, pcoding_ip_down = average_metagene_values(geneOrderIP, metageneIP, pcoding_genes)
pcoding_input_vals, pcoding_input_up, pcoding_input_down = average_metagene_values(geneOrderInput, metageneInput, pcoding_genes)

noncoding_ip_vals, noncoding_ip_up, noncoding_ip_down = average_metagene_values(geneOrderIP, metageneIP, noncoding_genes)
noncoding_input_vals, noncoding_input_up, noncoding_input_down = average_metagene_values(geneOrderInput, metageneInput, noncoding_genes)

# Normalize
pcoding_gene_vals = pcoding_ip_vals / pcoding_input_vals
pcoding_up_vals = pcoding_ip_up / pcoding_input_up
pcoding_down_vals = pcoding_ip_down / pcoding_input_down

noncoding_gene_vals = noncoding_ip_vals / noncoding_input_vals
noncoding_up_vals = noncoding_ip_up / noncoding_input_up
noncoding_down_vals = noncoding_ip_down / noncoding_input_down

# ----------------------------
# Plotting
# ----------------------------
fig, axes = plt.subplots(1,3,figsize=(12,6), gridspec_kw={'width_ratios':[1,2,1]}, sharey=True)
fig.subplots_adjust(wspace=0, hspace=0)

percentages = [i + 0.5 for i in range(100)]
percentages[0], percentages[-1] = 0, 100

# Upstream
axes[0].plot(range(-flankingSize,1), pcoding_up_vals, label='Pcoding')
axes[0].plot(range(-flankingSize,1), noncoding_up_vals, label='Noncoding', linestyle='dashed')
axes[0].set_xlim(-flankingSize,0)

# Gene body
axes[1].plot(percentages, pcoding_gene_vals, label='Pcoding')
axes[1].plot(percentages, noncoding_gene_vals, label='Noncoding', linestyle='dashed')
axes[1].set_xlim(0,100)
axes[1].xaxis.tick_top()

# Downstream
axes[2].plot(range(0,flankingSize+1), pcoding_down_vals, label='Pcoding')
axes[2].plot(range(0,flankingSize+1), noncoding_down_vals, label='Noncoding', linestyle='dashed')
axes[2].set_xlim(0,flankingSize)

# Labels
fig.text(0.25,0.04,'Distance from TSS (bp)', ha='center')
fig.text(0.5,0.08,'Position through gene (%)', ha='center')
fig.text(0.75,0.04,'Distance from TES (bp)', ha='center')
fig.text(0.06,0.5,'Normalized coverage depth', va='center', rotation='vertical')

axes[0].legend(loc='upper right')

plt.savefig(plotFilePath)
