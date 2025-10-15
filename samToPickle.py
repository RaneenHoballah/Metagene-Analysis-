"""
samToPickle.py

Converts a SAM alignment file into a pickle file containing base-wise coverage
for all genes in the genome, including flanking regions.

Inputs:
- SAM file
- GTF file of genes
- Output pickle file path
- Flanking size (bp)

Outputs:
- Pickle file containing gene coverages
- Optional intermediate debug pickle files (genesByChromosome.pickle, coverageDict_debug.pickle)

Dependencies:
- Python 3
- numpy, csv, pickle
"""

import csv
import numpy as np
import sys
import pickle
from collections import defaultdict
import os

# ----------------------------
# Arguments & Settings
# ----------------------------
if len(sys.argv) != 5:
    print("Usage: python3 samToPickle.py <input.sam> <genes.gtf> <output.pickle> <flankingSize>")
    sys.exit(1)

samFilePath = sys.argv[1]
gtfFilePath = sys.argv[2]
pickleFilePath = sys.argv[3]
flankingSize = int(sys.argv[4])

# Optional debug flag
debug = False

# Check input files exist
if not os.path.exists(samFilePath):
    raise FileNotFoundError(f"SAM file not found: {samFilePath}")
if not os.path.exists(gtfFilePath):
    raise FileNotFoundError(f"GTF file not found: {gtfFilePath}")

print(f"Processing SAM file: {samFilePath}")

# ----------------------------
# Read GTF File
# ----------------------------
gtfContents = list(csv.reader(open(gtfFilePath), delimiter='\t'))

# ----------------------------
# Read SAM File
# ----------------------------
samContents = list(csv.reader(open(samFilePath), delimiter='\t'))
samHeaders = samContents[:8]  # standard SAM header lines
samContents = samContents[8:]  # actual alignments

print(f"Total SAM entries: {len(samContents)}")
normalisingValue = len(samContents)

# ----------------------------
# Organize Genes by Chromosome
# ----------------------------
chromosomes = list(set([line[0] for line in gtfContents]))
genesByChromosome = {chromosome: defaultdict(tuple) for chromosome in chromosomes}

for line in gtfContents:
    if line[2] == 'transcript':
        chromosome = line[0]
        geneID = [entry for entry in line[8].split(';') if 'gene_id' in entry][0].split('"')[1]
        strand = line[6]
        start = int(line[3])
        end = int(line[4])
        genesByChromosome[chromosome][geneID] = (strand, start, end)

if debug:
    pickle.dump(genesByChromosome, open("genesByChromosome.pickle", "wb"))

# ----------------------------
# Determine chromosome lengths for coverage
# ----------------------------
maxChromosomeBasePositions = {
    chromosome: np.max([genesByChromosome[chromosome][geneID][2] for geneID in genesByChromosome[chromosome].keys()]) + 20000
    for chromosome in chromosomes
}

# ----------------------------
# Initialize coverage dictionary
# ----------------------------
coverageDict = {
    chromosome: {base: 0 for base in range(-20000, maxChromosomeBasePositions[chromosome])}
    for chromosome in chromosomes
}

# ----------------------------
# Populate coverage dictionary
# ----------------------------
for alignment in samContents:
    chromosome = alignment[2]
    start = int(alignment[3]) - 1
    size = len(alignment[9])
    end = start + size
    for position in range(start, end):
        coverageDict[chromosome][position] += 1 / normalisingValue

if debug:
    with open("coverageDict_debug.pickle", "wb") as f:
        pickle.dump(coverageDict, f)

# ----------------------------
# Compute gene coverages including flanking regions
# ----------------------------
geneCoverages = {
    chromosome: {geneID: None for geneID in genesByChromosome[chromosome].keys()}
    for chromosome in chromosomes
}

for chromosome in chromosomes:
    for geneID in geneCoverages[chromosome].keys():
        strand, start, end = genesByChromosome[chromosome][geneID]
        coordExtremes = [start - flankingSize, end + 1 + flankingSize]
        if strand == '-':
            coordExtremes.reverse()
            coordRange = list(range(coordExtremes[0], coordExtremes[1], -1))
        else:
            coordRange = list(range(coordExtremes[0], coordExtremes[1]))
        geneCoverages[chromosome][geneID] = [coverageDict[chromosome][base] for base in coordRange]

# ----------------------------
# Save gene coverages to pickle
# ----------------------------
pickle.dump(geneCoverages, open(pickleFilePath, 'wb'))
print(f"Pickle file created: {pickleFilePath}")
