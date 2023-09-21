import sys
import subprocess
import pysam
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import argparse

parser = argparse.ArgumentParser(description='Calculate complexity of scATAC-seq data')
parser.add_argument('-i', '--input_bam', help='Input bam file', required=True)
parser.add_argument('-c', '--cells', help='List of cells', required=True)
parser.add_argument('-t', '--threads', help='Number of threads', default=8)
parser.add_argument('-o', '--output', help='Output folder', default='complexity')
args = parser.parse_args()



subsample_fraction = [0.025,0.05,0.1,0.2, 0.3, 0.4, 0.5,0.6,0.7,0.8,0.9,1]
# subsample_fraction = [0.025]
threads = args.threads

input_bam_basename = os.path.basename(args.input_bam)
os.makedirs(args.output, exist_ok=True)
args.input_bam = os.path.abspath(args.input_bam)

bam_files       = [args.output + "/{}_subsampled_{}.bam".format(input_bam_basename.replace('.bam',''), fr) for fr in subsample_fraction]
fragments_files = [file.replace('.bam','_fragments.tsv') for file in bam_files]

cells_list = []
with open(args.cells,'r') as f:
    for line in f:
        line = line.rstrip()
        cells_list.append(line)


cmd_list = ["samtools view -@ {threads} --subsample {s} --subsample-seed 2002 -b -o {out} {input}".format(threads=threads, s=s, out=out, input=args.input_bam) for s,out in zip(subsample_fraction, bam_files)]
for cmd in cmd_list:
    sys.stderr.write('running: %s \n' % cmd)
    subprocess.run(cmd, shell=True)

# Index the bam files
index_cmd = ["samtools index {bam}".format(bam = file) for file in bam_files]
for i in index_cmd:
    sys.stderr.write('running: %s \n' % i)
    subprocess.run(i, shell=True)

sinto_cmd = ["sinto fragments -p {threads} -b {bam} -c {cells} -f {fragments}".format(threads=threads, bam=bam, cells=args.cells, fragments=fragments) for bam,fragments in zip(bam_files,fragments_files) ]

for i in sinto_cmd:
    sys.stderr.write('running: %s \n' % i)
    subprocess.run(i, shell=True)

def get_fragments(f):
    fragments = {bcd: 0 for bcd in cells_list}
    with open(f,'r') as f:
        for line in f:
            line = line.rstrip().split('\t')
            cell_barcode = line[3]
            fragments[cell_barcode] += 1
    return fragments

result = {}

# Summarize the results
for i,f in enumerate(subsample_fraction):
    nreads    = sum([int(l.split('\t')[2]) for l in pysam.idxstats(bam_files[i]).split('\n')[:-1]])
    fragments = get_fragments(fragments_files[i])
    result[f] = [nreads, np.median([int(i) for i in fragments.values()]), np.mean([int(i) for i in fragments.values()])]

df = pd.DataFrame(result).transpose()
df = df.rename(columns={0:'nreads',1:'median_fragments',2:'mean_fragments'})

plt.style.use('ggplot')
df.plot.scatter(x='nreads',y='median_fragments')

df.to_csv(args.output + '/complexity.tsv',sep='\t')
plt.savefig(args.output + '/complexity.png',format='png')

# Clenaup
for f in bam_files + fragments_files + [f + '.bai' for f in bam_files]:
    os.remove(f)