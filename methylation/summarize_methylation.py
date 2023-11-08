'''
summarize_methylation.py

python summarize_methylation.py -i INPUT.bedGraph

Purpose: Give summary of global methylation level given output from MethylDackel 

Arguments: 
	-i, --input: input .bedGraph file from MethylDackel 
	-t, --threshold: methylation level threshold, bases with methylation levels greater than this value will be marked as methylated, given as percentage (e.g, 1% should be 1)

Input file should be the output from MethyDackel extract command. The single cytosine methylation metric file or per CpG/CHG metric file can be used. 

Jessie Pelosi, 2023 
'''

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', dest = "input")
parser.add_argument('-t', '--threshold', dest = "thresh", default = 0)

args = parser.parse_args()

INFILE = args.input
THRESH = float(args.thresh) 
PREFIX = args.prefix 

import pandas as pd
bedGraph = pd.read_table(INFILE, header = None, names=["SeqName", "Start", "End", "MethylPercentage", "NumberMethylatedReads", "NumberUnmethylatedReads"])


total = bedGraph.shape[0]
methylated = bedGraph[bedGraph.MethylPercentage > THRESH].shape[0]

percent = (methylated/total) * 100 

print("Methylation level given input file:", INFILE, "is", percent,"%")
