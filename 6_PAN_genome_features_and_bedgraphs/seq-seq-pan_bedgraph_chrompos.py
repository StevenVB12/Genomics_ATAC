#!/usr/bin/env python

# Use as: se-seq-pan_bedgraph_chrompos.py -I start_end/H_e_dem_peaks.bed -g 1 -p genome_pos_H_e_dem.txt -o start_end/H_e_dem_peaks

# genome_pos_H_e_dem.txt
# Scaf		Length	CummStart
# Herato0101	22325789	0
# Herato0201	28654	22325789
# Herato0202	186807	22354443
# …		…	…

from __future__ import division
import argparse, re

parser = argparse.ArgumentParser()
parser.add_argument("-I", "--infile", help="bed file", required = True)
parser.add_argument("-g", "--genomeN", help="genome number in pan genome", required = True)
parser.add_argument("-p", "--scafgenomepos", help="table with cumulative genome positions of scaffolds", required = True)
parser.add_argument("-o", "--outfile", help="Output txt file name", action = "store")
args = parser.parse_args()

infile = args.infile
outfile = args.outfile
genomeNumber = args.genomeN
scafpos = args.scafgenomepos

s= open(outfile + '_start.txt',"a+")
s.write("\t".join((genomeNumber,"c","\n")))

e= open(outfile + '_end.txt',"a+")
e.write("\t".join((genomeNumber,"c","\n")))

posIn = open(scafpos, "r")
posLine = posIn.readline()

posDic = {}
while posLine:
    
    lineElem =  posLine.split('\t')
    posDic[lineElem[0]] = lineElem[2]
    posLine = posIn.readline()

fileIn = open(infile, "r")
fileLine = fileIn.readline()

n = 1
while fileLine:
    lineElem =  fileLine.split('\t')
    
    if(n ==1):
        start = int(lineElem[1]) + int(posDic[lineElem[0]]) + 1
    else:
        start = int(lineElem[1]) + int(posDic[lineElem[0]])
        
    end = int(lineElem[2]) + int(posDic[lineElem[0]])
    
    s.write("".join((str(start),"\n")))
    e.write("".join((str(end),"\n")))
    
    fileLine = fileIn.readline()
    
    n+=1
