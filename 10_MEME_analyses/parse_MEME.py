#!/usr/bin/env python

import argparse, re

parser = argparse.ArgumentParser()
parser.add_argument("-I", "--infile", help="meme.html", required = True)
args = parser.parse_args()

infile = args.infile
fileIn = open(infile, "r")

infileID = infile.split("/")

fileLine = fileIn.readline()
fileLine = fileLine.strip()

go = False
tom = False
tomList = []
memeLine = []
while fileLine:
    if "\"motifs\": [" in fileLine:
        go = True
    if go is True:
        if "\"consensus\"" in fileLine:
            cons = fileLine.split(":")
        if "\"evalue\"" in fileLine:
            pval = fileLine.split(":")
        if "\"tomtom_matches\"" in fileLine:
            tom = True
        while tom is True:
            if "\"id\":" in fileLine:
                fb = fileLine.split(":")
                tomList.append(fb[1].replace("\"", "").replace(",", ""))
            if "\"alt\":" in fileLine:
                alt = fileLine.split(":")
                tomList.append(alt[1].replace("\"", "").replace(",", ""))
                
                memeLine.append(infileID[1])
                memeLine.append(cons[1].replace("\"", "").replace(",", ""))
                memeLine.append(pval[1].replace("\"", "").replace(",", ""))
                memeLine = memeLine + tomList
                print(' '.join(memeLine))
                tomList = []
                memeLine = []
                
            if "\"prog\"" in fileLine:

                
                tom = False
                tomList = []
                memeLine = []
                
            fileLine = fileIn.readline()
            fileLine = fileLine.strip()
    if "\"groups\"" in fileLine:
        break
    fileLine = fileIn.readline()
    fileLine = fileLine.strip()