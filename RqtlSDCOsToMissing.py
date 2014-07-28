#!/usr/bin/python
# -*- coding: utf-8 -*-
#Set short double crossovers (tight double recombinations) that occur within a
#   designated distance (in cM) from R/qtl csvr formatted input file. Outputs
#   a file with short double crossovers set to missing with the suffix
#   "-SDCOsToMissingV02"
#
#Written by Ryan McCormick at Texas A&M University, July 2014,
#   as part of the the Truong and McCormick et al. (2014) manuscript.
#
#Provided as is, without warranty, without guarantee.
#
# *     This program is free software; you can redistribute it and/or
# *     modify it under the terms of the GNU General Public License,
# *     version 3, as published by the Free Software Foundation.
# *
# *     This program is distributed in the hope that it will be useful,
# *     but without any warranty; without even the implied warranty of
# *     merchantability or fitness for a particular purpose.  See the GNU
# *     General Public License, version 3, for more details.
# *
# *     A copy of the GNU General Public License, version 3, is available
# *     at http://www.r-project.org/Licenses/GPL-3
##############################################################################

import sys
import numpy
from collections import defaultdict
import scipy.stats as sp
import numpy as np

LILI_TABLE = []   #Global list of lists that will contain the genotypes and be modified

#####
#####
#### Converting short range DCOs to missing

if len(sys.argv) <= 1:
    print("RqtlSDCOsToMissing.py input.csv")
    sys.exit()
if sys.argv[1] == "--help" or sys.argv[1] == "-h":
    print("RqtlSDCOsToMissing.py input.csv")
    sys.exit()
try:
    LILI_TABLE = [line.strip() for line in open(sys.argv[1])]
    LILI_TABLE = [element.split(',') for element in LILI_TABLE]
except IOError:
    print("Cannot open target file. Please check your input:")
    print("\t$ python RqtlSDCOsToMissing.py input.csv")
    sys.exit()

bool_keepSDCO = False

if bool_keepSDCO == False:
    INDEX_CHR = 1
    INDEX_POS = 2
    INDEX_STARTSAMPLECOL = 3
    INDEX_STARTMARKERROW = 0
    FLOAT_DCOSIZE = 2.0
    currentChr = "NULL"
    for i in range(len(LILI_TABLE)):
        if LILI_TABLE[i][1] == "":
            INDEX_STARTMARKERROW = INDEX_STARTMARKERROW + 1
        else:
            break


    sys.stderr.write("\tConverting short DCOs to missing (DCO within %s cM)\n" % FLOAT_DCOSIZE)
    numSDCOgenotypes = 0
    for index_sample in range(INDEX_STARTSAMPLECOL, len(LILI_TABLE[INDEX_STARTMARKERROW])):
        #sys.stderr.write(LILI_TABLE[0][index_sample] + "\n")
        currentChr = "NULL"
        for index_marker in range(INDEX_STARTMARKERROW, len(LILI_TABLE)):
            if LILI_TABLE[index_marker][INDEX_CHR] != currentChr:
                #sys.stderr.write("Chr of current marker != current chr\n")
                currentChr = LILI_TABLE[index_marker][INDEX_CHR]
                continue

            #Get the previous non-missing genotype
            previousGenotype = LILI_TABLE[index_marker - 1][index_sample]
            i = 2  #Previous marker counter for missing genotypes; start at the next previous marker
            while previousGenotype == "-":
                if index_marker - i <= 0:   # Don't go out of range of the list
                    previousGenotype = "NULL"
                previousGenotype = LILI_TABLE[index_marker - i][index_sample]
                i = i + 1

            #Get the current non-missing genotype
            currentGenotype = LILI_TABLE[index_marker][index_sample]
            if currentGenotype == "-":
                continue

            #Compare the previous and current genotypes
            if previousGenotype != currentGenotype:   #Single crossover found
                #sys.stderr.write("No match; previous and current: " + previousGenotype + " " + currentGenotype + "\n")
                i = 1   #Forward marker counter
                markerDistance = 0.0
                bool_secondXO = False
                currentPos = float(LILI_TABLE[index_marker][INDEX_POS])
                while bool_secondXO == False and markerDistance <= FLOAT_DCOSIZE:   #While a second crossover isn't found and the markers are sufficiently close
                    if index_marker + i >= len(LILI_TABLE):
                        break
                    try:
                        nextGenotype = LILI_TABLE[index_marker + i][index_sample]
                    except IndexError:
                        sys.stdout.write("Index warning on sample: " + LILI_TABLE[0][index_sample] + " at row " + str(index_marker) + "\n")
                        sys.stdout.write("Length of table is: " + str(len(LILI_TABLE)) + " and i is: " + str(i) + "\n")
                    if LILI_TABLE[index_marker + i][INDEX_CHR] != currentChr:
                        #sys.stderr.write("Found new chromosome. Breaking\n")
                        break
                    while nextGenotype == "-":   #Get the next non-missing genotype
                        i = i + 1
                        try:
                            nextGenotype = LILI_TABLE[index_marker + i][index_sample]
                        except IndexError:   #The end of the genotypes have been reached; don't consider it a DCO.
                            sys.stdout.write("Index warning on sample: " + LILI_TABLE[0][index_sample] + " at row " + str(index_marker) + "\n")
                            sys.stdout.write("Length of table is: " + str(len(LILI_TABLE)) + " and i is: " + str(i) + "\n")
                            nextGenotype = currentGenotype
                    try:
                        nextPos = float(LILI_TABLE[index_marker + i][INDEX_POS])
                        markerDistance = nextPos - currentPos
                    except IndexError:
                        markerDistance = 0  #The end of the genotypes have been reached; set marker distance to 0
                    #sys.stderr.write("Current and Next Genotype: " + currentGenotype + " " + nextGenotype + "\n")
                    if currentGenotype != nextGenotype:
                        #sys.stderr.write("Found Double Crossover Within Size Range; setting current to missing\n")
                        LILI_TABLE[index_marker][index_sample] = "-"
                        numSDCOgenotypes = numSDCOgenotypes + 1
                        bool_secondXO = True
                    i = i + 1

str_outFile = sys.argv[1] + "-SDCOsToMissingV02"
file_outFile = open(str_outFile, 'w')

for i in range(len(LILI_TABLE)):
    file_outFile.write(",".join(LILI_TABLE[i]))
    file_outFile.write("\n")

file_outFile.close()
