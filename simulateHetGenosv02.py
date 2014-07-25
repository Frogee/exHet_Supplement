#!/usr/bin/python
# -*- coding: utf-8 -*-
#Simulate genotypes for a recombinant inbred population propagated by selfing.
#   The generation interval, proportion heterozygosity maintained per generation,
#   number of individuals, and map size (in cM) can be altered. Currently uses
#   the Haldane mapping function.
#
#Written by Ryan McCormick at Texas A&M University, July 2014,
#   for the Truong and McCormick et al. (2014) manuscript.
#   Some of the code (and variable names) are ports of the C code that
#   we added on top of a forked R/qtl code base, found here:
#   https://github.com/Frogee/qtl/blob/hetexp/src/hmm_bcsft.c
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
import numpy as np
from math import pow, sqrt, log, exp
import random

NUMIND = 1000
NUMMARK = 1000
GEN = 7
HET = 0.6373027617
MAPSIZE = 200
ERRORRATE = 0.01
MISSINGRATE = 0.05

def generateCountmat(rf, generation, het):

    transpr = [None] * 10
    countmat = [None] * 9
    r = rf
    t = generation
    B_11 = B_12 = B_14 = B_22 = B_23 = 0.0
    u = d = 0.0
    h = het
    hpowt = h2 = 0.0
    r2 = r3 = r4 = r5 = u2 = u3 = u4 = d2 = 0.0
    i = 0

    if ((r > 0.4999999) and (r < 0.5000001)):
        r = 0.4999999
    if (r < 0.0000001):
        r = 0.0000001

    hpowt = pow(h, t)
    h2 = pow(h, 2.0)
    r2 = pow(r, 2.0)
    r3 = pow(r, 3.0)
    r4 = pow(r, 4.0)
    r5 = pow(r, 5.0)

    u =  -(2.0*h*r - r + sqrt((r2 - 2.0*h*r + h)*(2.0*h*r - 2.0*r - h + r2 + 1.0)) - 2.0*h*r2 + r2)/(h + 2.0*r - 2.0*h*r + 2.0*h*r2 - 2.0*r2 - 1.0)
    u2 = pow(u, 2.0)
    u3 = pow(u, 3.0)
    u4 = pow(u, 4.0)

    d = 2.0*(pow((1.0-r),2))+8.0*u*r*(1-r)+ 2.0*(r2)+ 2.0*(u2)*((pow((1.0-r),2.0))+(r2))
    d2 = pow(d, 2.0)

    B_11 = (1.0/2.0)*((2.0*(- 8.0*r3*u3 + 8.0*r3*u2 + 12.0*r2*u3 - 12.0*r2*u2 - 2.0*d*r2*u + d*r2 - 4.0*r*u3 + 8.0*r*u2 + 2.0*d*r*u - 2.0*d*r - 2.0*u2 + d))/((d + 4.0*r*u2 - 2.0*u2)*(- 4.0*r2*u2 + 4.0*r*u2 - 2.0*u2 + d)) - (d*pow((-(4.0*r*u2 - 2.0*u2)/d),t))/(2.0*(d*u2 + 4.0*r*u4 - 2.0*u4)) + (4.0*hpowt*(- u*r2 + u*r))/(h*(4.0*r2*u2 - 4.0*r*u2 + 2.0*u2 - d*h)) + (d*pow(((4.0*r2*u2 - 4.0*r*u2 + 2.0*u2)/d),t)*(16.0*r2*u2 - 16.0*r3*u2 + 8.0*r4*u2 - d*h - 8.0*r*u2 + 2.0*u2 - 4.0*d*r2*u + 2.0*d*h*r + 4.0*d*r*u - 2.0*d*h*r2 - 4.0*d*h*r*u + 4.0*d*h*r2*u))/(2.0*(2.0*r2*u2 - 2.0*r*u2 + u2)*(32.0*r2*u4 - 32.0*r3*u4 + 16.0*r4*u4 + d2*h - 2.0*d*u2 - 16.0*r*u4 + 4.0*u4 + 4.0*d*r*u2 - 4.0*d*r2*u2 - 2.0*d*h*u2 + 4.0*d*h*r*u2 - 4.0*d*h*r2*u2)))

    B_14 = (1.0/2.0)*((d*pow((-(4.0*r*u2 - 2.0*u2)/d),t))/(2.0*(d*u2 + 4.0*r*u4 - 2.0*u4)) + (2.0*(- 8.0*r3*u3 + 12.0*r2*u3 - 2.0*d*r2*u + d*r2 - 4.0*r*u3 + 2.0*d*r*u))/((d + 4.0*r*u2 - 2.0*u2)*(- 4.0*r2*u2 + 4.0*r*u2 - 2.0*u2 + d)) + (4.0*hpowt*(- u*r2 + u*r))/(h*(4.0*r2*u2 - 4.0*r*u2 + 2.0*u2 - d*h)) + (d*pow(((4.0*r2*u2 - 4.0*r*u2 + 2.0*u2)/d),t)*(16.0*r2*u2 - 16.0*r3*u2 + 8.0*r4*u2 - d*h - 8.0*r*u2 + 2.0*u2 - 4.0*d*r2*u + 2.0*d*h*r + 4.0*d*r*u - 2.0*d*h*r2 - 4.0*d*h*r*u + 4.0*d*h*r2*u))/(2.0*(2.0*r2*u2 - 2.0*r*u2 + u2)*(32.0*r2*u4 - 32.0*r3*u4 + 16.0*r4*u4 + d2*h - 2.0*d*u2 - 16.0*r*u4 + 4.0*u4 + 4.0*d*r*u2 - 4.0*d*r2*u2 - 2.0*d*h*u2 + 4.0*d*h*r*u2 - 4.0*d*h*r2*u2)))

    B_12 = (1.0/4.0)*((4.0*d*(- r2 + r)*pow(((4.0*r2*u2 - 4.0*r*u2 + 2.0*u2)/d),t))/((2.0*u*r2 - 2.0*u*r + u)*(4.0*r2*u2 - 4.0*r*u2 + 2.0*u2 - d*h)) - (8.0*hpowt*(- u*r2 + u*r))/(h*(4.0*r2*u2 - 4.0*r*u2 + 2.0*u2 - d*h)))

    B_22 = (d*pow(((4.0*r2*u2 - 4.0*r*u2 + 2.0*u2)/d),t))/(4.0*(2.0*r2*u2 - 2.0*r*u2 + u2)) - (d*pow((-(4.0*r*u2 - 2.0*u2)/d),t))/(4.0*(2.0*r*u2 - u2))

    B_23 = (d*pow(((4.0*r2*u2 - 4.0*r*u2 + 2.0*u2)/d),t))/(4.0*(2.0*r2*u2 - 2.0*r*u2 + u2)) + (d*pow((-(4.0*r*u2 - 2.0*u2)/d),t))/(4.0*(2.0*r*u2 - u2))

    transpr[0] = B_11   #   /* AABB */
    transpr[1] = B_12   #   /* AABb aaBb*/
    transpr[2] = B_14   #   /* AAbb */
    transpr[3] = B_22   #   /* AaBb */
    transpr[4] = B_23   #   /* AabB */
    transpr[5] = B_11   #   /* aabb */
    transpr[6] = B_12   #   /* AaBB Aabb*/

    #//sprintf(text, "%s\t%f\t%f\t%f\n", "Marginal probabilities 7, 8, 9: ", transpr[7], transpr[8], transpr[9]);
    #//Rprintf(text);
    #/* marginal probabilities for a single marker from the joint probability function*/
    transpr[7] = transpr[0] + transpr[1] + transpr[2]   #    /* AA */
    transpr[7] = log(transpr[7])   #
    transpr[9] = transpr[7]   #				/* aa */
    transpr[8] = transpr[1] + transpr[3] + transpr[4] + transpr[1]
    transpr[8] = log(transpr[8])

    countmat[0] = transpr[0] * 100.0   # aa to aa
    countmat[1] = transpr[1] * 100.0 * 2.0   # aa to ab or bb to ab
    countmat[2] = transpr[3] * 100.0 + transpr[4] * 100.0   # ab to ab
    countmat[3] = transpr[2] * 100.0 * 2.0   # aa to bb or bb to aa
    countmat[4] = transpr[6] * 100.0 * 2.0   # ab to aa or ab to bb
    countmat[5] = transpr[5] * 100.0   # bb to bb

    countmat[0] = transpr[0] * 1.0   # aa to aa
    countmat[1] = transpr[1] * 1.0   # aa to ab
    countmat[2] = transpr[1] * 1.0   # bb to ab
    countmat[3] = transpr[3] * 1.0 + transpr[4] * 1.0   # ab to ab
    countmat[4] = transpr[2] * 1.0   # aa to bb
    countmat[5] = transpr[2] * 1.0   # bb to aa
    countmat[6] = transpr[6] * 1.0   # ab to aa
    countmat[7] = transpr[6] * 1.0   # ab to bb
    countmat[8] = transpr[5] * 1.0   # bb to bb

    return countmat

def drawSample(countmat):
    #Inspired by http://stackoverflow.com/questions/4437250/choose-list-variable-given-probability-of-each-variable
    #print("Within drawSample()")
    r = random.random()
    index = 0
    while(r >= 0 and index < len(countmat)):
        r -= countmat[index]
        index += 1
    return index - 1

def countmatIndexToGenoPair(index):
    #print("Within countmatIndexToGenoPair()")
    if index == 0:
        return ["aa", "aa"]
    elif index == 1:
        return ["aa", "ab"]
    elif index == 2:
        return ["bb", "ab"]
    elif index == 3:
        return ["ab", "ab"]
    elif index == 4:
        return ["aa", "bb"]
    elif index == 5:
        return ["bb", "aa"]
    elif index == 6:
        return ["ab", "aa"]
    elif index == 7:
        return ["ab", "bb"]
    elif index == 8:
        return ["bb", "bb"]

def rfIndexToRF(index):
    #print("Within countmatIndexToGenoPair()")
    if index == 0:
        return random.uniform(0.00, 0.0025)
    elif index == 1:
        return random.uniform(0.0025, 0.005)
    elif index == 2:
        return random.uniform(0.005, 0.10)
    elif index == 3:
        return random.uniform(0.00, 0.0001)

def initializeGenotypes(rf, gen, het):
    countmat = generateCountmat(rf, gen, het)
    geneticDistance = -50.0 * log(1.0 - 2.0 * rf)   #Haldane mapping function
    genoList = []
    genoList.append(["ID", "", "", ""])
    genoList.append(["Marker1", "1", "0", "0"])
    genoList.append(["Marker2", "1", str(rf), str(geneticDistance) ])
    #print("Within initializeGenotypes()")
    for i in range(NUMIND):
        #print("On individual " + str(i))
        #print(genoList)
        countmatIndex = drawSample(countmat)
        genoPair = countmatIndexToGenoPair(countmatIndex)
        genoList[0].append(str(i + 1))
        genoList[1].append(genoPair[0])
        genoList[2].append(genoPair[1])
    return genoList

def returnNextGenotype(countmat, prevGeno):
    #countmat[0] = transpr[0] * 1.0   # aa to aa
    #countmat[1] = transpr[1] * 1.0   # aa to ab
    #countmat[2] = transpr[1] * 1.0   # bb to ab
    #countmat[3] = transpr[3] * 1.0 + transpr[4] * 100.0   # ab to ab
    #countmat[4] = transpr[2] * 1.0   # aa to bb
    #countmat[5] = transpr[2] * 1.0   # bb to aa
    #countmat[6] = transpr[6] * 1.0   # ab to aa
    #countmat[7] = transpr[6] * 1.0   # ab to bb
    #countmat[8] = transpr[5] * 1.0   # bb to bb
    genoModCountmat = [None] * 9
    nextGeno = "-"
    if prevGeno == "aa":
        sumClasses = countmat[0] + countmat[1] + countmat[4]
        genoModCountmat[0] = countmat[0] / sumClasses
        genoModCountmat[1] = countmat[1] / sumClasses
        genoModCountmat[2] = 0.0
        genoModCountmat[3] = 0.0
        genoModCountmat[4] = countmat[4] / sumClasses
        genoModCountmat[5] = 0.0
        genoModCountmat[6] = 0.0
        genoModCountmat[7] = 0.0
        genoModCountmat[8] = 0.0
        countmatIndex = drawSample(genoModCountmat)
        nextGeno = (countmatIndexToGenoPair(countmatIndex))[1]

    if prevGeno == "ab":
        sumClasses = countmat[3] + countmat[6] + countmat[7]
        genoModCountmat[0] = 0.0
        genoModCountmat[1] = 0.0
        genoModCountmat[2] = 0.0
        genoModCountmat[3] = countmat[3] / sumClasses
        genoModCountmat[4] = 0.0
        genoModCountmat[5] = 0.0
        genoModCountmat[6] = countmat[6] / sumClasses
        genoModCountmat[7] = countmat[7] / sumClasses
        genoModCountmat[8] = 0.0
        countmatIndex = drawSample(genoModCountmat)
        nextGeno = (countmatIndexToGenoPair(countmatIndex))[1]

    if prevGeno == "bb":
        sumClasses = countmat[2] + countmat[5] + countmat[8]
        genoModCountmat[0] = 0.0
        genoModCountmat[1] = 0.0
        genoModCountmat[2] = countmat[2] / sumClasses
        genoModCountmat[3] = 0.0
        genoModCountmat[4] = 0.0
        genoModCountmat[5] = countmat[5] / sumClasses
        genoModCountmat[6] = 0.0
        genoModCountmat[7] = 0.0
        genoModCountmat[8] = countmat[8] / sumClasses
        countmatIndex = drawSample(genoModCountmat)
        nextGeno = (countmatIndexToGenoPair(countmatIndex))[1]

    return nextGeno

def addNextMarker(genoList):
    rfClassProbs = [0.90, 0.09, 0.01]
    rfClassIndex = drawSample(rfClassProbs)
    rf = rfIndexToRF(rfClassIndex)
    geneticDistance = -50.0 * log(1.0 - 2.0 * rf)   #Haldane mapping function
    previousGeneticDistance = float(genoList[-1][3])
    cumulativeGeneticDistance = previousGeneticDistance + geneticDistance
    if (cumulativeGeneticDistance > (MAPSIZE - (0.05 * MAPSIZE))):
        rf = rfIndexToRF(3)   #As the map approaches MAPSIZE, ensure all markers are small rfs
        geneticDistance = -50.0 * log(1.0 - 2.0 * rf)   #Haldane mapping function
        previousGeneticDistance = float(genoList[-1][3])
        cumulativeGeneticDistance = previousGeneticDistance + geneticDistance
    if (len(genoList) == NUMMARK):
        rf = 0.5 * (1.0 - exp(-(MAPSIZE - previousGeneticDistance) / 50))   #Make the last marker sufficient to make the map MAPSIZE
        geneticDistance = -50.0 * log(1.0 - 2.0 * rf)   #Haldane mapping function
        previousGeneticDistance = float(genoList[-1][3])
        cumulativeGeneticDistance = previousGeneticDistance + geneticDistance
    countmat = generateCountmat(rf, GEN, HET)
    genoList.append(["Marker_" + str(len(genoList)), "1", str(rf), str(cumulativeGeneticDistance)])
    for i in range(4, len(genoList[0])):
        prevGeno = genoList[-2][i]
        genoList[-1].append(returnNextGenotype(countmat, prevGeno))


#############
#############
#Begin main
initialRF = 0.05
genoList = initializeGenotypes(initialRF, GEN, HET)
for i in range(NUMMARK - 2):
    addNextMarker(genoList)


outFile = open("simulated_genotypes.csv", 'w')
cumulativeGeneticDistance = 0.0
for i in range(len(genoList)):
    for j in range(len(genoList[i])):
        if (j == len(genoList[i]) - 1):
            outFile.write(genoList[i][j] + "\n")
        else:
            outFile.write(genoList[i][j] + ",")
outFile.close()

outFile2 = open("simulated_genotypes_rqtl.csv", 'w')
#Column 0 is marker name
#Column 1 is the chromosome
#Column 2 is the rf and is to 0 to imitate genetic distance
#Column 3 is ignored
for i in range(len(genoList)):
    for j in range(len(genoList[i])):
        if (j == len(genoList[i]) - 1):
            outFile2.write(genoList[i][j] + "\n")
        elif ( j == 3):
            continue
        elif ( i == 0 and j == 2):
            outFile2.write(",")
        elif ( i != 0 and j == 2):
            outFile2.write("0" + ",")
        else:
            outFile2.write(genoList[i][j] + ",")
outFile2.close()

def generateError(genotype):
    altTypeProbs = [0.5, 0.5]
    altIndex = drawSample(altTypeProbs)
    if genotype == "aa":
        if altIndex == 0:
            return "ab"
        if altIndex == 1:
            return "bb"
    if genotype == "ab":
        if altIndex == 0:
            return "aa"
        if altIndex == 1:
            return "bb"
    if genotype == "bb":
        if altIndex == 0:
            return "aa"
        if altIndex == 1:
            return "bb"

outFile3 = open("simulated_genotypes_errors_and_missing_rqtl.csv", 'w')
#Column 0 is marker name
#Column 1 is the chromosome
#Column 2 is the rf and is to 0 to imitate genetic distance
#Column 3 is ignored
for i in range(len(genoList)):
    for j in range(len(genoList[i])):
        if (j == len(genoList[i]) - 1):
            outFile3.write(genoList[i][j] + "\n")
        elif ( j == 3):
            continue
        elif ( i == 0 and j == 2):
            outFile3.write(",")
        elif ( i != 0 and j == 2):
            outFile3.write("0" + ",")
        elif ( i != 0 and j > 3 ):
            genotype = genoList[i][j]
            errorProbs = [(1.0 - ERRORRATE), ERRORRATE]
            missingProbs = [(1.0 - MISSINGRATE), MISSINGRATE]
            #0 is correct, 1 is error or missing; error overrides missing if both occur
            isError = drawSample(errorProbs)
            isMissing = drawSample(missingProbs)
            if (isError == 1):
                genotype = generateError(genotype)
                outFile3.write(genotype + ",")
                #print("Error generated at row " + str(i) + " col " + str(j))
            elif (isMissing == 1):
                genotype = "-"
                outFile3.write(genotype + ",")
            else:
                outFile3.write(genotype + ",")
        else:
            outFile3.write(genoList[i][j] + ",")
outFile3.close()
