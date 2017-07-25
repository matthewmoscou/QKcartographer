#! /usr/bin/python

"""
%prog prefix population_type genotypes
Assesses segregation distortion within a population

prefix denotes the population identifier
population_type must adher to QTL Cartographer experimental design codes

Author: Matthew Moscou <matthew.moscou@tsl.ac.uk>
This performs the following:
   1. Imports genotypic data
   2. Generates R plot of segregation distortion

Improvements from existing script set (13 June 2017):

Future improvements to include:
"""

# modules
import commands

import math
from math import *

import optparse
from optparse import OptionParser

import os

import string


# global variables
paletteAHB = ["#E69F00", "#939393", "#56B4E9", "#D55E00"]
paletteAH = ["#E69F00", "#939393", "#D55E00"]
paletteBH = ["#E69F00", "#939393", "#D55E00"]
paletteAB = ["#E69F00", "#56B4E9", "#D55E00"]
color_name = ['Orange', 'Grey', 'Sky blue', 'Vermillion']
allele_name = {'A':'A', 'H':'H', 'B':'B', '-':'Missing'}


# import arguments and options
usage = "usage: %prog prefix population_type genotypes"
parser = OptionParser(usage=usage)
(options, args) = parser.parse_args()


# import population experimental design code
population_type = args[1]


# import genetic map
genetic_map_file = open(args[2], 'r')

mapped_markers = []
genetic_map = []
marker_distance = {}
chromosome_markers = {}
chromosomes = []

truth = False

for line in genetic_map_file.readlines():
	line = string.replace(line, '\n', '')
	sline = string.split(line, '\t')
	
	if truth:
		mapped_markers.append(sline[0])
		marker_distance[sline[0]] = float(sline[2])

		if sline[1] not in chromosome_markers.keys():
			chromosome_markers[sline[1]] = []
			chromosomes.append(sline[1])

		chromosome_markers[sline[1]].append(sline[0])
		genetic_map.append([sline[0], sline[3:]])

		if population_size != len(sline[3:]):
			print 'Unequal progeny in population at marker ' + sline[0]

		population_size = len(sline[3:])
	else:
		population_size = len(sline[3:])

	truth = True

genetic_map_file.close()


# set population type for analysis
if args[1] in ['RI0', 'RI1']:
	alleles = ['A', 'B', '-']
	palette = paletteAB
elif args[1] in ['B1']:
	alleles = ['A', 'H', '-']
	palette = paletteAH
elif args[1] in ['B2']:
	alleles = ['B', 'H', '-']
	palette = paletteBH
elif args[1][0] in ['S']:
	alleles = ['A', 'H', 'B', '-']
	palette = paletteAHB

# generate R plot
R_input_file = open(args[0] + '_segregation_distortion.txt', 'w')

R_input_file.write('Chromosome' + '\t' + 'cM' + '\t' + 'Allele' + '\t' + 'Frequency' + '\n')

for chromosome in chromosomes:
	for marker in chromosome_markers[chromosome]:
		for allele in alleles:
			R_input_file.write(chromosome + '\t' + str(marker_distance[marker]) + '\t' + allele_name[allele] + '\t' + str((genetic_map[mapped_markers.index(marker)][1].count(allele) * 1.0) / population_size) + '\n')

R_input_file.close()

named_alleles = []

for allele in alleles:
	named_alleles.append(allele_name[allele])

dir_path = os.path.dirname(os.path.realpath(__file__))

R_analysis_file = open(args[0] + '_segregation_distortion.R', 'w')

R_analysis_file.write('library(ggplot2)' + '\n')
R_analysis_file.write('setwd("' + dir_path + '")\n')
R_analysis_file.write('data = read.table(file="' + args[0] + '_segregation_distortion.txt", header=T)' + '\n')
R_analysis_file.write('data = data.frame(data)' + '\n')
R_analysis_file.write('data$Allele <- factor(data$Allele, levels=c("' + '","'.join(named_alleles) + '"))' + '\n')

#R_analysis_file.write('postscript(file="' args[0] + '_segregation_distortion.ps", width=9, height=3)' + '\n')
R_analysis_file.write('png(file="' + args[0] + '_segregation_distortion.png", width=1500, height=500)' + '\n')
R_analysis_file.write('ggplot(data, aes(cM, group=Chromosome)) + geom_line(aes(x=cM, y=Frequency, group=Allele, color=Allele)) + geom_point(aes(x=cM, y=Frequency, color=Allele)) + facet_grid (.~ Chromosome, scales = "free_x", space="free_x") + scale_colour_manual(values = c("' + '","'.join(palette) + '"))' + '\n')

R_analysis_file.write('dev.off()' + '\n')

R_analysis_file.close()

commands.getstatusoutput('R --vanilla < ' + args[0] + '_segregation_distortion.R > temp.txt')
commands.getstatusoutput('rm temp.txt')
