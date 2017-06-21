#! /usr/bin/python

"""
%prog prefix population_type phenotypes
Generates informative plots for phenotypic data

prefix denotes the population identifier
population_type must adher to QTL Cartographer experimental design codes

Author: Matthew Moscou <matthew.moscou@tsl.ac.uk>
This performs the following:
   1. Imports phenotypic data
   2. Generates R plot of the distribution of phenotypic data (histograms) and pairwise plots

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


# functions
def powersetOfSize(set, size):
	if size > 1:
		for i in range(len(set) - (size-1)):
			for subset in powersetOfSize(set[i+1:], size-1):
				yield [set[i]] + subset
	elif size == 1:
		for i in set:
			yield [i]
	else:
		yield []


# global variables
paletteAHB = ["#E69F00", "#939393", "#56B4E9", "#D55E00"]
paletteAH = ["#E69F00", "#939393", "#D55E00"]
paletteAB = ["#E69F00", "#56B4E9", "#D55E00"]
color_name = ['Orange', 'Grey', 'Sky blue', 'Vermillion']
allele_name = {'A':'A', 'H':'H', 'B':'B', '-':'Missing'}


# import arguments and options
usage = "usage: %prog prefix population_type phenotypes"
parser = OptionParser(usage=usage)
(options, args) = parser.parse_args()


# import population experimental design code
population_type = args[1]


# import genetic map
phenotype_file = open(args[2], 'r')

phenotypes = {}
phenotype_order = []

for line in phenotype_file.readlines():
	line = string.replace(line, '\n', '')
	sline = string.split(line, '\t')

	phenotypes[sline[0]] = sline[1:]
	phenotype_order.append(sline[0])

phenotype_file.close()


# export data for R plot
R_input_file = open(args[0] + '_phenotypes_ggplot.txt', 'w')

R_input_file.write('Line')

for phenotype in phenotype_order:
	R_input_file.write('\t' + phenotype)

R_input_file.write('\n')

for index in range(len(phenotypes[phenotype_order[0]])):
	R_input_file.write(str(index))
	
	for phenotype in phenotype_order:
		if phenotypes[phenotype][index] != '-':
			R_input_file.write('\t' + phenotypes[phenotype][index])
		else:
			R_input_file.write('\t' + 'NA')

	R_input_file.write('\n')

R_input_file.close()

# generate R plot code
# histograms for every phenotype
# all pairwise combinations?
R_analysis_file = open(args[0] + '_phenotypes_ggplot.R', 'w')

dir_path = os.path.dirname(os.path.realpath(__file__))

R_analysis_file.write('library(ggplot2)' + '\n')
R_analysis_file.write('setwd("' + dir_path + '")\n')
R_analysis_file.write('data = read.table(file="' + args[0] + '_phenotypes_ggplot.txt", header=T)' + '\n')
R_analysis_file.write('data = data.frame(data)' + '\n')

for phenotype in phenotypes.keys():
	R_analysis_file.write('png(file="' + args[0] + '_histogram_' + phenotype + '.png", width=500, height=500)' + '\n')
	R_analysis_file.write('ggplot(data, aes(' + phenotype + ')) + geom_histogram()' + '\n')
	R_analysis_file.write('dev.off()' + '\n')

for phenotypic_pair in powersetOfSize(phenotypes.keys(), 2):
	R_analysis_file.write('png(file="' + args[0] + '_pairwise_' + phenotypic_pair[0] + '_' + phenotypic_pair[1] + '.png", width=500, height=500)' + '\n')
	R_analysis_file.write('ggplot(data, aes(' + phenotypic_pair[0] + ',' + phenotypic_pair[1] + ')) + geom_point()' + '\n')
	R_analysis_file.write('dev.off()' + '\n')

R_analysis_file.close()

commands.getstatusoutput('R --vanilla < ' + args[0] + '_phenotypes_ggplot.R > temp.txt')
commands.getstatusoutput('rm temp.txt')

