#! /usr/bin/python

"""
%prog
Generates pre- and post-processed QTL tables and scripts for epistasis analysis

Author: Matthew Moscou <matthew.moscou@tsl.ac.uk>
This performs the following:
   1. Imports QTL experiment-wide thresholds and significant QTL
   2. Exports publication quality tables for significant QTLs
   3. Optional: Reads in quality controlled QTLs and generates an R script for epistasis analysis

Improvements from existing script set:

Future improvements to include:
	User selection of hypothesis tests to include
		H0H1, H0H2, H0H3, H1H3, H2H3
	Initial work in this area is limited, as I need to adapt previous scripts to handle only a single test output
	Need to implement earlier scripts to handle epistasis in backcross populations

"""


## import modules
import optparse
from optparse import OptionParser

import os

import string


## global variables
LRTS_LOD = 4.605186364


## functions
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


# import arguments and options
usage = "usage: %prog population_type alpha"
parser = OptionParser(usage=usage)
(options, args) = parser.parse_args()


# establish population type
# SFX, RI0, RI1
if len(args[0]) >= 2:
	if args[0][:2] == 'SF':
		hypothesis_tests = ['H0H1', 'H0H3']
	elif args[0][:2] == 'RI':
		hypothesis_tests = ['H0H1']
	elif args[0][0] == 'B':
		hypothesis_tests = ['H2H3']


### import permutations results
permutation_file = open('qtlcart.maximum_likelihood_alpha_' + args[1], 'r')

# {trait_name:[H0:H1 EWT, H0:H3 EWT], ...}
trait_EWT = {}

for line in permutation_file.readlines():
	sline = string.split(line)

	if len(hypothesis_tests) > 1:
		trait_EWT[sline[0]] = [sline[4], sline[1]]
	else:
		trait_EWT[sline[0]] = [sline[1]]

permutation_file.close()


### import map file
map_file = open('analysis/qtlcart.map', 'r')

# {chromosome:chromosome_name, ...}
chr_ID = {}

# {chromosome:{marker_index:marker_name, ...}, ...}
chr_marker_ID = {}
truth_marker = False
truth_chromosome = False

for line in map_file.readlines():
	sline = string.split(line)

	if len(sline) > 1:
		if sline[1] == 'MarkerNames':
			if not truth_marker:
				truth_marker = True
			else:
				truth_marker = False
		if sline[1] == 'ChromosomeNames':
			if not truth_chromosome:
				truth_chromosome = True
			else:
				truth_chromosome = False

	
	if truth_marker:
		if len(sline) > 2:
			if sline[0] not in chr_marker_ID.keys():
				chr_marker_ID[sline[0]] = {}

			chr_marker_ID[sline[0]][sline[1]] = sline[2]

	if truth_chromosome:
		chr_ID[sline[0]] = sline[1]

map_file.close()


### import Eqtl file and output QTL table

#	0	1	2	3	4		5	6	7	8	9	10		11
#	Index	Trait	Chr	cM	Peak Marker	EWT	LRTS	LOD	AEE	PVE
#	Index	Trait	Chr	cM	Peak Marker	EWT	LRTS	LOD	AEE	DEE	DEE / AEE	PVE
traits = []
hypothesis_trait_data = {}

for test in hypothesis_tests:
	hypothesis_trait_data[test] = {}

eqtl_file = open('analysis/qtlcart_H0H1.eqt', 'r')

truth = True
index = 1
trait = ''

for line in eqtl_file.readlines():
	sline = string.split(line)

	if len(sline) > 0:
		if sline[0] == '#End':
			truth = False

	if truth:
		if len(sline) > 2:
			if len(trait) > 0:
				if sline[0] != '#':
					if float(sline[4]) > float(trait_EWT[trait][0]):
						hypothesis_trait_data['H0H1'][trait].append([index, trait, chr_ID[sline[1]], sline[3], chr_marker_ID[sline[1]][sline[2]], trait_EWT[trait][0], sline[4], str(float(sline[4]) / LRTS_LOD), sline[5], sline[7]])
						index += 1

			if sline[2] == 'trait':
				trait = sline[4]
				traits.append(sline[4])
				hypothesis_trait_data['H0H1'][trait] = []

eqtl_file.close()

QTL_table = open('QTL_analysis_table_H0H1_draft.txt', 'w')

QTL_table.write('Index' + '\t' + 'Trait' + '\t' + 'Chr' + '\t' + 'cM' + '\t' + 'Peak Marker' + '\t' + 'EWT' + '\t' + 'LRTS' + '\t' + 'LOD' + '\t' + 'AEE' + '\t' + 'PVE' + '\t' + 'Notes' + '\n')

for trait in hypothesis_trait_data['H0H1'].keys():
	for element in hypothesis_trait_data['H0H1'][trait]:
		line = str(element[0])

		for item in element[1:]:
			line += ('\t' + str(item))

		QTL_table.write(line + '\n')

QTL_table.close()

if len(hypothesis_tests) > 1:
	eqtl_file = open('analysis/qtlcart_H0H3.eqt', 'r')
	
	truth = True
	index = 1
	trait = ''
	
	for line in eqtl_file.readlines():
		sline = string.split(line)
	
		if len(sline) > 0:
			if sline[0] == '#End':
				truth = False
	
		if truth:
			if len(sline) > 2:
				if len(trait) > 0:
					if sline[0] != '#':
						if float(sline[4]) > float(trait_EWT[trait][1]):
							hypothesis_trait_data['H0H3'][trait].append([index, trait, chr_ID[sline[1]], sline[3], chr_marker_ID[sline[1]][sline[2]], trait_EWT[trait][1], sline[4], str(float(sline[4]) / LRTS_LOD), sline[5], sline[6], str(float(sline[6]) / float(sline[5])), sline[7]])
							index += 1
	
				if sline[2] == 'trait':
					trait = sline[4]
					hypothesis_trait_data['H0H3'][trait] = []
	
	eqtl_file.close()


	QTL_table = open('QTL_analysis_table_H0H3_draft.txt', 'w')

	QTL_table.write('Index' + '\t' + 'Trait' + '\t' + 'Chr' + '\t' + 'cM' + '\t' + 'Peak Marker' + '\t' + 'EWT' + '\t' + 'LRTS' + '\t' + 'LOD' + '\t' + 'AEE' + '\t' + 'DEE' + '\t' + 'DEE / AEE' + '\t' + 'PVE' + '\t' + 'Notes' + '\n')

	for trait in hypothesis_trait_data['H0H3'].keys():
		for element in hypothesis_trait_data['H0H3'][trait]:
			line = str(element[0])

			for item in element[1:]:
				line += ('\t' + str(item))

			QTL_table.write(line + '\n')

QTL_table.close()

raw_input('Press enter if QTL table has been finalized.')


### import finalized QTL table
hypothesis_trait_QTL = {}

for test in hypothesis_tests:
	hypothesis_trait_QTL[test] = {}

QTL_table = open('QTL_analysis_table_H0H1_finished.txt', 'r')

line = QTL_table.readline()
line = QTL_table.readline()


while line:
	sline = string.split(line)

	if sline[1] not in hypothesis_trait_QTL['H0H1'].keys():
		hypothesis_trait_QTL['H0H1'][sline[1]] = []
	
	hypothesis_trait_QTL['H0H1'][sline[1]].append([sline[2], sline[3]])

	line = QTL_table.readline()

QTL_table.close()

if len(hypothesis_tests) > 1:
	QTL_table = open('QTL_analysis_table_H0H3_finished.txt', 'r')
	
	line = QTL_table.readline()
	line = QTL_table.readline()
	
	while line:
		sline = string.split(line)
	
		if sline[1] not in hypothesis_trait_QTL['H0H3'].keys():
			hypothesis_trait_QTL['H0H3'][sline[1]] = []
		
		hypothesis_trait_QTL['H0H3'][sline[1]].append([sline[2], sline[3]])
	
		line = QTL_table.readline()
	
	QTL_table.close()


### export Rqtl epistasis analysis base script
# start with base information
# general models

R_input = open('epistasis_QTL_analysis.R', 'w')

R_input.write('# Set working directory' + '\n')
R_input.write('setwd("' + os.getcwd() + '")' + '\n')
R_input.write('\n')
R_input.write('# Read in Rqtl package' + '\n')
R_input.write('library(qtl)' + '\n')
R_input.write('\n')
R_input.write('QKdata = read.cross(format="qtlcart", file="analysis/qtlcart.cro", mapfile="analysis/qtlcart.map")' + '\n')
R_input.write('QKdata = calc.genoprob(QKdata, step=2, error.prob=0.001)' + '\n')
R_input.write('QKdata = sim.geno(QKdata, step=2, n.draws=128, err=0.001)' + '\n')
R_input.write('\n')

for trait in traits:
	R_input.write('########## ' + trait + '\n')
	R_input.write('### H0:H1' + '\n')

	if trait in hypothesis_trait_QTL['H0H1'].keys():
		if len(hypothesis_trait_QTL['H0H1'][trait]) > 1:
			chromosome = ''
			cM = ''
			QTL = ''

			for index in range(len(hypothesis_trait_QTL['H0H1'][trait])):
				if index == 0:
					chromosome += ('"' + hypothesis_trait_QTL['H0H1'][trait][index][0] + '"')
					cM += (hypothesis_trait_QTL['H0H1'][trait][index][1])
					QTL += ('Q' + str(index + 1))
				else:
					chromosome += (',"' + hypothesis_trait_QTL['H0H1'][trait][index][0] + '"')
					cM += (',' + hypothesis_trait_QTL['H0H1'][trait][index][1])
					QTL += (' + Q' + str(index + 1))

			R_input.write('chr = c(' + chromosome + ')' + '\n')
			R_input.write('pos = c(' + cM + ')' + '\n')

			for qtl_interaction in powersetOfSize(range(len(hypothesis_trait_QTL['H0H1'][trait])), 2):
				QTL += (' + Q' + str(qtl_interaction[0] + 1) + ':Q' + str(qtl_interaction[1] + 1))
		
			R_input.write('phenotype.formula = y ~ ' + QTL + '\n')
			R_input.write('\n')
			R_input.write('qtl = makeqtl(QKdata, chr=chr, pos=pos)' + '\n')
			R_input.write('phenotype.fq = fitqtl(QKdata, pheno.col=' + str(traits.index(trait) + 1) + ', qtl=qtl, formula=phenotype.formula)' + '\n')
			R_input.write('summary(phenotype.fq)' + '\n')
			R_input.write('\n')

	if len(hypothesis_tests) > 1:
		R_input.write('### H0:H3' + '\n')
	
		if trait in hypothesis_trait_QTL['H0H3'].keys():
			if len(hypothesis_trait_QTL['H0H3'][trait]) > 1:
				chromosome = ''
				cM = ''
				QTL = ''
	
				for index in range(len(hypothesis_trait_QTL['H0H3'][trait])):
					if index == 0:
						chromosome += ('"' + hypothesis_trait_QTL['H0H3'][trait][index][0] + '"')
						cM += hypothesis_trait_QTL['H0H3'][trait][index][1]
						QTL += ('Q' + str(index + 1))
					else:
						chromosome += (',"' + hypothesis_trait_QTL['H0H3'][trait][index][0] + '"')
						cM += (',' + hypothesis_trait_QTL['H0H3'][trait][index][1])
						QTL += (' + Q' + str(index + 1))
	
				R_input.write('chr = c(' + chromosome + ')' + '\n')
				R_input.write('pos = c(' + cM + ')' + '\n')
	
				for qtl_interaction in powersetOfSize(range(len(hypothesis_trait_QTL['H0H3'][trait])), 2):
					QTL += (' + Q' + str(qtl_interaction[0] + 1) + ':Q' + str(qtl_interaction[1] + 1))
			
				R_input.write('phenotype.formula = y ~ ' + QTL + '\n')
				R_input.write('\n')
				R_input.write('qtl = makeqtl(QKdata, chr=chr, pos=pos)' + '\n')
				R_input.write('phenotype.fq = fitqtl(QKdata, pheno.col=' + str(traits.index(trait) + 1) + ', qtl=qtl, formula=phenotype.formula)' + '\n')
				R_input.write('summary(phenotype.fq)' + '\n')
				R_input.write('\n')
				R_input.write('\n')
			else:
				R_input.write('\n')
		else:
			R_input.write('\n')

R_input.close()
