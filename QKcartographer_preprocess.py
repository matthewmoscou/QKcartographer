#! /usr/bin/python

"""
%prog prefix population_type phenotypes genotypes
Converts genotypic and phenotypic data into QTL Cartographer input files

prefix denotes the population identifier
population_type must adher to QTL Cartographer experimental design codes

Author: Matthew Moscou <matthew.moscou@tsl.ac.uk>
This performs the following:
   1. Imports genotypic and phenotypic data
   2. Exports data into QTL Cartographer readable format

Improvements from existing script set (13 January 2017):

Future improvements to include:
"""

# modules
import math
from math import *

import optparse
from optparse import OptionParser

import string


# import arguments and options
parser = OptionParser()
(options, args) = parser.parse_args()


# import population experimental design code
population_type = args[1]


# import phenotypic trait data
trait_file = open(args[2], 'r')

number_individuals = -1
traits = {}

for line in trait_file.readlines():
	sline = string.split(line)

	if number_individuals >= 0:
		if number_individuals != len(sline[1:]):
			print 'Unequal individuals observed in trait:', sline[0]

	traits[sline[0]] = sline[1:]

	number_individuals = len(sline[1:])

trait_file.close()


# import genetic map
genetic_map_file = open(args[3], 'r')


mapped_markers = []
genetic_map = []
marker_distance = {}
chromosome_markers = {}

truth = False

for line in genetic_map_file.readlines():
	sline = string.split(line)
	
	if truth:
		mapped_markers.append(sline[0])
		marker_distance[sline[0]] = float(sline[2])

		if sline[1] not in chromosome_markers.keys():
			chromosome_markers[i] = []

		chromosome_markers[sline[1]].append(sline[0])
		genetic_map.append([sline[0], sline[3:]])

		if population_size != len(sline[3:]):
			print 'Unequal progeny in population at marker ' + sline[0]

		population_size = len(sline[3:])

	truth = True

genetic_map_file.close()


# identify the maximum number of markers and chromosome identifiers
maximum_markers = 0
chromosomes = chromosome_markers.keys()
chromosomes.sort()

for chromosome in chromosome_markers.keys():
	if len(chromosome_markers[chromosome]) > maximum_markers:
		maximum_markers = len(chromosome_markers[chromosome])
	


# generate Rmap file
rmap_output = open(args[0] + '_Rmap.inp', 'w')

rmap_output.write('# 6588597 bychromosome -filetype map.inp' + '\n')
rmap_output.write(' -type positions' + '\n')
rmap_output.write(' -function 2' + '\n')
rmap_output.write(' -Units cM' + '\n')
rmap_output.write(' -chromosomes ' + str(len(chromosomes)) + '\n')
rmap_output.write(' -maximum ' + str(maximum_markers) + '\n') # may need to adjust or generate estimator
rmap_output.write(' -named yes' + '\n')
rmap_output.write(' -start' + '\n')

for chromosome in chromosomes:
	rmap_output.write(' -Chromosome ' + chromosome + '\n')
	for marker in chromosome_markers[chromosome]:
		rmap_output.write('   ' + marker + '   ' + str(marker_distance[marker]) + '\n')

rmap_output.write(' -stop' + '\n')
rmap_output.write(' -end' + '\n')

rmap_output.close()


# generate Rcross file
rcross_output = open(args[0] + '_Rcross.inp', 'w')

rcross_output.write('# 6588597 -filetype cross.inp' + '\n')
rcross_output.write('-SampleSize ' + str(population_size) + '\n')
rcross_output.write('-Cross ' args[1] + '\n')
rcross_output.write('-traits ' + str(len(traits.keys())) + '\n')
rcross_output.write('-missingtrait -' + '\n')
rcross_output.write('-case yes' + '\n')
rcross_output.write('-TranslationTable' + '\n')
rcross_output.write(' AA    0     0' + '\n')
rcross_output.write(' Aa    1     1' + '\n')
rcross_output.write(' aa    2     2' + '\n')
rcross_output.write(' A-    12    12' + '\n')
rcross_output.write(' a-    10    10' + '\n')
rcross_output.write(' --    -1    -1' + '\n')
rcross_output.write('-start markers' + '\n')

for chromosome in chromosomes:
	for marker in chromosome_markers[chromosome]:
		rcross_output.write(marker)

		for element in genetic_map:
			if element[0] == marker:
				genotypes = element[1]

		for genetic_element in genotypes:
			if genetic_element == 'A':
				rcross_output.write(' 2')
			elif genetic_element == 'H':
				rcross_output.write(' 1')
			elif genetic_element == 'B':
				rcross_output.write(' 0')
			elif genetic_element == 'C':
				rcross_output.write(' 10')
			elif genetic_element == 'D':
				rcross_output.write(' 12')
			else:
				rcross_output.write(' -1')

		rcross_output.write('\n')
			
rcross_output.write('-stop markers' + '\n')
rcross_output.write('\n')
rcross_output.write('-missingtrait -' + '\n')
	

rcross_output.write('-start traits' + '\n')

for phenotype in traits.keys():
	rcross_output.write(phenotype)

	for trait_element in traits[phenotype]:
		rcross_output.write(' ' + trait_element)

	rcross_output.write('\n')
	
rcross_output.write('-stop traits' + '\n')
rcross_output.write('-quit' + '\n')

rcross_output.close()

