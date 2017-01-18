#! /usr/bin/python

"""
%prog prefix population_type alpha
Generates figures from QTL analyses using QTL Cartographer.

Author: Matthew Moscou <matthew.moscou@tsl.ac.uk>
This performs the following:
   1. Imports analyses from QTL Cartographer
   2. Imports results from permutation analysis
   3. Generates R script and runs R to generate plots of QTL analysis

Improvements from existing script set:

Future improvements to include:
	Implement the execution of scripts automatically, including permutations
	Add the ability to plot chromosomal maps of the following:
		Segregation distortion
		Bonferroni correction for multiple tests (0.05 / # of tests)
		Histogram of raw data from the population (ggplot2)
"""

# modules
import math
from math import log

import optparse
from optparse import OptionParser

import os

import sets

import string


## global variables
# hypothesis test information
hypothesis_tests_short = ['30', '31', '32', '10', '20']
hypothesis_tests = ['H0:H3', 'H1:H3', 'H2:H3', 'H0:H1', 'H0:H2']

# color palettes
palette = ["#D55E00", "#F0E442", "#009E73", "#E69F00", "#56B4E9"]
color_name = ['Vermillion', 'Yellow', 'Bluish green', 'Orange', 'Sky blue']
#chr_color = ["#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"]
#color_name = ['Orange', 'Sky blue', 'Bluish green', 'Yellow', 'Blue', 'Vermillion', 'Reddish Purple']
#palette = ['black', 'gold3', 'cyan2', 'red', 'blue']


# import arguments and options
usage = "usage: %prog prefix population_type alpha"
parser = OptionParser(usage=usage)
#parser.add_option("-x", "--execute", action="store_true", dest="execute", default=False, help="Execute analysis and permutations")
parser.add_option("-p", "--postscript", action="store_true", dest="postscript", default=False, help="Export postscript instead of PNG files")
(options, args) = parser.parse_args()


# functions
def str_vector(float_vector):
	string_vector = '('

	for element in float_vector:
		if len(string_vector) == 1:
			string_vector += str(element)
		else:
			string_vector += (', ' + str(element))
	
	string_vector += ')'

	return string_vector


# set working directory
filename = 'QTL_' + args[0] + '_' + args[1]


## IMPORT - maximum likelihood cutoffs
# example format
# trait		H0:H3		H1:H3		H2:H3		H0:H1		H0:H2
# T1		16.530298645	17.0694366	13.232320515	12.97881691	13.44939775
# T2		16.15350163	17.19938026	12.65990407	12.46655091	13.056576225
# ...

traits = []
trait_hypothesis_likelihood = {}

max_likelihood_file = open(os.getcwd() + '/qtlcart.maximum_likelihood_alpha_' + args[2], 'r')

truth = False

for line in max_likelihood_file.readlines():
	if truth:
		sline = string.split(line)

		traits.append(sline[0])

		trait_hypothesis_likelihood[sline[0]] = [float(sline[1]), float(sline[2]), float(sline[3]), float(sline[4]), float(sline[5])]
	
	truth = True

max_likelihood_file.close()

print ''

## VISUALIZATION - generate and execute R files for the visualization of QTL information

# import Zmapqtl output files
# storage of data in the 'trait_Zmapqtl'
# 0		1	2	3	4	5	6	7	8	9	10
# chr.bin	cM	H0:H3	H1:H3	H2:H3	H0:H1	H0:H2	add-H1	add-H3	dom-H2	dom-H3

# dictionary contains all information from Z files
# dataset -> trait ->
#   [[chr.bin, ..], [cM, ..], [LOD (H0:H3), ..], [LOD (H1:H3), ..], [LOD (H2:H3), ..], [LOD (H0:H1), ..], [LOD (H0:H2), ..], [add-H1, ..], [add-H3, ..], [dom-H2, ..], [dom-H3, ..]]
trait_Zmapqtl = {}
cM_positions = []
chr_cM_positions = {}

dataset = open(os.getcwd() + '/analysis/qtlcart.z', 'r')

line = dataset.readline()

start = False
initialization = True
old_chr_cM = 0.0

while line:
	if len(line) > 0:
		# capture trait ID
		if line in ['-Model           6      Model number\n', '-Model           3      Model number\n']:
			line = dataset.readline()
			sline = string.split(line)

			trait = sline[4][1:len(sline[4])-1]
			trait_Zmapqtl[trait] = [[], [], [], [], [], [], [], [], [], [], []]

		# start reading data for given trait
		if line == '-s\n':
			start = True
			line = dataset.readline()
	
		# end reading data for given trait
		if line == '-e\n':
			start = False
			initialization = False

		# if start has been initialized, store data for given trait
		if start:
			line = string.replace(line, ' inf ', ' 0.0 ')
			line = string.replace(line, ' -inf ', ' 0.0 ')
			line = string.replace(line, ' -666.0000000 ', ' 0.0 ')
			line = string.replace(line, ' 666.0000000 ', ' 0.0 ')
			line = string.replace(line, ' 783.3540119 ', ' 0.0 ')
			sline = string.split(line)

			if float(sline[1]) < 10:
				trait_Zmapqtl[trait][0].append(sline[0] + '.0' + sline[1])
			else:
				trait_Zmapqtl[trait][0].append(sline[0] + '.' + sline[1])

			#print line
			#print sline
			trait_Zmapqtl[trait][1].append(float(sline[2]) * 100)

			if float(sline[3]) < 0:
				trait_Zmapqtl[trait][2].append(0.0)
			else:
				trait_Zmapqtl[trait][2].append(float(sline[3]))

			if float(sline[4]) < 0:
				trait_Zmapqtl[trait][3].append(0.0)
			else:
				trait_Zmapqtl[trait][3].append(float(sline[4]))

			if float(sline[5]) < 0:
				trait_Zmapqtl[trait][4].append(0.0)
			else:
				trait_Zmapqtl[trait][4].append(float(sline[5]))

			if float(sline[10]) < 0:
				trait_Zmapqtl[trait][5].append(0.0)
			else:
				trait_Zmapqtl[trait][5].append(float(sline[10]))

			if float(sline[11]) < 0:
				trait_Zmapqtl[trait][6].append(0.0)
			else:
				trait_Zmapqtl[trait][6].append(float(sline[11]))

			trait_Zmapqtl[trait][7].append(float(sline[6]))
			trait_Zmapqtl[trait][8].append(float(sline[7]))
			trait_Zmapqtl[trait][9].append(float(sline[8]))
			trait_Zmapqtl[trait][10].append(float(sline[9]))

			# initialize cM positions for plotting
			if initialization:
				if int(sline[0]) not in chr_cM_positions.keys():
					if len(cM_positions) == 0:
						old_chr_cM = 0
					else:
						old_chr_cM = cM_positions[len(cM_positions) - 1] + 1.0

				cM_positions.append(float(sline[2]) * 100.0 + old_chr_cM)

				if int(sline[0]) not in chr_cM_positions.keys():
					chr_cM_positions[int(sline[0])] = []

				chr_cM_positions[int(sline[0])].append(float(sline[2]) * 100.0)

	line = dataset.readline()

dataset.close()

# in the future, add the following text based on a user flag
# this will add marker points below a plot
"""
# import marker positions (generation of tick marks on plot)
marker_positions_file = open(filename[:len(filename)-4] + '.IDs', 'r')

marker_positions = []

for line in marker_positions_file.readlines():
	sline = string.split(line)

	marker_positions.append(float(sline[2]))

marker_positions_file.close()
"""


os.popen('mkdir ' + os.getcwd() + '/figures')

analysis_traits = []

for suffix in traits:
	if suffix in trait_hypothesis_likelihood.keys():
		analysis_traits.append(suffix)

R_plot_graph = open(filename + '_figures.R', 'w')

# exporting scripts to generate graphs
for suffix in analysis_traits:
	trait = suffix

	print '\t', trait

	maximum_LOD = trait_hypothesis_likelihood[trait][0]

	for hypothesis_index in range(2, len(hypothesis_tests) + 2):
		for data_point in trait_Zmapqtl[trait][hypothesis_index]:
			if data_point != 0.0:
				if (data_point / (float(trait_hypothesis_likelihood[trait][hypothesis_index - 2]) / float(trait_hypothesis_likelihood[trait][0]))) > maximum_LOD:
					maximum_LOD = data_point / (float(trait_hypothesis_likelihood[trait][hypothesis_index - 2]) / float(trait_hypothesis_likelihood[trait][0]))

	maximum_additivity = max([abs(min(trait_Zmapqtl[trait][7])), abs(min(trait_Zmapqtl[trait][8])), max(trait_Zmapqtl[trait][7]), max(trait_Zmapqtl[trait][8])])
	maximum_dominance = max([abs(min(trait_Zmapqtl[trait][9])), abs(min(trait_Zmapqtl[trait][10])), max(trait_Zmapqtl[trait][10]), max(trait_Zmapqtl[trait][9])])

	############################################################################
	# plot graphs

	# open png or postscript file for LRTS/LOD plot
	if options.postscript:
		R_plot_graph.write('postscript(file="' + os.getcwd() + '/figures/' + filename + '_' + trait + '_LRTS.ps", width=9, height=3)' + '\n')
	else:
		R_plot_graph.write('png(file="' + os.getcwd() + '/figures/' + filename + '_' + trait + '_LRTS.png", width=1500, height=500)' + '\n')

	# plot LRTS/LOD curve of H0:H3
	R_plot_graph.write('plot(c' + str_vector(cM_positions) + ', c' + str_vector(trait_Zmapqtl[trait][2]) + ', type="l", xlab="", ylab="", main="", axes=F, xlim=c(' + str(max(cM_positions) * 0.038) + ', ' + str(max(cM_positions) * 0.962) + '), ylim=c(0, ' + str(maximum_LOD) + '), col="' + palette[0] + '", lwd=3)' + '\n')

	# plot experiment-wise threshold (EWT) across chromosomes
	R_plot_graph.write('lines(c(0, ' + str(max(cM_positions)) + '), c(' + str(trait_hypothesis_likelihood[trait][0]) + ', ' + str(trait_hypothesis_likelihood[trait][0]) + '), col="blue", lwd=3)' + '\n')
	
	R_plot_graph.write('box(lwd=3)' + '\n')

	#R_plot_graph.write('axis(1, at = c' + str_vector(marker_positions) + ', labels=F)' + '\n')
	
	# for replicate datasets, normalize to the maximum_likelihood of the average
	for hypothesis_index in range(2, len(hypothesis_tests) + 2):
		data_normalized = []
		
		for data_point in trait_Zmapqtl[trait][hypothesis_index]:
			#print data_point, trait_hypothesis_likelihood[trait][hypothesis_index - 2], trait_hypothesis_likelihood[trait][0]
			if data_point == 0.0:
				data_normalized.append(0.0)
			else:
				data_normalized.append(data_point / (float(trait_hypothesis_likelihood[trait][hypothesis_index - 2]) / float(trait_hypothesis_likelihood[trait][0])))

		R_plot_graph.write('lines(c' + str_vector(cM_positions) + ', c' + str_vector(data_normalized) + ', col="' + palette[hypothesis_index - 2] + '", lwd=3)' + '\n')


	chr_old = string.split(trait_Zmapqtl[trait][0][0], '.')[0]
	
	# determine chromosome positions
	for position in trait_Zmapqtl[trait][0]:
		chr_current = string.split(position, '.')[0]
		
		if chr_current != chr_old:
			R_plot_graph.write('lines(c(' + str(cM_positions[trait_Zmapqtl[trait][0].index(position)]) + ', ' + str(cM_positions[trait_Zmapqtl[trait][0].index(position)]) + '), c(' + str(maximum_LOD * -0.04) + ', ' + str(maximum_LOD * 1.04) + '), col="black", lwd=3)' + '\n')

		chr_old = chr_current

	R_plot_graph.write('dev.off()' + '\n')

	############################################################################

	# open PNG file for additivity plot
	if options.postscript:
		R_plot_graph.write('postscript(file="' + os.getcwd() + '/figures/' + filename + '_' + trait + '_additivity.ps", width=9, height=3)' + '\n')
	else:
		R_plot_graph.write('png(file="' + os.getcwd() + '/figures/' + filename + '_' + trait + '_additivity.png", width=1500, height=500)' + '\n')

	# plot LRTS/LOD curve of H0:H3
	R_plot_graph.write('plot(c' + str_vector(cM_positions) + ', c' + str_vector(trait_Zmapqtl[trait][8]) + ', type="l", xlab="", ylab="", main="", axes=F, xlim=c(' + str(max(cM_positions) * 0.038) + ', ' + str(max(cM_positions) * 0.962) + '), ylim=c(' + str(-1 * maximum_additivity) +', ' + str(maximum_additivity) + '), col="' + palette[4] + '", lwd=3)' + '\n')

	R_plot_graph.write('box(lwd=3)' + '\n')
	#R_plot_graph.write('axis(1, at = c' + str_vector(marker_positions) + ', labels=F)' + '\n')

	R_plot_graph.write('lines(c(0, ' + str(max(cM_positions)) + '), c(0, 0), col="black", lwd=3)' + '\n')
	R_plot_graph.write('lines(c' + str_vector(cM_positions) + ', c' + str_vector(trait_Zmapqtl[trait][7]) + ', col="' + palette[3] + '", lwd=3)' + '\n')

	chr_old = string.split(trait_Zmapqtl[trait][0][0], '.')[0]
	
	# determine chromosome positions
	for position in trait_Zmapqtl[trait][0]:
		chr_current = string.split(position, '.')[0]
		
		if chr_current != chr_old:
			R_plot_graph.write('lines(c(' + str(cM_positions[trait_Zmapqtl[trait][0].index(position)]) + ', ' + str(cM_positions[trait_Zmapqtl[trait][0].index(position)]) + '), c(' + str(maximum_additivity * -1.08) + ', ' + str(maximum_additivity * 1.08) + '), col="black", lwd=3)' + '\n')

		chr_old = chr_current

	R_plot_graph.write('dev.off()' + '\n')

	############################################################################

	# open PNG file for dominance plot
	if options.postscript:
		R_plot_graph.write('postscript(file="' + os.getcwd() + '/figures/' + filename + '_' + trait + '_dominance.ps", width=9, height=3)' + '\n')
	else:
		R_plot_graph.write('png(file="' + os.getcwd() + '/figures/' + filename + '_' + trait + '_dominance.png", width=1500, height=500)' + '\n')

	# plot LRTS/LOD curve of H0:H3
	R_plot_graph.write('plot(c' + str_vector(cM_positions) + ', c' + str_vector(trait_Zmapqtl[trait][10]) + ', type="l", xlab="", ylab="", main="", axes=F, xlim=c(' + str(max(cM_positions) * 0.038) + ', ' + str(max(cM_positions) * 0.962) + '), ylim=c(' + str(-1 * maximum_dominance) +', ' + str(maximum_dominance) + '), col="' + palette[4] + '", lwd=3)' + '\n')

	R_plot_graph.write('box(lwd=3)' + '\n')
	#R_plot_graph.write('axis(1, at = c' + str_vector(marker_positions) + ', labels=F)' + '\n')

	R_plot_graph.write('lines(c(0, ' + str(max(cM_positions)) + '), c(0, 0), col="black", lwd=3)' + '\n')
	R_plot_graph.write('lines(c' + str_vector(cM_positions) + ', c' + str_vector(trait_Zmapqtl[trait][9]) + ', col="' + palette[1] + '", lwd=3)' + '\n')

	chr_old = string.split(trait_Zmapqtl[trait][0][0], '.')[0]
	
	# determine chromosome positions
	for position in trait_Zmapqtl[trait][0]:
		chr_current = string.split(position, '.')[0]
		
		if chr_current != chr_old:
			R_plot_graph.write('lines(c(' + str(cM_positions[trait_Zmapqtl[trait][0].index(position)]) + ', ' + str(cM_positions[trait_Zmapqtl[trait][0].index(position)]) + '), c(' + str(maximum_dominance * -1.08) + ', ' + str(maximum_dominance * 1.08) + '), col="black", lwd=3)' + '\n')

		chr_old = chr_current

	R_plot_graph.write('dev.off()' + '\n')

	############################################################################

# run R script go generate figures

R_plot_graph.close()

os.popen('R --vanilla < ' + filename + '_figures.R > ' + filename + '.Routput' )
