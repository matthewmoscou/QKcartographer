#! /usr/bin/python

"""
%prog alpha number_of_permutations
Identifies experiment-wide threshold based on alpha from permuted analyses

Author: Matthew Moscou <matthew.moscou@tsl.ac.uk>
This performs the following:
   1. Imports permuted data
   2. Exports all peak QTL for additional analysis
   3. Exports experiment-wise thresholds based on alpha

Improvements from existing script set:

Future improvements to include:
"""

# modules
import math
from math import modf

import optparse
from optparse import OptionParser

import string


# functions
def quantile(data, k, presorted=False):
	__abstract__  = 'miscellaneous statistical functions'
	__copyright__ = 'Copyright (c) 2007-2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
	__license__   = 'See GLU license for terms by running: glu license'
	__revision__  = '$Id$'

	'''
	Return the k-th weighted quantile of values in data.  If the data are
	already sorted, then the presorted parameter should be set to avoid
	resorting the data.  Otherwise, a copy will be made and then sorted.

	@param      data: sequence of data elements that can be meaningfully ordered
	@type       data: if presorted, sequence.  If not presorted, iterable
	@param         k: the perventile value in the range [0..1]
	@type          k: float
	@param presorted: Indicator if the input data sequence is supplied in sorted order
	(default=False)
	@type  presorted: bool
	@return         : k-th weighted quantile of data

	Quantiles of data=0..16 from 0..1 every 0.5
	>>> observed = [ quantile(range(17), i/20.) for i in range(21) ]
	>>> expected = [0, 0.8, 1.6, 2.4, 3.2, 4, 4.8, 5.6, 6.4, 7.2, 8,
	...             8.8, 9.6, 10.4, 11.2, 12, 12.8, 13.6, 14.4, 15.2, 16]
	>>> # Verify that the sum of the errors is less than 10^-20
	>>> assert sum( abs(o-e) for o,e in zip(observed,expected) ) < 10e-20
	'''

	if not (0 <= k <= 1):
		raise ValueError('Quantile value %f out of range [0..1]' % k)

	if not presorted:
		data = sorted(data)

	if not data:
		raise ValueError('Input sequence may not be empty')

	n = len(data) - 1
	f,w = modf(n*k)
	w = int(w)

	if f < 1e-10:
		q = data[w]
	else:
		q = (1-f)*data[w] + f*data[w+1]

	return q


# import arguments and options
usage = "usage: %prog prefix alpha number_of_permutations"
parser = OptionParser(usage=usage)
(options, args) = parser.parse_args()


# read data from individual permutations
stem = 'qtlcart'

trait_test_maximum_value = {}
hypothesis_tests = ['H0:H3', 'H1:H3', 'H2:H3', 'H0:H1', 'H0:H2']
traits = []
trait = ''

for i in range(0,int(args[1])):
	permutation = open('permutations/' + stem + '_' + str(i) + '.z', 'r')
	line = permutation.readline()

	start = False
	
	while line:
		if len(line) > 0:
			if line in ['-Model           6      Model number\n', '-Model           3      Model number\n']:
				line = permutation.readline()
				sline = string.split(line)
			
				if i == 0 and len(trait) > 1:
					traits.append(trait)
					trait_test_maximum_value[trait] = {}
					trait_test_maximum_value[trait]['H0:H3'] = [maximum_value[0]]
					trait_test_maximum_value[trait]['H1:H3'] = [maximum_value[1]]
					trait_test_maximum_value[trait]['H2:H3'] = [maximum_value[2]]
					trait_test_maximum_value[trait]['H0:H1'] = [maximum_value[3]]
					trait_test_maximum_value[trait]['H0:H2'] = [maximum_value[4]]
				elif len(trait) > 1:
					trait_test_maximum_value[trait]['H0:H3'].append(maximum_value[0])
					trait_test_maximum_value[trait]['H1:H3'].append(maximum_value[1])
					trait_test_maximum_value[trait]['H2:H3'].append(maximum_value[2])
					trait_test_maximum_value[trait]['H0:H1'].append(maximum_value[3])
					trait_test_maximum_value[trait]['H0:H2'].append(maximum_value[4])
				
				trait = sline[4][1:len(sline[4])-1]
				maximum_value = [0.0, 0.0, 0.0, 0.0, 0.0]

			if line == '-s\n':
				start = True
				line = permutation.readline()
			if line == '-e\n':
				start = False
				
			if start:
				sline = string.split(line)

				if float(sline[3]) > maximum_value[0]:
					maximum_value[0] = float(sline[3])
				if float(sline[4]) > maximum_value[1]:
					maximum_value[1] = float(sline[4])
				if float(sline[5]) > maximum_value[2]:
					maximum_value[2] = float(sline[5])
				if float(sline[10]) > maximum_value[3]:
					maximum_value[3] = float(sline[10])
				if float(sline[11]) > maximum_value[4]:
					maximum_value[4] = float(sline[11])

		line = permutation.readline()
			
	if i == 0:
		traits.append(trait)
		trait_test_maximum_value[trait] = {}
		trait_test_maximum_value[trait]['H0:H3'] = [maximum_value[0]]
		trait_test_maximum_value[trait]['H1:H3'] = [maximum_value[1]]
		trait_test_maximum_value[trait]['H2:H3'] = [maximum_value[2]]
		trait_test_maximum_value[trait]['H0:H1'] = [maximum_value[3]]
		trait_test_maximum_value[trait]['H0:H2'] = [maximum_value[4]]

# export most significant QTL for every permutation and trait
maximum_likelihood = open(stem + '.maximum_likelihood_1000', 'w')

for trait in trait_test_maximum_value.keys():
	for test in trait_test_maximum_value[trait].keys():
		maximum_likelihood.write(trait + '\t' + test)

		for maximum_value in trait_test_maximum_value[trait][test]:
			maximum_likelihood.write('\t' + str(maximum_value))
		
		maximum_likelihood.write('\n')

maximum_likelihood.close()

# export experiment-wise treshold for every trait
maximum_likelihood_quantile_95 = open(stem + '.maximum_likelihood_alpha_' + args[0], 'w')

maximum_likelihood_quantile_95.write('trait')

for test in hypothesis_tests:
	maximum_likelihood_quantile_95.write('\t' + test)

maximum_likelihood_quantile_95.write('\n')

for trait in traits:
	maximum_likelihood_quantile_95.write(trait)

	for test in hypothesis_tests:
		maximum_likelihood_quantile_95.write('\t' + str(quantile(trait_test_maximum_value[trait][test], float(args[0]))))

	maximum_likelihood_quantile_95.write('\n')

maximum_likelihood_quantile_95.close()
