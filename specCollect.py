##############################################################################
# importing packages
##############################################################################

import numpy as np
import os
import sys
import shutil
import glob
import re

##############################################################################
'''
	
'''
##############################################################################

##############################################################################
# define useful sorting functions
##############################################################################

def convert( text ):
	if text.isdigit():
		return int(text)
	else:
		return text

def alphanum_key( key ):
	return [ convert(c) for c in re.split('([0-9]+)',key) ]

def niceSort( l ):
	return sorted(l,key = alphanum_key)

##############################################################################
# get list of calculation directories
##############################################################################

pwd = os.getcwd()
calcDirList = niceSort(glob.glob(pwd + '/W*'))

##############################################################################
# create specData.txt file from qtable or gammatable entries in each 
# calculation directory
##############################################################################

'''
Check to see whether qtable or gammatable entries exist
'''
if os.path.isfile(calcDirList[0] + '/qtable'):
	calcType = 'pw'
	print('Plane wave calculation detected in : ' + calcDirList[0])
elif os.path.isfile(calcDirList[0] + '/gammatable'):
	calcType = 'eb'
	print('Electron beam calculation detected in : ' + calcDirList[0])
else:
	print('Calculation results were not found. Exiting program')
	sys.exit()

pArg = input('Continue? (y/n) ')
if pArg is 'y':
	pass
else:
	print('Exiting program.')
	sys.exit()

'''
Copy data from data tables. For the two calculation types:

	pw 	:	aData will contain along its columns the excitation energy (eV)
			and the extinction, absorption, and scattering cross-sections 
			(nm^2), respectively.

	eb	:	aData will contain along its columns the electron energy-loss
			energy (eV), electron energy-loss probability density (eV^{-1}),
			and cathodoluminescence power (?), respectively.

'''

print('Collecting data.')
data = []
# Open new data file
# with open('specData.txt','w+') as dataFile:


if calcType is 'pw':
	# Iterate through calculation directories
	for d in calcDirList:
		with open(d + '/qtable','r') as table:
			lines = table.readlines()
			# Find the data header in qtable and read the data beneath
			for i in range(len(lines)):
				if lines[i].strip()[0:4] == 'aeff':
					data = re.split(' +',lines[i+1].strip())
					data = np.array([float(element) for element in data])
					break
				else:
					pass
			# Place data in an array
			if d == calcDirList[0]:
				aData = data
			else:
				aData = np.vstack( (aData,data) )
		table.close()
	# Format data
	aData = np.delete(aData,[0,5,6,7,8],axis = 1)
	aData[:,0] = 1239.841984 / 1000 / aData[:,0]

	# Get effective radius form parameter file
	pArg = input('Obtain effective radius from ddscat.par? (y/n) ')
	if pArg is 'y':
		dirPar = 'ddscat.par'
	else:
		pArg2 = input('Please enter name of parameter file. ')
		dirPar = pArg2
	try:
		with open(dirPar,'r') as parFile:
			lines = parFile.readlines()
			for i in range(len(lines)):
				if lines[i][0:15] == "'**** Effective":
					radEff = float(re.split(' +',lines[i+1])[0]) * 1000
					break
				else:
					pass
	except FileNotFoundError:
		print('Parameter file not found. Exiting program.')
		sys.exit()
	# Convert efficiencies output by DDA to cross-sections
	aData[:,1:] = aData[:,1:] * np.pi * radEff**2

	hdr = ('Energy (eV)  sigma_ext (nm^2)  sigma_abs (nm^2)' + 
		'  sigma_sca (nm^2)')


elif calcType is 'eb':
	# Iterate through calculation directories
	for d in calcDirList:
		with open(d + '/gammatable','r') as table:
			lines = table.readlines()
			# Find the data header in gammatable and read the data beneath
			for i in range(len(lines)):
				if lines[i].strip()[0:3] == 'eq.':
					data = re.split(' +',lines[i+1].strip())
					data = np.array([float(element) for element in data])
				else:
					pass
			# Place data in an array
			if d == calcDirList[0]:
				aData = data
			else:
				aData = np.vstack( (aData,data) )

		table.close()

	aData = np.delete(aData,[0],axis = 1)

	hdr = 'Energy (eV)  Gamma_EELS (eV^{-1})  CL Power (?)'

print('Writing data to file.')
np.savetxt('specData.txt',aData,fmt = '%-7.3f',delimiter = ' ',
	header = hdr,comments = '#')
print('Data successfully written.')

# dataFile.close()
