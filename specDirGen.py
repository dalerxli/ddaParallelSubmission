##############################################################################
# importing packages
##############################################################################

import numpy as np
import os as os
import sys as sys
import shutil as shutil

##############################################################################
# check to see if proper files exist
##############################################################################

##############################################################################
'''
The necessary files to create simulation directories properly are:
	
	fPar	:	The parameter file with all but the desired calculation energy
				correctly input. Usually named "ddscat.par"
	fExec	:	The executable file for DDA. Usually named "ddscat"
	fSh		:	The shape file to be used. Will be renamed as "shape.dat" in 
				the simulation directories

These files should be located in the same directory as specDirGen.py	
'''
##############################################################################

fPar = "ddscat.par"
fExec = 'ddscat'
fSh = 'shape.dat'

pwd = os.getcwd()

print('---')
print('---')

if os.path.isfile(pwd + "/" + fPar):
	print("DDSCAT parameter file found : " + fPar)
else:
	print("DDSCAT parameter file not found. Exiting program.")
	sys.exit()

if os.path.isfile(pwd + "/" + fExec):
	print("DDSCAT executable file found : " + fExec)
else:
	print("DDSCAT executable file not found. Exiting program.")
	sys.exit()

if os.path.isfile(pwd + "/" + fSh):
	print("Shape file found : " + fSh)
else:
	print("Shape file not found. Exiting program.")
	sys.exit()

##############################################################################
# initialize calculation parameters
##############################################################################

'''Initial and final energies (in eVs) of light to be used in 
	simulation'''
enInit = 3.0
enFin = 4.0

'''Number of calculations to be performed and wavelength interval between
	successive calculations'''
nCalc = 10
enStep = (enFin - enInit)/nCalc

'''Array of desired wavelengths (in microns) for calculation.'''
wlList = 1239.841984 / np.arange(enInit,enFin,enStep) / 1000

##############################################################################
# create directories to run each simulation wavelength
##############################################################################

'''Create simulation directories. If they already exist, give user a prompt to
	delete them and recreate (or not).'''

try:
	for i in np.arange(1,nCalc + 1,1):
		os.mkdir( "W" + str(i) )
	print("Directories created.")

except OSError:
	pArg = input("Directories already exist. Delete and recreate? (y/n) ")

	if pArg is 'y':
		print("Deleting directories...")
		for j in np.arange(1,1000000,1):
			if os.path.isdir( "W" + str(j) ):
				shutil.rmtree( "W" + str(j) )
			else:
				print("Directories deleted.")
				break

		for i in np.arange(1,nCalc + 1,1):
			os.mkdir( "W" + str(i) )
		print("Directories created.")

	else:
		print("Directories were not created. Exiting program.")
		sys.exit()

##############################################################################
# copy proper files and edit the parameter files
##############################################################################

'''Copy ddscat.par files, create temporary ones with proper edits, and then
	delete temporary files.'''

wlFlag = 0

for i in np.arange(1,nCalc + 1,1):
	shutil.copyfile( fPar,"W" + str(i) + "/" + fPar )
	with open("W" + str(i) + "/" + fPar,"r") as oldFile:
		with open("W" + str(i) + "/" + fPar + "_temp","w+") as newFile:
			for line in oldFile.readlines():

				if wlFlag == 1:
					newFile.write( 
						"{0:.4f}".format(wlList[i-1]) + " " + 
						"{0:.4f}".format(wlList[i-1]) + " " + 
						str(1) + " " + "'INV' = wavelengths" + 
						" (first,last,how many,how=LIN,INV,LOG)\n"
						)
					wlFlag = 0
				else:
					newFile.write(line)

				if line == "'**** Wavelengths (micron) ****'\n":
					wlFlag = 1
				else:
					pass
		newFile.close()
	oldFile.close()

	os.remove("W" + str(i) + "/" + fPar)
	os.rename("W" + str(i) + "/" + fPar + "_temp", "W" + str(i) + "/" + fPar)

print('Parameter file copied.')

##############################################################################
# copy executable file
##############################################################################

for i in np.arange(1,nCalc + 1,1):
	shutil.copyfile( fExec,"W" + str(i) + "/" + fExec )
	shutil.copystat( fExec,"W" + str(i) + "/" + fExec )

print('Executable file copied.')

##############################################################################
# copy shape file
##############################################################################

for i in np.arange(1,nCalc + 1,1):
	shutil.copyfile( fSh,"W" + str(i) + "/shape.dat" )

print('Shape file copied.')
