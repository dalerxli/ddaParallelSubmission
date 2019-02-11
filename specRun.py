##############################################################################
# importing packages
##############################################################################

import numpy as np
import os as os
import sys as sys
import shutil as shutil
import glob as glob
import subprocess

##############################################################################
# check to see if proper files and folders exist
##############################################################################

##############################################################################
'''
The necessary files and folders to submit simulation jobs are:
	
	fSl		:	The dummy .slurm file to be copied to each simulation 
				directory, edited, and run. Usually named 
				'launch_each_XXX.slurm', where 'XXX' is the name of the server
				partition to which the simulation jobs are being submitted.
	wlDlist	:	List of the directories in which jobs will be run for each 
				frequency requested. Usually named WN, where N is a number 
				indexing the wavelength to be run in that directory. Each 
				directory in the list is described with an absolute path.

These files should be located in the same directory as specRun.py	
'''
##############################################################################

fSl = "sampleCHEM.slurm"

pwd = os.getcwd()

wlDlist = glob.glob(pwd + '/W*/')

if os.path.isfile(pwd + "/" + fSl):
	print("SLURM submission example file found : " + fSl)
else:
	print("SLURM submission example file not found. Exiting program.")
	sys.exit()

if len(wlDlist) != 0:
	print("Simulation directories found : " + str(len(wlDlist)))
else:
	print("Simulation directories not found. Exiting program.")
	sys.exit()

##############################################################################
# copy files
##############################################################################

'''The jobName variable contains the job name header for the jobs to be 
	submitted. Each job will also include an index.'''
jobName = 'dftJob'

print('Job name header is : ' + jobName)

pArg = input('Continue? (y/n) ')
if pArg is 'y':
	pass
else:
	print('Exiting program.')
	sys.exit()

for i in range(len(wlDlist)):
	shutil.copyfile(pwd + '/' + fSl, pwd + '/W' + str(i+1) + '/' + fSl)

	with open(pwd + '/W' + str(i+1) + '/' + fSl,'r') as oldFile:
		with open(pwd + '/W' + str(i+1) + '/' + fSl + '_temp','w+') as newFile:

			for line in oldFile.readlines():
				if line[0:18] == '#SBATCH --job-name':
					newFile.write('#SBATCH --job-name=' + jobName + str(i+1) + 
					'\n' )
				elif line[0:17] == '#SBATCH --workdir':
					newFile.write('#SBATCH --workdir=' + pwd + '/W' + str(i+1)+
					 '\n')
				else:
					newFile.write(line)

		newFile.close()
	oldFile.close()

	os.remove(pwd + "/W" + str(i+1) + "/" + fSl)
	os.rename(pwd + "/W" + str(i+1) + "/" + fSl + "_temp",
		pwd + "/W" + str(i+1) + "/" + "W" + str(i+1) + '_' + fSl)

print('SLURM script copied.')

##############################################################################
# submit jobs
##############################################################################

print('Submitting jobs...')
for i in range(len(wlDlist)):
	subprocess.call(['sbatch',pwd + "/W" + str(i+1) + "/" + "W" +
		str(i+1) + '_' + fSl])
