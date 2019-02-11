#! /usr/bin/env python
##############################################################################
'''
This program generates a ddscat.par file for running a DDA simulation using 
parameters defined in a simplified input parameter file and a given shape 
file. 

The input files can be named whatever the user desires, as long as the program
is run from the command line as:

'python ddscatInputGen.py [mode flag] [parameter file] [shape file]'

This program currently recognizes the following arguments from the simplified 
parameter file:

in_wl       :   The initial wavelength in microns to run in a DDA spectrum 
                calculation
dc_dir      :   The absolute path to the dielectric data for material 1
dc_dir_sub  :   The absolute path to the dielectric data for material 2
dc_dir_list :   The absolute path to a directory containing ONLY the
                dielectric data files to be used in simulation. The files will
                be read in alphanumeric order
n_back      :   The value to be used for the background dielectric constant

ebeam_pos   :   The value for the electron beam center position
ebeam_ke    :   The value for the electron beam kinetic energy

All ebeam parameters are only used if this program is run in eddscat (-e) 
mode. The syntax of the parameter files requires the parameter names to be 
spelled exactly as given above, each followed by a colon, a single space, and
a numeric or alphabetic parameter value. All lines preceded by '#' will be
treated as comments.

The shape file is needed to obtain the total number of polarizable 
points to be used in the calculation and thus calculate the effective radius
needed by the DDA algorithm.
'''
##############################################################################

##############################################################################
# importing packages
##############################################################################

import re
import sys
import math
import os

##############################################################################
# defining a function to import parameter input files and organize the data
##############################################################################

def main():
    """
    Define command line arguments for parameter input and organization.
    """
    from argparse import ArgumentParser, RawDescriptionHelpFormatter
    from textwrap import dedent
    parser = ArgumentParser(description=dedent(main.__doc__),
                            formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    parser.add_argument('-d', '--ddscat', action='store_const',const='ddscat',
        dest='mode', help='Generate the ddscat.par file.')
    parser.add_argument('-e', '--eddscat', action='store_const', 
        const='eddscat',dest='mode', 
        help='Generate the ddscat.par file for eDDA.')
    parser.add_argument('parIn', help='The parameter.input file.')
    parser.add_argument('sTemp', help='The shape_temp.dat file.')
    args = parser.parse_args()

    print('---')
    print('---')
    if args.mode == 'ddscat':
        print('Generating parameter file ddscat.par for PLANE WAVE' + 
            ' excitation source.')
    elif args.mode == 'eddscat':
        print('Generating parameter file ddscat.par for ELECTRON BEAM' + 
            ' excitation source.')    
    else:
        print('Please choose the proper flag to indicate excitation source.')
        print('Exiting program.')
        sys.exit()
    
    # Obtain all the parameters
    par = []
    val = []
    with open(args.parIn) as file:
         for line in file:
             if "#" not in line and len(line.split()) is not 0:
                par.append(line.split()[0][:-1])
                if 'ebeam_pos' in line:
                   val.append(line.split()[1:4])
                else:
                   val.append(line.split()[1])
    
    # Template files to be used
    shape_temp = args.sTemp

    # ddscatpar
    if args.mode == 'ddscat':
        ddscat_out = generate_ddscat(par, val, shape_temp)
        # write to file
        with open('ddscat.par', 'w') as file:
            file.writelines(ddscat_out)
            file.close()
        print('---')
        print('---')
        print('Parameter file generated as ddscat.par')

    # ddscatpar eDDA
    if args.mode == 'eddscat':
        eddscat_out = generate_eddscat(par, val, shape_temp)
        # write to file
        with open('ddscat.par', 'w') as file:
            file.writelines(eddscat_out)
        print('---')
        print('---')
        print('Parameter file generated as ddscat.par')

##############################################################################
# defining an alphanumeric sorting algorithm
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
# defining a function to generate DDA parameter file
##############################################################################

# Generate "ddscat.par" based on "shape_temp.dat" 
def generate_ddscat(par, val, shape_temp):
    in_wl   = val[par.index("in_wl")]
    dip_space = val[par.index("dip_space")]
    n_back    = val[par.index("n_back")]

    '''
    If dc_dir_list is defined, ddscatInputGen.py will copy the filenames (in 
    alphanumeric order) in dc_dir_list. Otherwise, ddscatInputGen will copy
    filenames given by dc_dir and, if defined, dc_dir_sub. 
    '''
    try:
        dc_dir_list = val[par.index("dc_dir_list")]
        print('---')
        print('---')
        print('Acquiring list of materials from directory.')
        print('Chosen materials directory is :',dc_dir_list)
        fileList = niceSort(os.listdir(dc_dir_list))
        print('Sorted list of dielectric files is :',fileList)
        ncomp = len(fileList)
    except ValueError:
        dc_dir = val[par.index("dc_dir")]
        print('---')
        print('---')
        print('Acquiring individual particle and substrate materials' + 
            ' from directory.')
        print('Chosen material 1 directory is :',dc_dir)
        try:
            dc_dir_sub = val[par.index("dc_dir_sub")]   
            ncomp = 2 # currently we only consider 2 dfferent materials
            print('Chosen material 2 directory is :',dc_dir_sub)
        except ValueError:
           ncomp = 1

    print('---')
    print('---')
    print('Chosen initial wavelength is :',in_wl,'(um)')
    print('Chosen background dielectric constant is :',n_back)
    print('Chosen dipole spacing is :',dip_space,'(nm)')
    pArg = input('Continue? (y/n) ')
    if pArg is 'y':
        pass
    else:
        print('Exiting program.')
        sys.exit()

    # memory allocation 
    x = []
    y = []
    z = []
    with open(shape_temp) as file:
         data = file.readlines ()
    for line in data[7:]:
        line = line.split()
        x.append(int(line[1]))
        y.append(int(line[2]))
        z.append(int(line[3]))
    mem_allo_x = max(x) - min(x) + 10
    mem_allo_y = max(y) - min(y) + 10
    mem_allo_z = max(z) - min(z) + 10

    # effective radius
    effR = (3 * len(x) / (4 * math.pi))**(1 / 3.0) * float(dip_space)*10**(-3)
    effR = "{0:.4f}".format(effR)
    
    out =  ("' ========== Parameter file for v7.1.0 ==================='\n")
    out += ("'**** Preliminaries ****'\n")
    out += ("'NOTORQ' = CMTORQ*6 (DOTORQ, NOTORQ) -- either do or skip " + 
        "torque calculations\n")
    out += ("'PBCGS2' = CMDSOL*6 (PBCGS2, PBCGST, PETRKP) -- select " + 
        "solution method\n")
    out += ("'GPFAFT' = CMETHD*6 (GPFAFT, FFTWJ, CONVEX)\n")
    out += ("'LATTDR' = CALPHA*6 (GKDLDR, LATTDR, SCLDR)\n")
    out += ("'NOTBIN' = CBINFLAG (ALLBIN, ORIBIN, NOTBIN)\n")
    out += ("'**** Initial Memory Allocation ****'\n")
    out += ("%r %r %r = dimensioning allowance for target generation\n" 
        % (mem_allo_x, mem_allo_y, mem_allo_z))
    out += ("'**** Target Geometry and Composition ****'\n")
    out += ("'FROM_FILE' = CSHAPE*9 shape directive\n")
    out += ("no SHPAR parameters needed\n")
    out += ("%r     = NBACK = background refractive index\n" % float(n_back))
    out += ("%r         = NCOMP = number of dielectric materials\n" 
        % int(ncomp))
    if 'dc_dir_list' in locals():
        for i in range( len(fileList) ):
            out += (("'%s' = file with refractive index for the chosen" + 
            " material " + str(i+1) + "\n")
             % (dc_dir_list + "/" + fileList[i]))
    else:
        out += (("'%s' = file with refractive index for the chosen" + 
            " material \n") % dc_dir)
        if ncomp == 2:
            out += ("'%s' = file with refractive index for the substrate \n" 
                % dc_dir_sub)
        else:
            pass
    out += ("'**** Error Tolerance ****'\n")
    out += ("1.00e-5 = TOL = MAX ALLOWED (NORM OF |G>=AC|E>-ACA|X>)" + 
        "/(NORM OF AC|E>)\n")
    out += ("'**** Integration limiter for PBC calculations ****'\n")
    out += ("1.00e-2 = GAMMA (1e-2 is normal, 3e-3 for greater accuracy)\n")
    out += ("'**** Angular resolution for calculation of <cos>, etc. ****'\n")
    out += ("0.5     = ETASCA (number of angles is proportional to " + 
        "[(3+x)/ETASCA]^2 )\n")
    out += ("'**** Wavelengths (micron) ****'\n")
    out += (("%r %r 1 'INV' = wavelengths (first,last,how many," + 
        "how=LIN,INV,LOG)\n") % (float(in_wl), float(in_wl)))
    out += ("'**** Effective Radii (micron) **** '\n")
    out += ("%r %r 1 'LIN' = aeff (first,last,how many,how=LIN,INV,LOG)\n" 
        % (float(effR), float(effR)))
    out += ("'**** Define Incident Polarizations ****'\n")
    out += ("(0,0) (1.,0.) (0.,0.) = Polarization state e01 " + 
        "(k along x axis)\n")
    out += ("1 = IORTH  (=1 to do only pol. state e01; =2 to also do orth. "+ 
        "pol. state)\n")
    out += ("'**** Specify which output files to write ****'\n")
    out += ('0 = IWRKSC (=0 to suppress, =1 to write ".sca" file for ' + 
        'each target orient.\n')
    out += ('1 = IWRPOL (=0 to suppress, =1 to write ".pol" file for each ' + 
        '(BETA,THETA)\n')
    out += ("'**** Specify Target Rotations ****'\n")
    out += ("0.    0.   1  = BETAMI, BETAMX, NBETA  " + 
        "(beta=rotation around a1)\n")
    out += ("0.    0.   1  = THETMI, THETMX, NTHETA " + 
        "(theta=angle between a1 and k)\n")
    out += ("0.    0.   1  = PHIMIN, PHIMAX, NPHI " + 
        "(phi=rotation angle of a1 around k)\n")
    out += ("'**** Specify first IWAV, IRAD, IORI (normally 0 0 0) ****'\n")
    out += ("0   0   0    = first IWAV, first IRAD, first IORI " + 
        "(0 0 0 to begin fresh)\n")
    out += ("'**** Select Elements of S_ij Matrix to Print ****'\n")
    out += ("6       = NSMELTS = number of elements of S_ij to print " + 
        "(not more than 9)\n")
    out += ("11 12 21 22 31 41       = indices ij of elements to print\n")
    out += ("'**** Specify Scattered Directions ****'\n")
    out += ("'TFRAME' = CMDFRM (LFRAME, TFRAME for Lab Frame or " + 
        "Target Frame)\n")
    out += ("1 = NPLANES = number of scattering planes\n")
    out += ("0.  0. 180. 1 = phi, theta_min, theta_max (deg) for plane A\n")
    out += ("50. 0. 180. 1 = phi, theta_min, theta_max (deg) for plane B\n")

    return out

##############################################################################
# defining a function to generate eDDA parameter file
##############################################################################

# Generate "ddscat.par" for eDDA based on "shape_temp.dat" 
def generate_eddscat(par, val, shape_temp):
    in_wl   = val[par.index("in_wl")]
    dip_space = val[par.index("dip_space")]
    n_back    = val[par.index("n_back")]
    ebeam_pos = val[par.index("ebeam_pos")]
    ebeam_ke  = val[par.index("ebeam_ke")]
    '''
    If dc_dir_list is defined, ddscatInputGen.py will copy the filenames (in 
    alphanumeric order) in dc_dir_list. Otherwise, ddscatInputGen will copy
    filenames given by dc_dir and, if defined, dc_dir_sub. 
    '''
    try:
        dc_dir_list = val[par.index("dc_dir_list")]
        print('---')
        print('---')
        print('Acquiring list of materials from directory.')
        print('Chosen materials directory is :',dc_dir_list)
        fileList = os.listdir(dc_dir_list)
        ncomp = len(fileList)
    except ValueError:
        dc_dir = val[par.index("dc_dir")]
        print('---')
        print('---')
        print('Acquiring individual particle and substrate materials' + 
            ' from directory.')
        print('Chosen material 1 directory is :',dc_dir)
        try:
            dc_dir_sub = val[par.index("dc_dir_sub")]   
            ncomp = 2 # currently we only consider 2 dfferent materials
            print('Chosen material 2 directory is :',dc_dir_sub)
        except ValueError:
           ncomp = 1

    print('---')
    print('---')
    print('Chosen initial wavelength is :',in_wl,'(um)')
    print('Chosen background dielectric constant is :',n_back)
    print('Chosen dipole spacing is :',dip_space,'(nm)')
    pArg = input('Continue? (y/n) ')
    if pArg is 'y':
        pass
    else:
        print('Exiting program.')
        sys.exit()
  
    # memory allocation 
    x = []
    y = []
    z = []
    with open(shape_temp) as file:
         data = file.readlines()
    for line in data[7:]:
        line = line.split()
        x.append(int(line[1]))
        y.append(int(line[2]))
        z.append(int(line[3]))
    mem_allo_x = max(x) - min(x) + 10
    mem_allo_y = max(y) - min(y) + 10
    mem_allo_z = max(z) - min(z) + 10

    # effective radius
    effR = (3 * len(x) / (4 * math.pi))**(1 / 3.0) * float(dip_space)*10**(-3)
    effR = "{0:.4f}".format(effR)
    
    out =  ("' ========== Parameter file for v7.1.0 ==================='\n")
    out += ("'**** Preliminaries ****'\n")
    out += ("'NOTORQ' = CMTORQ*6 (DOTORQ, NOTORQ) -- either do or skip" + 
        " torque calculations\n")
    out += ("'PBCGS2' = CMDSOL*6 (PBCGS2, PBCGST, PETRKP) -- select " + 
        "solution method\n")
    out += ("'GPFAFT' = CMETHD*6 (GPFAFT, FFTWJ, CONVEX)\n")
    out += ("'LATTDR' = CALPHA*6 (GKDLDR, LATTDR, SCLDR)\n")
    out += ("'NOTBIN' = CBINFLAG (ALLBIN, ORIBIN, NOTBIN)\n")
    out += ("'**** Initial Memory Allocation ****'\n")
    out += ("%r %r %r = dimensioning allowance for target generation\n" 
        % (mem_allo_x, mem_allo_y, mem_allo_z))
    out += ("%r %r %r = x, y, z position of e-beam center\n" 
        % (float(ebeam_pos[0]), float(ebeam_pos[1]), float(ebeam_pos[2])))
    out += ("%r = initial KE of electron\n" % float(ebeam_ke))
    out += ("'**** Target Geometry and Composition ****'\n")
    out += ("'FROM_FILE' = CSHAPE*9 shape directive\n")
    out += ("no SHPAR parameters needed\n")
    out += ("%r         = NCOMP = number of dielectric materials\n" 
        % int(ncomp))
    if 'dc_dir_list' in locals():
        for i in range( len(fileList) ):
            out += (("'%s' = file with refractive index for the chosen" + 
            " material " + str(i+1) + "\n")
             % (dc_dir_list + "/" + fileList[i]))
    else:
        out += (("'%s' = file with refractive index for the chosen" + 
            " material \n") % dc_dir)
        if ncomp == 2:
            out += ("'%s' = file with refractive index for the substrate \n" 
                % dc_dir_sub)
        else:
            pass
    out += ("'%s' = file with refractive index for the chosen material \n" 
        % dc_dir)
    if ncomp == 2:
       out += ("'%s' = file with refractive index for the substrate \n" 
        % dc_dir_sub)
    out += ("'**** Error Tolerance ****'\n")
    out += ("1.00e-5 = TOL = MAX ALLOWED (NORM OF |G>=AC|E>-ACA|X>)" + 
        "/(NORM OF AC|E>)\n")
    out += ("'**** Integration limiter for PBC calculations ****'\n")
    out += ("1.00e-2 = GAMMA (1e-2 is normal, 3e-3 for greater accuracy)\n")
    out += ("'**** Angular resolution for calculation of <cos>, etc. ****'\n")
    out += ("0.5     = ETASCA (number of angles is proportional to " + 
        "[(3+x)/ETASCA]^2 )\n")
    out += ("'**** Wavelengths (micron) ****'\n")
    out += (("%r %r 1 'INV' = wavelengths (first,last,how many," + 
        "how=LIN,INV,LOG)\n") % (float(in_wl), float(in_wl)))
    out += ("'**** Effective Radii (micron) **** '\n")
    out += ("%r %r 1 'LIN' = aeff (first,last,how many,how=LIN,INV,LOG)\n" 
        % (float(effR), float(effR)))
    out += ("'**** Define Incident Polarizations ****'\n")
    out += ("(0,0) (1.,0.) (0.,0.) = Polarization state e01 (k along " + 
        "x axis)\n")
    out += ("1 = IORTH  (=1 to do only pol. state e01; =2 to also do orth. "+ 
        "pol. state)\n")
    out += ("'**** Specify which output files to write ****'\n")
    out += ('0 = IWRKSC (=0 to suppress, =1 to write ".sca" file for each ' + 
        'target orient.\n')
    out += ('1 = IWRPOL (=0 to suppress, =1 to write ".pol" file for each ' + 
        '(BETA,THETA)\n')
    out += ("'**** Specify Target Rotations ****'\n")
    out += ("0.    0.   1  = BETAMI, BETAMX, NBETA  (beta=rotation around " + 
        "a1)\n")
    out += ("0.    0.   1  = THETMI, THETMX, NTHETA (theta=angle between " + 
        "a1 and k)\n")
    out += ("0.    0.   1  = PHIMIN, PHIMAX, NPHI (phi=rotation angle of " + 
        "a1 around k)\n")
    out += ("'**** Specify first IWAV, IRAD, IORI (normally 0 0 0) ****'\n")
    out += ("0   0   0    = first IWAV, first IRAD, first IORI (0 0 0 to " + 
        "begin fresh)\n")
    out += ("'**** Select Elements of S_ij Matrix to Print ****'\n")
    out += ("6       = NSMELTS = number of elements of S_ij to print (not " + 
        "more than 9)\n")
    out += ("11 12 21 22 31 41       = indices ij of elements to print\n")
    out += ("'**** Specify Scattered Directions ****'\n")
    out += ("'TFRAME' = CMDFRM (LFRAME, TFRAME for Lab Frame or Target " + 
        "Frame)\n")
    out += ("1 = NPLANES = number of scattering planes\n")
    out += ("0.  0. 180. 1 = phi, theta_min, theta_max (deg) for plane A\n")
    out += ("50. 0. 180. 1 = phi, theta_min, theta_max (deg) for plane B\n")

    return out

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.exit(1)
