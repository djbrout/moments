#! /usr/global/paper/bin/python
from math import pi
import numpy as np

#from matplotlib.pyplot import *

import healpy as hp
import pyfits as pf
import os

pid = os.getpid()

def get_mem():
    lines = open('/proc/%d/status' % pid).readlines()
    print '\n'.join([L for L in lines if L.find('VmSize') != -1])

indir = 'bcc_v1.0_truth_orig/'
outdir = 'bcc_v1.0_hpix_photoz/'
#outdir = 'temp_output/'
form = '%.8e'

joe = '/data2/home/clampitt/bcc_v1.0/'
bcc_orig = '/data3/scratch/bcc_v1/aardvark_v1.0/truth/'
local = '/home/dbrout/bccml/corrected_healpix_v1.c/'
nside = 8

# Read command line argument
#if len(sys.argv) != 2:
#    sys.exit('Must provide one value.')
#tpix = int(sys.argv[1]) - 1

#ind = int(os.environ['JOB_ID']) - 1    
#print os.environ['SGE_TASK_ID']

ind = int(os.environ['SGE_TASK_ID']) - 1
#ind = 0
#tpix = np.loadtxt('potential_pix.txt')[ind]
print ind
import sys
#sys.exit()
tpix = np.loadtxt(joe+'des_pix.txt')[ind]
#tpix = np.loadtxt('potential_pix.txt')[0]

spixfile = 'source_files_into_pix%d.txt' % (tpix)
spix = np.loadtxt(joe+'source_pix_lists/' +spixfile)
print 'source pixels = ', spix

outfile = local+'aardvark_v1.0c_hpix_truth.%d.fit' % (tpix)

#key = ['ra', 'dec', 'm200', 'central', 'amag', 'tmag', 'z', 'photoz_gaussian']
#form = ['E', 'E', 'E', 'I', '5E', '5E', 'E', 'E']
key = ['ra', 'dec','z', 'photoz',
       'tmag','omag','amag','gamma1','gamma2','kappa']
form = ['E', 'E', 'E', 'E',
        'E', 'E', 'E', 'E', 'E', 'E']

pixlist = ['62',

           ]
esign = [-1,

         ]

knowsigns = False

ct = 0
new_data = {}
a = open(local+'report_v1.c_pix'+str(tpix)+'.txt','w')    
for j in range(len(spix)):

    if (spix[j] == 5000): continue

    esign_mult = 1
    if knowsigns:
        indx = -1
        for p in pixlist:
            indx += 1
            if spix[j] == p:
                esign_mult = esign[indx]
            
    a.write(outfile+'\n')
    infile = bcc_orig+'Aardvark_v1.0c_truth_des_rotated.%d.fit' % (spix[j])
    a.write(infile+'\n')
    
    hdulist = pf.open(infile)
    #print hdulist[1].columns, '\n\n'
    ra = hdulist[1].data.field('ra')
    dec = hdulist[1].data.field('dec')
    z = hdulist[1].data.field('z')

    amag_r = hdulist[1].data.field('AMAG')[:,1]
    m_r = hdulist[1].data.field('TMAG')[:,1]
    m_g = hdulist[1].data.field('TMAG')[:,0]
    gamma1 = hdulist[1].data.field('GAMMA1')*esign_mult
    gamma2 = hdulist[1].data.field('GAMMA2')*esign_mult
    omag_r = hdulist[1].data.field('OMAG')[:,1]
    kappa = hdulist[1].data.field('KAPPA')

    # Obtain gaussian photoz with appropriate scatter for LRG
    zgauss = np.random.normal(z, 0.03*(1.+z), len(z))

    # Shift RA to continuous interval
    con = (ra > 200.)
    ra[con] = ra[con] - 360.

    theta = 90. - dec
    loc = hp.ang2pix(nside, theta * pi/180., ra * pi/180.)
    pcut = (loc == tpix)

    ct = ct + len(ra[pcut])
    print 'j, num = ', j, '\t', ct
    print get_mem()

    for k in range(len(key)):
        if (key[k] == 'ra'):
            if (j == 0): new_data[key[k]] = ra[pcut]
            else: new_data[key[k]] = np.hstack((new_data[key[k]], ra[pcut]))
        elif (key[k] == 'abs_r'):
            if (j == 0): new_data[key[k]] = abs_r[pcut]
            else: new_data[key[k]] = np.hstack((new_data[key[k]], abs_r[pcut]))
        elif (key[k] == 'm_r'):
            if (j == 0): new_data[key[k]] = m_r[pcut]
            else: new_data[key[k]] = np.hstack((new_data[key[k]], m_r[pcut]))
        elif (key[k] == 'm_g'):
            if (j == 0): new_data[key[k]] = m_g[pcut]
            else: new_data[key[k]] = np.hstack((new_data[key[k]], m_g[pcut]))
        elif (key[k] == 'photoz'):
            if (j == 0): new_data[key[k]] = zgauss[pcut]
            else: new_data[key[k]] = np.hstack((new_data[key[k]], zgauss[pcut]))
        elif (key[k] == 'tmag'):
            if (j == 0): new_data[key[k]] = m_r[pcut]
            else: new_data[key[k]] = np.hstack((new_data[key[k]], m_r[pcut]))
        elif (key[k] == 'omag'):
            if (j == 0): new_data[key[k]] = omag_r[pcut]
            else: new_data[key[k]] = np.hstack((new_data[key[k]], omag_r[pcut]))
        elif (key[k] == 'amag'):
            if (j == 0): new_data[key[k]] = amag_r[pcut]
            else: new_data[key[k]] = np.hstack((new_data[key[k]], amag_r[pcut]))
        elif (key[k] == 'gamma1'):
            if (j == 0): new_data[key[k]] = gamma1[pcut]
            else: new_data[key[k]] = np.hstack((new_data[key[k]], gamma1[pcut]))
        elif (key[k] == 'gamma2'):
            if (j == 0): new_data[key[k]] = gamma2[pcut]
            else: new_data[key[k]] = np.hstack((new_data[key[k]], gamma2[pcut]))
        elif (key[k] == 'kappa'):
            if (j == 0): new_data[key[k]] = kappa[pcut]
            else: new_data[key[k]] = np.hstack((new_data[key[k]], kappa[pcut]))
        else:
            #print key[k]
            if (j == 0): new_data[key[k]] = hdulist[1].data.field(key[k])[pcut]
            else: new_data[key[k]] = np.hstack((new_data[key[k]], hdulist[1].data.field(key[k])[pcut]))
                                                
    hdulist.close()
a.close()
    
tmpcols = []
for i in range(len(key)):
    tmpcols.append(pf.Column(name=key[i], format=form[i],
                             array=new_data[key[i]]))

hdu = pf.PrimaryHDU()
tbhdu = pf.new_table(tmpcols)
thdulist = pf.HDUList([hdu, tbhdu])

thdulist.writeto(outfile)
thdulist.close()
