import numpy as np
import scipy as sp
import matplotlib.pylab as plt
import os
import rdfits as r
import mytools
import sys
import pyfits as pf

def getsquareregion(file_dir,topleft_secondrow,topright_firstrow,bottomleft_thirdrow):
	tl_sr = pf.open(file_dir+topleft_secondrow)
	tr_fr = pf.open(file_dir+topright_firstrow)
	bl_tr = pf.open(file_dir+bottomleft_thirdrow)
	tl_sr_cols = tl_sr[1].data
	tr_fr_cols = tr_fr[1].data
	bl_tr_cols = bl_tr[1].data
	minra = min(np.asarray(tl_sr_cols['RA']))
	maxra = max(np.asarray(tr_fr_cols['RA']))
	mindec = min(np.asarray(bl_tr_cols['DEC']))
	maxdec = max(np.asarray(tl_sr_cols['DEC']))

	return [minra,maxra,mindec,maxdec]


#OPTIONS!
#z_lens_min = [.1,.1,.1,.1,.1,.1]
z_lens_min = [.1]
#z_lens_max = [.5,.5,.5,.5,.5,.5]
z_lens_max = [.5]
#z_src_min = [.5,.5,.5,.5,.5,.5]
z_src_min = [.5]
#z_src_max = [.6,.6,.6,1.2,1.2,1.2]
z_src_max = [1.2]
#mag_cut = [23.0,24.0,99,23.0,24.0,99]
mag_cut = [24.0]
#file_root = ['_magcut_23_src_.5to.6','_magcut_24_src_.5to.6','_magcut_99_src_.5to.6','_magcut_23_src_.5to1.2','_magcut_24_src_.5to1.2','_magcut_99_src_.5to1.2']
file_root = ['_magcut_23_src_.5to1.2']
#for j in range(len(z_lens_min)):


file_dir = "/home/dbrout/bccml/corrected_healpix_v1.c/"
OUTDIR = "/home/dbrout/moments/catalogs"

if not os.path.exists('figures'):
	os.makedirs('figures')
	
	title_in = ""
	
if not os.path.exists(OUTDIR):
	os.makedirs(OUTDIR)
	title = str(title_in)
newOUTDIR = OUTDIR+"/"

#500 Contiguous Sq Deg Area #1
fits = ['aardvark_v1.0c_hpix_truth.428.fit',
	'aardvark_v1.0c_hpix_truth.429.fit',
	'aardvark_v1.0c_hpix_truth.430.fit',
	'aardvark_v1.0c_hpix_truth.431.fit',
	'aardvark_v1.0c_hpix_truth.460.fit',
	'aardvark_v1.0c_hpix_truth.461.fit',
	'aardvark_v1.0c_hpix_truth.462.fit',
	'aardvark_v1.0c_hpix_truth.463.fit',
	'aardvark_v1.0c_hpix_truth.492.fit',
	'aardvark_v1.0c_hpix_truth.493.fit',
	'aardvark_v1.0c_hpix_truth.494.fit',
	'aardvark_v1.0c_hpix_truth.495.fit',
	'aardvark_v1.0c_hpix_truth.524.fit',
	'aardvark_v1.0c_hpix_truth.525.fit',
	'aardvark_v1.0c_hpix_truth.526.fit',
	'aardvark_v1.0c_hpix_truth.527.fit'
	]

esign = [1,
	 1,
	 1,
	 1,
	 1,
	 1,
	 1,
	 1,
	 1,
	 1,
	 1,
	 1,
	 1,
	 1,
	 1,
	 1
	 ]


topleft_secondrow = 'aardvark_v1.0c_hpix_truth.460.fit'
topright_firstrow = 'aardvark_v1.0c_hpix_truth.431.fit'
bottomleft_thirdrow = 'aardvark_v1.0c_hpix_truth.493.fit'

#500 Contiguous Sq Deg Area #2
#fits = ['aardvark_v1.0_hpix_truth.400.fit',
#	'aardvark_v1.0_hpix_truth.401.fit',
#	'aardvark_v1.0_hpix_truth.402.fit',
#	'aardvark_v1.0_hpix_truth.403.fit',
#	'aardvark_v1.0_hpix_truth.432.fit',
#	'aardvark_v1.0_hpix_truth.433.fit',
#	'aardvark_v1.0_hpix_truth.434.fit',
#	'aardvark_v1.0_hpix_truth.435.fit',
#	'aardvark_v1.0_hpix_truth.464.fit',
#	'aardvark_v1.0_hpix_truth.465.fit',
#	'aardvark_v1.0_hpix_truth.466.fit',
#	'aardvark_v1.0_hpix_truth.467.fit',
#	'aardvark_v1.0_hpix_truth.496.fit',
#	'aardvark_v1.0_hpix_truth.497.fit',
#	'aardvark_v1.0_hpix_truth.498.fit',
#	'aardvark_v1.0_hpix_truth.499.fit'
#	]

#topleft_secondrow = 'aardvark_v1.0_hpix_truth.432.fit'
#topright_firstrow = 'aardvark_v1.0_hpix_truth.403.fit'
#bottomleft_thirdrow = 'aardvark_v1.0_hpix_truth.464.fit'

print 'getting square region'
where = getsquareregion(file_dir,topleft_secondrow,topright_firstrow,bottomleft_thirdrow)
minra = where[0]
maxra = where[1]
mindec = where[2]
maxdec = where[3]



tables = []
indx=0
for fit in fits:
	indx+=1
	print indx
	tables.append(pf.open(file_dir+fit))
	print file_dir+fit

	

z = []
photoz = []
TMAGr = []
AMAGr = []
OMAGr = []
RA = []
DEC = []
GAMMA1 = []
GAMMA2 = []
K = []
#e1 = []
#e2 = [] 

print 'catting tables'
cols = []
index = -1
for table in tables:
	cc =table[1].columns
	index += 1
	col = table[1].data
	z.extend(np.asarray(col.field('z')))
	photoz.extend(np.asarray(col.field("PHOTOZ")))
	TMAGr.extend(np.asarray(col.field("TMAG")))
	AMAGr.extend(np.asarray(col.field('AMAG')))
	OMAGr.extend(np.asarray(col.field('OMAG')))
	RA.extend(np.asarray(col.field("RA")))
	DEC.extend(np.asarray(col.field("DEC")))
	GAMMA1.extend(np.asarray(col.field("GAMMA1")).astype(float)*esign[index])
	GAMMA2.extend(np.asarray(col["GAMMA2"]).astype(float)*esign[index])
	K.extend(np.asarray(col.field("KAPPA")))
	

z = np.asarray(z)
photoz = np.asarray(photoz)
TMAGr = np.asarray(TMAGr)
AMAGr = np.asarray(AMAGr)
OMAGr = np.asarray(OMAGr)
RA = np.asarray(RA)
DEC = np.asarray(DEC)
GAMMA1 = np.asarray(GAMMA1)
GAMMA2 = np.asarray(GAMMA2)
K = np.asarray(K) 
#e1 = np.asarray(e1)
#e2 = np.asarray(e2)

for j in range(len(z_lens_min)):
	print 'option '+str(j)
	ww = [(RA > minra) & (RA < maxra) & (DEC > mindec) & (DEC < maxdec) & (TMAGr < mag_cut[j])]
	zn = z[ww]
	photozn = photoz[ww]
	TMAGrn = TMAGr[ww]
	AMAGrn = AMAGr[ww]
	OMAGrn = OMAGr[ww]
	RAn = RA[ww]
	DECn = DEC[ww]
	GAMMA1n = GAMMA1[ww]
	GAMMA2n = GAMMA2[ww]
	Kn = K[ww]


	print 'step 1'
	bg = [(zn >z_src_min[j]) & (zn < z_src_max[j])] # background galaxies are where z > zcut = .5
	RAbg = RAn[bg]
	DECbg = DECn[bg]
	GAMMA1bg = GAMMA1n[bg]
	GAMMA2bg = GAMMA2n[bg]
	weightsbg = np.ones(GAMMA1n[bg].shape)
	zbg = zn[bg]
	TMAGrbg = TMAGrn[bg]
	AMAGrbg = AMAGrn[bg]
	OMAGrbg = OMAGrn[bg]
	

	print 'step 2'
	fg = [(zn < z_lens_max[j]) & (zn > z_lens_min[j])] # foreground galaxies are where z < zcut = .5
	RAfg = RAn[fg]
	DECfg = DECn[fg]
	GAMMA1fg = GAMMA1n[fg]
	GAMMA2fg = GAMMA2n[fg]
	weightsfg = np.ones(GAMMA1n[fg].shape)
	zfg = zn[fg]
	Kfg = Kn[fg]
	TMAGrfg = TMAGrn[fg]
	AMAGrfg = AMAGrn[fg]
	OMAGrfg = OMAGrn[fg]

	print 'plotting'
	fig = plt.figure()
	plt.hist(zn,30, normed=0)
	plt.xlabel("redshift")
	plt.ylabel("Counts")
	plt.title("z")
	fig.savefig("./figures/entire_distribution"+file_root[j]+".png")
	print 'saving'

	if os.path.exists(newOUTDIR+'foreground'+file_root[j]+'.fits'):
		os.remove(newOUTDIR+'foreground'+file_root[j]+'.fits')
	mytools.write_fits_table(newOUTDIR+'foreground'+file_root[j]+'.fits', ['z','RA','DEC','tmagr','omagr','amagr','kappa'], [zfg,RAfg,DECfg,TMAGrfg, OMAGrfg, AMAGrfg,Kfg])
	if os.path.exists(newOUTDIR+'background'+file_root[j]+'.fits'):
		os.remove(newOUTDIR+'background'+file_root[j]+'.fits')
	mytools.write_fits_table(newOUTDIR+'background'+file_root[j]+'.fits', ['RA','DEC','S1','S2','W','z','omagr'],
				 [RAbg,DECbg,GAMMA1bg,GAMMA2bg,weightsbg,zbg,OMAGrbg])
	if os.path.exists(newOUTDIR+'catalogue'+file_root[j]+'.fits'):
		os.remove(newOUTDIR+'catalogue'+file_root[j]+'.fits')
	mytools.write_fits_table(newOUTDIR+'catalogue'+file_root[j]+'.fits',
				 ['RA','DEC','z','S1','S2','TMAGr','OMAGr','AMAGr','KAPPA','PHOTOZ'],
				 [RA, DEC, z, GAMMA1, GAMMA2, TMAGr, OMAGr, AMAGr, K, photoz])

