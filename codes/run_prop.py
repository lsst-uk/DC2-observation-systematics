import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits

#import rubin_sim.maf as maf

import healpy as hp

import os

import sys

sys.path.insert(0, '../obs_eff_codes/')

import snr_mag_lsst as sml

import measure_properties_with_systematics as mp

#dump save
import pickle


#arguements:
import argparse

parser = argparse.ArgumentParser(description = 'Measure correlation between magnitudes/redshifts and systematics: ')
#input
parser.add_argument('-imaglim', type=float, nargs='+', default=[24,25], help='magnitude limits imposed on the sample selected')
parser.add_argument('-ibandsysonly', type=int, default=0, help='set 0 to use systematic maps in each band, set 1 to use systematic maps in i-band only')
parser.add_argument('-outfile', default='', help='directory to store the output plot')
parser.add_argument('-sysmap', default='Median_FWHMeff', help='systematic map to load')
parser.add_argument('-dered', type=bool, default=False, help='whether use de-reddened magnitudes')

args = parser.parse_args()

###---functions---###

def dump_save(stuff,filename):
    '''This saves the dictionary and loads it from appropriate files'''
    with open(filename,'wb') as fout:
        pickle.dump(stuff,fout,pickle.HIGHEST_PROTOCOL)
        #json.dump(self.impute, fout, sort_keys=True, indent=3)
    #print('written impute ditionary:',filename)
    return 0

def dump_load(filename):
    with open(filename,'rb') as fin:
        stuff=pickle.load(fin)
        #self.impute = json.load(fin)
    #print('loaded impute ditionary:',filename)
    return stuff

def find_percentiles(x, percent):
    sortind = np.argsort(x)
    xvalue_percent = np.zeros(len(percent))
    for ii in range(len(percent)):
        ind = int(len(x)*percent[ii])
        xvalue_percent[ii] = x[sortind[ind]]
    return xvalue_percent

###---functions---###



###---load systematic maps---###

#mask for Y5
savedir='/global/cscratch1/sd/qhang/DESC_DC2_obs-dr6/'
fname = savedir+'obj_footprint_mask-nside-128.fits'
mask = hp.read_map(fname)

band_list = ['u', 'g', 'r', 'i', 'z', 'y']

sysmap_group = ''

#maf maps
metric_name_list = ['Median_fiveSigmaDepth',
    'Median_filtSkyBrightness', 'Median_FWHMeff', 
    'CoaddM5']
#sup maps
proplist = ['nexp_sum', 'psf_size_wmean']
#other maps
otherlist = ['ebv', 'stellar']

sysmap_store = {}

if args.sysmap in metric_name_list:
    sysmap_group = 'maf'
    Opsimdir = '/global/cscratch1/sd/qhang/minion_1016/MAF-5year/desc_maf/'
    MAF=args.sysmap
    for ii in range(len(band_list)):
        band=band_list[ii]
        tag='_%s_%s_and_nightlt1825_HEAL.fits'%(MAF, band)
        title='%s-band %s'%(band, MAF)
            
        fname=Opsimdir+'minion_1016_dc2'+tag
        sysmap_store[title] = (hp.read_map(fname))*mask
    print('Loaded: %s map.'%MAF)
    
if args.sysmap in proplist:
    sysmap_group = 'sup'
    supmapdir='/global/cscratch1/sd/qhang/minion_1016/MAF-5year/supreme_map/'
    prop = args.sysmap
    for ii in range(len(band_list)):
        band=band_list[ii]
        fname = supmapdir + 'supreme_dc2_dr6d_v3_%s_%s-nside-128.fits'%(band, prop)
        title = '%s-band %s'%(band, prop)
        sysmap_store[title] = (hp.read_map(fname))*mask
    
if args.sysmap in otherlist:
    sysmap_group = 'other'
    if args.sysmap == 'ebv':
        othermapdir = '/global/cscratch1/sd/qhang/other_systematic_maps/'
        fname = othermapdir + 'ebv_ring_rot_512.fits'
        mapin = hp.read_map(fname)
        mapin = hp.ud_grade(mapin, 128)
        sysmap_store['ebv'] = mapin*mask
    elif args.sysmap == 'stellar':
        fname = othermapdir + 'allwise_total_rot_512.fits'
        mapin = hp.read_map(fname)
        mapin = hp.ud_grade(mapin, 128)
        sysmap_store['stellar'] = mapin*mask

###---load systematic maps---###


###---get pixels for percentiles in systematic maps---###

percent = [0.2, 0.4, 0.6, 0.8]

percentile={}

if sysmap_group != 'other':
    for ii in range(len(band_list)):
        band=band_list[ii]
        key='%s-band %s'%(band,args.sysmap)
        usemap = (sysmap_store[key])[mask.astype(bool)]
        percentile[key] = find_percentiles(usemap, percent)
        print(key, percentile[key])
elif sysmap_group == 'other':
    usemap = (sysmap_store[args.sysmap])[mask.astype(bool)]
    percentile[args.sysmap] = find_percentiles(usemap, percent)
    print(args.sysmap, percentile[args.sysmap])

###---get pixels for percentiles in systematic maps---###


###---binning---###

i_lim = np.array(args.imaglim,dtype='double')

mag_bins = {
    'u': np.linspace(21,28,8),
    'g': np.linspace(21,29,8),
    'r': np.linspace(21,29,8),
    'i': np.linspace(i_lim[0],i_lim[1],8),
    'z': np.linspace(20,27,8),
    'y': np.linspace(20,27,8),
}

#Y10 requirement
pz_bins = np.linspace(0.2, 1.2, 11)

mag_bins_cens = {}

for key in mag_bins.keys():
    mag_bins_cens[key] = ((mag_bins[key] + np.roll(mag_bins[key],1))*0.5)[1:]

pz_bins_cens = ((pz_bins + np.roll(pz_bins,1))*0.5)[1:]

###---binning---###


###===run the algorithm===###

indir = '/global/cscratch1/sd/qhang/DESC_DC2_obs-dr6/DC2_obj_with_pz'

#select the pixels:
selected_pix = {}

if sysmap_group != 'other':

    if args.ibandsysonly == 0:
        for ii in range(len(band_list)):
            band=band_list[ii]
            key='%s-band %s'%(band,args.sysmap)
            selected_pix[key] = mp.select_pixels_from_sysmap(sysmap_store[key], mask, percentile[key])

    if args.ibandsysonly == 1:
        for ii in range(len(band_list)):
            band=band_list[ii]
            key='%s-band %s'%(band,args.sysmap)
            dummy_key = 'i-band %s'%(args.sysmap)
            selected_pix[key] = mp.select_pixels_from_sysmap(sysmap_store[dummy_key],mask, percentile[dummy_key])
            
elif sysmap_group == 'other':
    for ii in range(len(band_list)):
        band=band_list[ii]
        key='%s-band %s'%(band,args.sysmap)
        selected_pix[key] = mp.select_pixels_from_sysmap(sysmap_store[args.sysmap], mask, percentile[args.sysmap])

mean, std = mp.scan_over_tracts(selected_pix, mag_bins, pz_bins, indir, dered=args.dered, imaglim = i_lim)

###===run the algorithm===###



###---generate plot and save---###

fig,axarr=plt.subplots(3,6,figsize=[18,9])

label=['q1', 'q2','q3','q4','q5']
ref_ind = 2

for ii, bb in enumerate(band_list):
    plt.sca(axarr[0, ii])
    key = '%s-band %s'%(bb, args.sysmap)
    subkey = 'magerr_%s_cModel'%bb
    data = mean[key][subkey]
    error = std[key][subkey]
    for kk in range(data.shape[0]):
        plt.errorbar(mag_bins_cens[bb], data[kk,:]/data[ref_ind,:], 
                     yerr=error[kk,:]/data[ref_ind,:], label=label[kk])
    plt.title(key)
    plt.ylabel(subkey+'/ref')
    plt.xlabel('mag_%s_cModel'%bb)
    plt.ylim([0.90, 1.10])
    plt.legend()
    
    plt.sca(axarr[1, ii])
    key = '%s-band %s'%(bb, args.sysmap)
    subkey = 'mag_%s_cModel - mag_%s_truth'%(bb,bb)
    data = mean[key][subkey]
    error = std[key][subkey]
    for kk in range(data.shape[0]):
        plt.errorbar(mag_bins_cens[bb], data[kk,:], yerr=error[kk,:], label=label[kk])
    plt.title(key)
    plt.ylabel(subkey)
    plt.xlabel('mag_%s_cModel'%bb)
    
    plt.sca(axarr[2, ii])
    key = '%s-band %s'%(bb, args.sysmap)
    subkey = 'redshift_truth-%s'%bb
    data = mean[key][subkey]
    error = std[key][subkey]
    for kk in range(data.shape[0]):
        plt.errorbar(pz_bins_cens, (data[kk,:]-data[ref_ind,:])/(1+data[ref_ind,:]), 
                     yerr=error[kk,:]/(1+data[ref_ind,:]), label=label[kk])
    plt.title(key)
    plt.ylabel(subkey+'/(1+ref)')
    plt.xlabel('photoz_mode')
    plt.ylim([-0.01,0.01])
    plt.fill_between(pz_bins_cens, -0.002, 0.002, 
                     alpha=0.05, color='k')

plt.tight_layout()


#tag:
tag = 'DC2-dr6-prop-w-'
tag = tag + args.sysmap
tag = tag + '-imaglim-%d-%d'%(args.imaglim[0],args.imaglim[1])
if args.ibandsysonly == 1:
    tag = tag + '-ibandsysonly'
if args.dered == True:
    tag = tag + '-dered'

outname = args.outfile + tag + '.png'

plt.savefig(outname, bbox_inches='tight')
print('Saved: %s'%outname)
    
###---generate plot and save---###