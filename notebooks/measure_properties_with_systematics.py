import numpy as np
from astropy.io import fits
import healpy as hp
import os
import pickle

#other dependences
import sys
sys.path.insert(0, '../obs_eff_codes/')
import snr_mag_lsst as sml

#dependences for dust:
import dustmaps
from dustmaps.sfd import SFDQuery
from astropy.coordinates import SkyCoord
from dustmaps.config import config
config['data_dir'] = '/global/cfs/cdirs/lsst/groups/PZ/PhotoZDC2/run2.2i_dr6_test/TESTDUST/mapdata' #update this path when dustmaps are copied to a more stable location!


#functions:
def scan_over_tracts(required_pixels, mag_bins, pz_bins, indir, nside=128, dered=False):
    """
    ===INPUTS===
    
    required pixels: a dictionary of pixels to extract; it is a list of length nobs containing the pixels 
    
    mag_bins: bins in magnitude in each band
    
    pz_bins: bins in photo-z
    
    indir: where the DC2 catalogue sits; in future this should be replaced by the GCR catalogue codes
    
    ===FUNCTION===
    
    For each set of the required pixels, compute:
    - median and mean mag_err_cModel, in bins of magnitude
    - median and mean mag_cModel - mag_truth, in bins of magnitude
    - mean true redshift, in bins of photo-z
    
    ###Note: if use ebv map need to look at de-reddened magnitudes
    
    ===OUTPUTS===
    
    A dictionary of all these properties
    
    (Option: save to a pickle file.)
    
    #return a list of selected pixels in each tract
    
    
    """
    
    systematics_keys = list(required_pixels.keys())
    nquantiles = len(required_pixels[systematics_keys[0]])
    print('%s sets of pixels.'%nquantiles)
        
    DATAHOLDER={}
    for syskey in systematics_keys:
        band = syskey[0]
        props = [
                    'magerr_%s_cModel'%band, 
                    'mag_%s_cModel - mag_%s_truth'%(band,band),
                    'redshift_truth',
                ]
        DATAHOLDER[syskey]={}
        for propkey in props:
            if 'mag' in propkey:
                nbins = len(mag_bins[band])-1
                outkey = propkey
            if 'redshift' in propkey:
                nbins = len(pz_bins)-1
                outkey = propkey+'-%s'%band
            DATAHOLDER[syskey][outkey]=np.zeros((nquantiles, nbins, 3))
            
    #here error checks:
    if(0):
        print('initialised dataholder.')
        print('level 0 keys: ', DATAHOLDER.keys())
        print('level 1 keys: ')
        for key in DATAHOLDER[syskey].keys():
            print(key, DATAHOLDER[syskey][key].shape) 
        return 0
    
    tract_info=np.loadtxt('/global/cscratch1/sd/qhang/DESC_DC2_obs-dr6/tract_info.txt')
    
    tot=0
    
    for ii, tract in enumerate(tract_info[:,0]):
        
        if tot%50==0:
            print('Done %s tracts.'%tot)
        
        #catdir='/global/cscratch1/sd/qhang/DESC_DC2_obs-dr6/DC2_obj_with_pz-tract-%s.fits'%int(tract)
        catdir = indir + '-tract-%s.fits'%int(tract)
        
        #check if directory exist:
        if os.path.isfile(catdir) == True:
            
            fin=fits.open(catdir)
            
            ra=fin[1].data['ra']
            dec=fin[1].data['dec']
            cat_pix=hp.ang2pix(nside, ra, dec, lonlat=True)
            #anymore selection to pass on? e.g. i magnitude > 24
            imag=fin[1].data['mag_i_cModel']
            selmag = np.where(imag>24)[0] 
            
            #if deredden needed, compute the dereddened magnitudes for the catalogue:
            #if dered==True:
                #mag_dered_matrix = deredden_magnitudes(fin)
            
            for syskey in systematics_keys:
                band = syskey[0]
                props = [
                    'magerr_%s_cModel'%band, 
                    'mag_%s_cModel - mag_%s_truth'%(band,band),
                    'redshift_truth',
                ]
                
                for ll in range(nquantiles):
                    pix2use = required_pixels[syskey][ll]
                    sel=np.in1d(cat_pix, pix2use)
                    num_ind=np.arange(len(sel))
                    sel=num_ind[sel]
                    sel=np.intersect1d(sel, selmag)

                    if len(sel)>0:
                        catalog_prop_dist=get_mean_meansq_len(fin, sel, props, mag_bins, pz_bins, band=band, dered=dered)
                        #for each key, this will return a M*3 array where M is the number of binning in mag_bins or pz_bins
                        for key in catalog_prop_dist.keys():
                            #now prop_dist_alltracts[key] is nquantiles*M*3 array
                            DATAHOLDER[syskey][key][ll, :, :]+=catalog_prop_dist[key]
                
            tot+=1
    
    print('Computing mean and std')
    #now process these:
    MEANOUT={}
    STDOUT={}
    for syskey in systematics_keys:
        MEANOUT[syskey] = {}
        STDOUT[syskey] = {}
        for key in DATAHOLDER[syskey].keys():
            ndata = DATAHOLDER[syskey][key][:,:,2].astype(float)
            mean = DATAHOLDER[syskey][key][:,:,0]/ndata
            meansq = DATAHOLDER[syskey][key][:,:,1]/ndata
            std = np.sqrt(abs(meansq-mean**2))/np.sqrt(ndata)
            
            #print(mean.shape, std.shape)
            #print(mean[0,:], meansq[0,:], ndata[0,:])
            
            MEANOUT[syskey][key] = mean
            STDOUT[syskey][key] = std
        
    return MEANOUT, STDOUT


def get_mean_meansq_len(fin, sel, props, mag_bins, pz_bins, band = 'u', dered = False):
    
    #preset props
    #binned in magnitudes, or photo-z
    
    catalog_prop_dist={}
    
    for ii, key in enumerate(props):
        if '-' not in key:
            #check which quantity to bin:
            if 'magerr' in key:
                bins = mag_bins[band]
                outkey = key
                #now check if need to redden:
                binkey = 'mag_%s_cModel'%band
                
                if dered == False:
                    digi = np.digitize(fin[1].data[binkey][sel], bins)
                elif dered == True:
                    usemg = deredden_magnitudes(fin, sel, band=band)
                    digi = np.digitize(usemg, bins)
                                                
            if 'redshift' in key:
                bins = pz_bins
                binkey = 'photoz_mode'
                outkey = key+'-%s'%band
            
                digi = np.digitize(fin[1].data[binkey][sel], bins)
        
            #calculate mean, meansq, and ndata for each magnitude band:
            catalog_prop_dist[outkey] = np.zeros(((len(bins)-1), 3))
            for kk in range(len(bins)-1):
                sel2 = digi==kk+1
                #need to make sure that the sum is not nan or inf
                usedata = fin[1].data[key][sel][sel2]
                ind_nan = np.isnan(usedata)
                ind_inf = np.isinf(usedata)
                usedata = usedata[(~ind_nan)&(~ind_inf)]
                catalog_prop_dist[outkey][kk, 0] = np.sum(usedata)
                catalog_prop_dist[outkey][kk, 1] = np.sum(usedata**2)
                catalog_prop_dist[outkey][kk, 2] = len(usedata)
            
        elif '-' in key:
            item=key.split(' - ')
            bins = mag_bins[band]
            binkey = 'mag_%s_cModel'%band
            
            if dered == False:
                meas=fin[1].data[item[0]][sel]
            elif dered == True:
                meas = deredden_magnitudes(fin, sel, band=band)
            
            true=fin[1].data[item[1]][sel]
            diff=meas-true
            
            if dered == False:
                digi = np.digitize(fin[1].data[binkey][sel], bins)
            elif dered == True:
                usemg = deredden_magnitudes(fin, sel, band=band)
                digi = np.digitize(usemg, bins)
            
            catalog_prop_dist[key] = np.zeros(((len(bins)-1), 3))
            for kk in range(len(bins)-1):
                sel2 = digi==kk+1
                usedata = diff[sel2]
                ind_nan = np.isnan(usedata)
                ind_inf = np.isinf(usedata)
                usedata = usedata[(~ind_nan)&(~ind_inf)]
                catalog_prop_dist[key][kk, 0] = np.sum(usedata)
                catalog_prop_dist[key][kk, 1] = np.sum(usedata**2)
                catalog_prop_dist[key][kk, 2] = len(usedata)
    
    return catalog_prop_dist


#de-redden routine here:
def deredden_magnitudes(df, sel, band='u'):
    # set the A_lamba/E(B-V) values for the six LSST filters 
    maglist = np.array(['u','g','r','i','z','y'])
    band_a_ebv = np.array([4.81,3.64,2.70,2.06,1.58,1.31])
    
    coords = SkyCoord(df[1].data['ra'][sel], df[1].data['dec'][sel], unit = 'deg',frame='fk5')
    sfd = SFDQuery()
    ebvvec = sfd(coords)
    #df['ebv'] = ebvvec
    
    ind = np.where(maglist == band)
    #for i,band in enumerate(['u','g','r','i','z','y']):
    mag_dered= df[1].data[f'mag_{band}_cModel'][sel]-ebvvec*band_a_ebv[ind]
    
    return mag_dered



def select_pixels_from_sysmap(mapin, mask, map_bins, added_range=True):
    """
    map_bins: array of bin edges for the pixel selection, ascending
    added_range: whether include values < map_bins[0] or > map_bins[1]
    """
    
    pix = np.arange(len(mapin))
    usemap = mapin[mask.astype(bool)]
    usepix = pix[mask.astype(bool)]
    
    pixels=[]
    
    if added_range==True:
        ind = usemap < map_bins[0]
        pixels.append(usepix[ind])
    
    for ii in range(len(map_bins)-1):
        ind = (usemap >= map_bins[ii])&(usemap < map_bins[ii+1])
        pixels.append(usepix[ind])
        
    if added_range==True:
        ind = usemap >= map_bins[-1]
        pixels.append(usepix[ind])
    
    return pixels




#density

def scan_over_tracts_galmap(nside):
    """
    required pixels: nside=128
    bins: dictionary with matching keys to props, containing bin edges for the histogram
    """
    
    galmap=np.zeros(int(12*nside**2))
        
    tot=0
    for ii, tract in enumerate(tract_info[:,0]):
        
        if tot%50==0:
            print('Done %s tracts.'%tot)
        
        catdir='/global/cscratch1/sd/qhang/DESC_DC2_obs-dr6/DC2_obj_with_pz-tract-%s.fits'%int(tract)
        
        #check if directory exist:
        if os.path.isfile(catdir) == True:
            
            fin=fits.open(catdir)
            
            sel=select_objects(fin)
            
            if len(sel)>0:
            
                ra=fin[1].data['ra'][sel]
                dec=fin[1].data['dec'][sel]
                #check if inside the pixels input
                cat_pix=hp.ang2pix(nside, ra, dec, lonlat=True)

                cc=np.histogram(cat_pix, bins=np.arange(int(12*nside**2)+1))

                galmap+=cc[0]

            tot+=1
    
    return galmap


def binned_corr(densmap, sysmap, mask, ranges):
    
    #densmap
    Ndens = len(densmap.keys())
    Nbins = 9
    
    binedge=np.linspace(ranges[0], ranges[1], Nbins+1)
    bincen = ((binedge+np.roll(binedge,1))*0.5)[1:]

    bin_index = np.digitize(sysmap[mask.astype(bool)],binedge)

    #print(np.unique(bin_index))
    galnorm_z_stats=np.zeros((Nbins,2,Ndens))

    for jj in range(Ndens):
        for ii in range(Nbins):
            ind = np.where(bin_index==(ii+1))[0]
            #use this to compute the scatter around mean in the gal pixels
            galnorm_z_stats[ii,0,jj]=np.mean(densmap[jj][mask.astype(bool)][ind])
            npix = len(ind)
            galnorm_z_stats[ii,1,jj]=np.std(densmap[jj][mask.astype(bool)][ind])/np.sqrt(npix)

    return bincen, galnorm_z_stats


