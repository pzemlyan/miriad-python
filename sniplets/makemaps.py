        #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 15:08:15 2017

@author: pete
"""
def getExtent(recsf,recst):
    recs = recsf
    recs2 = recst
    from astropy import wcs
    w = wcs.WCS(recs[0].header)
    w2 = wcs.WCS(recs2[0].header)
    
    ix = arange(0,recs[0].header['NAXIS2'])
    iy = arange(0,recs[0].header['NAXIS1'])
    
    coords1 = zeros((len(iy),4))
    coords2 = zeros((len(ix),4))
    coords1[:,0] = iy
    coords2[:,1] = ix
    ra = w.all_pix2world(coords1,1,ra_dec_order=True)[:,0]
    dec = w.all_pix2world(coords2,1,ra_dec_order=True)[:,1]
    
    coords2[:,0] = ra[0]
    coords2[:,1] = dec
    
    coords1[:,0] = ra
    coords1[:,1] = dec[0]
    
    
    tox = w2.all_world2pix(coords2,1)[:,1]
    toy = w2.all_world2pix(coords1,1)[:,0]
    
    
    ix = arange(0,recs2[0].header['NAXIS2'])
    iy = arange(0,recs2[0].header['NAXIS1'])
    
    coords1 = zeros((len(iy),4))
    coords2 = zeros((len(ix),4))
    
    coords1[:,0] = iy
    coords2[:,1] = ix
    
    ra = w2.all_pix2world(coords1,1,ra_dec_order=True)[:,0]
    dec = w2.all_pix2world(coords2,1,ra_dec_order=True)[:,1]
    
    coords2[:,0] = ra[0]
    coords2[:,1] = dec
    
    coords1[:,0] = ra
    coords1[:,1] = dec[0]
    
    
    frx = w.all_world2pix(coords2,1)[:,1]
    fry = w.all_world2pix(coords1,1)[:,0]
    
    return [fry[0],fry[-1],frx[0],frx[-1]]
    

def getExtent2(recsf,recst):
    recs = recsf
    recs2 = recst
    from astropy import wcs
    w = wcs.WCS(recs[0].header)
    w2 = wcs.WCS(recs2[0].header)
    
    ix = arange(0,recs[0].header['NAXIS1'])
    iy = arange(0,recs[0].header['NAXIS2'])
    
    coords1 = zeros((len(iy),4))
    coords2 = zeros((len(ix),4))
    coords1[:,0] = iy
    coords2[:,1] = ix
    ra = w.all_pix2world(coords1,1,ra_dec_order=True)[:,0]
    dec = w.all_pix2world(coords2,1,ra_dec_order=True)[:,1]
    
    coords2[:,0] = ra[0]
    coords2[:,1] = dec
    
    coords1[:,0] = ra
    coords1[:,1] = dec[0]
    
    
    tox = w2.all_world2pix(coords2,1)[:,1]
    toy = w2.all_world2pix(coords1,1)[:,0]
    
    
    from scipy import interpolate
    
    fx = interpolate.interp1d(ix,tox)
    fy = interpolate.interp1d(iy,toy)
    
    return [fx,fy]

def getTX(recsf,recst):
    recs = recsf
    recs2 = recst
    from astropy import wcs
    w = wcs.WCS(recs[0].header)
    w2 = wcs.WCS(recs2[0].header)
    
    ix = arange(0,recs[0].header['NAXIS1'])
    iy = arange(0,recs[0].header['NAXIS2'])
    
    coords1 = zeros((len(iy),4))
    coords2 = zeros((len(ix),4))
    coords1[:,0] = iy
    coords2[:,1] = ix
    ra = w.all_pix2world(coords1,1,ra_dec_order=True)[:,0]
    dec = w.all_pix2world(coords2,1,ra_dec_order=True)[:,1]
    
    coords2[:,0] = ra[0]
    coords2[:,1] = dec
    
    coords1[:,0] = ra
    coords1[:,1] = dec[0]
    
    
    tox = w2.all_world2pix(coords2,1)[:,1]
    toy = w2.all_world2pix(coords1,1)[:,0]
    
    
    
    from scipy import interpolate
    
    fx = interpolate.interp1d(ix,tox)
    fy = interpolate.interp1d(iy,toy)
    
    return [fx,fy]


def getFromWorld(recs):
    from astropy import wcs
    w = wcs.WCS(recs[0].header)
    
    ix = arange(-100,recs[0].header['NAXIS1']+100)
    iy = arange(-100,recs[0].header['NAXIS2']+100)
    
    coords1 = zeros((len(iy),4))
    coords2 = zeros((len(ix),4))
    coords1[:,0] = iy
    coords2[:,1] = ix
    ra = w.all_pix2world(coords1,1,ra_dec_order=True)[:,0]
    dec = w.all_pix2world(coords2,1,ra_dec_order=True)[:,1]
    
    from scipy import interpolate
    
    fx = interpolate.interp1d(ra,ix)
    fy = interpolate.interp1d(dec,iy)
    
    return [fx,fy]

def getToWorld(recs):
    from astropy import wcs
    w = wcs.WCS(recs[0].header)
    
    ix = arange(-100,recs[0].header['NAXIS1']+100)
    iy = arange(-100,recs[0].header['NAXIS2']+100)
    
    coords1 = zeros((len(iy),4))
    coords2 = zeros((len(ix),4))
    coords1[:,0] = iy
    coords2[:,1] = ix
    ra = w.all_pix2world(coords1,1,ra_dec_order=True)[:,0]
    dec = w.all_pix2world(coords2,1,ra_dec_order=True)[:,1]
    
    from scipy import interpolate
    
    fx = interpolate.interp1d(ix,ra)
    fy = interpolate.interp1d(iy,dec)
    
    return [fx,fy]
