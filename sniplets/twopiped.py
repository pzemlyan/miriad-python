#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  6 18:38:45 2017

@author: pete
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 13:01:38 2016

MIRIAD mirpy3 pipeline

@author: pete
"""
from termcolor import colored
import miriad
import mirtask
from mirexec import *
from astropy import units as u;
import glob


def loadImg (url):
    import miriad
    import numpy as np
    v = miriad.ImData(url)
    xy = v.open('rw')    
    data = np.zeros((xy.axes[2],xy.axes[1],xy.axes[0]))
    for i in range(xy.axes[2]) :
        data[i] = xy.readPlane(axes=[i]).data
    return (data,xy)


def printout(t,s,cr= 'green',verbose=False):
    if (verbose):
        print(s)
    print(colored(t,cr))
    
    
freqs = loadtxt('co.cat',usecols=0)/1e3


import os
os.mkdir('fits')
os.mkdir('tmp')
    
sources = ('../nn/s255n_miriad/s255n-extended-usb.contsub','../nn/s255n_miriad/s255n-compact-usb.contsub','230.usb.nb','../co2-1/n_reg/mcm.uvm/')
sources2 = ('../nn/s255n_miriad/s255n-extended-lsb.contsub','../nn/s255n_miriad/s255n-compact-lsb.contsub','230.lsb.nb','../co2-1/n_reg/mcm.uvm/)



for uvr in [0,5,10,20,30,40]:
    try:
        printout(uvr,'','red')
        if freqs > 222:
            makemaps(sources,'CO',freqs,[0,200])
        else:
            makemaps(sources2,'ch3oh.229.758'.format(uvr),freqs,[0,200])
    except TaskFailError as e:
        printout('Wrong weights for {0},'.format(f),'red')


freqs = loadtxt('ocs.cat',usecols=0)/1e3


makemaps(('s255n_ch3oh2169.cal/',),'ch3oh_m_216',None,[0,200])
makemaps(('s255n_ch3oh2297.cal/',),'ch3oh_m_229',None,[0,200])
makemaps(('s255n_ch3oh2783.cal/',),'ch3oh_m_278',None,[0,200])


for f in freqs:
    try:
        if f > 222:
            makemaps(sources,'ocs_{0:0.4f}'.format(f),f,[0,200])
        else:
            makemaps(sources2,'ocs_{0:0.4f}'.format(f),f,[0,200])
    except TaskFailError as e:
        printout('Wrong weights for {0},'.format(f),'red')


def makemaps(src,mol,freq,uvr,verbose = False):
    r'''
    run pipeleine on specified files named by molecula
    
    Parameters
    -------------
    
    files : list
        calibrated visibility datasets
    names : list
        name for output files

    Returns
    ------------
    status of execution
    
    '''
    
    def find_nearest(array,value):
        idx = (np.abs(array-value)).argmin()
        return idx


    vdataset = ''
    for f in src:
        vdataset+=f+','
    vdataset=vdataset[:-1]

    if freq != None:    
        for f in src:
            t = TaskPutHead()
            t.in_ = f+'/restfreq'
            t.value = freq
            s = t.snarf()
        
        printout('frequnces set',s)

    
#    d = {}
#    for k in _keywords:
#        d[k] = ''
#    t.set(**d)
    m = 'tmp/'+mol
    t = TaskInvert()
    t.vis = vdataset
    t.map = m+".mp"
    t.beam = m+".bm"
    t.imsize =   [128,128]
#    t.select = 'uvrange({0},{1})'.format(uvr[0],uvr[1])
    t.cell = 0.6
    t.line = ['vel','240','-50','0.5','0.5']
    t.robust = 1
    t.mode = 'dft'
#    t.cell = 0.85955216805
#    t.imsize = 128,128
    t.options = 'double'
    s = t.snarf()

    printout('inverting finished [{0}]'.format(t.line),s,'red',True)
    t = TaskClean()
    t.map = m+".mp"
    t.beam = m+".bm"
    t.out = m+'.cl'
    t.niters = 300
    t.options = 'pad'
    t.phat=0.4
    t.cutoff = 0.3
    s = t.snarf()

    printout('Cleaning done',s)
    t = TaskRestore()
    t.map = m+".mp"
    t.beam = m+".bm"
    t.out = m+'.re'
    t.model = m+'.cl'
    s = t.snarf()
    printout('Restoring complete',s)
    
    t = TaskLinMos()
    t.in_ = m+".re"
    t.out = m+'.dt'
    s = t.snarf()
    printout('Restoring complete',s)

    t = TaskFits()
    t._in = m+'.dt'
    t.op = 'xyout'
    t.out = 'fits/'+mol+'.nosfc.fits'
    s = t.snarf()
    printout('1-st FITS file exported',s)
    
    


    
    
    
