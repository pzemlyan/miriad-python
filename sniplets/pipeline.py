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
    
    
files = glob.glob('../*.cal')


import os
os.mkdir('fits')
os.mkdir('uv')
    
source = ('../s255n_c18o.cal')

uvout = selfcal(source,'c18o',.5)    

for f in files:
    try:
        printout(f,'','red')
        makecalmap(f,f[9:f.index('.cal')],uvout)
    except TaskFailError as e:
        printout('Wrong weights for '+moleculas,'','red')

def selfcal(file,moleculas,cutoff=2.5,verbose = False):
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


    
    vdataset = miriad.VisData(file)
    m = 'tmp.'+moleculas
    
    
    t = TaskUVSpec()
    t.interval = 1000
    t.options = 'avall'
    t.axis = 'v,amp'
    t.log = 'temp'
    t.vis = vdataset
    t.snarf()
    
    d = loadtxt('temp')
    
#    d = {}
#    for k in _keywords:
#        d[k] = ''
#    t.set(**d)
    t = TaskInvert()
    t.vis = vdataset
    t.map = m+".mp"
    t.beam = m+".bm"
    t.imsize =   [128,128]

    t.cell = 0.4
    cw = diff(d[:,0]).max()
    if (cw > 3) : cw = abs(diff(d[:,0])[0])
    bw = d[:,0].max()-d[:,0].min()
    t.line = ['vel',int(bw/cw)-6,d[:,0].min()+3*cw,cw,cw]
    t.robust = 1
#    t.cell = 0.85955216805
#    t.imsize = 128,128
    t.options = 'double,systemp'
    s = t.snarf()

    printout('inverting finished [{0}]'.format(t.line),s)
    imdata = loadImg(m+".mp")[0]
    std = imdata.std()
    t = TaskClean()
    t.map = m+".mp"
    t.beam = m+".bm"
    t.out = m+'.cl'
    t.niters = 100000
    t.cutoff = cutoff
    s = t.snarf()    

    printout('Cleaning done',s)
    t = TaskRestore()
    t.map = m+".mp"
    t.beam = m+".bm"
    t.out = m+'.re'
    t.model = m+'.cl'
    s = t.snarf()
    printout('Restoring complete',s)
    t = TaskFits()
    t._in = m+'.re'
    t.op = 'xyout'
    t.out = 'fits/'+moleculas+'.nosfc.fits'
    s = t.snarf()
    printout('1-st FITS file exported',s)
    
    i = 0
    
    lastuv = file 
    t = TaskUVAver()
    t.vis = lastuv
    lastuv = 'uv/'+moleculas+'{0}.phc.uv'.format(i)

    t.out = lastuv
    s = t.snarf()
    printout('Aplying inital gains',s)
    
    t = TaskSelfCal()
    t.vis = lastuv
    t.interval = 20
    t.model = m+'.cl'
    t.options = 'phase'
    s = t.snarf()
    printout('Phase calibration done',s,verbose=True)
    
    t = TaskUVSpec()
    t.interval = 1000
    t.options = 'avall'
    t.axis = 'v,amp'
    t.log = 'temp'
    t.vis = lastuv
    t.snarf()
    
    d = loadtxt('temp')

    
    files = glob.glob('tmp*')
    import shutil
    for f in files:
        shutil.rmtree(f)

    t = TaskInvert()
    t.vis = lastuv
    t.map = m+".mp"
    t.beam = m+".bm"
    t.imsize =   [128,128]

    t.cell = 0.4
    t.line = ['vel',int(bw/cw)-6,d[:,0].min()+3*cw,cw,cw]
    t.robust = 1
#    t.cell = 0.85955216805
#    t.imsize = 128,128
    t.options = 'double,systemp'
    s = t.snarf()
    printout('Inverting phased ',s)
    
    imdata = loadImg(m+".mp")[0]
    std = imdata.std()
    t = TaskClean()
    t.map = m+".mp"
    t.beam = m+".bm"
    t.out = m+'.cl'
    t.niters = 100000
    t.cutoff = cutoff
    s = t.snarf()    
    printout('Cleaning phased',s)
    

    while(i >= 0):
        i+=1
        t = TaskUVAver()
        t.vis = lastuv
        lastuv = 'uv/'+moleculas+'{0}.amc.uv'.format(i)
        t.out = lastuv
        s = t.snarf()
        printout('Aplying phased gains',s)
        
        t = TaskSelfCal()
        t.vis = lastuv
        t.interval = 20
        t.model = m+'.cl'
        t.options = 'amp'
        s = t.snarf()
        printout('Amp calibration done',s,verbose=True)
        
        t = TaskUVSpec()
        t.interval = 1000
        t.options = 'avall'
        t.axis = 'v,amp'
        t.log = 'temp'
        t.vis = lastuv
        t.snarf()
        
        d = loadtxt('temp')
    
        
        files = glob.glob('tmp*')
        import shutil
        for f in files:
            shutil.rmtree(f)
    
        
        t = TaskInvert()
        t.vis = lastuv
        t.map = m+".mp"
        t.beam = m+".bm"
        t.imsize =   [128,128]
    
        t.cell = 0.4
        t.line = ['vel',int(bw/cw)-6,d[:,0].min()+3*cw,cw,cw]
        t.robust = 1
    #    t.cell = 0.85955216805
    #    t.imsize = 128,128
        t.options = 'double,systemp'
        s = t.snarf()
        printout('Inverting amp ',s)
        
        imdata = loadImg(m+".mp")[0]
        std = imdata.std()
        t = TaskClean()
        t.map = m+".mp"
        t.beam = m+".bm"
        t.out = m+'.cl'
        t.niters = 100000
        t.cutoff = cutoff
        s = t.snarf()    
        printout('Cleaning amp',s)
        
        t = TaskRestore()
        t.map = m+".mp"
        t.beam = m+".bm"
        t.out = m+'.re'
        t.model = m+'.cl'
        s = t.snarf()
        printout('Restoring amp complete',s)
        
        t = TaskFits()
        t._in = m+'.re'
        t.op = 'xyout'
        t.out = 'fits/'+moleculas+'{0}.amp.fits'.format(i)
        s = t.snarf()
        printout('sfc FITS file exported',s)

        s = input('Continue?[y]')
        if (len(s) > 0 and s[0] != 'y'):
            return lastuv
    
def makecalmap(file,moleculas,source,verbose = False):
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


    
    vdataset = miriad.VisData(file)
    m = 'tmp.'+moleculas
    
    
    t = TaskUVSpec()
    t.interval = 1000
    t.options = 'avall'
    t.axis = 'v,amp'
    t.log = 'temp'
    t.vis = vdataset
    t.snarf()
    
    d = loadtxt('temp')
    
#    d = {}
#    for k in _keywords:
#        d[k] = ''
#    t.set(**d)
    t = TaskInvert()
    t.vis = vdataset
    t.map = m+".mp"
    t.beam = m+".bm"
    t.imsize =   [128,128]

    t.cell = 0.4
    cw = diff(d[:,0]).max()
    if (cw > 3) : cw = abs(diff(d[:,0])[0])
    bw = d[:,0].max()-d[:,0].min()
    t.line = ['vel',int(bw/cw)-6,d[:,0].min()+3*cw,cw,cw]
    t.robust = 1
#    t.cell = 0.85955216805
#    t.imsize = 128,128
    t.options = 'double,systemp'
    s = t.snarf()

    printout('inverting finished [{0}]'.format(t.line),s)
    imdata = loadImg(m+".mp")[0]
    std = imdata.std()
    t = TaskClean()
    t.map = m+".mp"
    t.beam = m+".bm"
    t.out = m+'.cl'
    t.niters = 100000
    t.cutoff = d[:,1].std()*2
    s = t.snarf()    

    printout('Cleaning done',s)
    t = TaskRestore()
    t.map = m+".mp"
    t.beam = m+".bm"
    t.out = m+'.re'
    t.model = m+'.cl'
    s = t.snarf()
    printout('Restoring complete',s)
    t = TaskFits()
    t._in = m+'.re'
    t.op = 'xyout'
    t.out = 'fits/'+moleculas+'.nosfc.fits'
    s = t.snarf()
    printout('1-st FITS file exported',s)
    
    
    i = 0
    for sr in source:
        lastuv = file 
        t = TaskUVAver()
        t.vis = lastuv
        lastuv = "uv/"+moleculas+'{0}.phc.uv'.format(i)
        t.out = lastuv
        s = t.snarf()
        printout('Aplying inital gains',s)
    
        t = TaskGPCopy()
        t.vis = sr
        t.out = lastuv
        t.mode = 'apply'
        s = t.snarf()
        printout('Aplying selfcalibrated weights from ',s)
        
        i+=1
    
    
    t = TaskUVSpec()
    t.interval = 1000
    t.options = 'avall'
    t.axis = 'v,amp'
    t.log = 'temp'
    t.vis = lastuv
    t.snarf()
    
    d = loadtxt('temp')

    
    files = glob.glob('tmp*')
    import shutil
    for f in files:
        shutil.rmtree(f)

    t = TaskInvert()
    t.vis = lastuv
    t.map = m+".mp"
    t.beam = m+".bm"
    t.imsize =   [128,128]

    t.cell = 0.4
    t.line = ['vel',int(bw/cw)-6,d[:,0].min()+3*cw,cw,cw]
    t.robust = 1
#    t.cell = 0.85955216805
#    t.imsize = 128,128
    t.options = 'double,systemp'
    s = t.snarf()
    printout('Inverting selfcalibrated ',s)
    
    imdata = loadImg(m+".mp")[0]
    std = imdata.std()
    t = TaskClean()
    t.map = m+".mp"
    t.beam = m+".bm"
    t.out = m+'.cl'
    t.niters = 100000
    t.cutoff = d[:,1].std()*2
    s = t.snarf()    
    printout('Cleaning selfcalibrated',s)
    

            
    t = TaskRestore()
    t.map = m+".mp"
    t.beam = m+".bm"
    t.out = m+'.re'
    t.model = m+'.cl'
    s = t.snarf()
    printout('Restoring sfc complete',s)
    
    t = TaskFits()
    t._in = m+'.re'
    t.op = 'xyout'
    t.out = "fits/"+moleculas+'{0}.amp.fits'.format(i)
    s = t.snarf()
    printout('sfc FITS file exported',s)

    
    
    
