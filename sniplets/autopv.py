#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 18:50:56 2017

@author: pete
"""
import glob

# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 15:02:05 2014

@author: pete
"""
from astropy import units as u;
from astropy.io.fits.hdu.image import PrimaryHDU
from astropy.wcs import WCS
from scipy.constants import c,h,k,pi
import numpy as np
from lmfit import minimize, Parameters, Parameter, report_fit, report_errors
from math import log, atan, sin
from matplotlib.widgets import Button
from enum import Enum
from astropy.io import fits
import pickle

from matplotlib.gridspec import GridSpec
from matplotlib.pyplot import gcf, draw_if_interactive, delaxes

files = glob.glob('*.fits')


#1-2 pts:
pts = load('../../estimates/pv/n2hp/path.npy')

figure(3)
subplot(131)
makeAutoPv(files,pts,f1)
plotPV(open('../../ch3oh_masers.cat'),array(p1),array(p2),recs,'o','y',5)
plotPV(open('../../water_maser.cat'),array(p1),array(p2),recs,'s','y',5)
xlim(-37000, 40000)
ylim(-9,12)
ax = gca()
ax.annotate("SMA1-NE", xy=(7227, -0.3), xytext=(0, -5),arrowprops=dict(arrowstyle="->",color='g'),color='g')
ax.annotate("SMA1-SW", xy=(10217, 2.6), xytext=(13000, 7),arrowprops=dict(arrowstyle="->",color='g'),color='g')
figure(3)
subplot(132)
makeAutoPv(files,pts2,f1)
plotPV(open('../../ch3oh_masers.cat'),array(p1),array(p1+n2),recs,'o','y',5)
plotPV(open('../../water_maser.cat'),array(p1),array(p1+n2),recs,'s','y',5)
xlim(-20000, 40000)
ylim(-7,7)
figure(3)
subplot(133)
makeAutoPv(files,pts3,f1)
plotPV(open('../../ch3oh_masers.cat'),array(p1),array(p2),recsf,'o','y')
plotPV(open('../../water_maser.cat'),array(p1),array(p2),recsf,'s','y')
xlim(-46000, 68000)
ylim(-15,10)


os.makedirs('png/pv/ch3cn2',exist_ok=True)  
for ang in linspace(0,pi,50):    
#        print('\r running '+f+" at {0:0.1f}".format(ang/2/pi*360),end='')
    pts = [p1,(p1[0]+20*cos(ang),p1[1]+20*sin(ang))]
    subplot(111)
    makeAutoPv(files,pts,f1,subplots=True)
#        xlim(-15000, 21000)
#    ylim(-20, 20)

#    np.savetxt(f[:-5]+'pvd.np',pv)
#    pvds.append(pv)
    xlim(0e3,25e3)
    ylim(-10,10)
    subplot(121)
    contour(nh3bg,levels=linspace(3,6,10),extent=getExtent(recsch3,recsf))
    imshow(fbgch3,vmax=3,vmin=0.7)
    ylim(20,45)
    xlim(20,45)
    savefig('png/pv/ch3cn2/anim.{0}.pvd.png'.format(ang))
    print(" at {0:0.1f}".format(ang/2/pi*360),end='')
os.popen(r'convert -delay 20 -loop 0 png/pv/ch3cn2/anim*.png png/pv/ch3cn2.gif').read()        
os.popen(r'notify-send "anim done"').read()        


def makeAutoPv(files,pts,f1,subplots=False):
    cmsps = iter(['kk','r','g','b','m'])
    lvs = iter([2,2.2,1.3,2.3,1.5])
    lst = iter(['solid','dashed','solid','dashed','solid'])
    cnt = iter([7,5,6,5,5,4])
    
    figure(3,figsize=(15,7))
    cla()
    pvds = []
    
    fl = True
    recsf = fits.open(f1)
    
    for f in files:
        print(f)
        recs = fits.open(f)
        cube = recs[0].data[0]
        spec = cube[:,cube.shape[1]//2,cube.shape[2]//2]
        maks = spec>spec.mean()+spec.std()
        fbg = cube[where(maks),:,:][0].max(0)
        pvd = Positionerovelociter(recs,cube)
        fgr = None
        if subplots:
            fgr = figure(3,figsize=(50,50))
            subplot(121)
        else:
            fgr = figure(4,figsize=(10, 10))
        cla()
        pvd.loadImg(fbg)
        
        pvd.target(fgr)
        fx,fy = getExtent2(recsf,recs)
        pvd.points = []
        for p in pts:
            pvd.points.append([fx(p[0]),fy(p[1])])
        pv = pvd.pvd(1.5,1,15)
        
        pvd.plot()
        title(f )
    #    savefig('png/mp/1.'+f[:-5]+'.png')
        figure(3,figsize=(50, 50))
        if subplots:
            subplot(122)
    #    cla()
        if not fl:
            levels = logspace(log10(pv.std()*next(lvs)),log10(pv.max()),next(cnt))
            cm = next(cmsps)
            if cm == 'kk':
                cs = contourf(pv,extent=pvd.extent(recs[0]),cmap='cool',levels=levels)
            else:
                cs = contour(pv,extent=pvd.extent(recs[0]),colors=cm,levels=levels,linestyles=next(lst))
        
    #        clabel(cs, levels[1::2],  # label every second level
    #           inline=1,
    #           fmt='C$^{18}O$',
    #           fmt=f[:f.index('.')],
    #           fontsize=12)
        if fl:
            imshow(pv,extent=pvd.extent(recs[0]),vmin=0.2,vmax=pv.max(),cmap='gray')
            cb = colorbar()
            cb.set_label('F, Ян/луч')
            fl=not(fl)
        ax = gca()
        ax.set_aspect('auto')
        xlabel('v, м/c')
        ylabel('Положение, "')
    #    savefig('png/ppv/'+f[:-5]+'/anim.{0}.pvd.png'.format(ang))
    #    title(f )
    #    xlim(-4000, 20000)
    #    ylim(-12, 12)
    #    np.savetxt(f[:-5]+'pvd.np',pv)
    #    savefig('png/pv/1.'+f[:-5]+'.pvd.eps')
    #    pvds.append(pv)
    #    savefig('png/pv/1.'+f[:-5]+'.pvd.png')
        print(f)


def makeAutoPolyPv(files,pts,f1):
    cmsps = iter(['r','b','g','kk','m'])
    lvs = iter([1,1.7,1.5,2.3,1.5])
    lst = iter(['solid','dashed','solid','dashed','solid'])
    cnt = iter([7,5,6,5,5,4])
    
    figure(3)
    cla()
    pvds = []
    
    fl = True
    recsf = fits.open(f1)
    
    for f in files:
        print(f)
        recs = fits.open(f)
        cube = recs[0].data[0]
        spec = cube[:,cube.shape[1]//2,cube.shape[2]//2]
        maks = spec>spec.mean()+spec.std()
        fbg = cube[where(maks),:,:][0].max(0)
        pvd = PolyPositionerovelociter(cube)
        fgr = figure(4,figsize=(10, 10))
        cla()
        pvd.loadImg(fbg)
        
        pvd.target(fgr)
        
        
        fx,fy = getExtent2(recsf,recs)
        pvd.points = []
        for p in pts:
            pvd.points.append([fx(p[0]),fy(p[1])])
        pv = pvd.pvd(4,1)
        
        pvd.plot(recs[0],20)
        title(f )
    #    savefig('png/mp/1.'+f[:-5]+'.png')
        figure(3,figsize=(10, 10))
    #    cla()
        if not fl:
            levels = logspace(log10(pv.std()*next(lvs)),log10(pv.max()*0.9),next(cnt))
            print('{1} levels is {0}'.format(levels,f))
            cm = next(cmsps)
            if cm == 'kk':
                cs = contourf(pv,extent=pvd.extent(recs[0]),cmap='cool',levels=levels)
            else:
                cs = contour(pv,extent=pvd.extent(recs[0]),colors=cm,levels=levels,linestyles=next(lst))
        
    #        clabel(cs, levels[1::2],  # label every second level
    #           inline=1,
    #           fmt='C$^{18}O$',
    #           fmt=f[:f.index('.')],
    #           fontsize=12)
        if fl:
            imshow(pv,extent=pvd.extent(recs[0]),vmin=0,vmax=pv.max(),cmap='gist_stern')
            cb = colorbar()
            cb.set_label('F, Ян/луч')
            fl=not(fl)
        ax = gca()
        ax.set_aspect('auto')
        xlabel('v, м/c')
        ylabel('Положение, "')
    #    savefig('png/ppv/'+f[:-5]+'/anim.{0}.pvd.png'.format(ang))
    #    title(f )
    #    xlim(-4000, 20000)
    #    ylim(-12, 12)
    #    np.savetxt(f[:-5]+'pvd.np',pv)
    #    savefig('png/pv/1.'+f[:-5]+'.pvd.eps')
    #    pvds.append(pv)
    #    savefig('png/pv/1.'+f[:-5]+'.pvd.png')
        print(f)



pts1 = [(63.060228131656693, 36.801449872878422),
 (52.153885796742934, 43.134164777021908),
 (51.450250807393651, 57.030955816670087),
 (51.450250807393651, 70.927746856318265),
 (55.672060743489304, 82.713632927918624),
 (62.884319384319376, 86.75953411667696),
 (71.503848003847992, 84.824537895966458),
 (76.429292929292927, 78.140005497148337),
 (81.0029203600632, 68.640933140933129),
 (80.827011612725869, 63.187761973476242),
 (74.142479213907762, 57.734590806019369),
 (67.985673057101621, 49.994605923177332),
 (67.80976430976429, 38.736446093588938),
 (58.486600700886406, 34.162818662818651)]

pts= [(58.017419088847646, 17.43592386449528),
 (45.215934858791996, 39.719989005703276),
 (44.267676767676761, 44.461279461279446),
 (54.698515769944322, 55.208204493918764),
 (49.008967223252924, 81.127258984401834),
 (49.167010238438792, 92.980485123342248),
 (52.011784511784498, 98.670033670033661),
 (56.436988936988925, 100.72459286744999),
 (71.293032364460913, 88.081151652580218),
 (88.993850065278608, 83.023775166632305),
 (97.686215900501594, 71.802721088435362),
 (104.32402253830823, 59.475365903937323),
 (104.64010856867998, 39.719989005703276),
 (101.16316223459079, 37.507386793101062),
 (91.680581323438446, 49.834741977599109),
 (86.623204837490533, 49.518655947227359),
 (83.462344533773077, 40.036075036075026),
 (88.361678004535122, 26.918504775647619),
 (79.511269154126282, 18.384181955610515),
 (67.974129045557603, 20.438741153026854),
 (58.175462104033514, 17.277880849309408)]

figure(3)
subplot(111)
makeAutoPolyPv(files,pts,f1)
plotPPV(open('../../water_maser.cat'),pts,recs,'s','y',1)
plotPPV(open('../../ch3oh_masers.cat'),pts,recs,'o','y',1)
plotPPV(open('../../clumps.cat'),pts,recs,(5, 1),'w',1,draw_text=True)
plotPPV(open('../../metanol_probes.cat'),pts,recs,'o','g',1,draw_text=True)
xlim(0,15e3)

text(7.25e3,(137+160)/2,'1',color='k',size='x-large',horizontalalignment='center',verticalalignment='center')
text(7.25e3,7.5,'1',color='k',size='x-large',horizontalalignment='center',verticalalignment='center')
text(7.8e3,(10+45)/2,'2',color='k',size='x-large',horizontalalignment='center',verticalalignment='center')
text(8.3e3,(45+80)/2,'3',color='k',size='x-large',horizontalalignment='center',verticalalignment='center')
text(9.1e3,(70+105)/2,'4',color='k',size='x-large',horizontalalignment='center',verticalalignment='center')
text(7.8e3,(100+140)/2,'2',color='k',size='x-large',horizontalalignment='center',verticalalignment='center')
text(10e3,(109+132)/2,'5',color='k',size='x-large',horizontalalignment='center',verticalalignment='center')
#text(9e3,1.3e3,'4',color='k',size='x-large',horizontalalignment='center',verticalalignment='center')
#text(8.6e3,1.2e3,'3.1',color='k',size='x-large',horizontalalignment='center',verticalalignment='center')
text(10.3e3,(59+71)/2,'5',color='k',size='x-large',horizontalalignment='center',verticalalignment='center')
text(9.1e3,(7+23)/2,'4',color='k',size='x-large',horizontalalignment='center',verticalalignment='center')

plot([7.2e3,7.2e3],[137,160],'-.g',linewidth=2,zorder=2)
plot([7.2e3,7.2e3],[0,15],'-.g',linewidth=2,zorder=2)
plot([7.8e3,7.8e3],[10,45],'-.g',linewidth=2,zorder=2)
plot([8.3e3,8.3e3],[45,80],'-.g',linewidth=2,zorder=2)
plot([9.1e3,9.1e3],[70,105],'-.g',linewidth=2,zorder=2)
plot([7.8e3,7.8e3],[100,140],'-.g',linewidth=2,zorder=2)
plot([10.0e3,10.0e3],[109,132],'-.g',linewidth=2,zorder=2)
plot([10.29e3,10.29e3],[59,71],'-.g',linewidth=2,zorder=2)
plot([9.1e3,9.1e3],[7,23],'-.g',linewidth=2,zorder=2)
#plot([10.29e3,10.29e3],[95,115],'-.g',linewidth=1,zorder=2)
#plot([9.76e3,9.76e3],[120,135],'-.g',linewidth=1)
plot([9.76e3,9.76e3],[5,25],'-.g',linewidth=1,zorder=2)

plot([7.508e3,7.508e3],[110,140],'-.g',linewidth=2,zorder=2)
plot([7.508e3,7.508e3],[0,40],'-.g',linewidth=2,zorder=2)
plot([8.3e3,8.3e3],[20,60],'-.g',linewidth=2,zorder=2)
plot([9e3,9e3],[55,110],'-.g',linewidth=2,zorder=2)
plot([9.5e3,9.5e3],[93,130],'-.g',linewidth=2,zorder=2)
plot([9.5e3,9.5e3],[5,25],'-.g',linewidth=2,zorder=2)


clr = 'k'

text(7.708e3,122.5,'1',color=clr,size='x-large',horizontalalignment='center',verticalalignment='center',zorder=2)
text(7.708e3,15,'1',color=clr,size='x-large',horizontalalignment='center',verticalalignment='center',zorder=2)
text(8.3e3,40,'2',color=clr,size='x-large',horizontalalignment='center',verticalalignment='center',zorder=2)
text(9e3,77.5,'3',color=clr,size='x-large',horizontalalignment='center',verticalalignment='center',zorder=2)
text(9.5e3,(93+130)/2,'4',color=clr,size='x-large',horizontalalignment='center',verticalalignment='center',zorder=2)
text(9.5e3,(5+25)/2,'4',color=clr,size='x-large',horizontalalignment='center',verticalalignment='center',zorder=2)
#text(9.766e3,15,'6',color=clr,size='x-large',horizontalalignment='center',verticalalignment='center',zorder=2)
#text(9.766e3,127,'6',color=clr,size='x-large',horizontalalignment='center',verticalalignment='center')
ylim(0,135)

text(6.82e3,10,'1',color='r',size='x-large',horizontalalignment='center',verticalalignment='center')
text(6.82e3,20,'1',color='r',size='x-large',horizontalalignment='center',verticalalignment='center')
text(8e3,30,'2',color='r',size='x-large',horizontalalignment='center',verticalalignment='center')
text(8e3,40,'2',color='r',size='x-large',horizontalalignment='center',verticalalignment='center')
text(8.93e3,50,'3.1',color='r',size='x-large',horizontalalignment='center',verticalalignment='center')
text(8.93e3,60,'3',color='r',size='x-large',horizontalalignment='center',verticalalignment='center')
text(8.93e3,10,'3',color='r',size='x-large',horizontalalignment='center',verticalalignment='center')
#text(7.508e3,112,'4',color='k',size='x-large',horizontalalignment='center',verticalalignment='center')
#text(7.508e3,112,'4',color='r',size='x-large',horizontalalignment='center',verticalalignment='center')
text(7.508e3,20,'4',color='r',size='x-large',horizontalalignment='center',verticalalignment='center')
text(10.29e3,30,'5',color='r',size='x-large',horizontalalignment='center',verticalalignment='center')
#text(9.766e3,20,'6',color='k',size='x-large',horizontalalignment='center',verticalalignment='center')
text(9.766e3,40,'6',color='r',size='x-large',horizontalalignment='center',verticalalignment='center')

text(6.82e3,10,'1',color='r',size='x-large',horizontalalignment='center',verticalalignment='center')
text(6.82e3,20,'1',color='r',size='x-large',horizontalalignment='center',verticalalignment='center')


scatter(8e3,12,c='y',marker=(5, 1),zorder=3)
scatter(8.7e3,21,c='y',marker=(5, 1),zorder=3)
scatter(9e3,5,c='y',marker=(5, 1),zorder=3)
text(8e3,11.5,'SMA1',color='k',size='x-large',horizontalalignment='center',verticalalignment='top')
text(7.9e3,21,'SMA2',color='k',size='x-large',horizontalalignment='right',verticalalignment='center')
text(10e3,5,'SMA3',color='k',size='x-large',horizontalalignment='left',verticalalignment='center')




text(7.2e3,1.3e3,'1',color='k',size='x-large',horizontalalignment='center',verticalalignment='center')
text(7.8e3,1.3e3,'2',color='k',size='x-large',horizontalalignment='center',verticalalignment='center')
text(8.3e3,1.3e3,'3',color='k',size='x-large',horizontalalignment='center',verticalalignment='center')
text(9e3,1.3e3,'4',color='k',size='x-large',horizontalalignment='center',verticalalignment='center')
text(9.5e3,1.3e3,'5',color='k',size='x-large',horizontalalignment='center',verticalalignment='center')
plot([7.2e3,7.2e3],[123,140],'-.g',linewidth=2,zorder=2)
plot([7.2e3,7.2e3],[0,18],'-.g',linewidth=2,zorder=2)
plot([7.8e3,8.3e3],[18,45],'-.g',linewidth=2,zorder=2)
plot([8.3e3,8.3e3],[45,55],'-.g',linewidth=2,zorder=2)
plot([9e3,9e3],[55,75],'-.g',linewidth=2,zorder=2)
plot([9.5e3,9.5e3],[93,130],'-.g',linewidth=2,zorder=2)
plot([9.5e3,9.5e3],[5,25],'-.g',linewidth=2,zorder=2)





ax = gca()
lon = ax.coords[0]
lat = ax.coords[1]
lon.set_major_formatter('dd:mm:ss.s')
lon.set_major_formatter('hh:mm:ss.s')
lon.set_axislabel(r'$\alpha$')
lat.set_axislabel(r'$\delta$')
lon.set_ticks(spacing=4 * u.arcsec, color='white', exclude_overlapping=True)
lat.set_ticks(spacing=4 * u.arcsec, color='white', exclude_overlapping=True)
lon.display_minor_ticks(True)
lat.display_minor_ticks(True)
lat.set_minor_frequency(5)
lon.set_minor_frequency(5)
cla()
imshow(fbg,vmax=100,vmin=0,cmap='gray')
pvd.points
pvd.points = pvd.points[:-1]
pvd.armed = False
pvd.plot(recs[0],10)
#contour(c18bg,levels=    ,extent=getExtent(recsn,recsc),linestyles=['dotted','dotted','solid','solid','solid','solid'],colors='#20a0ff')
contour(nh3bg,levels=[5],extent=getExtent(recsn,recsnh),colors='r',linestyles='dashed')

figure(2)
cla()
fl = True
recsf = fits.open(f1)
for f in files:
    recs = fits.open(f)
    cube = recs[0].data[0]
    spec = cube[:,cube.shape[1]//2,cube.shape[2]//2]
    maks = spec>spec.mean()+spec.std()
    maks = spec>spec.mean()+spec.std()
    fbg = cube[where(maks),:,:][0].max(0)
    pvd = Positionerovelociter(cube)
    fgr = figure(1,figsize=(10, 10))
    cla()
    pvd.loadImg(fbg)
    
    pvd.target(fgr)
    
    
    fx,fy = getExtent2(recsf,recs)
    pvd.points = []
    for p in pts:
        pvd.points.append([fx(p[0]),fy(p[1])])
    pv = pvd.pvd(3,1,20)
    
    pvd.plot()
    title(f )
#    savefig('png/mp/1.'+f[:-5]+'.png')
    figure(3,figsize=(10, 10))
    subplot(222)
#    cla()
    if not fl:
        levels = logspace(log(0.2,10),log(0.4,10),5)
        cs = contour(pv,extent=pvd.extent(recs[0]),colors='#50ffa0',levels=levels)
    
#        clabel(cs, levels[1::2],  # label every second level #50ffa0
#           inline=1,
#           fmt='C$^{18}O$',
#           fmt=f[:f.index('.')],
#           fontsize=12)
    if fl:
        imshow(pv,extent=pvd.extent(recs[0]),vmin=pv.std()/2,cmap='gray')
        fl=not(fl)
    ax = gca()
    ax.set_aspect('auto')
    ylabel('v, м/c')
    ylabel('Положение, "')
    
     
for f in files:
    recs = fits.open(f)
    cube = recs[0].data[0]
    ps = where(cube.mean(0) == cube.mean(0).max())
    spec = cube[:,ps[0][0],ps[1][0]]
    maks = spec>spec.mean()+2*spec.std()
    nh3bg = cube[where(maks),:,:][0].max(0)
    tvel = zeros(cube.shape)
    from astropy import wcs
    w2 = wcs.WCS(recs[0].header)
    spec2 = spec
    chan2 = arange(len(spec2))+1
    fr = zeros((chan2.shape[0],4))
    fr[:,2] = chan2
    fr = w2.all_pix2world(fr,1)
    vel21 = fr[:,2]
    for (xi,yi),val in np.ndenumerate(cube[0]):
        tvel[:,xi,yi] = vel21
    cc = cube[maks]
    vd = (cc*tvel[maks]).sum(0)/cc.sum(0)

import os
from subprocess import call


for f in files:
    recs = fits.open(f)
    cube = recs[0].data[0]
    spec = cube[:,63,63]
    maks = spec>spec.mean()+spec.std()
    nh3bg = cube[where(maks),:,:][0].max(0)
    os.makedirs('png/pv/'+f[:-5],exist_ok=True)  
    for ang in linspace(0,pi,50):    
#        print('\r running '+f+" at {0:0.1f}".format(ang/2/pi*360),end='')
        pts = [(cube.shape[1]/2,cube.shape[2]/2),(cube.shape[1]/2+50*cos(ang),cube.shape[2]/2+50*sin(ang))]
        pvd = Positionerovelociter(cube)
        fgr = figure(11,figsize=(20, 10))
        
        pvd.loadImg(nh3bg)
        subplot(121)
        cla()
        pvd.target(fgr)
        pvd.points = pts
        pv = pvd.pvd(3,1,20)
        
        pvd.plot()
    #    ylim(0,128)
        title(f )
#        savefig('png/mp/1.'+f[:-5]+'.png')
        subplot(122)
        cla()
        imshow(pv,extent=pvd.extent(recs[0]),vmin=pv.std()/2)
        ax = gca()
        ax.set_aspect('auto')
        contour(pv,extent=pvd.extent(recs[0]),cmap='hot',levels=linspace(0.05,pv.max(),20))
        xlim(0,14000)
        ylim(-12,12)
        title(f+" at {0:0.1f}".format(ang/2/pi*360))
#        xlim(-15000, 21000)
    #    ylim(-20, 20)
    
    #    np.savetxt(f[:-5]+'pvd.np',pv)
    #    pvds.append(pv)
        
        savefig('png/pv/'+f[:-5]+'/anim.{0}.pvd.png'.format(ang))
        print('\r'+f+" at {0:0.1f}".format(ang/2/pi*360),end='')
    os.popen(r'convert -resize 960x476 -delay 20 -loop 0 png/pv/'+f[:-5]+'/anim*.png png/pv/'+f[:-5]+'.gif').read()        
    print(f+' done')


files = ['']

n1 = array((p2[0]-p1[0],p2[1]-p1[1]))
n1 = n1/norm(n1)
n2 = array((n1[1],-n1[0]))
point = array(p1)
po = []
for i in range(-number,number):
    co = n1*i*step+n1*step/2
    for (xi,yi),val in ndenumerate(nh3bg):
        vec = array((yi,xi))-point-co
        if abs(vdot(vec,n1)) < step and abs(vdot(vec,n2)) < width :
            if not inte.mask[xi,yi,0]:
                if par[xi][yi] != '' and par[xi][yi] != 'none':
                    p = par[xi][yi].params
                    po.append([p['dv1'].value,-i*step*recsnh[0].header['cdelt1']*u.deg.to(u.arcsec)])
    


for f in files:
    recs = fits.open(f)
    cube = recs[0].data[0]
    spec = cube[:,63,63]
    maks = spec>spec.mean()+spec.std()
    nh3bg = cube[where(maks),:,:][0].max(0)
    os.makedirs('png/ppv/'+f[:-5],exist_ok=True)  
    pts = [(cube.shape[1]/2,cube.shape[2]/2),(cube.shape[1]/2+50*cos(ang),cube.shape[2]/2+50*sin(ang))]
    pvd = Positionerovelociter(cube)
    fgr = figure(1,figsize=(20, 10))
    
    pvd.loadImg(nh3bg)
    subplot(121)
    cla()
    pvd.target(fgr)
    pvd.points = pts
    pv = pvd.pvd(4,1,30)
    
    pvd.plot()
#    ylim(0,128)
    title(f)
#        savefig('png/mp/1.'+f[:-5]+'.png')
    subplot(122)
    cla()
    imshow(pv,extent=pvd.extent(recs[0]),vmin=pv.std()/2)
    ax = gca()
    ax.set_aspect('auto')
    contour(pv,extent=pvd.extent(recs[0]),cmap='hot',levels=linspace(0.05,pv.max(),20))
    xlim(0,14000)
    ylim(-12,12)
    title(f+" at {0:0.1f}".format(ang/2/pi*360))
#        xlim(-15000, 21000)
    #    ylim(-20, 20)
    
    #    np.savetxt(f[:-5]+'pvd.np',pv)
    #    pvds.append(pv)
        
    savefig('png/pv/'+f[:-5]+'/anim.{0}.pvd.png'.format(ang))
    os.popen(r'convert -resize 960x476 -delay 20 -loop 0 png/pv/'+f[:-5]+'/anim*.png png/pv/'+f[:-5]+'.gif').read()        
    print(f+' done')
