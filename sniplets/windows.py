# -*- coding: utf-8 -*-
"""
Created on Wed Jun 11 13:57:40 2014

@author: pete
"""
def splitLines(vis,lines, intensity=-2.7):
    import mirexec
    import numpy as np
    import os
    visb=vis
    mirexec.addEnvironmentAutotools('/home/pete/miriad/exec')
    
    task = mirexec.TaskUVList()
    task.set(vis=visb) #TODO aaaa
    task.set(options='spectra')
    task.set(select='',line='',scale='',recnum='',log='')
    s = task.snarf()[0]
    s = s[7:]
    
    wn = []
    tmp = []
    for p in s:
        if(len(p)!=0):
            tmp.append([float(i) for i in p.decode().split(':')[1].split()])
        else:
            wn.append(tmp)
            tmp = []
    
    arr = np.array(wn)
    wth = arr.shape[2]
    hgt = arr.shape[1]
    arr2 = np.zeros((hgt,arr.shape[2]*arr.shape[0]))
    for i in range(arr.shape[0]):
        arr2[:,i*wth:(i+1)*wth] = arr[i,:,:]
        
    intr = np.zeros((arr2.shape[1],2))
    intr[:,0] = arr2[3]+arr2[2]*arr2[4]
    intr[:,1] = arr2[3]
    line = np.loadtxt(lines,usecols=(0,2))
    prs = line[np.where (line[:,1] > intensity)][:,0]/1000
    print('selected {} lines'.format(len(prs)))
    windows = []
    freqs = []
    fw = {}
    if intr[0,0]>intr[0,1]:
        intr = np.fliplr(intr)
    for frq in prs:
        for f1,f2 in intr:
            if  f1< frq < f2: 
                w = np.where(intr==f1)[0][0]
                wins = [w+1]
                windows.append(w+1)
                freqs.append(frq)
                if (intr[0,1]< intr[1,1]):
                    if (min(f1,f2)<frq<min(f1,f2)+abs(f2-f1)/2)  and w != 0:
                        print('left side')
                        windows.append(w)
                        wins.append(w)
                    elif f1 != intr[-1,0]:
                        print('right side')
                        windows.append(w+2)
                        wins.append(w+2)
                else :
                    if (min(f1,f2)<frq<min(f1,f2)+abs(f2-f1)/2):
                        print('left side')
                        windows.append(w+2)
                        wins.append(w+2)
                    elif f1 != intr[-1,0]   and w != 0:
                        print('right side')
                        windows.append(w)
                        wins.append(w)
                fw[frq] = wins

    
    if len(windows) == 0:
        print('nothing found!')
        return
    windows = np.unique(windows)
    freqs = np.unique(freqs)
    for frq in fw.keys(): 
        print('windows {0} at freq {1}'.format(','.join(map(str,fw[frq])),frq))
    
    newpath = r'windows' 
    if not os.path.exists(newpath):
        os.makedirs(newpath)
    
    task = mirexec.TaskUVCat()
    task.set(vis=visb,
             select='window({0})'.format(','.join(map(str,windows))), 
             out='windows/n{0}.uv'.format('_'.join(map(str, [np.where(prs == i)[0][0] for i in freqs])))
             ,stokes='',options='')
    task.snarf()

    for frq in fw.keys():
        ww = fw[frq]
                
        task = mirexec.TaskUVCat()
        task.set(vis=visb,
                 select='window({0})'.format(','.join(map(str,ww))), 
                 out='windows/f{0}.uv'.format(frq)
                 ,stokes='',options='')
        task.snarf()
        
        task = mirexec.TaskUvPutHead()
        task.set(vis='windows/f{0}.uv'.format(frq),hdvar='restfreq',type='r',varval=frq,out='windows/f{0}.corr.uv'.format(frq),length='',table='',time0='')
        task.snarf()
        
    return (intr,windows,freqs) 