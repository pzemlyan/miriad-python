# -*- coding: utf-8 -*-
"""
Created on Mon May 30 16:05:29 2016

@author: pete
"""


def filterPar(par, key, rang, indexes=['1','2']):
    r'''
    Filter params 2D list
    
    Parameters
    -------------
    
    par : list
        2d list constst of Parameters to filter
    key : str
        key to filter
    rang : turple
         two floats, range
    indexes : turple
         two str, indexes in Parameter keys 

    Returns
    ------------
    parameter 2d list

    '''
    fil = [['' for z in range(len(par[0]))] for y in range(len(par))]
    for xi in range(len(par)):
        for yi in range(len(par[0])):
            p = par[xi][yi]
            if p in ['','none']: continue
            
            for i in indexes:
                if key+i in p.params:
                    if rang[0] < p.params[key+i] < rang[1]:
                        if fil[xi][yi] == '':
                            fil[xi][yi] = lmfit.minimizer.MinimizerResult()
                            fil[xi][yi].params = Parameters()
                            
                        ti = '1'
                        for ti in indexes:
                            if key+ti not in fil[xi][yi].params: break
                        for kz in p.params.keys():
                            if kz[-1] == i:
                                fil[xi][yi].params[kz[:-1]+ti] = Parameter(name=kz[:-1]+ti,value=p.params[kz].value)
                    
    return fil