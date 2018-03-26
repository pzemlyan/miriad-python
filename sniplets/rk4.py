# -*- coding: utf-8 -*-
"""
Created on Mon May  5 12:39:09 2014

@author: pete
"""
from pylab import *
from numpy.linalg.linalg import norm as leng
def RK4(f):
    return lambda t, y, dt: (
            lambda dy1: (
            lambda dy2: (
            lambda dy3: (
            lambda dy4: (dy1 + 2*dy2 + 2*dy3 + dy4)/6
            )( dt * f( t + dt  , y + dy3   ) )
	    )( dt * f( t + dt/2, y + dy2/2 ) )
	    )( dt * f( t + dt/2, y + dy1/2 ) )
	    )( dt * f( t       , y         ) )
 
def theory(t,x): return (1-x[1]*t*t)
 
from math import sqrt 
 

#dr = RK4(lambda t,x: array([x[0],-sign(x[1])*1/x[1]**2]))
dr = RK4(lambda t,x: vstack(
(
(1/(leng(x[1])**2))
*(-x[1]/leng(x[1]))
,x[0])
))
t, r, dt = 0.,array(([0.1,-1,0.1],[1,0.1,0.1])), .1
at,ax=[t],[r]
while t <= 1000:
    if abs(round(t) - t) < 1e-5:
        print("r(%2.1f)\t= %4.6f" % ( t, leng(r[1])))
    t, r = t + dt, r+dr(t,r,dt)
    at.append(t),ax.append(r)
    
