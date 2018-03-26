# -*- coding: utf-8 -*-
"""
Created on Wed Mar  2 11:39:21 2016

@author: pete
"""

nu1 = nh31.restfreq
nu2 = nh32r.restfreq

Tex = lambda B,nu: h*nu/k/log(h*nu/B/k+1)
Ftex = lambda T,nu: 1/(exp(h*nu/k/T)-1)
Nu = lambda Tex,dv,tau,nu: 1.6e13*Ftex(Tex,nu)*dv*tau
N = lambda Nu, Tex: Nu*(1+exp(h*nu/k/Tex))

B11 = p1['B1']
B12 = p1['B2']

B21 = p2['B1']
B22 = p2['B2']

N11 = N(Nu(Tex(B11),1.3,nh31.tau(velx,p1.params).max()),Tex(B11))
N12 = N(Nu(Tex(B12),1.3,nh32.tau(velx,p1.params).max()),Tex(B12))

N21 = N(Nu(Tex(B21),1.3,nh32r.tau(velx,p2.params).max()),Tex(B21))
N22 = N(Nu(Tex(B22),1.3,nh32b.tau(velx,p2.params).max()),Tex(B22))

Tkinr = -41.5*k/k/log()