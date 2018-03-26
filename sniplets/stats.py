#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 13 18:39:32 2017

j1-j2<=j<=j1+j2

{abc}
{def}
@author: pete
"""

def s6j(a,b,c,d,e,f):
    from math import factorial as fc
    from math import sqrt
    
    trangle = lambda a,b,c : a-b<=c<=a+b
    delta = lambda a,b,c : sqrt(fc(a+b-c)*fc(a-b+c)*fc(-a+b+c)/fc(a+b+c+1))
    if (not trangle(a,b,c)) or (not trangle(c,d,e)) or (not trangle(a,e,f)) or (not trangle(b,c,f)):
        return 0
    
    summ = 0
    n = 0
    while n < a+b+c+d+e+f:
        if      (n-a-b-c) < 0 or\
                (n-c-d-e) < 0 or\
                (n-a-e-f) < 0 or\
                (n-b-d-f) < 0 or\
                (a+b+d+e-n) < 0 or\
                (a+c+d+f-n) < 0 or\
                (b+c+e+f-n) < 0: 
            n+=1
            continue
        summ+=(-1)**n*fc(n+1)/(fc(n-a-b-c)*fc(n-c-d-e)*fc(n-a-e-f)*fc(n-b-d-f)\
                              *fc(a+b+d+e-n)*fc(a+c+d+f-n)*fc(b+c+e+f-n))
        n+=1
        print(n)
    summ*=delta(a,b,c)*delta(c,d,e)*delta(a,e,f)*delta(b,d,f)
    return summ

I = 1
Ja = cat[:,0]
F1a = cat[:,1]
Fa = cat[:,2]
J = cat[:,3]
F1 = cat[:,4]
F = cat[:,5]
R = zeros(Ja.shape)
for i in range(len(Ja)):
    R[i] = (2*F1a[i]+1)*(2*F1[i]+1)*(2*Fa[i]+1)*(2*F[i]+1)/(2*I+1)/(2*I+1)
    R[i]*=(s6j(I,F1a[i],Ja[i],1,J[i],F1[i])**2)*s6j(I,Fa[i],F1a[i],1,F1[i],F[i])**2