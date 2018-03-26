# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 11:08:03 2016

@author: pete
"""
def second(par,cube,nh3bg,ix,iy):
    double = 0
    red = 0
    blue = 0
    none = 0         
    badfit = 0

    def res1(p,x,y,c):
        return (y - nh32r.fitfunc(x,p))**2
    
    def res2(p,x,y,c):
        return (y - nh32r.fitfunc(x,p)-nh32b.fitfunc(x,p))**2


    par2 = [['' for z in range( 0,nh3bg.shape[1])] for y in range(0,nh3bg.shape[0])]
    for xi in ix:
        for yi in iy:
            spec = cube[:,xi,yi]
            if par[xi][yi] != '' and par[xi][yi] != 'none':
#                print('\r({0},{1}-in),       '.format(xi,yi),end="")
                if xi == 28 and yi == 26:
                    print(r'\n\n\naaaaa\n\n\n\n')
                
                if 'dv2' in par[xi][yi].params:
                    params = Parameters()
                    for k in par[xi][yi].params: 
                        params.add(k,par[xi][yi].params[k].value,min=par[xi][yi].params[k].min,max=par[xi][yi].params[k].max)
                    params['dv1'].min = params['dv1'].value-5e2
                    params['dv1'].max = params['dv1'].value+5e2
                    params['dv2'].min = params['dv2'].value-5e2
                    params['dv2'].max = params['dv2'].value+5e2
                    params['sigma1'].min = params['sigma1'].value-5e1
                    params['sigma1'].max = params['sigma1'].value+5e1
                    params['sigma2'].min = params['sigma2'].value-5e1
                    params['sigma2'].max = params['sigma2'].value+5e1
                    mini = lmfit.Minimizer(res2,params,fcn_args=(x,spec,[]))
                    par2hned = mini.minimize(method='nelder')
                    par2[xi][yi] = minimize(res2,par2hned.params,args=(x,spec,[]),)
                    double+=1
                    mak = (nh32r.fitfunc(x,par2[xi][yi].params)>2e-3)+(nh32b.fitfunc(x,par2[xi][yi].params)>2e-3)
                    st = spec[where(logical_not(mak))].std()

                    maa = sqrt(res2(par2[xi][yi].params,velx,spec,[])).max()
                    if maa > st*2.2:
                        print('\r({0},{1}-badfit),'.format(xi,yi),end="")
                        par2[xi][yi] = 'none'
                        red+=1
                    else: 
                        print('\r({0},{1}+ (d:{2},r:{3},b:{4},bf:{5})),'.format(xi,yi,double,red,blue,badfit),end="")
                else:
                    params = Parameters()
                    for k in par[xi][yi].params: 
                        params.add(k,par[xi][yi].params[k].value,min=par[xi][yi].params[k].min,max=par[xi][yi].params[k].max)
                    params['dv1'].min = params['dv1'].value-5e2
                    params['dv1'].max = params['dv1'].value+5e2
                    params['sigma1'].min = params['sigma1'].value-5e1
                    params['sigma1'].max = params['sigma1'].value+5e1
                    params['Nu1'].value = params['Nu1'].value
                    params['B1'].value = 6
                    
                    mini = lmfit.Minimizer(res1,params,fcn_args=(x,spec,[]))
                    par2hned = mini.minimize(method='nelder')
                    par2[xi][yi] = minimize(res1,par2hned.params,args=(x,spec,[]),)
                    double+=1
                    mak = (nh32r.fitfunc(x,par2[xi][yi].params)>2e-3)
                    st = spec[where(logical_not(mak))].std()
                    maa = sqrt(res1(par2[xi][yi].params,velx,spec,[])).max()
#                    if maa > st*2.2:
#                        print('\r({0},{1}-badfit),'.format(xi,yi),end="")
#                        par2[xi][yi] = 'none'
#                        red+=1
#                    else: 
                    print('\r({0},{1}+ (d:{2},r:{3},b:{4},bf:{5})),'.format(xi,yi,double,red,blue,badfit),end="")
            
    
            else:
                print('\r({0},{1}-),'.format(xi,yi),end="")
                par2[xi][yi] = 'none'
            if xi == 28 and yi == 26:
                print(r'\n\n\naaaaa\n\n\n\n')

    return par2

# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 11:08:03 2016

@author: pete
"""
