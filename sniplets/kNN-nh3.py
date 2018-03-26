#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 14 00:03:59 2017
"""

n_grid = 50

B = linspace(0,12,n_grid)
sigma = linspace(200,1500,n_grid)
Nu = logspace(27,28.47,n_grid)
dv = linspace(6e3,12e3,n_grid)


data = zeros((x.shape[0],n_grid,n_grid,n_grid,n_grid),dtype='float16')
y = zeros((1,n_grid,n_grid,n_grid,n_grid))

data = []
y = []
for xi in ix:
    for yi in iy:
        if not mask[xi,yi,0]:
            spec = zeros(x.shape)
            p1 = Parameters()
            p1.add('sigma1',sigma[xi,yi,0])
            p1.add('dv1',df1[xi,yi,0])
            p1.add('B1',B1[xi,yi,0])
            p1.add('Nu1',Nu1[xi,yi,0])
            spec+=n2h1.fitfunc(x,p1)
            if not mask[xi,yi,1]:
                p1.add('sigma2',sigma[xi,yi,1])
                p1.add('dv2',df1[xi,yi,1])
                p1.add('B2',B1[xi,yi,1])
                p1.add('Nu2',Nu1[xi,yi,1])
                spec+=n2h2.fitfunc(x,p1)
                y.append(df1[xi,yi])
            else: 
                y.append([df1[xi,yi,0],0])
            data.append(spec)

           

spec = zeros(x.shape)
data.append(spec)
y.append([0,0])
X = data
#data.append(spec)

            
data = array(data)
y = array(y)


for (iB,iS,iN,iV),z in ndenumerate(y[0]):
    p1 = Parameters()
    p1.add('sigma1',sigma[iS])
    p1.add('dv1',dv[iV])
    p1.add('B1',B[iB])
    p1.add('Nu1',Nu[iN])

    data[:,iB,iS,iN,iV] = nh32r.fitfunc(x,p1)
    y[0,iB,iS,iN,iV] = dv[iV]
    
    
data.shape =(data.shape[0],-1)
data = data.swapaxes(0,1)
y.shape = (1,-1)
y = y.swapaxes(0,1)

print('sample dataset is generated')


X = data
X[isnan(X)]=0
y[isnan(y)]=0

Xp = cube2[:].reshape(cube2.shape[0],-1).swapaxes(0,1)

from sklearn.neighbors import KNeighborsRegressor
estimator = KNeighborsRegressor(weights='distance')
estimator.fit(X, y)
y_pred = estimator.predict(Xp)


dvf = ma.array(y_pred).reshape(128,128,1)
dvf.mask=dvf < 5e3
figure()
subplot(121)
imshow(dvf[:,:,0],cmap='seismic',vmin=dvf.min(),vmax=dvf.max())
contour(cube.max(0),levels=linspace(0.2,40,12))
subplot(122)
imshow(dvf[:,:,1],cmap='seismic',vmin=dvf.min(),vmax=dvf.max())
contour(cube.max(0),levels=linspace(0.2,40,12))
y = y.reshape(y.shape)

B2 = linspace(0,37,n_grid)
sigma2 = linspace(100,1000,n_grid/2)    
Nu2 = logspace(0,6.3,n_grid/2)
dv2 = linspace(8.4e3,12e3,n_grid)

subplot(151)
step(x,spec,'-k')
plot(x,X[s[0,0]])
xlim(-5e3,20e3)
subplot(152)
step(x,spec,'-k')
plot(x,X[s[0,1]])
xlim(-5e3,20e3)
subplot(153)
step(x,spec,'-k')
plot(x,X[s[0,2]])
xlim(-5e3,20e3)
subplot(154)
step(x,spec,'-k')
plot(x,X[s[0,3]])
xlim(-5e3,20e3)
subplot(155)
step(x,spec,'-k')
plot(x,X[s[0,4]])
xlim(-5e3,20e3)


data = zeros((x.shape[0],n_grid,n_grid//2,n_grid//2,n_grid,n_grid,n_grid//2,n_grid//2,n_grid),dtype='float16')
y = zeros((2,n_grid,n_grid//2,n_grid//2,n_grid,n_grid,n_grid//2,n_grid//2,n_grid))

for (iB,iS,iN,iV,iB2,iS2,iN2,iV2),z in ndenumerate(y[0]):
    p1 = Parameters()
    p1.add('sigma1',sigma[iS])
    p1.add('dv1',dv[iV])
    p1.add('B1',B[iB])
    p1.add('N1',dv[iV])

    p1.add('sigma2',sigma2[iS2])
    p1.add('dv2',dv2[iV2])
    p1.add('B2',B2[iB2])
    p1.add('N2',dv2[iV2])


    data[:,iB,iS,iN,iV,iB2,iS2,iN2,iV2] = n2h1.fitfunc(x,p1)+n2h2.fitfunc(x,p1)
    y[0,iB,iS,iN,iV,iB2,iS2,iN2,iV2] = dv[iV]
    y[1,iB,iS,iN,iV,iB2,iS2,iN2,iV2] = dv[iV2]



n2hbg = cube.max(0)

p1 = Parameters()
p1.add('sigma1',500,min=2e2,  max=2e3)
p1.add('dv1',8331,   min=6076,max=13252,vary=True)
p1.add('B1',31,   min=0,  max=100)
p1.add('Nu1',8e5,min=2e3, max=7e6)

ps = [p1]

p2 = Parameters()

p2.add('sigma1',874,min=4e2,  max=10e3)
p2.add('dv1',7960,   min=6076,max=12.6e3,vary=True)
p2.add('B1',31,   min=0,  max=100)
p2.add('Nu1',5.2e5,min=2e3, max=2e6)

p2.add('sigma2',432,min=2e2,  max=2e3)
p2.add('dv2',10035,   min=8.0e3,max=12252,vary=True)
p2.add('B2',11,   min=0,  max=100)
p2.add('Nu2',8.1e5,min=2e3, max=2e6)

ps2 = [p2]


def bestFit(ps,spec):
    std = (spec-nh31(vel21,ps[0])).std()
    pr = ps[0]
    for p in ps:
        if (spec-nh31(vel21,p)).std() < std:
            std = (spec-nh31(vel21,p)).std()
            pr = p
    return pr
        
def bestFit2(ps,spec):
    std = (spec-nh31(vel21,ps[0])-nh32(vel21,ps[0])).std()
    pr = ps[0]
    for p in ps:
        if (spec-nh31(vel21,p)-nh32(vel21,p)).std() < std:
            std = (spec-nh31(vel21,p)-nh32(vel21,p)).std()
            pr = p
    return pr



def res1(p,x,y,c):
    return (y - nh31.fitfunc(x,p))

def res2(p,x,y,c):
    return (y - nh31.fitfunc(x,p)-nh32.fitfunc(x,p))
    
def res3(p,x,y,c):
    return (y - nh31.fitfunc(x,p)-nh32.fitfunc(x,p)-nh33.fitfunc(x,p))

parnh = [['' for z in range( 0,n2hbg.shape[1])] for y in range(0,n2hbg.shape[0])]
figure(1)
import time
cla()
ion()
show()
imshow(n2hbg, vmax=cube2[:,six.min():six.max(),siy.min():siy.max()].max())
xlim(siy.min(),siy.max())
ylim(six.min(),six.max())
for xi in ix:
    for yi in iy:
        if not ( parnh[xi][yi] == '' or parnh[xi][yi] == 'none'):
            if not (hasattr(parnh[xi][yi],'params') and 'B2' in parnh[xi][yi].params):
                continue
        spec = cube[:,xi,yi]
        if (n2hbg[xi,yi]>spec[where(logical_not(maksh))].std()*factor):# and parnh[xi][yi] == '':
            p = Parameters()
            r,s =  estimator.kneighbors(spec.reshape(1,-1),1)
            ind = unravel_index(s[0][0],(1,n_grid,n_grid,n_grid,n_grid))
            p.add('sigma1',sigma[ind[2]],min=2e2,  max=2e3)
            p.add('dv1',y[s[0][0]],   min=6076,max=13252,vary=True)
            p.add('B1',B[ind[1]],   min=0,  max=100)
            p.add('Nu1',Nu[ind[3]],min=2e3, max=7e6)

#            p = bestFit(ps,spec)
 #           v0 = vel21[where(spec == spec.max())][0]
 #           p['dv1'].value = v0
            mini = lmfit.Minimizer(res1,p,fcn_args=(vel21,spec,[]))


            parnh[xi][yi] = mini.minimize(method='leastsq')
            print('\r({0},{1}-knn ),'.format(xi,yi),end="")
            plot(yi,xi,'o',color=(0,0,(parnh[xi][yi].params['dv1']-6000)/6000))
            pyplot.draw()
            plt.pause(0.01)
            time.sleep(0.05)
            
#            if (spec-n2h1.fitfunc(vel21,p))[where(maksh)].std() > 3*spec[where(logical_not(mak))].std():
#                plot(yi,xi,'og')
#                ps.append(par2hned.params)
            
            
            
            mak = (nh31.fitfunc(vel21,parnh[xi][yi].params)>2e-1)

            pr1 = parnh[xi][yi] 
            st1 = (spec-nh31.fitfunc(vel21,parnh[xi][yi].params)).std()

            st = spec.std()
#            b2 = par[xi][yi].params['B2']
            if st1 > st:
                plot(yi,xi,'xw')
                print('\r({0},{1}-),'.format(xi,yi),end="")
                parnh[xi][yi] = 'none'
#                continue
            
            if (spec-nh31.fitfunc(vel21,pr1.params))[where(maksh)].std() > 1.5*spec[where(logical_not(mak))].std():
                p = bestFit2(ps2,spec)
                
                mini = lmfit.Minimizer(res2,p,fcn_args=(vel21,spec,[]))
                par2hned = mini.minimize(method='nedler')
    
    
                parnh[xi][yi] = par2hned
    
                mak = (nh31.fitfunc(vel21,parnh[xi][yi].params)>2e-3)+(nh32.fitfunc(vel21,parnh[xi][yi].params)>2e-3)

                st = (spec-nh31.fitfunc(vel21,parnh[xi][yi].params)-(nh32.fitfunc(vel21,parnh[xi][yi].params))).std()
                if (st*1.3 > st1):
                    parnh[xi][yi] = pr1 
                    plot(yi,xi,'>b')
                else:
                    plot(yi,xi,'>r')
                    if (spec-nh31.fitfunc(vel21,p)-nh32.fitfunc(vel21,p)).std() > 3*(spec-nh31.fitfunc(vel21,par2hned.params)-nh32.fitfunc(vel21,par2hned.params)).std():
                        plot(yi,xi,'xr')
                        ps2.append(par2hned.params)

                print('\r({0},{1}y),'.format(xi,yi),end="")
    #                if par[xi][yi].params['df1'].stderr > 100000 or А памятник растопленному снегу не par[xi][yi].params['df2'].stderr > 100000 :
    #                par[xi][yi] = 'none'
    #                print('\r({0},{1}x),'.format(xi,yi),end="")
#                if st*1.5<maa :
#                    plot(yi,xi,'oy')
#                    print('\r({0},{1}-badfit),'.format(xi,yi),end="")
#                    parnh[xi][yi] = 'none'
#                    badfit+=1

        else:
            plot(yi,xi,'xk')
            print('\r({0},{1}-),'.format(xi,yi),end="")
            parnh[xi][yi] = 'none'



for xi in ix:
    for yi in iy:
        if not ( parnh[xi][yi] == '' or parnh[xi][yi] == 'none'):
            if not (hasattr(parnh[xi][yi],'params') and 'B2' in parnh[xi][yi].params):
                continue
            spec = cube[:,xi,yi]
            mini = lmfit.Minimizer(res2,parnh[xi][yi].params,fcn_args=(vel21,spec,[]))
            par2hned = mini.minimize(method='leastsq')
        
        
            parnh[xi][yi] = par2hned
        
            plot(yi,xi,'xr')
           
ax = gca()
lon = ax.coords[0]
lat = ax.coords[1]

lon.set_major_formatter('dd:mm:ss.s')
lat.set_major_formatter('dd:mm:ss.s')
lon.set_axislabel(r'$\alpha$')
lat.set_axislabel(r'$\delta$')
lat.set_ticks(spacing=15. * u.arcsec,color='w')
lon.set_ticks(spacing=15. * u.arcsec,color='w')
lon.display_minor_ticks(True)
lat.display_minor_ticks(True)
lat.set_minor_frequency(15)
lon.set_minor_frequency(15)




text(6.82e3,2e3,'1',color='k',size='x-large',horizontalalignment='center',verticalalignment='center')
text(7.508e3,2e3,'4',color='k',size='x-large',horizontalalignment='center',verticalalignment='center')
text(8e3,2e3,'2',color='k',size='x-large',horizontalalignment='center',verticalalignment='center')
text(8.93e3,2e3,'3',color='k',size='x-large',horizontalalignment='center',verticalalignment='center')
text(8.6e3,1.2e3,'3.1',color='k',size='x-large',horizontalalignment='center',verticalalignment='center')
text(10.29e3,2e3,'5',color='k',size='x-large',horizontalalignment='center',verticalalignment='center')
text(9.766e3,2e3,'6',color='k',size='x-large',horizontalalignment='center',verticalalignment='center')


plot([6.82e3,6.82e3],[1.5e3,1.7e2],'-k')
plot([7.508e3,7.508e3],[1.5e3,1.5e2],'-k')
plot([8.93e3,8.93e3],[1.5e3,1.5e2],'-k')
plot([8.6e3,8.6e3],[0.9e3,1.4e2],'-k')
plot([10.29e3,10.29e3],[1.5e3,.8e2],'-k')
plot([9.766e3,9.766e3],[1.5e3,9],'-k')
text(7.54e3,2e3,'4',color='k',size='x-large',horizontalalignment='center',verticalalignment='center')
text(8e3,2e3,'2',color='k',size='x-large',horizontalalignment='center',verticalalignment='center')
text(8.9e3,2e3,'3',color='k',size='x-large',horizontalalignment='center',verticalalignment='center')
text(8.6e3,1.2e3,'3.1',color='k',size='x-large',horizontalalignment='center',verticalalignment='center')
text(10.2e3,2e3,'5',color='k',size='x-large',horizontalalignment='center',verticalalignment='center')
text(9.7e3,2e3,'6',color='k',size='x-large',horizontalalignment='center',verticalalignment='center')
