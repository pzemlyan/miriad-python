#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 14 00:03:59 2017
"""

n_grid = 20

B = uniform(2,37,n_grid)
sigma = uniform(600,2000,n_grid)
Nu = logspace(0,6.3,n_grid)
dv = uniform(5e3,11e3,5*n_grid)


    y[0,i] = dv[i]
    y[1,i] = B[i]
    y[2,i] = Nu[i]
    y[3,i] = sigma1[i]


data = zeros((x.shape[0],n_grid,n_grid,n_grid,5*n_grid),dtype='float16')
y = zeros((1,n_grid,n_grid,n_grid,5*n_grid))

chisqr = zeros(n2hbg.shape)

for (ix,iy),v in ndenumerate(chisqr):
    if not ( parnh[ix][iy] == '' or parnh[ix][iy] == 'none'):
        chisqr[ix,iy] = par[ix][iy].residual.std()

data = []
y = []
for (xi,yi),v in ndenumerate(n2hbg):
    if chisqr[xi,yi] == 0 or chisqr[xi,yi] > 2.5 : continue
    if not mask[xi,yi,0]:
        if not mask[xi,yi,1]:
            y.append(array([
                    df1[xi,yi,0],
                    df1[xi,yi,1],
                    ]))
        else: 
            y.append(array([
                    df1[xi,yi,0],
                    -1e7,
                    ]))
        data.append(cube[:,xi,yi])

           

for i in range(10):
    spec = uniform(low=-1,high=1,size=x.shape)*0.04
    data.append(spec)
    y.append([-1e7,-1e7])

data = array(data)
y = array(y)
X = data

#data.append(spec)

cat = loadtxt('../../../cats/n2h.cat')
n2hr = Analyse.N2HModel(cat,'1')           
n2hb = Analyse.N2HModel(cat,'2')           


n_grid = 10000

B = uniform(0.5,45,n_grid)
sigma1 = uniform(200,2000,n_grid)
Nu = uniform(1e3,1e6,n_grid)
dv = uniform(6.5e3,12e3,n_grid)

data = zeros((x.shape[0],n_grid+7),dtype='float16')
y = zeros((8,n_grid+7))

for i in range(n_grid):
    p1 = Parameters()
    p1.add('sigma1',sigma1[i])
    p1.add('dv1',dv[i])
    p1.add('B1',B[i])
    p1.add('N1',Nu[i])
    data[:,i] = n2hr.fitfunc(x,p1)

    y[0,i] = dv[i]
    y[1,i] = B[i]
    y[2,i] = Nu[i]
    y[3,i] = sigma1[i]
    if i%1000==0:
        print(i*100/n_grid)

    

B2 = uniform(3,20,n_grid)
sigma2 = uniform(200,600,n_grid)
Nu2 = uniform(1e20,4e22,n_grid)
ddv = uniform(0,4e3,n_grid)
data = zeros((x.shape[0],n_grid+7),dtype='float16')
y = zeros((8,n_grid+7))

for i in range(n_grid):
    p1 = Parameters()
    p1.add('sigma1',sigma1[i])
    p1.add('dv1',dv[i])
    p1.add('B1',B[i])
    p1.add('Nu1',Nu[i])

    p1.add('sigma2',sigma2[i])
    p1.add('dv2',ddv[i]+dv[i])
    p1.add('B2',B2[i])
    p1.add('Nu2',Nu2[i])

    data[:,i] = nh31.fitfunc(x,p1)+nh32.fitfunc(x,p1)

    y[0,i] = dv[i]
    y[1,i] = B[i]
    y[2,i] = Nu[i]
    y[3,i] = sigma1[i]
    y[4,i] = ddv[i]+dv[i]
    y[5,i] = B2[i]
    y[6,i] = Nu2[i]
    y[7,i] = sigma2[i]
    if i%1000==0:
        print(i*100/n_grid)



for i in range(7):
    data[:,n_grid+i] = zeros(len(spec))

    y[0,n_grid+i] = -1e7
    y[1,n_grid+i] = -1
    y[2,n_grid+i] = -5
    y[3,n_grid+i] = -100
    y[4,n_grid+i] = -1e7
    y[5,n_grid+i] = -1
    y[6,n_grid+i] = -5
    y[7,n_grid+i] = -100
            

for (iB,iS,iN,iV),z in ndenumerate(y[0]):
    p1 = Parameters()
    p1.add('sigma1',sigma[iS])
    p1.add('dv1',dv[iV])
    p1.add('B1',B[iB])
    p1.add('N1',Nu[iN])

    data[:,iB,iS,iN,iV] = n2hr.fitfunc(x,p1)
    y[0,iB,iS,iN,iV] = dv[iV]

data = []
tk = []
for (xi,yi),val in ndenumerate(nh3bg):
    if not tkinimg.mask[xi,yi,0]:
        data.append(np.append(cube[:,xi,yi],cube2[:,xi,yi]))
        tk.append(tkinimg[xi,yi,0])
    
    
data = data.swapaxes(0,1)
y = y.swapaxes(0,1)



print('sample dataset is generated')


X = data
X[isnan(X)]=0

Xp = cube[:].reshape(cube.shape[0],-1).swapaxes(0,1)
#Xpp = concatenate((Xp,cube[:].reshape(cube.shape[0],-1).swapaxes(0,1)),1)
from sklearn.neighbors import KNeighborsRegressor
from sklearn.linear_model import RANSACRegressor, TheilSenRegressor
from sklearn.neural_network import MLPRegressor
estimator = MLPRegressor()

estimator = KNeighborsRegressor(n_neighbors=3,weights='distance')
estimator.fit(X, y)
y_pred = estimator.predict(Xp)
dvf = ma.array(y_pred).reshape(256,256,2)

p = Parameters()
r =  estimator.predict(spec.reshape(1,-1))[0]
p.add('dv1',r[0],   min=6076,max=13252,vary=True)
p.add('B1',r[1],   min=0,  max=100)
p.add('N1',r[2],min=2e3, max=7e6)
p.add('sigma1',r[3],min=2e2,  max=2e3)

p.add('dv2',r[4],   min=6076,max=13252,vary=True)
p.add('B2',r[5],   min=0,  max=100)
p.add('N2',r[6],min=2e3, max=7e6)
p.add('sigma2',r[7],min=2e2,  max=2e3)




dvf = ma.array(y_pred).reshape(128,128,2)
dvf.mask=dvf<6e3
dvf.mask[:,:,0]|=cube.max(0) < 5
dvf.mask[:,:,1]|=cube.max(0) < 5
dvf.mask[:,:,2]=cube.max(0) < 0.5
dvf.mask[:,:,3]=cube.max(0) < 0.5
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
p1.add('N1',8,min=1e4, max=1e6)

ps = [p1]

p2 = Parameters()

p2.add('sigma1',874,min=6e2,  max=1e3)
p2.add('dv1',7960,   min=6076,max=9.5e3,vary=True)
p2.add('B1',37,   min=0.5,  max=40)
p2.add('N1',5.2e4,min=1e5, max=5e6)

p2.add('sigma2',332,min=3e2,  max=8e2)
p2.add('dv2',12035,   min=8.5e3,max=14e3,vary=True)
p2.add('B2',7,   min=0,  max=20)
p2.add('N2',1.1e3,min=1e4, max=1e6)

ps2 = [p2]


def bestFit(ps,spec):
    std = (spec-n2h1(vel21,ps[0])).std()
    pr = ps[0]
    for p in ps:
        if (spec-n2h1(vel21,p)).std() < std:
            std = (spec-n2h1(vel21,p)).std()
            pr = p
    return pr
        
def bestFit2(ps,spec):
    std = (spec-n2h1(vel21,ps[0])-n2h2(vel21,ps[0])).std()
    pr = ps[0]
    for p in ps:
        if (spec-n2h1(vel21,p)-n2h2(vel21,p)).std() < std:
            std = (spec-n2h1(vel21,p)-n2h2(vel21,p)).std()
            pr = p
    return pr





def res1(p,x,y,c):
    return (y - n2h1.fitfunc(x,p))

def res2(p,x,y,c):
    return (y - n2hr.fitfunc(x,p)-n2hb.fitfunc(x,p))
    
def res3(p,x,y,c):
    return (y - nh31.fitfunc(x,p)-nh32.fitfunc(x,p)-nh33.fitfunc(x,p))

parnh = [['' for z in range( 0,n2hbg.shape[1])] for y in range(0,n2hbg.shape[0])]
figure(1)
import time
cla()
ion()
show()
imshow(n2hbg, vmax=cube2[:,six.min():six.max(),siy.min():siy.max()].max())
figure(1)
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
            r =  estimator.predict(spec.reshape(1,-1))[0]
            p.add('sigma1',r[3],min=2e2,  max=2e3)
            p.add('dv1',r[0],   min=6076,max=13252,vary=True)
            p.add('B1',r[1],   min=0,  max=35)
            p.add('N1',r[2],min=2e5, max=7e6)

#            p = bestFit(ps,spec)
 #           v0 = vel21[where(spec == spec.max())][0]
 #           p['dv1'].value = v0
            mini = lmfit.Minimizer(res1,p,fcn_args=(vel21,spec,[]))


            pr1 = mini.minimize(method='leastsq')
            print('\r({0},{1}st-knn ),'.format(xi,yi),end="")
            plot(yi,xi,'oy')
            pyplot.draw()
            plt.pause(0.01)
            time.sleep(0.05)
            
#            if (spec-n2h1.fitfunc(vel21,p))[where(maksh)].std() > 3*spec[where(logical_not(mak))].std():
#                plot(yi,xi,'og')
#                ps.append(par2hned.params)
            
            
            
            mak = (n2h1.fitfunc(vel21,pr1.params)>2e-1)

            
            st1 = (spec-n2h1.fitfunc(vel21,pr1.params)).std()

            st = spec.std()
#            b2 = par[xi][yi].params['B2']
            if st1 > st:
                plot(yi,xi,'xw')
                print('\r({0},{1}-),'.format(xi,yi),end="")
                parnh[xi][yi] = 'none'
                continue
            
            if (spec-n2h1.fitfunc(vel21,pr1.params))[where(maksh)].std() > 1.5*spec[where(logical_not(mak))].std():
                p = bestFit2(ps2,spec)
                
                mini = lmfit.Minimizer(res2,p,fcn_args=(vel21,spec,[]))
                par2hned = mini.minimize(method='nedler')
                par2hned = mini.minimize(method='leastsq',params=par2hned.params)
    
    

    
                if (pr1.chisqr*0.6 < par2hned.chisqr):
                    parnh[xi][yi] = pr1 
                    plot(yi,xi,'>b')
                else:
                    plot(yi,xi,'>r')
                    parnh[xi][yi] = par2hned
                    if (spec-n2h1.fitfunc(vel21,p)-n2h2.fitfunc(vel21,p)).std() > 1.5*(spec-n2h1.fitfunc(vel21,par2hned.params)-n2h2.fitfunc(vel21,par2hned.params)).std():
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
lon.set_major_formatter('hh:mm:ss.s')
lon.set_axislabel(r'$\alpha$')
lat.set_axislabel(r'$\delta$')
lon.set_ticks(spacing=20. * u.arcsec,color='w')
lat.set_ticks(spacing=20. * u.arcsec,color='w')
lon.display_minor_ticks(True)
lat.display_minor_ticks(True)
lat.set_minor_frequency(5)
lon.set_minor_frequency(5)
lat.set_ticklabel_position('r')
lat.set_axislabel_position('r')

masers = open('../../clumps.cat')
masers = open('../../metanol_probes.cat')
masers = masers.readlines()

fx,fy = getFromWorld(recs)
for m in masers:
    if len(m)<7:
        continue
    m = m.split()
    c = SkyCoord(m[1],m[2],unit=(u.hourangle, u.deg),frame='fk5')
    text(fx(c.ra)+1,fy(c.dec),m[0][-4:],color='k',size='large',horizontalalignment='left',verticalalignment='center')
    scatter([fx(c.ra)],[fy(c.dec)],c='y',marker='*',zorder=3)


masers = open('../../clumps.cat')
masers = masers.readlines()

fx,fy = getFromWorld(recs)
for m in masers:
    if len(m)<7:
        continue
    m = m.split()
    c = SkyCoord(m[1],m[2],unit=(u.hourangle, u.deg),frame='fk5')
    text(fx(c.ra)+1,fy(c.dec),m[0][-4:],color='c',size='large',horizontalalignment='left',verticalalignment='center')
    scatter([fx(c.ra)],[fy(c.dec)],c='y',marker=(5, 1),zorder=3)


text(6.82e3,2e3,'1',color='k',size='x-large',horizontalalignment='center',verticalalignment='center')
text(7.808e3,2e3,'4',color='k',size='x-large',horizontalalignment='center',verticalalignment='center')
text(8.3e3,2e3,'2',color='k',size='x-large',horizontalalignment='center',verticalalignment='center')
text(9.23e3,2e3,'3',color='k',size='x-large',horizontalalignment='center',verticalalignment='center')
#text(8.6e3,1.2e3,'3.1',color='k',size='x-large',horizontalalignment='center',verticalalignment='center')
text(10.6e3,2e3,'5',color='k',size='x-large',horizontalalignment='center',verticalalignment='center')
text(9.766e3,2e3,'6',color='k',size='x-large',horizontalalignment='center',verticalalignment='center')


plot([7.2e3,7.2e3],[1.2e3,2e2],'--k',linewidth=0.5)
plot([7.808e3,7.808e3],[1.2e3,2e2],'--k',linewidth=0.5)
plot([8.3e3,8.3e3],[1.2e3,2e2],'--k',linewidth=0.5)
plot([9e3,9e3],[1.2e3,2e2],'--k',linewidth=0.5)
plot([9.5e3,9.5e3],[1.2e3,2e2],'--k',linewidth=0.5)
#plot([10.5e3,10.5e3],[9e2,2e2],'-k')


text(7.2e3,1.3e3,'1',color='k',size='x-large',horizontalalignment='center',verticalalignment='center')
text(7.8e3,1.3e3,'2',color='k',size='x-large',horizontalalignment='center',verticalalignment='center')
text(8.3e3,1.3e3,'3',color='k',size='x-large',horizontalalignment='center',verticalalignment='center')
text(9e3,1.3e3,'4',color='k',size='x-large',horizontalalignment='center',verticalalignment='center')
#text(8.6e3,1.2e3,'3.1',color='k',size='x-large',horizontalalignment='center',verticalalignment='center')
text(9.5e3,1.3e3,'5',color='k',size='x-large',horizontalalignment='center',verticalalignment='center')
#text(10.5e3,2e3,'6',color='k',size='x-large',horizontalalignment='center',verticalalignment='center')


masers = open('../../clumps.cat')
masers = masers.readlines()

fx,fy = getFromWorld(recs)
for m in masers:
    if len(m)<7:
        continue
    m = m.split()
    scatter(float(m[3])*1e3,1246,c='k',marker=(5, 1),zorder=3)
    
    




def getImages(par,n2h1,n2h2,nh3bg):
    tau = ma.zeros((nh3bg.shape[0],nh3bg.shape[1],2))
    inte =ma.zeros((nh3bg.shape[0],nh3bg.shape[1],2))
    df1 =ma.zeros((nh3bg.shape[0],nh3bg.shape[1],2))
    sigma1 =ma.zeros((nh3bg.shape[0],nh3bg.shape[1],2))
    B1 =ma.zeros((nh3bg.shape[0],nh3bg.shape[1],2))
    Nu1 =ma.zeros((nh3bg.shape[0],nh3bg.shape[1],2))
    
    edf1 =ma.zeros((nh3bg.shape[0],nh3bg.shape[1],2))
    esigma1 =ma.zeros((nh3bg.shape[0],nh3bg.shape[1],2))
    eB1 =ma.zeros((nh3bg.shape[0],nh3bg.shape[1],2))
    eNu1 =ma.zeros((nh3bg.shape[0],nh3bg.shape[1],2))
    
    tau.mask=False
    inte.mask=False
    df1.mask=False
    sigma1.mask=False
    B1.mask=False
    Nu1.mask=False
    
    edf1.mask=False
    esigma1.mask=False
    eB1.mask=False
    eNu1.mask=False
    
    xi,yi=0,0
    for xi in range(nh3bg.shape[0]):
        for yi in range(nh3bg.shape[0]):
            p = par[xi][yi]
            if (p == '') or (p == 'none'):
                tau.mask[xi,yi]|= True
                inte.mask[xi,yi]|= True
                Nu1.mask[xi,yi]|= True
                sigma1.mask[xi,yi]|= True
                df1.mask[xi,yi]|= True
                B1.mask[xi,yi]|= True
    
                eNu1.mask[xi,yi]|= True
                esigma1.mask[xi,yi]|= True
                edf1.mask[xi,yi]|= True
                eB1.mask[xi,yi]|= True
    
    
            else :
                if hasattr(n2h1,'tau') :
                    tau[xi,yi,0]=n2h1.tau(x,p.params).max()
                inte[xi,yi,0] = n2h1.fitfunc(x,p.params).max() 
                Nu1[xi,yi,0] = p.params['N1']
                edf1[xi,yi,0] = p.params['dv1'].stderr
                sigma1[xi,yi,0] = p.params['sigma1']
                esigma1[xi,yi,0] = p.params['sigma1'].stderr
                df1[xi,yi,0] = p.params['dv1']
                B1[xi,yi,0] = p.params['B1']
                
                if 'B2' in p.params:
                    B1[xi,yi,1] = p.params['B2']
#                    tau[xi,yi,1]=g2.tau(x,p.params).max()
                    inte[xi,yi,1] = n2h2.fitfunc(x,p.params).max() 
                    edf1[xi,yi,1] = p.params['dv2'].stderr
                    esigma1[xi,yi,1] = p.params['sigma2'].stderr

                    Nu1[xi,yi,1] = p.params['N2']
                    sigma1[xi,yi,1] = p.params['sigma2']
                    df1[xi,yi,1] = p.params['dv2']
    
    tau[where(tau==0)]=ma.masked
    inte[where(inte==0)]=ma.masked
    edf1[where(inte==0)]=ma.masked

    return (inte,tau,Nu1,sigma1,df1,B1,edf1,esigma1)

xi = 60
shift = 0
for yi in range(27,120,10):
        spec = cube[:,xi,yi]
        specco = cubeco[:,int(fx(xi)),int(fx(yi))]
        specco/=specco.max()/spec.max()
        p = Parameters()
        r =  estimator.predict(spec.reshape(1,-1))[0]
        p.add('sigma1',r[3],min=2e2,  max=2e3)
        p.add('dv1',r[0],   min=6076,max=13252,vary=True)
        p.add('B1',r[1],   min=0,  max=100)
        p.add('N1',r[2],min=2e3, max=7e6)

#            p = bestFit(ps,spec)
 #           v0 = vel21[where(spec == spec.max())][0]
 #           p['dv1'].value = v0
        mini = lmfit.Minimizer(res1,p,fcn_args=(vel21,spec,[]))


        pr1= mini.minimize(method='leastsq')
        p = pr1.params
        step(x,shift+spec/spec.max(),'-k')
        step(xco,shift+specco/spec.max(),'-k',linewidth=0.5)
        plot(x,shift+n2h1.fitfunc(vel21,p)/spec.max())
        step(x,shift+(spec-n2h1.fitfunc(vel21,p))/spec.max(),'-g',linewidth=0.5)
        shift+=1.1+0.1
        
        
        
xs = []
ys = []
fx,fy = getExtent2(recsco,recsnh)
for (xi,yi),v in ndenumerate(n2hbg):
    if not dvf.mask[xi,yi,0]:
        if not dvf.mask[xi,yi,1]:
            xs.append(dvf[xi,yi,0])
            xs.append(dvf[xi,yi,1])
            ys.append(co[int(fx(xi)),int(fy(yi))])
            ys.append(co[int(fx(xi)),int(fy(yi))])
        else: 
            xs.append(dvf[xi,yi,0])
            ys.append(co[int(fx(xi)),int(fy(yi))])


from sklearn.cluster import SpectralClustering
scl = SpectralClustering(
        n_clusters=30, eigen_solver='arpack',
        affinity="nearest_neighbors")

from sklearn.cluster import DBSCAN
scl = DBSCAN()


tkinimg = np.ma.zeros((nh3bg.shape[0],nh3bg.shape[1],2))
for (xi,yi),val in np.ndenumerate(nh3bg):
    p1 = par1[xi][yi]
    p2 = par2[xi][yi]
    spec = cube[:,xi,yi]
    if (spec.max()>spec[where(logical_not(maks))].std()*factor):

        if p1 != '' and p1 != 'none' and p2 != '' and p2 != 'none':
            if ('B1' in p1.params):
                tkinimg[xi,yi,0] = nh31.Tex(nh32r,par1[xi][yi].params,par2[xi][yi].params,velx)
            else :         tkinimg[xi,yi,0] = np.ma.masked
    
            if ('B2' in p1.params):
                tkinimg[xi,yi,1] = nh32.Tex(nh32b,par1[xi][yi].params,par2[xi][yi].params,velx)
            else :         tkinimg[xi,yi,1] = np.ma.masked
        else:
            tkinimg[xi,yi] = np.ma.masked
    else:
        tkinimg[xi,yi] = np.ma.masked

tkinimg = np.ma.zeros(cube[0].shape)

for (xi,yi),vl in ndenumerate(cube[1]):
    if (not dvf.mask[xi,yi,0]) and (not dvf2.mask[xi,yi,0]):
        p1 = Parameters()
        p1.add('sigma1',dvf[xi,yi,3])
        p1.add('dv1',dvf[xi,yi,0])
        p1.add('B1',dvf[xi,yi,1])
        p1.add('Nu1',dvf[xi,yi,2])
    
        p2 = Parameters()
        p2.add('sigma1',dvf2[xi,yi,3])
        p2.add('dv1',dvf2[xi,yi,0])
        p2.add('B1',dvf2[xi,yi,1])
        p2.add('Nu1',dvf2[xi,yi,2])
        tkinimg[xi,yi] = nh31.Tex(nh312,p1,p2,velx)
    else:
        tkinimg[xi,yi] = np.ma.masked


y[0,i] = dv[i]
y[1,i] = B[i]
y[2,i] = Nu[i]
y[3,i] = sigma1[i]
for i in ids[0]:
    p = Parameters()
    p.add('dv1',dv[i],   min=6076,max=13252,vary=True)
    p.add('B1',B[i],   min=0,  max=100)
    p.add('N1',Nu[i])
    p.add('sigma1',sigma1[i])
    subplot(1,7,where(ids[0]==i)[0]+1)
    plot(x,n2h1.fitfunc(x,p),'-r')
    ax = gca()
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    xlim(-3e3,20e3)
    
    
parnhv2 = [['' for z in range( 0,n2hbg.shape[1])] for y in range(0,n2hbg.shape[0])]
for (xi,yi),vl in ndenumerate(cube[1]):
    if (not dvf.mask[xi,yi,0]):
        spec = cube[:,xi,yi]
        if (not dvf.mask[xi,yi,1]):
            p = Parameters()
            r =  estimator.predict(spec.reshape(1,-1))[0]
            p.add('sigma1',r[3],min=2e2,  max=2e3)
            p.add('dv1',r[0],   min=6076,max=13252,vary=False)
            p.add('B1',r[1],   min=0,  max=35)
            p.add('N1',r[2],min=2e5, max=7e6)
            p.add('sigma2',r[7],min=2e2,  max=2e3)
            p.add('dv2',r[4],   min=6076,max=13252,vary=True)
            p.add('B2',r[5],   min=0,  max=35)
            p.add('N2',r[6],min=2e5, max=7e6)
            mini = lmfit.Minimizer(res2,p,fcn_args=(v,spec,[]))


            pr1 = mini.minimize(method='leastsq')
            print('\r({0},{1}++),'.format(xi,yi),end="")
            parnhv2[xi][yi] = pr1 
        p = Parameters()
        r = estimator.predict(spec.reshape(1,-1))[0]
        p.add('sigma1',r[3],min=2e2,  max=2e3)
        p.add('dv1',r[0],   min=6076,max=13252,vary=False)
        p.add('B1',r[1],   min=0,  max=35)
        p.add('N1',r[2],min=2e5, max=7e6)
        mini = lmfit.Minimizer(res1,p,fcn_args=(v,spec,[]))


        pr1 = mini.minimize(method='leastsq')
        print('\r({0},{1}+),'.format(xi,yi),end="")
        parnhv2[xi][yi] = pr1 
    else: print('\r({0},{1}-),'.format(xi,yi),end="")

