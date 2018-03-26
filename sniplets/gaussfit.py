# -*- coding: utf-8 -*-
"""
Created on Thu Sep  8 16:47:54 2016

@author: pete
"""

cs1 = Analyse.GaussModel('1',restfreq=recs[0].header['restfreq'])
cs2 = Analyse.GaussModel('2',restfreq=recs[0].header['restfreq'])

spec = cubes[:,64,64]


ws = wcs.WCS(recss[0].header)

chan = arange(len(spec))+1

fr = zeros((chan.shape[0],4))
fr[:,2] = chan
fr = ws.all_pix2world(fr,1)
x = fr[:,2]
xx = linspace(x[0],x[-1],len(x)*10)
vel = xx#(xx*u.Hz).to(u.km/u.s, equivalencies=u.doppler_radio(recs[0].header['RESTFREQ']*u.Hz)).value
velx = x#(x*u.Hz).to(u.km/u.s, equivalencies=u.doppler_radio(recs[0].header['RESTFREQ']*u.Hz)).value

params = Parameters()
params.add('sigma2',450,min=2e1,  max=8e2)
params.add('dv2',8.1e3, min=x.min(),max=x.max(),vary=True)
params.add('B2',8.3,  min=.1,  max=10)

params.add('sigma1',450,min=2e1,  max=6e2)
params.add('dv1',10e3,  min=x.min(),max=x.max(),vary=True)
params.add('B1',8.462,   min=0.1,  max=10)
params.add('a',0,min=-10,max=10,vary=True)
params.add('b',0,min=-1,max=1,vary=True)
params.add('c',0.1,min=-1,max=1,vary=True)



def residual(p,x,y,c):
    return (y - cs1.fitfunc(x,p)-cs2.fitfunc(x,p)-p['a'].value*x**2-p['b'].value*x-p['c'].value)**2


ppar = minimize(residual,params,args=(x,spec,[]),maxfev=10000)
figure()
step(velx+(velx[1]-velx[0])/2,spec,'k')
p = ppar.params
plot(xx,cs1.fitfunc(xx,ppar.params)-(-p['a'].value*xx**2-p['b'].value*xx-p['c'].value),'r')
plot(xx,cs2.fitfunc(xx,ppar.params)-(-p['a'].value*xx**2-p['b'].value*xx-p['c'].value),'b')
params = ppar.params
params['a'].vary = False
params['b'].vary = False
params['c'].vary = False
params['a'].value = 0
params['b'].value = 0
params['c'].value = 0
maks = (cs1.fitfunc(x,params)>2e-3)+cs2.fitfunc(x,params)>2e-3
bg = cubes[where(maks),:,:][0].max(0)


p1 = Parameters()
p1.add('sigma1',552,min=2e2,  max=6e2)
p1.add('dv1',7846,   min=3376,max=13252,vary=True)
p1.add('B1',.10,   min=0,  max=100)


p2 = Parameters()

p2.add('sigma1',552,min=2e2,  max=8e2)
p2.add('dv1',10240,   min=3376,max=12252,vary=True)
p2.add('B1',.4,   min=0,  max=100)

p2.add('sigma2',591,min=2e2,  max=8e2)
p2.add('dv2',7130,   min=3376,max=12252,vary=True)
p2.add('B2',.4,   min=0,  max=100)



p1 = Parameters()
p1.add('sigma1',552,min=2e2,  max=6e2)
p1.add('dv1',1710,   min=-1e4,max=1e4,vary=True)
p1.add('B1',.10,   min=.1,  max=100)


p2 = Parameters()

p2.add('sigma1',52,min=2e2,  max=8e3)
p2.add('dv1',1710,   min=-1e4,max=1e4,vary=True)
p2.add('B1',.1,   min=.1,  max=100)

p2.add('sigma2',552,min=2e2,  max=8e2)
p2.add('dv2',1810,   min=-1e4,max=1e4,vary=True)
p2.add('B2',.10,   min=.1,  max=100)

def res1(p,x,y,c):
    return (y - cs1.fitfunc(x,p))**2

def res2(p,x,y,c):
    return (y - cs1.fitfunc(x,p)-cs2.fitfunc(x,p))**2

pars = [['' for z in range( 0,bg.shape[1])] for y in range(0,bg.shape[0])]
figure(1)
import time
cla()
ion()
show()
imshow(bg, vmax=bg[six.min():six.max(),siy.min():siy.max()].max())
xlim(siy.min(),siy.max())
ylim(six.min(),six.max())
for xi in ix:
    for yi in iy:
        
        spec = cubes[:,xi,yi]
        if (bg[xi,yi]>spec[where(logical_not(maks))].std()*factor) and pars[xi][yi] == '':

            p = p1
            v0 = velx[where(spec == spec.max())][0]
            p['dv1'].value = v0
            mini = lmfit.Minimizer(res1,p,fcn_args=(x,spec,[]))
            parhned = mini.minimize(method='nelder')
            print('\r({0},{1}+           '.format(xi,yi),end="")
            plot(yi,xi,'o')
            pyplot.draw()
            plt.pause(0.01)
            time.sleep(0.05)
#            double+=1
            mak = (cs1.fitfunc(x,parhned.params)>2e-1)

            st1 = spec[where(logical_not(maks))].std()
            maa1 = (sqrt(res1(parhned.params,velx,spec,[])))
            maa1 = maa1[maks].std()

            b1 = parhned.params['B1']
#            b2 = pars[xi][yi].params['B2']
            if cs1.fitfunc(x,parhned.params).max() < st1:
                plot(yi,xi,'xw')
                print('\r({0},{1}-),'.format(xi,yi),end="")
                pars[xi][yi] = 'none'
#                continue
            
            
            else :
 #           if :
                p = p2
                
                mini = lmfit.Minimizer(res2,p,fcn_args=(x,spec,[]))
                p2hned = mini.minimize(method='nelder')

                mak = (cs1.fitfunc(x,p2hned.params)>2e-3)+(cs2.fitfunc(x,p2hned.params)>2e-3)

                st = spec[where(logical_not(mak))].std()
                maa = (sqrt(res2(p2hned.params,velx,spec,[])))
                maa = maa[mak].std()

                print('\r({0},{1}y),'.format(xi,yi),end="")
#            if pars[xi][yi].params['df1'].stderr > 100000 or pars[xi][yi].params['df2'].stderr > 100000 :
#                pars[xi][yi] = 'none'
#                print('\r({0},{1}x),'.format(xi,yi),end="")

                
                if(p2hned.chisqr*1.5<parhned.chisqr):
                    plot(yi,xi,'>r')
                    pars[xi][yi] = p2hned
                else:
                    pars[xi][yi] = parhned
                    plot(yi,xi,'>g')
                    
                if st*1.5<maa :
                    plot(yi,xi,'oy')
                    print('\r({0},{1}-badfit),'.format(xi,yi),end="")
                    pars[xi][yi] = 'none'
                
#                    badfit+=1

        else:
            plot(yi,xi,'xk')
            print('\r({0},{1}-),'.format(xi,yi),end="")
#x            pars[xi][yi] = 'none'
