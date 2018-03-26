from numpy import *
from matplotlib.pyplot import *
cube = zeros((127,246,246))
cube = zeros((119,82,82))
BMAJ=7.7541E-04*u.deg
BMIN=  6.8027E-04*u.deg
dat =recs[0].data[0]

def loadCube(rec,crop=[0,0],step=1,aver=3):
    dat =rec.data[0]
    cube = np.zeros((dat.shape[0]-abs(crop[0])-abs(crop[1]),dat.shape[1]/step,dat.shape[2]/step))
    for (x,y),val in ndenumerate(cube[0]):
        if crop[1] == 0 :
            cube[:,x,y] = dat[crop[0]:,x*step:x*step+aver,y*step:y*step+aver].mean(1).mean(1)
        else:            
            cube[:,x,y] = dat[crop[0]:crop[1],x*step:x*step+aver,y*step:y*step+aver].mean(1).mean(1)
    return cube

def loadCubeConvolve(rec,std=1):
    from astropy.convolution import convolve, convolve_fft
    from astropy.convolution import Gaussian2DKernel
    g = Gaussian2DKernel(stddev=std)
    dat =rec.data[0]
    cube = np.zeros(dat.shape)
    for i in range(dat.shape[0]):
        cube[i] =  convolve_fft(dat[i],g)
    return cube
    

    
def moment0(cube,v,mask):
    vv = zeros((cube.shape[1],cube.shape[2]))
    dv = abs(v[0]-v[1])
    for (ix,iy),val in ndenumerate(vv):
        for ic in range(cube.shape[0]):
            if mask[ic]:
                vv[ix,iy] += cube[ic,ix,iy]*dv

    return vv

def moment(cube,v,mask):
    vv = zeros((cube.shape[1],cube.shape[2]))
    for (ix,iy),val in ndenumerate(vv):
        for ic in range(cube.shape[0]):
            if mask[ic]:
                vv[ix,iy] += cube[ic,ix,iy]*v[ic]
        vv[ix,iy]/=cube[mask,ix,iy].sum()
    return vv

def moment2(cube,v,mask):
    m1 = moment(cube,v,mask)
    vv = zeros((cube.shape[1],cube.shape[2]))
    for (ix,iy),val in ndenumerate(vv):
        for ic in range(cube.shape[0]):
            if mask[ic]:
                vv[ix,iy] += cube[ic,ix,iy]*(v[ic]-m1[ix,iy])**2
        vv[ix,iy]/=cube[mask,ix,iy].sum()
        vv[ix,iy]=sqrt(vv[ix,iy])
    return vv


cube[:4]=0
cube[-4:]=0

spec = cube[:,150,132]

crpix = recs[0].header['crpix1']
cdelt3  = recs[0].header['deltav']
crvav  = recs[0].header['velo-lsr']
nax  = recs[0].header['naxis1']

crpix = recs[0].header['crpix3']
cdelt3  = recs[0].header['cdelt3']
crvav  = recs[0].header['crval3']
nax  = recs[0].header['naxis3']

from astropy import wcs
w = wcs.WCS(recs[0].header)

chan = arange(len(spec))+1

fr = zeros((chan.shape[0],4))
fr[:,2] = chan
fr = w.all_pix2world(fr,1)
x = fr[:,2]
xx = linspace(x[0],x[-1],len(x)*10)
vel = xx#(xx*u.Hz).to(u.km/u.s, equivalencies=u.doppler_radio(recs[0].header['RESTFREQ']*u.Hz)).value
velx = x#(x*u.Hz).to(u.km/u.s, equivalencies=u.doppler_radio(recs[0].header['RESTFREQ']*u.Hz)).value
#cube[-5:,:,:]=ma.masked
#cube[:4,:,:]=ma.masked
#spec = cube[:,150,132]
nh31 = Analyse.NH3Model(cat,ind='1') 
nh32 = Analyse.NH3Model(cat,ind='2')

#nh31 = Analyse.NH3Model(cat,ind='1',restfreq=23722633.335e3)
#nh32 = Analyse.NH3Model(cat,ind='2',restfreq=23722633.335e3)
params = Parameters()
params.add('sigma2',150,min=2e2,  max=8e2)
params.add('dv2',8.1e3, min=7625,max=12252,vary=True)
params.add('B2',8.3,  min=nh31.Bbg,  max=10)
params.add('Nu2',4.8e28,  min=1e20, max=1e29)

params.add('sigma1',150,min=2e2,  max=6e2)
params.add('dv1',10.4e3,   min=8790,max=12252,vary=True)
params.add('B1',8.462,   min=nh31.Bbg,  max=10)
params.add('Nu1',4.9e28,min=2e20, max=7e29)
params.add('a',0,min=-10,max=10)
params.add('b',0,min=-1,max=1)
params.add('c',0,min=-1,max=1)



def residual(p,x,y,c):
    return (y - nh31.fitfunc(x,p)-nh32.fitfunc(x,p)-p['a'].value*x**2-p['b'].value*x-p['c'].value)


par = minimize(residual,params,args=(x,spec,[]),maxfev=10000)
figure()
step(velx+(velx[1]-velx[0])/2,spec,'k')
plot(xx,nh31.fitfunc(xx,par.params)+par.params['a'].value*xx**2+par.params['b'].value*xx+par.params['c'].value,'r')
plot(xx,nh32.fitfunc(xx,par.params)+par.params['a'].value*xx**2+par.params['b'].value*xx+par.params['c'].value,'b')
params = par.params
params['a'].vary = False
params['b'].vary = False
params['c'].vary = False
params['a'].value = 0
params['b'].value = 0
params['c'].value = 0
maks = (nh31.fitfunc(x,params)>2e-3)+nh32.fitfunc(x,params)>2e-3
nh3bg = cube[where(maks),:,:][0].max(0)

bases = zeros(cube.shape)   
base = lambda p,x: p['a'].value*(x)**2+p['b'].value*(x)+p['c'].value+p['d'].value*(x)**3
def resbase(p,x,y,c):
    return ((y-base(p,x))**2)*c


p = Parameters()
p.add('a',0.0,min=-1,max=1)
p.add('b',0.00,min=-3,max=3)
p.add('c',0.00,min=-10,max=10)
p.add('d',0.00,min=-10,max=10)
ix = arange(0,len(x),1)

for xi in range(0,cube.shape[1]):
    for yi in range(0,cube.shape[2]):  
        spec = cube[:,xi,yi]
        
#        spec = spec.reshape(spec.shape)
        ppar = minimize(resbase,p,args=(ix,spec,logical_not(maks)))
        bases[:,xi,yi] = base(ppar.params,ix)
        print('\r({0},{1}+),'.format(xi,yi),end="")

cube -= bases
double = 0
red = 0
blue = 0
none = 0         
badfit = 0

def res1(p,x,y,c):
    return (y - nh31.fitfunc(x,p))**2

def res2(p,x,y,c):
    return (y - nh31.fitfunc(x,p)-nh32.fitfunc(x,p))**2
    
def res3(p,x,y,c):
    return (y - nh31.fitfunc(x,p)-nh32.fitfunc(x,p)-nh33.fitfunc(x,p))**2

ix = arange(109,144)
iy = arange(112,144 )
factor = 3
par = [['' for z in range( 0,nh3bg.shape[1])] for y in range(0,nh3bg.shape[0])]
figure(1)
import time
cla()
ion()
show()
imshow(nh3bg, vmax=cube[:,six.min():six.max(),siy.min():siy.max()].max())
for xi in ix:
    for yi in iy:
        
        spec = cube[:,xi,yi]
        if (nh3bg[xi,yi]>spec[where(logical_not(maks))].std()*factor):
            plot(yi,xi,'og')
            par[xi][yi] = minimize(residual,params,args=(x,spec,[]),)
            print('\r({0},{1}+ (d:{2},r:{3},b:{4},bf:{5})),'.format(xi,yi,double,red,blue,badfit),end="")
            pyplot.draw()
            plt.pause(0.01)
            time.sleep(0.05)
            double+=1
            b1 = par[xi][yi].params['B1']
            b2 = par[xi][yi].params['B2']
            if spec[where(logical_not(maks))].std()*3 > nh31.fitfunc(x,par[xi][yi].params).max():
 #           if :
                p = Parameters()
                for k in params: 
                    if k.endswith('2'):
                        p.add(k,par[xi][yi].params[k].value,min=par[xi][yi].params[k].min,max=par[xi][yi].params[k].max)
                
                p['dv2'].min = 0                
                p['dv2'].max = 2e4
                p.add('a',0,vary=False)
                p.add('b',0,vary=False)
                p.add('c',0,vary=False)

                par[xi][yi] = minimize(res2,p,args=(x,spec,[]),)

                for k in params: 
                    if k.endswith('1'):
                        par[xi][yi].params.add(k,params[k].value,min=params[k].min,max=params[k].max)
                        
                par[xi][yi].params['B1'].value = nh31.Bbg
                par[xi][yi].params['B1'].vary = False

                red+=1
                plot(yi,xi,'ob')
                print('\r({0},{1}y),'.format(xi,yi),end="")
            if spec[where(logical_not(maks))]. std()*3 > nh32.fitfunc(x,par[xi][yi].params).max():
#            if dv2.stderr/dv2.value > 1:
                p = Parameters()
                for k in params: 
                    if k.endswith('1'):
                        p.add(k,par[xi][yi].params[k].value,min=params[k].min,max=params[k].max)
                
                p['dv1'].max = 0          
                p['dv1'].max = 2e4
                p.add('a',0,vary=False)
                p.add('b',0,vary=False)
                p.add('c',0,vary=False)
                par[xi][yi] = minimize(res1,p,args=(x,spec,[]),)

                for k in params: 
                    if k.endswith('2'):
                        par[xi][yi].params.add(k,params[k].value)
                par[xi][yi].params['B2'].value = nh31.Bbg
                par[xi][yi].params['B2'].vary = False
                blue+=1
                plot(yi,xi,'or')
                print('\r({0},{1}z),'.format(xi,yi),end="")
                
            
#            if par[xi][yi].params['df1'].stderr > 100000 or par[xi][yi].params['df2'].stderr > 100000 :
#                par[xi][yi] = 'none'
#                print('\r({0},{1}x),'.format(xi,yi),end="")
            st = spec[where(logical_not(maks))].std()
            maa = (sqrt(residual(par[xi][yi].params,velx,spec,[])))
            maa = maa[maks].std()
            if st*1.5<maa :
                plot(yi,xi,'oy')
                print('\r({0},{1}-badfit),'.format(xi,yi),end="")
                par[xi][yi] = 'none'
                badfit+=1

        else:
            plot(yi,xi,'ok')
            print('\r({0},{1}-),'.format(xi,yi),end="")
            par[xi][yi] = 'none'
                
            
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
for xi in ix:
    for yi in iy:
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
            tau[xi,yi,0]=nh31.tau(x,p.params).max()
            inte[xi,yi,0] = nh31.fitfunc(x,p.params).max() 
            Nu1[xi,yi,0] = p.params['Nu1']
            sigma1[xi,yi,0] = p.params['sigma1']
            df1[xi,yi,0] = p.params['dv1']
            B1[xi,yi,0] = p.params['B1']
            edf1[xi,yi,0] = p.params['dv1'].stderr

            if 'B2' in p.params:
                B1[xi,yi,1] = p.params['B2']
                tau[xi,yi,1]=nh32.tau(x,p.params).max()
                inte[xi,yi,1] = nh32.fitfunc(x,p.params).max() 
                Nu1[xi,yi,1] = p.params['Nu2']
                sigma1[xi,yi,1] = p.params['sigma2']
                df1[xi,yi,1] = p.params['dv2']

tau[where(tau==0)]=ma.masked
inte[where(inte==0)]=ma.masked

    
def plotImg(img,par,ix,iy):
    nptsx = len(ix)
    nptsy = len(iy)
    #x=linspace(nh3bg.shape[0]-1,0,nptsx,dtype=int)  
    #y=linspace(0,nh3bg.shape[1]-1,nptsy,dtype=int)  
    six=linspace(ix[-1]-1,ix[0],nptsx,dtype=int)  
    siy=linspace(iy[0],iy[-1]-1,nptsy,dtype=int)  
    figure(1)
    cla()
    imshow(img, vmax=img[six.min():six.max(),siy.min():siy.max()].max())
    i=0
    for xi in six:
        for yi in siy:
    #        text(yi,xi,"{0}".format(i), horizontalalignment='center',     verticalalignment='center')
            if par[xi][yi] != '' and par[xi][yi] != 'none':
                plot(yi,xi,'or')
            else :
                plot(yi,xi,'xk')
            i+=1


def plotSpecMap(cube,par,ix,iy,velx,vel,mol1,mol2,par2):
    nptsx = 10
    nptsy = 10
    #x=linspace(nh3bg.shape[0]-1,0,nptsx,dtype=int)  
    #y=linspace(0,nh3bg.shape[1]-1,nptsy,dtype=int)  
    six=linspace(ix[-1]-1,ix[0],nptsx,dtype=int)  
    siy=linspace(iy[0],iy[-1]-1,nptsy,dtype=int)  
    if ix != '':
        six=ix[0:-1:2][::-1]
        siy=iy[0:-1:2]
        
    f = figure(2)
    f.canvas.mpl_connect('button_press_event', click)
    cla()
    import matplotlib.gridspec as gridspec
    gs = gridspec.GridSpec(len(six), len(siy))
    gs.update(wspace=0.00,hspace=0)
    axes = [[0 for yi in six] for xi in siy]
    for xi in six:
        for yi in siy:
            ind = (where(six==xi)[0][0]),(where(siy==yi)[0][0])
            ax = subplot(gs[ind])
            axes[ind[1]][ind[0]] = ax
#            fx = int(tox[xi*3])
#            fy = int(toy[yi*3])
            ax.text(0,0,"({0}:{1}):{2}".format(xi,yi,ind),verticalalignment='top')
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)
    #        ax.set_ylim(cube[:,ix.min():ix.max(),iy.min():iy.max()].min()/2,cube[:,ix.min():ix.max(),iy.min():iy.max()].max())
            ax.set_ylim(cube[:,ix.min():ix.max(),iy.min():iy.max()].min()/4,cube[:,ix.min():ix.max(),iy.min():iy.max()].max())
            ax.set_xlim(-13000,37000)
            spec = cube[:,xi,yi]
            step(velx,spec,'k',linewidth=2)
#            spec2 = cube2[:,fx,fy]
#            step(vel21,spec2,'k')
            
            plot([nh31.v0[10]+params['dv2'].value,nh31.v0[10]+params['dv2'].value],[0,1],'-b',linewidth=0.3)
            plot([nh31.v0[10]+params['dv1'].value,nh31.v0[10]+params['dv1'].value],[0,1],'-r',linewidth=0.3)
    #        err = spec[where(logical_not(maks))]. std()*3
    #        plot([velx[0],velx[-1]],[err,err],'-y',linewidth=0.3)
            #        plot(xxx+(xxx[1]-xxx[0])/2,bases[:,xi,yi],'-g')
    #        ax.fill_between(velx+(velx[1]-velx[0])/2,maks*cube.max()/2,facecolor='red',alpha=0.3)
            if (not par == '') and par[xi][yi] != '' and par[xi][yi] != 'none':
#                plot(vel,nh31.fitfunc(xx,par[xi][yi].params),'r',linewidth=2)
#                plot(vel,nh32.fitfunc(xx,par[xi][yi].params),'b',linewidth=2)
                plot(vel,mol1.fitfunc(vel,par[xi][yi].params),'--r')
                plot(vel,mol2.fitfunc(vel,par[xi][yi].params),'--b')
                if par2 != '':
                    plot(vel,mol1.fitfunc(vel,par2[xi][yi].params),'r')
                    plot(vel,mol2.fitfunc(vel,par2[xi][yi].params),'b')
                
def click(ev):
    print('click!')
    for xi in six:
        for yi in siy:
            ind = where(six==xi)[0][0],where(siy==yi)[0][0]
            if ev.inaxes==axes[ind[1]][ind[0]]:
                fx = int(tox[xi*3])
                fy = int(toy[yi*3])
                print('{0},{1}'.format(xi,yi))
                figure(3)
                cla()
                spec = cube[:,xi,yi]
                step(velx+(velx[1]-velx[0])/2,spec,'k')
                spec2 = cube2[:,frx[xi],fry[yi]]
                step(vel21,spec2,'k')
                plot([velx[0],velx[-1]],[err,err],'-y',linewidth=0.3)
#                plot(xxx+(xxx[1]-xxx[0])/2,bases[:,xi,yi],'-g')
#                fill_between(velx+(velx[1]-velx[0])/2,maks*spec.max()/2,facecolor='red',alpha=0.3)
                if par[xi][yi] != '' and par[xi][yi] != 'none':
                    plot(vel,nh31.fitfunc(xx,par[xi][yi].params),'r')
                    plot(vel,nh32.fitfunc(xx,par[xi][yi].params),'b')
                    plot(vel,n2h1.fitfunc(xx,par2hned[fx][fy].params),'r')
                    plot(vel,n2h2.fitfunc(xx,par2hned[fx][fy].params),'b')
                    plot(vel,n2h1.fitfunc(xx,par2h[fx][fy].params),'--r')
                    plot(vel,n2h2.fitfunc(xx,par2h[fx][fy].params),'--b')
                
                plot()
                show()


step(vel1+(vel1[1]-vel1[0])/2,spec1,'k',linewidth=1)
step(vel2+(vel2[1]-vel2[0])/2,spec2,'k',linewidth=3)
plot(vel,nh31.fitfunc(vel,p1.params),'-r',linewidth=1)
plot(vel,nh32.fitfunc(vel,p1.params),'-b',linewidth=1)
plot(vel,nh32r.fitfunc(vel,p2.params),'--r',linewidth=3)
plot(vel,nh32b.fitfunc(vel,p2.params),'--b',linewidth=3)

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
    



for (xi,yi),val in np.ndenumerate(df1[:,:,0]):   
    if not (df1.mask[xi,yi,0]):
        ax.scatter(xi,yi,df1[xi,yi,0],s=inte[xi,yi,0]^2*20,c='r')
    if not (df1.mask[xi,yi,1]):
        ax.scatter(xi,yi,df1[xi,yi,1],s=inte[xi,yi,1]^2*20,c='b')
    
    
    
from astropy import wcs
w = wcs.WCS(recs[0].header)
w2 = wcs.WCS(recs2[0].header)

ix = arange(0,cube1.shape[1])
iy = arange(0,cube1.shape[1])

fig = plt.figure()
#ax = fig.add_axes([0.1,0.1,0.9,0.9],projection=w, slices=['x','y',0,0])
#tr = ax.get_transform(w2)

coords1 = zeros((len(iy),4))
coords2 = zeros((len(ix),4))
coords1[:,0] = iy
coords2[:,1] = ix
ra = w.all_pix2world(coords1,1,ra_dec_order=True)[:,0]
dec = w.all_pix2world(coords2,1,ra_dec_order=True)[:,1]

coords2[:,0] = ra[0]
coords2[:,1] = dec

coords1[:,0] = ra
coords1[:,1] = dec[0]


tox = w2.all_world2pix(coords2,1)[:,1]
toy = w2.all_world2pix(coords1,1)[:,0]


ix = arange(0,cube2.shape[1])
iy = arange(0,cube2.shape[1])

coords1 = zeros((len(iy),4))
coords2 = zeros((len(ix),4))

coords1[:,0] = iy
coords2[:,1] = ix

ra = w2.all_pix2world(coords1,1,ra_dec_order=True)[:,0]
dec = w2.all_pix2world(coords2,1,ra_dec_order=True)[:,1]

coords2[:,0] = ra[0]
coords2[:,1] = dec

coords1[:,0] = ra
coords1[:,1] = dec[0]


frx = w.all_world2pix(coords2,1)[:,1]
fry = w.all_world2pix(coords1,1)[:,0]



tr = ax.get_transform(w2)
ax.imshow(nh3bg,vmax=10)
ax.contour(cube2.mean(0),transform=tr,cmap='viridis',linewidth=10*ones((10)))
ax.plot(ra[0],dec[0],'or',transform=ax.get_transform('world'))
ax.plot(ra[-1],dec[-1],'or',transform=ax.get_transform('world'))
ax.plot(toy[0],tox[0],'xk',transform=tr)
ax.plot(toy[-1],tox[-1],'xk',transform=tr)

for i in ix:
    for j in iy:
        



n2h = Analyse.N2HModel(catnh,ind='1')
crpix = recs2[0].header['crpix3']
cdelt3  = recs2[0].header['cdelt3']
crvav  = recs2[0].header['crval3']
nax  = recs2[0].header['naxis3']


nx = arange(-crpix*cdelt3+crvav,(nax-crpix)*cdelt3+crvav,cdelt3)
nxx = arange(-crpix*cdelt3+crvav,(nax-crpix)*cdelt3+crvav,cdelt3/5)
nvel = xx#(xx*u.Hz).to(u.km/u.s, equivalencies=u.doppler_radio(recs[0].header['RESTFREQ']*u.Hz)).value
nvelx = x#(x*u.Hz).to(u.km/u.s, equivalencies=u.doppler_radio(recs[0].header['RESTFREQ']*u.Hz)).value


spec2 = cube2[:,tox[48*3],toy[43*3]]
chan2 = arange(len(spec2))+1
fr = zeros((chan2.shape[0],4))
fr[:,2] = chan2
fr = w2.all_pix2world(fr,1)
vel21 = fr[:,2]

def res2h(p,x,y,c):
    return (n2h1.fitfunc(x,p)+ n2h2.fitfunc(x,p)-y)**2
    



cube[88:93 ] = 0
cube[93:96 ] = 0
cube[83] = 0
cube[82] = 0

tiy = arange(toy[iy[0]*3],toy[iy[-1]*3],dtype=int)
tix = arange(tox[ix[0]*3],tox[ix[-1]*3],dtype=int)


par2h = [['' for z in range( 0,n2hbg.shape[1])] for y in range(0,n2hbg.shape[0])]
par2hned = [['' for z in range( 0,n2hbg.shape[1])] for y in range(0,n2hbg.shape[0])]

def res2h(p,x,y,c):
    return (n2h1.fitfunc(x,p)+ n2h2.fitfunc(x,p)-y)**2

def res1h(p,x,y,c):
    return (n2h1.fitfunc(x,p)-y)**2


for xi in ix:
    for yi in iy:
        fx = int(frx[xi])
        fy = int(fry[yi])
        if par[fx][fy] != '' and par[fx][fy] != 'none':
            spec = cube4[:,xi,yi]
            p = par[fx][fy]
            params = Parameters()
            params.add('sigma1',424,min=2e2,  max=8e2)
            params.add('dv1',p.params['dv1'].value,   min=3376,max=12252,vary=True)
            params.add('B1',6.5,   min=nh31.Bbg+1,  max=20)
            params.add('Ntot1',1.7e15,min=2e6, max=7e15)
            if 'dv2' in p.params:
                
                params.add('sigma2',736,min=2e2,  max=10e2)
                params.add('dv2',p.params['dv2'].value,   min=3376,max=12252,vary=True)
                params.add('B2',6.9,   min=nh31.Bbg+1,  max=20)
                params.add('Ntot2',5.4e16,min=2e6, max=7e16)
                par2h[xi][yi] = minimize(res2h,params,args=(vel21,spec,[]),)
                mini = lmfit.Minimizer(res2h,params,fcn_args=(vel21,spec,[]))
                par2hned[xi][yi] = par2h[xi][yi]
            else :
                par2h[xi][yi] = minimize(res1h,params,args=(vel21,spec,[]),)
                mini = lmfit.Minimizer(res1h,params,fcn_args=(vel21,spec,[]))
                par2hned[xi][yi] = par2h[xi][yi]
            ax.plot(yi,xi,'or')
            print('\r({0},{1}-),'.format(xi,yi),end="")
            pyplot.draw()
            plt.pause(0.01)
            time.sleep(0.05)
            

for xi in ix:
    for yi in iy:
        if par2h[xi][yi] != '' and par2h[xi][yi] != 'none':
            fx = int(frx[xi]/3)
            fy = int(fry[yi]/3)
            plot(fx,fy,'xr')
            
            
figure(7)
tiy = arange(toy[iy[0]],toy[iy[-1]],2,dtype=int)
tix = arange(tox[ix[0]],tox[ix[-1]],2,dtype=int)



imshow(n2hbg, vmax=n2hbg[tix.min():tix.max(),tiy.min():tiy.max()].max())

      
for xi in tix:
    for yi in tiy:
        fx = int(frx[xi])
        fy = int(fry[yi])
        p = par[fx][fy]
        if p != '' and p != 'none':
            plot(yi,xi,'og')
            spec = cube2[:,xi,yi]
            params = Parameters()
            params.add('sigma1',600,min=2e2,  max=10e2)
            params.add('dv1',p.params['dv1'],   min=3376,max=12252,vary=True)
            params.add('B1',71,   min=nh31.Bbg,  max=100)
            params.add('Ntot1',1e20,min=2e10, max=7e29)
            
            params.add('sigma2',422,min=2e2,  max=6e2)
            params.add('dv2',p.params['dv2'],   min=3376,max=12252,vary=True)
            params.add('B2',58,   min=nh31.Bbg,  max=100)
            params.add('Ntot2',5.29e19,min=2e10, max=7e29)
            mini = lmfit.Minimizer(res2h,params,fcn_args=(vel21,spec,[]))
            par2hned[xi][yi] = mini.minimize(method='nelder')
            par2h[xi][yi] = minimize(res2h,par2hned[xi][yi].params,args=(vel21,spec,[]),)
            plot(yi,xi,'or')
            print('\r({0},{1}+),'.format(xi,yi),end="")
        else: plot(yi,xi,'ob')
        
for xi in ix:
    for yi in iy:
        fx = int(frx[xi])
        fy = int(fry[yi])
        if par[fx][fy] != '' and par[fx][fy] != 'none' and par2h[xi][yi] == '':
            spec = cube2[:,xi,yi]
            p = par[fx][fy]
            v1 = p.params['dv1']
            v2 = p.params['dv2']
            params = Parameters()
            params.add('sigma1',150,min=2e2,  max=1e4)
            params.add('dv1',p.params['dv1'].value,   min=(v1+v2)/2,max=14000,vary=True)
            params.add('B1',6,   min=nh31.Bbg,  max=100)
            params.add('Ntot1',9e13,min=2e10, max=7e29)
            
            params.add('sigma2',150,min=2e2,  max=1e4)
            params.add('dv2',p.params['dv2'].value,   min=3376,max=(v1+v2)/2,vary=True)
            params.add('B2',6,   min=nh31.Bbg,  max=100)
            params.add('Ntot2',4.9e13,min=2e10, max=7e29)
            mini = lmfit.Minimizer(res2h,params,fcn_args=(vel3,spec,[]))
            par2hned[xi][yi] = mini.minimize(method='nelder')
            par2h[xi][yi] = minimize(res2h,par2hned[xi][yi].params,args=(vel3,spec,[]),)
            p = par2h[xi][yi].params
            pn = par2hned[xi][yi].params
            
            if res2h(p,vel21,spec,[]).mean() > res2h(pn,vel21,spec,[]).mean():
                par2h[xi][yi] = par2hned[xi][yi]
                p = pn
                
            plot(yi,xi,'o',color=(n2h1.fitfunc(vel21,p).max()/cube2.max(),0,n2h2.fitfunc(vel21,p).max()/cube2.max()))
            print('\r({0},{1}-),'.format(xi,yi),end="")
            pyplot.draw()
            plt.pause(0.01)
            time.sleep(0.05)
        else: plot(yi,xi,'xk')
        if par2h[xi][yi] != '':
            p = par2h[xi][yi].params
            plot(yi,xi,'o',color=(n2h1.fitfunc(vel21,p).max()/cube2.max(),0,n2h2.fitfunc(vel21,p).max()/cube2.max()))

 

def saveParsTxt(path,par):
    i = 0
    j = 0            
    for ps in par:
        for p in ps:
            f = open(path+'/p{0}_{1}.txt'.format(i,j),'w+')
            if p != '' and p != 'none':
                p.params.dump(f)
            else: 
                f.write('no data')
            f.close()
            i+=1
        i=0
        j+=1

def loadParsTxt(path,par):
    i = 0
    j = 0            
    for ps in par:
        for p in ps:
            f = open(path+'/p{0}_{1}.txt'.format(i,j),'w+')
            if f.readlines() 
                p = Parameters
                p.params.dump(f)
            else: 
                f.write('no data')
            f.close()
            i+=1
        i=0
        j+=1
            
        



p1 = Parameters()
p1.add('sigma1',400,min=2e2,  max=2e3)
p1.add('dv1',9000,   min=3376,max=13252,vary=True)
p1.add('B1',10,   min=nh31.Bbg,  max=100)
p1.add('Nu1',4.9e28,min=2e10, max=7e29)


p2 = Parameters()

p2.add('sigma1',374,min=2e2,  max=2e3)
p2.add('dv1',10035,   min=3376,max=12252,vary=True)
p2.add('B1',10,   min=nh31.Bbg,  max=100)
p2.add('Nu1',2.2e28,min=2e10, max=7e29)

p2.add('sigma2',271,min=2e2,  max=6e2)
p2.add('dv2',8023,   min=3376,max=12252,vary=True)
p2.add('B2',10,   min=nh31.Bbg,  max=100)
p2.add('Nu2',1.63e28,min=2e10, max=7e29)

def res1(p,x,y,c):
    return (y - nh31.fitfunc(x,p))**2

def res2(p,x,y,c):
    return (y - nh31.fitfunc(x,p)-nh32.fitfunc(x,p))**2
    
def res3(p,x,y,c):
    return (y - nh31.fitfunc(x,p)-nh32.fitfunc(x,p)-nh33.fitfunc(x,p))**2

par = [['' for z in range( 0,nh3bg.shape[1])] for y in range(0,nh3bg.shape[0])]
figure(1)
import time
cla()
ion()
show()
imshow(nh3bg)
xlim(siy.min(),siy.max())
ylim(six.min(),six.max())
for xi in ix:
    for yi in iy:
        
        spec = cube[:,xi,yi]
        if (nh3bg[xi,yi]>spec[where(logical_not(maks))].std()*factor and spec.sum() > 1e-5) and par[xi][yi] == '':

            p = p1
            v0 = velx[where(spec == spec.max())][0]
            p['dv1'].value = v0
            mini = lmfit.Minimizer(res1,p,fcn_args=(x,spec,[]))
            par2hned = mini.minimize(method='nelder')

            par[xi][yi] = minimize(res1,par2hned.params,args=(x,spec,[]),)
            print('\r({0},{1}+ (d:{2},r:{3},b:{4},bf:{5})),'.format(xi,yi,double,red,blue,badfit),end="")
            plot(yi,xi,'o',color=(0,0,(par[xi][yi].params['dv1']-3000)/14000))
            pyplot.draw()
            plt.pause(0.01)
            time.sleep(0.05)
            double+=1
            mak = (nh31.fitfunc(x,par[xi][yi].params)>2e-1)

            st = spec[where(logical_not(maks))].std()
            maa = (sqrt(res1(par[xi][yi].params,velx,spec,[])))
            maa = maa[maks].std()

            b1 = par[xi][yi].params['B1']
#            b2 = par[xi][yi].params['B2']
            if nh31.fitfunc(x,par[xi][yi].params).max() < st*2:
                plot(yi,xi,'xw')
                print('\r({0},{1}-),'.format(xi,yi),end="")
                par[xi][yi] = 'none'x = np.random.randn(n_samples)

#                continue
            
            
            else :
                if st < maa:
 #           if :
                    p = p2
                    
                    mini = lmfit.Minimizer(res2,p,fcn_args=(x,spec,[]))
                    par2hned = mini.minimize(method='nelder')
        
        
                    par[xi][yi] = minimize(res2,par2hned.params,args=(x,spec,[]),)
    
                    mak = (nh31.fitfunc(x,par[xi][yi].params)>2e-3)+(nh32.fitfunc(x,par[xi][yi].params)>2e-3)
    
                    red+=1
                    st = spec[where(logical_not(mak))].std()
                    maa = (sqrt(res2(par[xi][yi].params,velx,spec,[])))
                    maa = maa[mak].std()
    
                    plot(yi,xi,'>r')
                    print('\r({0},{1}y),'.format(xi,yi),end="")
    #            if par[xi][yi].params['df1'].stderr > 100000 or par[xi][yi].params['df2'].stderr > 100000 :
    #                par[xi][yi] = 'none'
    #                print('\r({0},{1}x),'.format(xi,yi),end="")
                if st*1.5<maa :
                    plot(yi,xi,'oy')
                    print('\r({0},{1}-badfit),'.format(xi,yi),end="")
                    par[xi][yi] = 'none'
                    badfit+=1

        else:
            plot(yi,xi,'xk')
            print('\r({0},{1}-),'.format(xi,yi),end="")
#x            par[xi][yi] = 'none'

                 
syncube = zeros(cube.shape)
for xi in ix:
    for yi in iy:
        p = par[xi][yi]
        if (p == '') or (p == 'none'):
            pass
        else:
            syncube[:,xi,yi] = nh31.fitfunc(x,p.params)
            if 'B2' in p.params:
                syncube[:,xi,yi] += nh32.fitfunc(x,p.params)
            

                         
                 
def newwcs(w,factor=1):
    wint = wcs.WCS(naxis=4)
    wint.wcs.crpix = w.wcs.crpix/factor
    wint.wcs.cdelt = w.wcs.cdelt*factor
    wint.wcs.radesys = w.wcs.radesys
    wint.wcs.mjdobs = w.wcs.mjdobs
    wint.wcs.dateobs = w.wcs.dateobs
    wint.wcs.crval = w.wcs.crval
    wint.wcs.ctype = w.wcs.ctype
    wint._naxis1 = w._naxis1/factor
    wint._naxis2 = w._naxis2/factor
    wint.cunit = w.wcs.cunit
    return wint 
    
    
for ps in par:
    for p in ps:
        if p.params['dv1'].value <  p.params['dv2'].value:
            dv1 = p.params['dv1'].value          
            
            
            
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
                if 'N1' in p.params:
                    Nu1[xi,yi,0] = p.params['N1']
                else: Nu1[xi,yi,0] = p.params['Nu1']
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

                    if 'N2' in p.params:
                        Nu1[xi,yi,1] = p.params['N2']
                    else: Nu1[xi,yi,0] = p.params['Nu1']
                    sigma1[xi,yi,1] = p.params['sigma2']
                    df1[xi,yi,1] = p.params['dv2']
    
    tau[where(tau==0)]=ma.masked
    inte[where(inte==0)]=ma.masked
    edf1[where(inte==0)]=ma.masked

    return (inte,tau,Nu1,sigma1,df1,B1,edf1,esigma1)

def plotOnlySpecMap(cube,cube2,xy,shape):
    #x=linspace(nh3bg.shape[0]-1,0,nptsx,dtype=int)  
    #y=linspace(0,nh3bg.shape[1]-1,nptsy,dtype=int)  
    six=arange(xy[0]-shape[0],xy[0]+shape[0])
    siy=arange(xy[1]-shape[1],xy[1]+shape[1])
    f = figure()
    import matplotlib.gridspec as gridspec
    gs = gridspec.GridSpec(len(six), len(siy))
    gs.update(wspace=0.00,hspace=0)
    axes = [[0 for yi in six] for xi in siy]
    for xi in six:
        for yi in siy:
            ind = (where(six==xi)[0][0]),(where(siy==yi)[0][0])
            ax = subplot(gs[ind])
            axes[ind[1]][ind[0]] = ax
#            fx = int(tox[xi*3])
#            fy = int(toy[yi*3])
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)
    #        ax.set_ylim(cube[:,ix.min():ix.max(),iy.min():iy.max()].min()/2,cube[:,ix.min():ix.max(),iy.min():iy.max()].max())
            ax.set_ylim(cube[:,ix.min():ix.max(),iy.min():iy.max()].min()/4,cube[:,ix.min():ix.max(),iy.min():iy.max()].max())
            ax.set_xlim(-10000,20000)
            spec = cube[:,xi,yi]
            step(vel21,spec,'k',linewidth=1)

            spec2 = cube2[:,xi,yi]
            step(velx,3*spec2,'k',linewidth=2)
#            spec2 = cube2[:,fx,fy]
#            step(vel21,spec2,'k')
            

def plotSubSpecMap(cube,cube2,par,ix,iy,velx,vel,mol1,mol2):
    six=ix
    siy=iy
        
    figure()
    cla()
    import matplotlib.gridspec as gridspec
    gs = gridspec.GridSpec(len(six), len(siy))
    gs.update(wspace=0.00,hspace=0)
    axes = [[0 for yi in six] for xi in siy]
    for xi in six:
        for yi in siy:
            ind = (where(six==xi)[0][0]),(where(siy==yi)[0][0])
            ax = subplot(gs[ind])
            axes[ind[1]][ind[0]] = ax
#            fx = int(tox[xi*3])
#            fy = int(toy[yi*3])
#            ax.text(0,0,"({0}:{1}):{2}".format(xi,yi,ind),verticalalignment='top')
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)
    #        ax.set_ylim(cube[:,ix.min():ix.max(),iy.min():iy.max()].min()/2,cube[:,ix.min():ix.max(),iy.min():iy.max()].max())
            ax.set_ylim(cube[:,ix.min():ix.max(),iy.min():iy.max()].min()/4,cube[:,ix.min():ix.max(),iy.min():iy.max()].max())
            ax.set_xlim(-13000,37000)
            spec = cube[:,xi,yi]
            step(velx,spec,'k',linewidth=2)
#            spec2 = cube2[:,fx,fy]
#            step(vel21,spec2,'k')
            
            plot([mol1.v0[10]+params['dv2'].value,mol1.v0[10]+params['dv2'].value],[0,5],'-b',linewidth=0.3)
            plot([mol1.v0[10]+params['dv1'].value,mol1.v0[10]+params['dv1'].value],[0,5],'-r',linewidth=0.3)
    #        err = spec[where(logical_not(maks))]. std()*3
    #        plot([velx[0],velx[-1]],[err,err],'-y',linewidth=0.3)
            #        plot(xxx+(xxx[1]-xxx[0])/2,bases[:,xi,yi],'-g')
    #        ax.fill_between(velx+(velx[1]-velx[0])/2,maks*cube.max()/2,facecolor='red',alpha=0.3)
            if (not par == '') and par[xi][yi] != '' and par[xi][yi] != 'none':
#                plot(vel,nh31.fitfunc(xx,par[xi][yi].params),'r',linewidth=2)
#                plot(vel,nh32.fitfunc(xx,par[xi][yi].params),'b',linewidth=2)
                plot(vel,mol1.fitfunc(vel,par[xi][yi].params),'--r')
                plot(vel,mol2.fitfunc(vel,par[xi][yi].params),'--b')
                if par2 != '':
                    plot(vel,mol1.fitfunc(vel,par2[xi][yi].params),'r')
                    plot(vel,mol2.fitfunc(vel,par2[xi][yi].params),'b')


for pn in pos:
    xi = pn[0]
    yi = pn[1]
    p = par2h[xi][yi].params
    spec = cube2[:,xi,yi]
    fx = int(frx[xi])
    fy = int(fry[yi])
    
    figure()
    step(x,cube1[:,fx,fy]*3,'k',linewidth=2)
    step(vel21,cube2[:,xi,yi],'k',linewidth=1)
    plot(vel21,n2h2.fitfunc(vel21,p),'-b')
    plot(vel21,n2h1.fitfunc(vel21,p),'-r')
    plot(vel21,n2h2.fitfunc(vel21,p)+n2h1.fitfunc(vel21,p),'-g')
    plot(vel21,(nh31.fitfunc(vel21,par[fx][fy].params)+nh32.fitfunc(vel21,par[fx][fy].params))*3,'--g',linewidth=2)
    legend(['$NH_3$','N_2H^+','"red" fit','"blue" fit', '$N_2H^+ fit$','$NH_3$ fit'])
    xlabel('$V,m/s$')
    ylabel("$T_b,K$")
    rc = w2.all_pix2world([[128,128,0,0]],1,ra_dec_order=True)
    off = (rc - w2.all_pix2world([[yi,xi,0,0]],1,ra_dec_order=True))*u.deg.to(u.arcsec)
    title('Fit at RA:{0:0.1f}, DEC:{1:0.1f}"'.format(off[0][0],off[0][1]))
    xlim(-12756,29406)


pos = []
def onclick(event):
    xc = int(event.ydata)
    yc = int(event.xdata)
    pos.append([xc,yc])

for xi in ix:
    for yi in iy:
        if par[xi][yi] == '' or par[xi][yi] == 'none':
            continue
        if nh31.fitfunc(x,par[xi][yi].params).max() < st*2.5:
            plot(yi,xi,'xw')
            print('\r({0},{1}-),'.format(xi,yi),end="")
            par[xi][yi] = 'none'

for xi in ix:
    for yi in iy:
        if par[xi][yi] == '' or par[xi][yi] == 'none':
            continue
        p = par[xi][yi].params
        spec = cube1[:,xi,yi]
        st = spec.std()
        if 'dv2' in p:
            if res2(p,x,spec,[]).std() > st*1.5: par[xi][yi] = 'none'
        else:
            if res1(p,x,spec,[]).std() > st*1.5: par[xi][yi] = 'none'


f = interp2d(frx, fry, img, kind='cubic')
for xi in ix:
    for yi in iy:
        plot(inte[xi,yi,0],f(xi,yi),'xk')
#        plot(cubenh[:,xi,yi].max(),cubenh2[:,xi,yi].max(),'xk')
        
for xi in ix:
    for yi in iy:
        plot(cubenh[:,xi,yi].max(),f(xi,yi),'x  ')

def createBeamArtist(rec,tr=None):
    from matplotlib.patches import Ellipse

    bmin = abs(rec.header['bmin']/rec.header['cdelt1'])
    bmaj = abs(rec.header['bmaj']/rec.header['cdelt1'])
    bpa = rec.header['bpa']
    if tr != None:
        circle1 = Ellipse(xy=(bmaj,bmin), width=bmaj, height=bmin, angle=bpa, transform=tr)
    else : circle1 = Ellipse(xy=(bmaj,bmin), width=bmaj, height=bmin, angle=bpa)
    return(circle1)



sincube = zeros((len(xx),nh3bg.shape[0],nh3bg.shape[1]))

for xi in ix:
    for yi in iy:
        if par[xi][yi] == '' or par[xi][yi] == 'none':
            continue
        pr = par[xi][yi].params
        sincube[:,xi,yi] =  nh31.fitfunc(xx,pr)