from astropy.io.fits.hdu.image import PrimaryHDU
from scipy.constants import c,h,k,pi
import numpy as np
from math import log, atan, sin, sqrt
from enum import Enum
from astropy.io import fits
from astropy import units as u;
from astropy.wcs import WCS
from lmfit import minimize, Parameters, Parameter, report_fit, report_errors

def loadFromFiles(url, mod=''):
    from os import listdir as ls
    """
    Load files from dir
    :arg url: url
    :arg mod: suffix (234232123[mod].fits)
    :returns: arry of fits datas
    """
    files = ls(url)
    filt = []
    for f in files:
        if f.endswith(mod+'.fits'): filt.append(f)
    files = filt
    recs = []
    for f in files:
        recs.append(f)
    return recs
    
def select(data, key, val):
    """
    filter fits array by given key and value
    :arg data: arraay to filter
    :arg key: to filter
    :arg val: tto filter
    :returns: result
    """
    key = key.upper()
    res = []
    for r in data:
        r = fits.open(r)
        if r[0].header[key].startswith(val): res.append(r)
        else: r.close()
    return res
    
def select2(data, key, val):
    """
    filter fits array by given key and value
    :arg data: arraay to filter
    :arg key: to filter
    :arg val: text to filter
    :returns: result
    """
    key = key.upper()
    res = []
    for r in data:
        v = r[0].header[key]
        if type(v) == str:
            if v.startswith(val): res.append(r)
            else: r.close()
        elif v == val: res.append(r)
        else: r.close()
    return res

def fold(data,df=20*u.MHz, orient=True):
    """
    fold spectra at data file
    :arg data: fits file to fold
    :arg df: frequency throw
    :arg orient: spike orientation, True - pos/neg, False - neg/pos 
    :returns: folded spectra in array
    """
    cs = abs(df.to(u.Hz)/(data[0].header['CDELT1']*u.Hz))
    return -data[0].data[:,:,cs:]+data[0].data[:,:,:-cs]
    
def plotPositions(data):
    for r in data:
        w = WCS(r[0].header)
        coord = w.wcs_pix2world(np.array([[0,0,0]]),0,ra_dec_order=True)
        plot(coord[0][0],coord[0][1],'xr')
        
def plotSpectralMap(data):
    import matplotlib.pyplot as plt
    fig = plt.figure();
    ra,deg = [],[]
    for r in data:
        w = WCS(r[0].header)
        coord = w.wcs_pix2world(np.array([[0,0,0]]),0,ra_dec_order=True)
        ra.append(coord[0][0])
        deg.append(coord[0][1])
    ra0,deg0 = min(ra),min(deg)
    size = (max(ra)-min(ra),max(deg)-min(deg))
    axr = np.unique(np.array(ra))
    axd = np.unique(np.array(deg))
    axr.sort()
    axd.sort()
    deltar = size[0]/(len(axr)-1)
    deltad = size[1]/(len(axd)-1)
    arr = [[[] for i in range(len(axd))] for i in range(len(axr))]
    for i in range(len(ra)):
        ax = plt.subplot2grid((len(axr),len(axd)),(int(round((ra[i]-ra0)/deltar)),int(round((deg[i]-deg0)/deltad))))
        ax.get_xaxis().set_ticks([])
        ax.get_yaxis().set_ticks([])
        ax.plot(fold(data[i])[0,0,:],alpha = 0.5)
        arr[np.where(axr == ra[i])[0][0]][np.where(axd == deg[i])[0][0]].append(data[i])
    fig.subplots_adjust(wspace=0)
    fig.subplots_adjust(hspace=0)
    return np.array(arr)

def constructCube(arr,spec):
    import numpy.ma as ma
    cube = ma.array(np.zeros((arr.shape[0],arr.shape[1],spec)))
    systemp = []
    
    for (x,y),val in np.ndenumerate(arr):
        if len(val)== 0: cube[x,y] = ma.masked
        else:
            d,w = [],[]
            for dat in val: 
                d.append(fold(dat)[0][0])
                w.append(dat[0].header['tsys']*dat[0].header['inttime'])
                systemp.append(w[-1])
            d=np.array(d)
            w=np.array(w)
            cube[x,y,:d.shape[1]] = np.average(d,axis=0,weights=w)
            
    rec = 0
    rec2= 1
    for (x,y),val in np.ndenumerate(arr):
        if len(val)!= 0:  
            rec = val[0][0]
            rec2 = arr[x+1,y+1][0][0]
            break
        
    w = WCS(rec.header)
    w2 = WCS(rec2.header)
    coord2 = w2.wcs_pix2world(np.array([[0,0,0]]),0,ra_dec_order=True)
    coord = w.wcs_pix2world(np.array([[0,0,0]]),0,ra_dec_order=True)
    
    hdu = fits.PrimaryHDU(np.array(cube,dtype=np.float32))
    hdulist = fits.HDUList([hdu])
    head = hdu.header
    head['CTYPE2'] = 'RA'
    head['CRVAL2'] = coord[0][0]
    head['CDELT2'] = coord2[0][0]-coord[0][0]
    head['CRPIX2'] = x+1
    head['CTYPE3'] = 'DEC'
    head['CRVAL3'] = coord[0][1]
    head['CDELT3'] = coord2[0][1]-coord[0][1]
    head['CRPIX3'] = y+1
    skip = ['CTYPE2','CRVAL2','CDELT2','CRPIX2','CTYPE3','CRVAL3','CDELT3','CRPIX3','SIMPLE','BITPIX','NAXIS','NAXIS1','NAXIS2','NAXIS3','EXTEND','COMMENT','HISTORY']
    for k in rec.header:
        if k not in skip:
            head.append(k,rec.header[k])
            head[k]= rec.header[k]
            head.comments[k] = rec.header.comments[k]
    return hdulist
    

def baseline(xy, func, step):
    
    '''
    Iteratically subtract baseline from 1-D array
    
    param xy 1-D array to fit
    '''    
    xy = xy.copy()
    from math import floor
    model = np.zeros(xy.shape)    
    
    if len(xy.shape)!=1:
        raise Exception('xy not 1-D array')
    
    def residual(params, x, data, error):
        c = []
        for key in params:
            c.append(params[key].value)
        model = func(x,c)
        return (data-model)
    params = Parameters()
    params.add('a', value=1)
    params.add('b', value=0)
    params.add('c', value=0)
    params.add('d', value=0)
    params.add('e', value=0)
    params.add('f', value=0)
    for st in np.arange(0,xy.shape[0],step):
        end = st+step
#        print('iterating from [{0} {1}]'.format(st,end))
        if (st > 0):
            data = xy[st-step/2:end+step/2]
        data = xy[st:end]
        x = np.linspace(0,len(data)-1,len(data))
        minn = minimize(residual, params, args=(x,data,[]))
#        report_fit(minn.params)
        c = []
        for key in minn.params:
            c.append(minn.params[key].value)
        model[st:end]=func(x,c)
    return model
        
    
    
    
def fitfunc(x,c):
    
    y = np.zeros(x.shape)
    c = c[:4]
    for power,par in enumerate(c):
         y+=par*x**power
    return y
    
def fitgaus(params,x,constants):
    w0 = params['w0'].value
    I  = params['I'].value
    sigma = params['sigma'].value
    k = params['k'].value
    b = params['b'].value
    d = -0.5*(((x-w0)/sigma)**2)
    aa = np.exp(d)*I+k*(x-w0)+b
    return aa
    
def resgaus(params, xr, yd, constants):
    model = fitgaus(params,xr,constants)
    return (yd-model)

def removeBaselines(cube):
    def plfunc(w0,I,sigma,k,b,x):
        d = -0.5*(((x-w0)/sigma)**2)
        aa = np.exp(d)*I+k*(x-w0)+b
        return aa
    
    constants = np.array([0])
    for (x,y),val in np.ndenumerate(cube[:,:,1]):
        spectra = cube[x,y]
        
        frist = spectra-baseline(spectra,fitfunc,600)
        sigma = np.std(frist)
        arrs = []
        ca = []
        pos = np.where(abs(frist)>1.5*sigma)[0]
        for i in range(len(pos)-1):
            ca.append(pos[i])
            if pos[i+1] - pos[i]!=1:
                if len(ca) > 1:
                    arrs.append(ca)
                    ca = []
        arrs = np.array(arrs)
        ids = []
        for i in range(len(arrs)):
            if len(arrs[i]) < 3:
                ids.append(i)
        arrs = np.delete(arrs,ids)

        w = []
        s = []
        I = []
        b = []
        k = []
        for ran in arrs:
            f = round(np.array(ran).mean())
            if (abs(frist[f])>1.5*sigma) and f > 200 and f < len(spectra)-200:
                params = Parameters()
                params.add('w0',f,min=f-10,max=f+10)
                params.add('I',spectra.max(),min = -10, max = 10)
                params.add('sigma', 1, min=0.01, max = 10)
                params.add('k',0)
                params.add('b',1)
                xr = np.arange(f-100,f+100,1,dtype=int)
                minimize(resgaus, params, args=(xr, spectra[f-100:f+100],constants))
#                print('({0},{1}):{2:.0f}'.format(x,y,x))
                report_fit(params)
                if params['w0'].value not in w:
                    w.append(params['w0'].value)
                    s.append(params['sigma'].value)
                    I.append(params['I'].value)
                    k.append(params['k'].value)
                    b.append(params['b'].value)
        spectra2 = spectra.copy()
#        figure()
#        plot(spectra,'-r', alpha = 0.2)
        for i in range(len(w)):
            xr = np.arange(w[i]-s[i]*3,w[i]+s[i]*3,1,dtype=int)
            model = plfunc(w[i],I[i],s[i],0,0,xr)
            if (spectra2[xr]-model).std() < spectra2[xr].std():
                plot(xr,plfunc(w[i],I[i],s[i],k[i],b[i],xr),'-g')
                spectra2[xr]-=model
#        plot(spectra2,'-b',alpha = 0.2)
        cube[x,y]-=baseline(spectra2,fitfunc,320)

def removeBaselinesWave(cube):
    from statsmodels.robust import stand_mad
    from scipy import stats
    import numpy as np
    import pywt
    
    def plfunc(w0,I,sigma,k,b,x):
        d = -0.5*(((x-w0)/sigma)**2)
        aa = np.exp(d)*I+k*(x-w0)+b
        return aa
    
    constants = np.array([0])
    for (x,y),val in np.ndenumerate(cube[1,:,:]):
        print(x,y)
        spectra = cube[:,x,y]
        
        frist = spectra-baseline(spectra,fitfunc,600)
        sigma = np.std(frist)
        arrs = []
        ca = []
        pos = np.where(abs(frist)>1.5*sigma)[0]
        for i in range(len(pos)-1):
            ca.append(pos[i])
            if pos[i+1] - pos[i]!=1:
                if len(ca) > 1:
                    arrs.append(ca)
                    ca = []
        arrs = np.array(arrs)
        ids = []
        for i in range(len(arrs)):
            if len(arrs[i]) < 3:
                ids.append(i)
        arrs = np.delete(arrs,ids)

        w = []
        s = []
        I = []
        b = []
        k = []
        for ran in arrs:
            f = round(np.array(ran).mean())
            if (abs(frist[f])>1.5*sigma) and f > 200 and f < len(spectra)-200:
                params = Parameters()
                params.add('w0',f,min=f-10,max=f+10)
                params.add('I',spectra.max(),min = -10, max = 10)
                params.add('sigma', 1, min=0.01, max = 10)
                params.add('k',0)
                params.add('b',1)
                xr = np.arange(f-100,f+100,1,dtype=int)
                minn = minimize(resgaus, params, args=(xr, spectra[f-100:f+100],constants))
#                print('({0},{1}):{2:.0f}'.format(x,y,x))
#                report_fit(minn.params)
                if minn.params['w0'].value not in w:
                    w.append(minn.params['w0'].value)
                    s.append(minn.params['sigma'].value)
                    I.append(minn.params['I'].value)
                    k.append(minn.params['k'].value)
                    b.append(minn.params['b'].value)
        spectra2 = spectra.copy()
        for i in range(len(w)):
            xr = np.arange(w[i]-s[i]*3,w[i]+s[i]*3,1,dtype=int)
            model = plfunc(w[i],I[i],s[i],0,0,xr)
            if (spectra2[xr]-model).std() < spectra2[xr].std():
#                plot(xr,plfunc(w[i],I[i],s[i],k[i],b[i],xr),'-g',linewidth=2)
                spectra2[xr]-=model

        noisy_coefs = pywt.wavedec(spectra2, 'db8', level=11, mode='per')
        sigma = stand_mad(noisy_coefs[-1])
        uthresh = 2*sigma*np.sqrt(2*np.log(len(spectra2)))
        denoised = noisy_coefs[:]
        
        denoised[1:] = (pywt.threshold(i,uthresh) for i in denoised[1:])
        
        signal = pywt.waverec(denoised, 'db8', mode='per')
        cube[:,x,y] = spectra-signal
                
                
def readCat(catalog):
    data = np.genfromtxt(catalog,usecols=(0,2,4))
    cat_freqs = (data[:,0]*u.MHz).to(u.Hz).value
    cat_SM2 = data[:,1]
    cat_En = data[:,2]
    arr = open(catalog).read().split('\n')
    cat_J = np.zeros(cat_SM2.shape)
    cat_K = np.zeros(cat_SM2.shape)
    i = 0
    for s in arr:
        if (len(s) == 0): break
        cat_J[i] = int(s[61:63].strip())
        cat_K[i] = int(s[63:65].strip())
        i+=1
    return(cat_freqs,cat_En,cat_J,cat_K)
                
def rot(integ, K, J, Ed,df,rms):
    from scipy.constants import c,h,k,pi,e
    from scipy.constants import k as ka
    from math import log,exp
    gk = 1 if K == 0 else 2
    gi =  0.5 if K % 3 == 0 else 0.25
    print(gi,gk,Ed,df,K,J,integ)
    Eu = df*h/ka+Ed*h/ka*c*100
    val = log(integ/gi/gk/(J**2-K**2))
    err = rms*gi*gk*(J**2-K**2)/integ
    return(val,Eu,err)
    
def resCH3C2H(params,x,y,freqs):
    model = fitCH3C2H(params,x,freqs)
    return (y-model)


def fitCH3C2H(params,x,freqs):
    w0 = params['w0'].value
    I = []
    for ch in ['0','1','2','3','4','5','6','7','8','9']:
        i = params.get('I'+ch)
        if (i == None): break
        else: I.append(i.value)
        
    
    sigma = params['sigma'].value
    y = 0
    for i in range(len(freqs)):
        f = freqs[i]
        ex = -0.5*(((x-f-w0)/sigma)**2)
        y+=I[i]*np.exp(ex)
    return y

def diag(subf,subs,cat_freqs,cat_K,cat_J,cat_En):
    '''
    vdopler shift!
    '''
    cat_freqs = cat_freqs.copy()
    xy = []
    constants = []
    params = Parameters()
    params.add('w0',0.0105,min=-0.01,max=0.04)
    params.add('I0',0.050,min = 0.01, max = 10)
    params.add('I1',0.050,min = 0.01, max = 10)
    params.add('I2',0.050,min = 0.01, max = 10)
    params.add('I3',0.050,min = 0.01, max = 10)
    params.add('I4',0.050,min = 0.01, max = 10)
    params.add('sigma', 0.0001, min=0.00001, max = 0.01)
    minm = minimize(resCH3C2H, params, args=(subf, subs,cat_freqs/1e9))
    lines = fitCH3C2H(minm.params,subf,cat_freqs/1e9)
    report_fit(minm.params)
    I = []
    for i in range(len(cat_freqs)):
        I.append(minm.params.get('I{0}'.format(i)).value)
    for i in range(len(cat_freqs)):
        f = cat_freqs[i]/1e9
        par = minm.params.copy()
        for j in range(len(cat_freqs)):
            par['I{0}'.format(j)].value = 0
        par['I{0}'.format(i)].value = I[i]
        integ = ((fitCH3C2H(minm.params,subf,cat_freqs/1e9)>0)*fitCH3C2H(minm.params,subf,cat_freqs/1e9))
        intl,eu,error = rot(integ.sum(),cat_K[i],cat_J[i],cat_En[i],cat_freqs[i],(subs-lines).std())
        xy.append((intl,eu,error))
    xy = np.array(xy)
    xy = np.flipud(xy)
    fitfunc = lambda x,k,b: -x/k+b        
    def residual(params, x, yd, error):
        k = params['k'].value
        b = params['b'].value
        model = fitfunc(x,k,b)
        return (yd-model)
    
    params = Parameters()
    params.add('k', value=40)
    params.add('b', value=-0.8)
    minm = minimize(residual, params, args=(xy[:-1,1],xy[:-1,0],[]))
    report_fit(minm.params)
    x = xy[:,1]
    y = xy[:,0]
    dy = xy[:,2]
    Dx = len(x)*((x**2).sum())+(x.sum())**2
    su = 0
    for i in range(len(dy)):
       su+=(dy[i]*(x[i]-x.sum())/Dx )**2
    dB = math.sqrt(su)
    
    return (minm.params['k'].value.dB,xy)
    
'''    
recs = loadFromFiles('.','m')
recs = select(recs,'object','77')
recs = select2(recs,'line','h13')



arr = plotSpectralMap(recs)
cube = constructCube(arr,32768)
spectra = cube.mean(axis=0).mean(axis=0)
'''

def loadSingle(sou,line):
    recs = loadFromFiles('.','s')
    recs = select(recs,'object',sou)
    recs = select2(recs,'line',line)
    recs = select2(recs,'OBSMODE','DSW SIG')
    ra,deg = [],[]
    for r in recs:
        w = WCS(r[0].header)
        coord = w.wcs_pix2world(np.array([[0,0,0]]),0,ra_dec_order=True)
        ra.append(coord[0][0])
        deg.append(coord[0][1])
    
    axr = np.unique(np.array(ra))
    axd = np.unique(np.array(deg))
    axr.sort()
    axd.sort()
    arr = [[[] for i in range(len(axd))] for i in range(len(axr))]
    for i in range(len(ra)):
        arr[np.where(axr == ra[i])[0][0]][np.where(axd == deg[i])[0][0]].append(recs[i])
    
    arr = np.array(arr)
    import numpy.ma as ma
    spec = 32768
    cube = ma.array(np.zeros((arr.shape[0],arr.shape[1],spec)))
    systemp = []
    
    spec = 32768
    systemp = []
    
    x = 0
    y = 0
    val = arr[x,y]
    d,w = [],[]
    
    d,w = [],[]
    for dat in val: 
        d.append(dat[0].data)
        w.append(dat[0].header['tsys']*dat[0].header['inttime'])
        systemp.append(w[-1])
    
    d=np.array(d)
    w=np.array(w)
    cube[x,y] = np.average(d,axis=0,weights=w)
    spec = cube[0,0]
    return spec


rec = True
fl = False
ran = []
xs = []
def onclick(event):
    if rec:
        ran.append(event.xdata)
        plot([event.xdata,event.xdata],[spec.min(),spec.max()],'-g')
        if len(ran) == 2:
            xs.append(ran.copy())
            ran.clear()
            


cid = fig.canvas.mpl_connect('button_press_event', onclick)

86552000000.0