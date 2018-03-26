
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 15:02:05 2014

@author: pete
"""
from astropy import units as u;
from astropy.io.fits.hdu.image import PrimaryHDU
from scipy.constants import c,h,k,pi
import numpy as np
from lmfit import minimize, Parameters, Parameter, report_fit, report_errors
from math import log, atan, sin
from matplotlib.widgets import Button
from enum import Enum
from astropy.io import fits

from matplotlib.gridspec import GridSpec
from matplotlib.pyplot import gcf, draw_if_interactive, delaxes

def getVel(head):
    chan = np.arange(0,head['NAXIS1'],1)
    return head['VELO-LSR']+(-head['CRPIX1']+chan)*head['DELTAV']

def mysubplot2grid(shape, loc, rowspan=1, colspan=1, **kwargs):
    fig = gcf()
    s1, s2 = shape
    subplotspec = GridSpec(s1, s2).new_subplotspec(loc,
    rowspan=rowspan,
    colspan=colspan)
    a = fig.add_subplot(subplotspec, **kwargs)
    bbox = a.bbox
    byebye = []
    for other in fig.axes:
        if other==a: continue
        if bbox.fully_overlaps(other.bbox):
            byebye.append(other)
    for ax in byebye: delaxes(ax)
    
    draw_if_interactive()
    return a


class MyButton(Button) :
    def set_label(self, label):
        self.label.remove()
        self.label = self.ax.text(0.5, 0.5, label,
                     verticalalignment='center',
                     horizontalalignment='center',
                     transform=self.ax.transAxes)
        self.ax.plot()
        

def rangexy(x,y):
    return ((x[0]-y[0])**2+(x[1]-y[1])**2)**0.5

def interpolate(image,axes,r):
    from scipy.ndimage import zoom
    x = np.linspace(0,axes[0].shape[0],axes[0].shape[0])
    xr = np.linspace(0,axes[0].shape[0],axes[0].shape[0]*r)
    y = np.linspace(0,axes[1].shape[0],axes[1].shape[0])
    yr = np.linspace(0,axes[1].shape[0],axes[1].shape[0]*r)
    return zoom(image,r),(np.interp(xr,x,axes[0]),np.interp(yr,y,axes[1]))


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
    params.add('a', value=0)
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
        minimize(residual, params, args=(x,data,[]))
#        report_fit(params)
        c = []
        for key in params:
            c.append(params[key].value)
        model[st:end]=func(x,c)
    return model
        

def fitfunc(x,c):
    
    y = np.zeros(x.shape)
    c = c[:4]
    for power,par in enumerate(c):
         y+=par*x**power
    return y
        
        
        


class Unit:
    def __init__(self, data):
        from astropy import wcs
        if (isinstance(data, PrimaryHDU)): #FITS file
            self.wcs = wcs.WCS(data.header)
            naxis = self.wcs.naxis
            self.slices = [0]*naxis
            self.slices[0] = 'x'
            self.slices[1] = 'y'
            if (data.header['ctype1'].startswith('FREQ')):
                print('ctype1')
                self.freq  = data.header['restfreq']
                self.crval = self._freq2vel(data.header['crval1'])
                self.cinc  = self._freq2vel(self.freq+data.header['cdelt1'])-self._freq2vel(self.freq)
                self.crpix = data.header['crpix1']
            elif (data.header['ctype3'].startswith('FREQ')):
                print('ctype3')
                self.freq  = data.header['crval3']
                self.crval = self._freq2vel(data.header['crval3'])
                self.cinc  = self._freq2vel(self.freq+data.header['cdelt3'])-self._freq2vel(self.freq)
                self.crpix = data.header['crpix3']
            
        else:
            self.crval = data.getArrayItem('crval3')[0]
            self.cinc = data.getArrayItem('cdelt3')[0]
            self.crpix = data.getArrayItem('crpix3')[0]
            self.freq = data.getArrayItem('restfreq')[0]*1e9
    
            naxis=data.getArrayItem('naxis')[0]
            self.slices = [0]*naxis
            self.slices[0] = 'y'
            self.slices[1] = 'x'
            w = wcs.WCS(naxis=naxis)
            crpix = []
            cdelt = []
            crval = []
            ctype = []
            for i in range(naxis):
                crpix.append(data.getArrayItem('crpix'+str(i+1))[0]*u.rad.to(u.deg))
                cdelt.append(data.getArrayItem('cdelt'+str(i+1))[0]*u.rad.to(u.deg))
                crval.append(data.getArrayItem('crval'+str(i+1))[0]*u.rad.to(u.deg))
                ctype.append(data.getArrayItem('ctype'+str(i+1)).decode())
            w.wcs.crpix = crpix
            w.wcs.cdelt = np.array(cdelt)
            w.wcs.crval = crval
            w.wcs.ctype = ctype
            self.wcs = w
            # Set up an "Airy's zenithal" projection
            # Vector properties may be set with Python lists, or Numpy arrays
                
        self.eq = [
            (u.chan,u.km/u.s,self._chanel2vel,self._vel2chanel),
            (u.chan,u.Hz,self._chanel2freq,self._freq2chanel),
        ]
        
    def register(self):
#        u.add_enabled_units([self.xpix,self.ypix,self.rradx,self.rrady,self.rsecx,self.rsecy])
        u.add_enabled_equivalencies(self.eq)        
        
    
    def yaxisImg(self):
        start = self.yval-self.yref*abs(self.ydelta)
        stop = start+self.ycount*self.ydelta
        return np.linspace(start,stop,self.ycount)
        
    def xaxisImg(self):
        start = (self.xval+self.xref*abs(self.xdelta))
        stop = start+self.xdelta*self.xcount
        return np.linspace(start,stop,self.xcount)

    def _chanel2vel (self,ch):
        return (ch*self.cinc+self.crval-self.crpix*self.cinc)
    
    def _chanel2freq (self,ch):
        from scipy.constants import c as lightspd
        cms = self._chanel2vel(ch)
        return self.freq*(1+cms/lightspd)
    
    def _freq2vel (self,f):
        from scipy.constants import c as lightspd
        return (1-(f/self.freq-1)*lightspd)*1e-3

    def _vel2freq (self,vel):
        from scipy.constants import c as lightspd
        return ((1-vel*1e3)/lightspd+1)*self.freq
    
    def _vel2chanel (self,vel):
        return (vel-self.crval+(self.crpix*self.cinc))/self.cinc
            
    
    def _freq2chanel (self,f):
        return self._vel2chanel( self._freq2vel(f))

class XUnit:
    
    xpix = u.def_unit('xpix')
    ypix = u.def_unit('ypix')    
    
    def __init__(self, data):
        self.ycount = data.getArrayItem("naxis2")[0]
        self.yref = data.getArrayItem("crpix2")[0]
        self.yval = data.getArrayItem("crval2")[0]
        self.ydelta = data.getArrayItem("cdelt2")[0]
        
        self.xcount = data.getArrayItem("naxis1")[0]
        self.xref = data.getArrayItem("crpix1")[0]
        self.xval = data.getArrayItem("crval1")[0]
        self.xdelta = data.getArrayItem("cdelt1")[0]

        xpix2sec = lambda x: self.xval-self.xref*self.xdelta+self.xdelta*x
        ypix2sec = lambda y: self.yval-self.yref*self.ydelta+self.ydelta*y

        sec2xpix = lambda x: (x-self.xval+self.xref*self.xdelta)/self.xdelta
        sec2ypix = lambda y: (y-self.yval+self.yref*self.ydelta)/self.ydelta

        scale = u.rad.to(u.arcsec)

        
        sec2relsecx = lambda y: (self.xval-y)*scale
        sec2relsecy = lambda y: (-self.yval+y)*scale


        pix2rsecx = lambda x: (x-self.xref)*(self.xdelta*scale)
        pix2rsecy = lambda y: (y-self.yref)*(self.ydelta*scale)

        rsecx2pix = lambda x: -x/self.xdelta/scale+self.xref
        rsecy2pix = lambda y: y/self.ydelta/scale+self.yref
        
        self.xpix = u.def_unit('Pixel X')
        self.ypix = u.def_unit('Pixel Y')

        self.rradx = u.def_unit('Relative Radians X')
        self.rrady = u.def_unit('Relative Radians Y')

        self.rsecx = u.def_unit(r'RA offset, "')
        self.rsecy = u.def_unit(r'DEC offset, "')

        self.eq = [
            (self.xpix,u.arcsec,xpix2sec,sec2xpix),
            (self.ypix,u.arcsec,ypix2sec,sec2ypix),
            (self.xpix,self.rsecx,pix2rsecx,rsecx2pix),
            (self.ypix,self.rsecy,pix2rsecy,rsecy2pix),
            (self.rsecx,u.arcsec,lambda x:x/scale+self.xval,sec2relsecx),
            (self.rsecy,u.arcsec,lambda y:y/scale+self.yval,sec2relsecy),
            (self.rradx,self.rsecx,lambda a:a/scale,lambda a:a*scale),
            (self.rrady,self.rsecy,lambda a:a/scale,lambda a:a*scale),
            (self.xpix,self.ypix,lambda a: a),
        ]
        

    def yaxisImg(self):
        start = self.yval-self.yref*abs(self.ydelta)
        stop = start+self.ycount*self.ydelta
        return np.linspace(start,stop,self.ycount)
        
    def xaxisImg(self):
        start = (self.xval+self.xref*abs(self.xdelta))
        stop = start+self.xdelta*self.xcount
        return np.linspace(start,stop,self.xcount)


    def register(self):
        for i in u.get_current_unit_registry().all_units:
            if ('x_pixels' in i.names):
                return

        u.add_enabled_units([self.xpix,self.ypix])
        u.add_enabled_equivalencies(self.eq)        
        
    
    def yaxisImg(self):
        return np.arange(self.yval-self.yref*self.ydelta,self.yval+(self.ycount-self.yref)*self.ydelta,self.ydelta)
        
    def xaxisImg(self):
        return np.arangce(self.xval-self.xref*self.xdelta,self.xval+(self.xcount-self.xref)*self.xdelta,self.xdelta)


def resWithMask(f,mask):
    return lambda p,x,y,c: mask*f(p,x,y,c)
    
import matplotlib.pyplot as plt

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx

h = 6.626e-27
c = 29979245800.0

k = 1.38e-16
#ba = BMAJ*BMIN*2*pi/2.355**2
#eq = u.brightness_temperature(ba,nh31.restfreq*u.Hz)

class Analyse:

    def toTb(cube, header):   
        r'''
        Convert Jy/beam scale to K
        
        Parameters
        -------------
        
        cube : array_like
            data to convert
        header : astropy.io.fits.header.Header
            with "bmaj" and "bmin" keys, in degrees.
        
        Returns
        ------------
        scaled data
        
        '''
        xar = header['bmaj']*u.deg
        yar = header['bmin']*u.deg
        beam_area = ((2*np.pi/2.355/2.355)*(xar*yar))
        return cube*(u.Jy).to(u.K, equivalencies=u.brightness_temperature(beam_area,header['restfreq']*u.Hz))
        

    
    def Jy_beam_toCGS_rad(data,Bmaj,Bmin):
        return ((data*u.Jy/(2*pi/2.355/2.355*Bmaj*Bmin)).cgs).value
    
    class MultiGaussModel:
        indexes = ['1','2','3','4','5','6','7','8','9','0','a','b','c','d','e','f']
        class GaussianSegment:
            def __init__(self,freq,eq,SM2,J,En,index='1'):
                self.freq=freq
                self.eq=eq
                self.index=index
                self.SM2 = SM2
                self.J = J
                self.En = En

                
            def __call__(self,x,p):
                return self.fitfunc(x,p)
            
            def fitfunc(self,p : Parameters,x):
#                print(p)
                i = self.index
                dv = p['dv'].value
                I  = p['I'+i].value
                sigma = p['sigma'].value
                
                b = -0.5*(((x-(self.freq*u.Hz).to(u.m/u.s, equivalencies=self.eq).value-dv)/sigma)**2)
                aa = np.exp(b)*I
                return aa


        def __init__(self,rec):
            self.segments = []
            self.bmaj = rec.header['bmaj']*u.deg.to(u.rad)
            self.bmin = rec.header['bmin']*u.deg.to(u.rad)
            self.par = Parameters()
            self.eq = u.doppler_radio(rec.header['restfreq']*u.Hz)
            

            
        def appendGauss(self, segment : GaussianSegment):
            self.segments.append(segment)
        
        def createGauss(self,catalog,index):
            i = index
            data = np.genfromtxt(catalog,usecols=(0,2,4))
            freq = (data[i,0]*u.MHz).to(u.Hz).value
            SM2 = data[i,1]
            En = data[i,2]
            
            arr = open(catalog).read().split('\n')
            arr.remove(arr[-1])
            J = int(arr[i][61:63].strip())
            s = Analyse.MultiGaussModel.GaussianSegment(freq,eq,SM2,J,En,Analyse.MultiGaussModel.indexes[i])
            self.segments.append(s)
            
        def fit(self,spec,v,par=None):
            if (par == None):
                par = Parameters()
                par.add('dv',8e3,min=7e3,max=11e3)
                par.add('sigma',0.8e3,min=0.3e3,max=1.5e4)
                for i in range(len(self.segments)):
                    par.add('I'+Analyse.MultiGaussModel.indexes[i],0.01,min=0,max=2)
                
            def residual(p,x,y):
                mspec =zeros(x.shape)
                for s in self.segments:
                    mspec+=s(p,x)
                return y-mspec
            
            #print(par)
            mini = Minimizer(residual,par,fcn_args=(v,spec))
            self.p1 = mini.minimize(method='nelder')
            self.std = self.p1.residual.std()    
            self.dv = abs(v[1]-v[0])
            return self.p1
            
        def calculate(self):
            ln = []
            energy = []
            freqs = []
            error = []
            v = linspace(-25e3,200e3,3000)
            for mol in self.segments:
                J = mol.J
                freq = mol.freq
                nu = freq
                down_energy = mol.En
                integral = mol(self.p1.params,v).sum()*(v[1]-v[0])
                N = ((integral*(c)**2)/(mol.SM2*self.bmaj*self.bmin*nu**3)*(3/16/2.5/1.13)*10e27)                
                print(N,J,freq,down_energy)
                if (N <= 0):
                    continue
                ln.append(log10(Nдщп1))
                freqs.append(freq)
                energy.append(freq*h/k+down_energy*h/k*c)
                error.append(self.std * (self.dv)**0.5/ln[-1])
            plt.plot(energy,ln,'or')
            plt.xlabel('$Energy, K$')
            plt.ylabel('$log(N_u/g_u)$')
#            for i in range(len(energy)):
#                plt.text(energy[i]-4,ln[i],"{0:0.3f}".format(freqs[i]/1e9), horizontalalignment='right',verticalalignment='top')
            arr = np.array((energy,ln,error)).T
            arr.shape
            self.arr = arr.copy()
            print(arr)
            sort = np.argsort(arr,axis=0)[:,0]
            sortt = np.zeros(arr.shape)
            sortt[:,0] = sort
            sortt[:,1] = sort
            sortt[:,2] = sort
            arr = arr[sort]        
            print(arr.shape)
            self.x = arr[:,0]
            self.y = arr[:,1]
            self.errors = arr[:,2]
            return self.x,self.y,self.errors
    
            
        def fitLine(self):
            fitfunc = lambda x,k,b: k*x+b        
            def residual(params, x, yd, error):
                k = params['k'].value
                b = params['b'].value
                model = fitfunc(x,k,b)
                return (yd-model)
            params = Parameters()
            params.add('k', value=-1)
            params.add('b', value=13)
            m = minimize(residual, params, args=(self.x,self.y,[]))
            params = m.params
            report_fit(params)
            plt.plot(self.x,self.y,'or')
            plt.errorbar(self.x,self.y,yerr=self.errors,fmt='or')
            xpoints = np.arange(min(self.x)-20,max(self.x)+20,3)
            plt.plot(xpoints,fitfunc(xpoints,params['k'],params['b']),'-b')
    #        print (1/(log(10)*k*params['k'].value))
            plt.text(xpoints.mean()+10,fitfunc(xpoints.mean(),params['k'],params['b']),
            'T$_{kin}$'+'={0:0.1f}'.format((-1/(log(10)*params['k'].value))))
        
    
    class N2HModel:
        
        def __init__(self,cat,ind='',restfreq=279511.760*1e6, J=3):
            """
            Create molecula model object
            
            Parameters
            ----------
            cat : array_like 
                data catalog from (Pagani et al, 2008)
            ind : str  
                to split paramters
            restfreq : float 
                frequency of transition. changed from 2.79512090e+11
                
                
            
            """            
            
            from numpy import exp

            self.cat_freqs = cat[:,6]*1e6
            self.cat_Ajk = cat[:,7]
            self.B0 = 46586.88e6
            self.Qrot = lambda Tex:k*Tex/h/self.B0+1/3 ###Mangum & Shirley 2015
            self.gu = (2*cat[:,2]+1)
            self.restfreq = restfreq
            self.v0 = (self.cat_freqs*u.Hz).to(u.m/u.s, equivalencies=u.doppler_radio(self.restfreq*u.Hz)).value
            self.index = ind
            self.Bbg = h*self.restfreq/k/(exp(h*self.restfreq/k/2.72)-1)
            self.B = lambda T: h*self.restfreq/k*(1/(exp(h*self.restfreq/k/T)-1))
            self.J = J
            self.Eu = h*self.B0*self.J*(self.J+1)
            
        def __call__(self,x,p):
            return self.fitfunc(x,p)

            
            
        def fitfunc(self,x,p : Parameters):
            r"""
            modeled specra in velocity axis
            
            :math:`T_{b}[K] = (B^p_{ind}+ \frac{h \nu} {k} \frac{1}{e^{h\nu / T_{bg}}-1})e^{-\tau_\nu}`
            
            :math:`\tau_\nu = \sum_{i=0}^{n_{hfs}} \frac{h\nu}{4\pi}\frac{A_{i_{ul}} N^p_{tot} Q_{rot}} {B^p_{ind} e^{Eu/kT}}\phi_i(v)`
            
            :math:`\phi_i(v) = \frac{1}{\sigma\sqrt(2\pi)}e^{-(x-v_{i}-v^p_0)^2/2\sigma^2}`

            :math:`T_{bg} = 2.72`
            
            :math:`B(T) = \frac{h \nu} {k} \frac{1} {e^{h \nu / kT}-1}`

            Parameters
            ----------
            x : array_like
                velocity axis
            p : lmfit.parameter.Parametrs
                Should contain: B,sigma,dv,Ntot followed by ``self.index`` (``'B1'``, for example)
                
                
            Returns:
            ----------
                modeled spectra
            """
            from numpy import exp

            return (p['B'+self.index].value-self.Bbg)*(1-exp(-self._tau(x,p)))
        
        def _tau(self,x,p):
            from numpy import exp

            sigma = p['sigma'+self.index].value
            dv = p['dv'+self.index].value
            B = p['B'+self.index].value
            N = p['N'+self.index].value
            tau = np.zeros(x.shape)
            for i in range(len(self.cat_freqs)):
                phi = exp(-(x-self.v0[i]-dv)**2/2/sigma**2)/(sigma*(2*pi)**0.5)
                tau += self.cat_Ajk[i]*N*phi
            return tau

    class N2HVaryModel:
        
        def __init__(self,cat,ind='',restfreq=2.79512090e+11, J=3):
            """
            Create molecula model object
            
            Parameters
            ----------
            cat : array_like 
                data catalog from (Pagani et al, 2008)
            ind : str  
                to split paramters
            restfreq : float 
                frequency of transition
            """            
            
            from numpy import exp

            self.cat_freqs = cat[:,6]*1e6
            self.cat_Ajk = cat[:,7]
            self.B0 = 46586.88e6
            self.Qrot = lambda Tex:k*Tex/h/self.B0+1/3 ###Mangum & Shirley 2015
            self.gu = (2*cat[:,2]+1)
            self.restfreq = restfreq
            self.v0 = (self.cat_freqs*u.Hz).to(u.m/u.s, equivalencies=u.doppler_radio(self.restfreq*u.Hz)).value
            self.index = ind
            self.Bbg = h*self.restfreq/k/(exp(h*self.restfreq/k/2.72)-1)
            self.B = lambda T: h*self.restfreq/k*(1/(exp(h*self.restfreq/k/T)-1))
            self.J = J
            self.Eu = h*self.B0*self.J*(self.J+1)
            
            
        def fitfunc(self,x,p : Parameters):
            r"""
            modeled specra in velocity axis
            
            :math:`T_{b}[K] = (B^p_{ind}+ \frac{h \nu} {k} \frac{1}{e^{h\nu / T_{bg}}-1})e^{-\tau_\nu}`
            
            :math:`\tau_\nu = \sum_{i=0}^{n_{hfs}} \frac{h\nu}{4\pi\nu_{i_{group}}}\frac{A_{i_{ul}} N^p_{tot} Q_{rot}} {B^p_{ind} e^{Eu/kT}}\phi_i(v)`
            
            :math:`\phi_i(v) = \frac{1}{\sigma\sqrt(2\pi)}e^{-(x-v_{i}-v^p_0)^2/2\sigma^2}`

            :math:`T_{bg} = 2.72`
            
            :math:`B(T) = \frac{h \nu} {k} \frac{1} {e^{h \nu / kT}-1}`

            Parameters
            ----------
            x : array_like
                velocity axis
            p : lmfit.parameter.Parametrs
                Should contain: B,sigma,dv,Ntot followed by ``self.index`` (``'B1'``, for example)
                
                
            Returns:
            ----------
                modeled spectra
            """
            from numpy import exp

            return (p['B'+self.index].value-self.Bbg)*(1-exp(-self._tau(x,p)))
        
        def _tau(self,x,p):
            from numpy import exp

            sigma = p['sigma'+self.index].value
            dv = p['dv'+self.index].value
            nu1 = p['nu1'+self.index].value
            nu2 = p['nu2'+self.index].value
            N = p['N'+self.index].value
            
            tau = np.zeros(x.shape)
            for i in range(len(self.cat_freqs)):
                nu = (nu1 if i<7 else nu2 if i > 30 else 1)
                phi = exp(-(x-self.v0[i]-dv)**2/2/sigma**2)/(sigma*(2*pi)**0.5)
                tau += self.cat_Ajk[i]*N/nu*phi
            return tau

#                Nu = Ntot*self.Qrot(self.B(B))/self.gu[i]*exp(self.Eu/k/self.B(B))
 #              tau += (h*self.restfreq/4/pi)*self.cat_Ajk[i]/B*Nu*phi


    class SymtopModel:
            
        def __call__(self,x,p):            
            return self.fitfunc(x,p)+1
            
                
        def fitfunc(self,x,p):
            from numpy import exp
            return (p['B'+self.index].value-self.Bbg)*(1-exp(-self.tau(x,p)))
        
        def tau(self,x,p):
            sigma = p['sigma'+self.index].value
            dv = p['dv'+self.index].value
            B = p['B'+self.index].value
            Nu = p['Nu'+self.index].value
            tau = np.zeros(x.shape)
            for i in range(len(self.cat_freqs)):
                phi = exp(-(x-self.v0[i]-dv)**2/2/sigma**2)/(sigma*(2*pi)**0.5)
                tau += (16*(pi**3)*(self.restfreq**4)*self.mu**2*self.S*Nu*self.cat_g[i]/(3*c**3*B))*phi
            return tau
        
            

    class NH3Model(SymtopModel):
        def __init__(self,cat,ind='',restfreq=23694506002.0):
            self.cat_freqs = cat[:,0]
            self.cat_g = cat[:,1]
            self.restfreq = restfreq
            self.v0 = (self.cat_freqs*u.Hz).to(u.m/u.s, equivalencies=u.doppler_radio(self.restfreq*u.Hz)).value
            self.S = 1/2
            self.mu = 1.468e-18
            self.index = ind
            self.Bbg = h*self.restfreq/k/(exp(h*self.restfreq/k/2.72)-1)
            
            
        def Tex(self,other,p1,p2,vx):
            '''
            Calculation of the kinetic temperature based on (1-1) and (2-2) transitions
            
            see http://adsabs.harvard.edu/abs/1992ApJ...388..467M
            
            '''
            from numpy.ma import log,exp
            from math import isnan
            i = self.index

            Tb1,Tb2 = 0,0
            if self.restfreq < other.restfreq:
                Tb1 = self.fitfunc(vx,p1).max()
                Tb2 = other.fitfunc(vx,p2).max()
                tau1 = self.tau(vx,p1).max()
#                tau2 = other.tau(vx,p2).max()
                s1 = p1['sigma'+i].value*2e-3
                s2 = p2['sigma'+i].value*2e-3
            else :
                Tb2 = self.fitfunc(vx,p1).max()
                Tb1 = other.fitfunc(vx,p2).max()
#                   tau2 = self.tau(vx,p1).max()
                tau1 = other.tau(vx,p2).max()
                s2 = p1['sigma'+i].value*2e-3
                s1 = p2['sigma'+i].value*2e-3

            Tex = -41.5*(log(-0.283*s2/tau1/s1*log(1-Tb2/Tb1*(1-exp(-tau1)))))**-1
            
            
            def err(p,Tex):
                Tk= p['Tk'].value
                return (Tex*(1+Tk/41.5*np.log(1+0.8*np.exp(-21.5/Tk)))-Tk)**2
            
            if isnan(Tex):
                return -1
            p = Parameters()
            p.add('Tk',20,min=0.1,max=4000)
            par = minimize(err,p,args={Tex})
            print('Tk:{0}'.format(par.params['Tk'].value))
            return  par.params['Tk'].value
            

    class CH3CNModel(SymtopModel):
        def __init__(self,cat,ind='',restfreq=220747261.7):
            self.cat_freqs = cat[:,0]
            self.cat_g = cat[:,1]
            self.restfreq = restfreq
            self.v0 = (self.cat_freqs*u.Hz).to(u.m/u.s, equivalencies=u.doppler_radio(self.restfreq*u.Hz)).value
            self.S = 1/2
            self.mu = 3.92197e-18
            self.index = ind
            self.Bbg = h*self.restfreq/k/(exp(h*self.restfreq/k/2.72)-1)
            

    
    
    def fitImage(self, xyfdata, func, xy, start_params,constants):
        def cp(par):
            pms = Parameters()
            for n in par: pms.add(n,par[n].value,min=par[n].min,max=par[n].max)
            return pms
            
        if (hasattr(self,'results')):
            print("pre") 
            self.pre = self.results
        else :
            print("none") 
            self.pre = np.zeros((len(start_params),xyfdata.shape[1],xyfdata.shape[2]))
        
        self.results = [[0 for x in range(xyfdata.shape[2])] for x in range(xyfdata.shape[1])]
        self.process = np.zeros((xyfdata.shape[1],xyfdata.shape[2]))
        self.process[xy] = 1
        self.results[xy[0]][xy[0]] = start_params
        i = 1
        while (self.process == 1).any():
            i+=1
#            if i == 4: return
            while (self.process == 1).any():
                for (xi,yi),value in np.ndenumerate(xyfdata[0]):
                    if (self.process[xi,yi] == 1):
                        spectra = xyfdata[:,xi,yi]
                        std = np.std((spectra-fit(self.pre,xi,yi,np.linspace(0,xyfdata.shape[0],xyfdata.shape[0])))[40:70])
                        if (np.std(spectra)*3>spectra.max()):
                            self.process[xi,yi] = -1
                        else:
                            print("empty {0}, processing {1} at {2} with std {3}".format((self.process == 0).sum()/self.process.size,(self.process == 1).sum()/self.process.size, (xi,yi),std))
                            pm = self.results[xi][yi]
                            for key in start_params:
                                start_params[key].value = pm[key].value
                            if((std < 0.5) and (self.pre[0,xi,yi]!=0)):
                                print('old')
                                res_param = Parameters()
                                for param in start_params:
                                    res_param.add(param,0)
                            else :
                                res_param = self.fitPoint(spectra,start_params,func, constants)
#                            print(res_param)
                            self.results[xi][yi] = cp(res_param)
                            self.process[xi,yi] = i
                        if (yi < xyfdata.shape[2]-1 and self.process[xi,yi+1] == 0):
                            self.process[xi,  yi+1] = 0.5
                            self.results[xi][yi+1] = cp(res_param)
                        if (yi > 0 and self.process[xi,yi-1] == 0):
                            self.process[xi,  yi-1] = 0.5
                            self.results[xi][yi-1] = cp(res_param)
                        if (xi < xyfdata.shape[1]-1 and self.process[xi+1,yi] == 0):
                            self.process[xi+1, yi]  = 0.5
                            self.results[xi+1][yi] = cp(res_param)
                        if (xi > 0 and self.process[xi-1,yi] == 0):
                            self.process[xi-1,  yi] = 0.5
                            self.results[xi-1][yi] = cp(res_param)
                        #plt.imshow(self.process)
            for (xi,yi),value in np.ndenumerate(self.process):
                    if value == 0.5: self.process[xi,yi] = 1
            
                    
        results = np.zeros((len(start_params),xyfdata.shape[1],xyfdata.shape[2]))
        self.results = np.array(self.results)
        for (xi,yi),value in np.ndenumerate(self.results):
            i = 0
            res = self.results[xi,yi]
            if res['w0'].value == 0:
                results[0:self.pre.shape[0],xi,yi] = self.pre[:,xi,yi]
            for param in res:
                results[i,xi,yi] = res[param].value
                i+=1
        
        self.results = results
                    
                    
    def fitPoint(self,spectra, params, fitfunc, constants):
        yd = spectra
        xr = np.linspace(0,yd.shape[0]-1,yd.shape[0])
#        plt.figure()
#        plt.plot(xr,yd)
        def residual(params, xr, yd, constants):
            
            model = fitfunc(params,xr,constants)
            return (yd-model)
        
#        print('input params:')
#        print([residual,params,xr,yd,constants])

        minimize(residual, params, args=(xr, yd,constants))
#        report_fit(params)
#        plt.plot(xr,fitfunc(params,xr,constants))
        

        return params        
        


                
        

            
            
class Shapes(Enum):
    rect = 1
    circle = 2
    ellipse = 3
    point = 4
    polygon = 5

class MouseInput(Enum):
    define_shape = 0
    view_spectra = 1

class Visual:
                
                        
        
        
    def __init__(self):
        self.v_unit = u.chan
        self.x_unit = u.arcsec
        self.y_unit = u.arcsec
        self.plot_ind = 0
        self.subplots = [0 for x in range(5)]
        self.limits_vel_en = -1
        self.points = ['red','green','blue','cyan','pink']
        self.plot_subplots  = [0]*5
        self.subplots_dots = []
        self.hist = []
        self.bgsimg = None
        self.draw_beam = True
        self.app = None
        self.mode = MouseInput.view_spectra
        self.patch = None
        self.d = Visual.Data(self)
        self.linewidths = 5
        self.contrs = []


    def scaleToJytB(self):
        xar = self.data.header['bmaj']*u.rad
        yar = self.data.header['bmin']*u.rad
        self.beam_area = ((2*np.pi/2.355/2.355)*(xar*yar))
        fr = self.units.freq*u.Hz
        self.d.xyfdata = (self.d.xyfdata*u.K).to(u.Jy, equivalencies=u.brightness_temperature(self.beam_area,fr)).value
        self.iunit = u.Jy/u.rad**2

    def scaleToJytPix(self):
        xar = self.data.getArrayItem('bmaj')*u.rad
        yar = self.data.getArrayItem('bmin')*u.rad
        scale = u.rad.to(u.arcsec)
        xdelta = abs(self.data.getArrayItem("cdelt1")[0])*scale
        ydelta = abs(self.data.getArrayItem("cdelt2")[0])*scale
        self.beam_area = ((2*np.pi/2.355/2.355)*(xar*yar)).to(u.arcsec**2).value
        fr = self.units.freq*u.Hz
        self.d.xyfdata = self.d.xyfdata/(self.beam_area/(xdelta*ydelta))
        self.iunit = u.Jy
        
    def calEnergy(self,points,data):
        
        self.calculetUpperEnergy(points.mean(1).mean(1).sum(),J,SnlMu2)        

    def calculetUpperEnergy(self,integral,J,SnlMu2):
        xar = ((self.data.getArrayItem('bmaj')*u.rad).to(u.arcsec)).value
        yar = ((self.data.getArrayItem('bmin')*u.rad).to(u.arcsec)).value
        nu = self.units.freq
        Ntog = ((integral*(c*100)**2)/(SnlMu2*xar*yar*nu**3)*(3/16/2.5/1.13)*10e27)
        print(Ntog,integral,SnlMu2)
        return Ntog

    def scaleToK(self):
        xar = self.data.getArrayItem('bmaj')*u.rad
        yar = self.data.getArrayItem('bmin')*u.rad
        self.beam_area = ((2*np.pi/2.355/2.355)*(xar*yar))
        self.fr = self.units.freq*u.Hz
        self.d.xyfdata = (self.d.xyfdata*u.Jy).to(u.K, equivalencies=u.brightness_temperature(self.beam_area,self.fr)).value
        self.d.xyfdata = self.d.xyfdata/(((self.f[-1]*u.chan).to(u.MHz)).value-((self.f[0]*u.chan).to(u.MHz)).value)
        self.iunit = u.K/u.MHz

    
    
    def setVelUnit(self, unit):
        self.v_unit = unit;
        
    def setFitsData(self,data):
        self.data = data
        self.d.xyfdata = data.data.copy()
        if self.d.xyfdata.shape[0] == 1: self.d.xyfdata = self.d.xyfdata[0]
        for (x,y,z), val in np.ndenumerate(self.d.xyfdata):
            if (val != val):
                self.d.xyfdata[x,y,z] = 0

#        self.d.xyfdata = np.swapaxes(self.d.xyfdata,2,0)
        
        self.f = np.arange(0,self.d.xyfdata.shape[0]);
        self.units = Unit(self.data)
#        self.units = Unit(self.data)self.subplots[self.plot_ind]
        
        
    def setData(self,data):
        self.data = data
        self.d.xyfdata = np.ones([data.axes[2],data.axes[1],data.axes[0]])
        for i in range(data.axes[2]) :
            self.d.xyfdata[i] = data.readPlane(axes=[i]).data
        self.f = np.arange(0,self.d.xyfdata.shape[0]);
        self.units = Unit(self.data)
        
#        self.units.register()
        
    def setBg(self, img, axes=None):
        self.bgimg = img
        if axes == None:
            axes = self.units.wcs

        self.sbglevs = np.linspace(img.min(),img.max(),5)
        self.bgimg_axes = axes

    def setSBg(self, img, axes=None):
        self.bgsimg = img
        if axes == None:
            axes = self.units.wcs

        self.bgsimg_axes = axes
        
    def on_move(self,event):
        if event.inaxes!=self.subp.axes: return
        xc = event.xdata
        yc = event.ydata
        xscale = 1
        yscale = 1
        if (self.mode == MouseInput.define_shape):
            if (self.d.shape == Shapes.point):
                pass
            elif (self.d.shape == Shapes.rect):
                if (len(self.d.point_of_interest) == 1):
                    self.patch._width = (-self.d.point_of_interest[0][0]+xc)*xscale
                    self.patch._height = (-self.d.point_of_interest[0][1]+yc)*yscale
            elif (self.d.shape == Shapes.circle):
                if (len(self.d.point_of_interest) == 1):
                    r = rangexy(self.d.point_of_interest[0],(xc,yc))*yscale
                    print ('CIRCLE Rad %f, x1,%f,y1%f, x2,%f,y2%f'%(r,xc,yc,self.d.point_of_interest[0][0],self.d.point_of_interest[0][1]))
                    self.patch.set_radius(r)
            elif (self.d.shape == Shapes.ellipse):
                if (len(self.d.point_of_interest) == 1):
                    r = rangexy(self.d.point_of_interest[0],(xc,yc))*yscale
                    print ('CIRCLE Rad %f, x1,%f,y1%f, x2,%f,y2%f'%(r,xc,yc,self.d.point_of_interest[0][0],self.d.point_of_interest[0][1]))
                    self.patch.set_radius(r)
                elif (len(self.d.point_of_interest) == 2):
                    xy1 = self.d.point_of_interest[0]
                    xy2 = self.d.point_of_interest[1]
                    xy3 = (xc,yc)
    
                    centr = ((xy1[0]*self.units.xpix).to(self.x_unit).value,(xy1[1]*self.units.ypix).to(self.y_unit).value)
                    long = ((xy2[0]*self.units.xpix).to(self.x_unit).value,(xy2[1]*self.units.ypix).to(self.y_unit).value)
                    short = ((xy3[0]*self.units.xpix).to(self.x_unit).value,(xy3[1]*self.units.ypix).to(self.y_unit).value)
                    A = centr
                    B = long
                    C = short
                    AB = (B[0]-A[0],B[1]-A[1])
                    AC = (C[0]-A[0],C[1]-A[1])
                    sinc = ((1-((AB[0]*AC[0]+AB[1]*AC[1])/rangexy(A,B)/rangexy(C,A))**2)**0.5)
                    print (sinc)
                    width = rangexy(C,A)*sinc
                    self.patch.height = width*2

            plt.draw()
                
        
        
    
    def onclick(self,event):
        self.click = event
            
        if event.inaxes!=self.subp.axes: return

        xc = int(event.ydata)
        yc = int(event.xdata)
        print ('button=%d, x=%i, y=%i, xdata=%f, ydata=%f'%(event.button, xc,yc, event.xdata, event.ydata))

        if (self.mode == MouseInput.define_shape):
            if (self.d.shape == Shapes.point):
                print('adding {0},{1} to shapes'.format(xc,yc))
                self.d.point_of_interest=((xc,yc),)
            elif (self.d.shape == Shapes.rect):
                if (len(self.d.point_of_interest) == 1):
                    self.d.point_of_interest+=((xc,yc),) 
                else :
                    self.d.point_of_interest=((xc,yc),) 
            elif (self.d.shape == Shapes.circle):
                if (len(self.d.point_of_interest) == 1):
                    self.d.point_of_interest+=((xc,yc),) 
                else :
                    self.d.point_of_interest=((xc,yc),) 
            elif (self.d.shape == Shapes.ellipse):
                if (len(self.d.point_of_interest) == 1 or len(self.d.point_of_interest) == 2):
                    self.d.point_of_interest+=((xc,yc),) 
                else :
                    self.d.point_of_interest=((xc,yc),) 
            self.plotAreaOfInterest()
            plt.draw()
        else :
            self.subplots[self.plot_ind] = plt.subplot2grid((5,3),(self.plot_ind,2)).axes
            
#            xax = self.bgimg_axes[0]
#            yax = self.bgimg_axes[1]
            x =  xc
            y =  yc
#            self.bgsimg = np.zeros(self.bgimg.shape)
#            self.bgsimg[x,y]  = 1
    #        plt.figure()
            self.results = np.zeros((self.bgimg.shape[0],self.bgimg.shape[1],6))
            self.subplot(x,y)
            self.hist.append((x,y))
    #        pl.title([pl[1]],['x = {1}, y = {2}'.format(round(event.xdata),round(event.ydata))])

     
            if len( self.subplots_dots) < 5 :
                self.subplots_dots.append(self.subp.plot(event.xdata,event.ydata,marker='o',color=self.points[self.plot_ind])[0])
                
            else : 
                self.subplots_dots[self.plot_ind].remove()
                self.subplots_dots[self.plot_ind] = (self.subp.plot(event.xdata,event.ydata,marker='o',color=self.points[self.plot_ind])[0])
    
        '''
        i = 0
        for (xi,yi),value in np.ndenumerate(self.bgsimg):
            for (xi,yi),value in np.ndenumerate(self.bgsimg):
                if (value == 0.5):
                    i+=1
                    self.bgsimg[xi,yi] = 1
                    self.subplot(xi,yi)
                    print (i)
                    if (i == 200):
                        return
'''

#        if event.inaxes!=self.subp.axes: return
                
        
    def subplot(self,x,y):
        if self.limits_vel_en == 0:
            xr = ((self.f*u.chan).to(self.v_unit)).value[self.limits_vel[0]:self.limits_vel[1]]
            yd = self.d.xyfdata[self.limits_vel[0]:self.limits_vel[1],x,y]
        else :
            xr = ((self.f*u.chan).to(self.v_unit)).value
            yd = self.d.xyfdata[:,x,y]
        mx = np.linspace(self.f.min(),self.f.max(),1000)
        if hasattr(self,'average'):
            yd = self.d.xyfdata[:,x-self.average:x+self.average,y-self.average:y+self.average].mean(1).mean(1)
            
            
        self.plot_ind=(self.plot_ind+1)%5
        ax1 = plt.gca()
        ax1.get_xaxis().get_major_formatter().set_useOffset(False)
        self.plot_subplots[self.plot_ind] = ax1
        ax1.step(xr,yd,color=self.points[self.plot_ind])
        ax1.set_xlabel(self.v_unit.to_string())
        ax1.set_ylabel(self.iunit.name)
        plt.text(1, 1,"({0:0.1f},{1:0.1f})".format(x,y),
                        horizontalalignment='right',
                        verticalalignment='top', size=10,
                        transform = self.subplots[self.plot_ind].transAxes)

#        ca = plt.gca().twiny()
#        ca.step(self.f,yd,color=self.points[self.plot_ind])
#        ca.set_xlabel("chan")
        xx = np.linspace(0,self.d.xyfdata.shape[0],self.d.xyfdata.shape[0])
        if (hasattr(self,'an')):
            aa = self.subfit(self.an.results,x,y,mx)        
            print(yd)
            kuku = self.subfit(
                            self.an.results,x,y,xx
                            ) 
            self.aa= yd-kuku
            ax1.plot(mx,aa,'-g')
            ax1.legend(['data','fit, std = {0}'.format(
                    np.std( self.aa[40:70])
                )])
#        ca.set_ylabel(self.iunit.name)
        plt.grid()
        

#        plt.title(self.data.getArrayItem('object').decode('utf-8')+" $(%0.2f,%0.2f),''$"%((x*self.units.xpix).to(self.units.rsecx).value,(y*self.units.ypix).to(self.units.rsecy).value))
        plt.plot()
#        self.optimize(x,y,ax1)
#        self.hist.append((x,y))
#        pl.title([pl[1]],['x = {1}, y = {2}'.format(round(event.xdata),round(event.ydata))])
        
 
 #       if len( self.subplots_dots) < 5 :
 #           self.subplots_dots.append(self.subp.plot(event.xdata,event.ydata,marker='o',color=self.points[self.plot_ind])[0])
            
 #       else : 
 #           self.subplots_dots[self.plot_ind].remove()
 #           self.subplots_dots[self.plot_ind] = (self.subp.plot(event.xdata,event.ydata,marker='o',color=self.points[self.plot_ind])[0])
            
    def thin2 (self, p, x, sp, ch):
        from numpy import size, exp
        sigma = p['sigma'].value
        w0 = p['w0'].value
        tamp = p['tamp'].value
        Ajk = p['Ajk'].value
        v1 = p['v1'].value
        v2 = p['v2'].value
        s=0
        
        s += Ajk/v1*(1-exp(-self.tau2(x,sp,ch,[tamp,0,0],w0,sigma)))
        s += Ajk*(1-exp(-self.tau2(x,sp,ch,[0,tamp,0],w0,sigma)))
        s += Ajk/v2*(1-exp(-self.tau2(x,sp,ch,[0,0,tamp],w0,sigma)))
#        print(s)
        return s
    
    def tau2(self,x,sp,ch,tamp,w0,sigma):
        from numpy import exp
        t = 0
        for i in range(np.size(sp)):
            ind = 0 if i<6 else 1 if i < 12  else 2
            am = tamp[ind]
            self.aa = np.array(-((x-w0-ch[i])**2)/2/sigma)
#            print(self.aa)
            t += (sp[i]*am*0.39894228/sigma)*exp(self.aa)
 #           print(t)
        return t
    
    
        
    def optimizeDoubleGaussImage(self):
        x = self.d.point_of_interest[0][0]
        y = self.d.point_of_interest[0][0]
        spectra = self.d.xyfdata[:,x,y]
        w0 = np.where(spectra == spectra.max())[0][0]
        params = Parameters()
        params.add('w0',w0, min = w0-20, max = w0+20)
        params.add('dw',w0, min = -20, max = +20)
        params.add('I',spectra.max(),min = 0, max = 1000)
        params.add('sigma', 1, min=0.01, max = 20)
        params.add('I2',spectra.max(),min = 0, max = 1000)
        params.add('sigma2', 1, min=0.01, max = 20)
        
        def fitfunc(p,x,constants):
            w0 = params['w0'].value
            I  = params['I'].value
            sigma = params['sigma'].value
            dw = params['dw'].value
            I2  = params['I2'].value
            sigma2 = params['sigma2'].value
            b = -0.5*(((x-w0)/sigma)**2)
            aa = np.exp(b)*I
            
            aa+= np.exp(-0.5*(((x-dw-w0)/sigma2)**2))*I2
            return aa
        
        if not hasattr(self,'an'):
            self.an = Analyse()
        self.an.fitImage(self.d.xyfdata,fitfunc,(x,y),params,np.array([0]))

    def optimizeGaussImage(self):
        x = self.d.point_of_interest[0][0]
        y = self.d.point_of_interest[0][1]
        spectra = self.d.xyfdata[:,x,y]
        w0 = np.where(spectra == spectra.max())[0][0]
        params = Parameters()
        params.add('w0',w0, min = w0-20, max = w0+20)
        params.add('I',spectra.max(),min = 0, max = 1000)
        params.add('sigma', 1, min=0.01, max = 20)
        
        def fitfunc(p,x,constants):
            w0 = params['w0'].value
            I  = params['I'].value
            sigma = params['sigma'].value
            b = -0.5*(((x-w0)/sigma)**2)
            aa = np.exp(b)*I
            return aa
            
        self.an = Analyse()
        self.an.fitImage(self.d.xyfdata,fitfunc,(x,y),params,np.array([0]))

            
    def optimize(self,cat):
        read = np.genfromtxt(cat)
        F = read[:,2]
        A = read[:,7]
        s = (2*F+1)*A
        w = read[:,6]*1e6
        ch = self.units._freq2chanel(w)
        constants = np.array([s,ch])
    #errfunc = lambda p, x, y: doublegauss(p, x) - y
        fitfunc = lambda params,xr,constants: self.thin2(params,xr,constants[0],constants[1])#+self.thin2(pd1,xr,s,ch)

        params = Parameters()
        params.add('sigma',   0.5,  min=0.01, max=2)
        params.add('w0', 14.2,  min=0, max=30)
        params.add('tamp', 51.44, min=0)
        params.add('Ajk', 18.2, min=0)
        params.add('var1', 1,  min=0.5, max=1.5)
        params.add('var2', 1,  min=0.5, max=1.5)
        params.add('const', value=0.1)
        an = Analyse()
        an.fitImage(self.d.xyfdata,fitfunc,(10,10),params,constants)
        
#        plt.plot((np.linspace(0,np.max(xr),5024)*u.chan).to(self.v_unit).value,fitfunc(p1,[1,1,0,0,1,1],np.linspace(0,np.max(xr),5024),s,ch)+params['const'].value,'-', color= [0,1,0,0.3])
#        plt.plot((np.linspace(0,np.max(xr),5024)*u.chan).to(self.v_unit).value,fitfunc([1,1,0,0,1,1],p2,np.linspace(0,np.max(xr),5024),s,ch)+params['const'].value,'-', color= [1,0,0,0.3] )
#        plt.plot((np.linspace(0,np.max(xr),5024)*u.chan).to(self.v_unit).value,fitfunc(p1,np.linspace(0,np.max(xr),5024),s,ch)+params['const'].value,'-k' )
#        plt.step(((xr+0.5)*u.chan).to(self.v_unit),yd,color=self.points[self.plot_ind])
#        self.ch = ch
#        self.w0 = params['w0'].value
#        w1 = params['w01'].value
#        vel = (((np.array(ch)+params['w0'].value)*u.chan).to(self.v_unit)).value
#        vel1 = (((np.array(ch)+w1)*u.chan).to(self.v_unit))
#        plt.bar(vel,s*3e2,0.01)
#        plt.bar(vel1,s*4e1,0.01)
#        plt.legend([r'$\tau={0}$'.format(self.tau2(xr,s,ch,[params['tamp'].value,params['tamp'].value,params['tamp'].value],params['w0'].value,params['sigma'].value)[ch[10]+params['w0'].value])])
#        print([r'$\tau={0}$'.format(self.tau2(xr,s,ch,[params['tamp'].value,params['tamp'].value,params['tamp'].value],params['w0'].value,params['sigma'].value)[ch[10]+params['w0'].value])])
#        print('finished!')
        
    def fitModel(self, model,y):
        ch = model.ch
        fitfunc = model.fitfunc
        resduial = lambda p,x,y,s,ch: y-fitfunc(p,x,s)
        p0 = model.startparams
        result = minimize(resduial, p0, args=(x, y,model.consts))
        return result
        
        
    def optimize2(self,cat):
        read = np.genfromtxt(cat)
        x = self.d.point_of_interest[0][0]
        y = self.d.point_of_interest[0][1]
#v.        ax1 = self.plot_subplots[self.plot_ind]
        
        
        print('optimizing %i,%i'%(x,y))
        
        xr = self.f
        yd = self.d.xyfdata[:,x,y]
        
        F = read[:,2]
        A = read[:,7]
        s = (2*F+1)*A
        w = read[:,6]*1e6
        ch = self.units._freq2chanel(w)
    #errfunc = lambda p, x, y: doublegauss(p, x) - y
        def fitfunc(params,xr,s,ch):
            pd0 = {}
            pd0['sigma'] = params['sigma']
            pd0['w0'] = params['w0']
            pd0['tamp'] = params['tamp']
            pd0['Ajk'] = params['Ajk']
            pd0['v1'] = params['var1']
            pd0['v2'] = params['var2']
            pd1 = {}
            pd1['sigma'] = params['sigma1']
            pd1['w0'] = params['w01']
            pd1['tamp'] = params['tamp1']
            pd1['Ajk'] = params['Ajk1']
            pd1['v1'] = params['var11']
            pd1['v2'] = params['var21']
            return self.thin2(pd0,xr,s,ch)+self.thin2(pd1,xr,s,ch)
        def residual(params, xr, yd, s, ch):
            
            model = fitfunc(params,xr,s,ch)+params['const'].value
#            model = fitfunc([p0,p1,p2,p3,p4,p5],xr,s,ch)+p6
        
            return (yd-model)
        #p0 = [  0.1,  0,  110000,  0.6,         1,  1]
        #p0 = [0.5,14, 1, 1, 1, 1]
        
        p0 = [1.4999999999999993, 14.226822571742712, 51.441131251068143, 18.20801956881299, 0.50003757072406307, 1.4949296569912287]

        if self.app != None:
            p0 = self.app
        if (self.results[x,y,0]!=0):
            p0 = self.results[x,y]
        print('input params:')
        print(p0)
        params = Parameters()
        params.add('sigma',   value= p0[0],  min=p0[0]-1, max=p0[0]+1)
        params.add('w0', value=p0[1],  min=p0[1]-10, max=p0[1]+5)
        params.add('tamp', value= p0[2], min=0)
        params.add('Ajk', value=p0[3], min=0)
        params.add('var1', value=p0[4],  min=0.5, max=1.5)
        params.add('var2', value=p0[5],  min=0.5, max=1.5)
        params.add('const', value=0.1)

        params.add('sigma1',   value= p0[0],  min=p0[0]-1, max=p0[0]+1)
        params.add('w01', value=p0[1]+4,  min=p0[1]-5, max=p0[1]+5)
        params.add('tamp1', value= p0[2], min=0)
        params.add('Ajk1', value=p0[3], min=0)
        params.add('var11', value=p0[4],  min=0.5, max=1.5)
        params.add('var21', value=p0[5],  min=0.5, max=1.5)
#        print(xr)
#        print(yd)
        print(s,ch)
        result = minimize(residual, params, args=(xr, yd,s,ch))
        
        # calculate final result
        final = yd + result.residual
        
        # write error report
 #       report_fit(params)
        
        
        p0 = params['sigma'].value
        p1 = params['w0'].value
        p2 = params['tamp'].value
        p3 = params['Ajk'].value
        p4 = params['var1'].value
        p5 = params['var2'].value
        p1 = [p0,p1,p2,p3,p4,p5]

        p01 = params['sigma1'].value
        p11 = params['w01'].value
        p21 = params['tamp1'].value
        p31 = params['Ajk1'].value
        p41 = params['var11'].value
        p51 = params['var21'].value
        p2 = [p01,p11,p21,p31,p41,p51]

        print('output params:')
        print(p1)
        self.app = [p1,p2]
        
        plt.plot((np.linspace(0,np.max(xr),5024)*u.chan).to(self.v_unit).value,fitfunc(p1,[1,1,0,0,1,1],np.linspace(0,np.max(xr),5024),s,ch)+params['const'].value,'-', color= [0,1,0,0.3])
        plt.plot((np.linspace(0,np.max(xr),5024)*u.chan).to(self.v_unit).value,fitfunc([1,1,0,0,1,1],p2,np.linspace(0,np.max(xr),5024),s,ch)+params['const'].value,'-', color= [1,0,0,0.3] )
        plt.plot((np.linspace(0,np.max(xr),5024)*u.chan).to(self.v_unit).value,fitfunc(p1,p2,np.linspace(0,np.max(xr),5024),s,ch)+params['const'].value,'-k' )
        plt.step(((xr+0.5)*u.chan).to(self.v_unit),yd,color=self.points[self.plot_ind])
        self.ch = ch
        self.w0 = params['w0'].value
#        w1 = params['w01'].value
        vel = (((np.array(ch)+params['w0'].value)*u.chan).to(self.v_unit)).value
#        vel1 = (((np.array(ch)+w1)*u.chan).to(self.v_unit))
        plt.bar(vel,s*3e2,0.01)
#        plt.bar(vel1,s*4e1,0.01)
        plt.legend([r'$\tau={0}$'.format(self.tau2(xr,s,ch,[params['tamp'].value,params['tamp'].value,params['tamp'].value],params['w0'].value,params['sigma'].value)[ch[10]+params['w0'].value])])
        print([r'$\tau={0}$'.format(self.tau2(xr,s,ch,[params['tamp'].value,params['tamp'].value,params['tamp'].value],params['w0'].value,params['sigma'].value)[ch[10]+params['w0'].value])])
        print('finished!')
        
        
        self.results[x,y] = p1
        
        xi = x
        yi = y
        return

        for (xi,yi),value in np.ndenumerate(self.bgsimg):
            if (value == 2):
                if (xi < self.bgimg.shape[0]-1 and self.bgsimg[xi+1,yi] == 0):
                    self.bgsimg[xi+1,yi] = 0.5
                if (xi > 0 and self.bgsimg[xi-1,yi] == 0):
                    self.bgsimg[xi-1,yi] = 0.5
                if (yi < self.bgimg.shape[0]-1 and self.bgsimg[xi,yi+1] == 0):
                    self.bgsimg[xi,yi+1] = 0.5
                if (xi > 0 and self.bgsimg[xi,yi-1] == 0):
                    self.bgsimg[xi,yi-1] = 0.5
            
       
        if (xi < self.bgimg.shape[0]-1 and self.bgsimg[xi+1,yi] != 2):
            if (self.results[xi+1,yi,0]!= 0 ):
                self.results[xi+1,yi] = (p1+self.results[xi+1,yi,0])/2
            else:
                self.results[xi+1,yi] = p1

        if (xi > 0 and self.bgsimg[xi-1,yi] != 2):
            if (self.results[xi-1,yi,0]!= 0 ):
                self.results[xi-1,yi] = (p1+self.results[xi-1,yi,0])/2
            else:
                self.results[xi-1,yi] = p1
        if (yi < self.bgimg.shape[0]-1 and self.bgsimg[xi,yi+1] != 2):
            if (self.results[xi,yi+1,0]!= 0 ):
                self.results[xi,yi+1] = (p1+self.results[xi,yi+1,0])/2
            else:
                self.results[xi,yi+1] = p1
        if (xi > 0 and self.bgsimg[xi,yi-1] != 2):
            if (self.results[xi,yi+1,0]!= 0 ):
                self.results[xi,yi+1] = (p1+self.results[xi-1,yi,0])/2
            else:
                self.results[xi,yi+1] = p1


    def appendContours(self,contours,xax,yax,cmap,clevs):
        '''
        arg contours
        arg xax
        arg yax
        arg cmap
        arg clevs
        '''
        self.contrs.append((contours,xax,yax,cmap,clevs))
    
    def plot(self):       
        self.fig = plt.figure()
        units = self.units
        self.subp = mysubplot2grid((5,3),(0,0),rowspan=5,colspan=2, projection=self.units.wcs,slices=self.units.slices)
        ax = self.subp
        lon = ax.coords[0]
        lat = ax.coords[1]
        
#        self.subp = plt.subplot2grid((6,3),(0,0),rowspan=6,colspan=3)
#        self.subp.axes.set_aspect('equal')
#        plt.grid()
#        plt.title(self.data.getArrayItem('object').decode('utf-8'))
        if (self.bgsimg != None):
            if self.bgsimg_axes != units.wcs:
                cp = self.subp.imshow(
                    self.bgsimg, cmap='gray', interpolation="nearest",  transform=ax.get_transform(self.bgsimg_axes),origin='lower')
            else: cp = self.subp.imshow(
                self.bgsimg, cmap='gray', interpolation="nearest",origin='lower')
#            cb = plt.colorbar() 
        for (cntr,wcs,cmap,clevs) in self.contrs:
            cp = self.subp.contour(cntr,clevs
                ,linewidths=self.linewidths)
            if self.clabel: plt.clabel(self.cp, inline=1, fontsize=10)
            self.subp.set_cmap(cmap)
            cb = plt.colorbar()
            
        if (self.bgimg != None):
            if self.bgimg_axes != units.wcs:
                self.cp = self.subp.contour(
                    self.bgimg,self.sbglevs, cmap='Blues'
                    ,linewidths=self.linewidths,
                    transform=ax.get_transform(self.bgimg_axes))
                
            else: self.cp = self.subp.contour(
                self.bgimg,self.sbglevs, cmap='Blues'
                ,linewidths=self.linewidths)
            if self.clabel: plt.clabel(self.cp, inline=1, fontsize=10)
            #self.subp.set_cmap('Blues')
#            cb = plt.colorbar()
        if self.draw_beam:
            from matplotlib.patches import Ellipse
            xd = self.data.getArrayItem('bmaj')[0]*u.rad.to(u.arcsec)
            yd = self.data.getArrayItem('bmaj')[0]*u.rad.to(u.arcsec)
            beam = Ellipse((xd*2,yd*2),xd,yd,0)
            self.subp.add_patch(beam)
            
#        if hasattr(self,'iunit'):
#            cb.set_label(str(self.iunit))
#        if hasattr(self,'d'):
#        if (not hasattr(self,'buttons')):
#            pl = plt.subplot2grid((6,3),(5,0))      
#            self.b1 = MyButton(pl,"Spectra")
#            self.b1.on_clicked(self.updateClickType)
#            p2 = plt.subplot2grid((6,3),(5,1))      
 #           self.b2 = MyButton(p2,"Freq limit")
        
#            self.b2.on_clicked(self.manageLimits)
            


        cid = self.fig.canvas.mpl_connect('button_press_event', self.onclick)
        crid = self.fig.canvas.mpl_connect('motion_notify_event', self.on_move)
        plt.show()
        return

    def manageLimits(self, event):
        self.limits_vel_en = 0
    
    def updateClickType(self, event):
        self.mode = MouseInput((self.mode.value+1)%2)
        self.b1.set_label(self.mode.name)
        
        
    def inEllipse(self,m):
        centr = self.d.point_of_interest[0]
        big = self.d.point_of_interest[1]
        small = self.d.point_of_interest[2]
        CB = (big[0]-centr[0],big[1]-centr[1])
        a = rangexy(centr,big)
        b = rangexy(centr,small)
        c = (a**2-b**2)**0.5
        r = rangexy(CB,(0,0))
        CB = (CB[0]/r*c,(CB[1]/r*c))
        F1 = (centr[0]+CB[0],centr[1]+CB[1])
        F2 = (centr[0]-CB[0],centr[1]-CB[1])
        
        return rangexy(F1,m)+rangexy(F2,m) <= 2*a

    def inCircle(self,m):
        centr = self.d.point_of_interest[0]
        r = self.d.point_of_interest[1]
        
        return rangexy(centr,m)<= rangexy(centr,r)

        
    def plotAreaOfInterest(self): 
        if (self.patch != None):
            self.patch.remove()
        from matplotlib.patches import Rectangle, Circle, Ellipse, Polygon
            
        currentAxis = self.subp


        if (self.d.shape == Shapes.point):
            xy1 = self.d.point_of_interest[0]
            print ("POINT: centr:{0}".format(xy1))
            self.patch = currentAxis.add_patch(Polygon(np.array(xy1), facecolor="grey"))
        elif (self.d.shape == Shapes.rect):
            xy1 = self.d.point_of_interest[0]
            if (len(self.d.point_of_interest) == 1):
                xy2 = xy1
            else :                
                xy2 = self.d.point_of_interest[1]
            xy1 = ((xy1[0]*self.units.xpix).to(self.x_unit).value,(xy1[1]*self.units.ypix).to(self.y_unit).value)
            xy2 = ((xy2[0]*self.units.xpix).to(self.x_unit).value,(xy2[1]*self.units.ypix).to(self.y_unit).value)
            print ("rect: 1:{0}, 2:{1}".format(xy1,xy2))
            self.patch = currentAxis.add_patch(Rectangle((xy1),xy2[0]-xy1[0],xy2[1]-xy1[1], facecolor=[0.5,0.5,0.5,0.5]))
        elif (self.d.shape == Shapes.circle):
            xy1 = self.d.point_of_interest[0]
            if (len(self.d.point_of_interest) == 1):
                xy2 = xy1
            else :                
                xy2 = self.d.point_of_interest[1]

            xy1 = ((xy1[0]*self.units.xpix).to(self.x_unit).value,(xy1[1]*self.units.ypix).to(self.y_unit).value)
            xy2 = ((xy2[0]*self.units.xpix).to(self.x_unit).value,(xy2[1]*self.units.ypix).to(self.y_unit).value)
            rad = rangexy(xy1,xy2)
            print ("CIRCLE: centr:{0}, rad:{1}".format(xy1,rad))
            self.patch = currentAxis.add_patch(Circle(xy1, rad, facecolor=[0.5,0.5,0.5,0.5]))
        elif (self.d.shape == Shapes.ellipse):
            xy1 = self.d.point_of_interest[0]
            if (len (self.d.point_of_interest) == 1):
                xy2 = xy1
                rad = 2
                xy1 = ((xy1[0]*self.units.xpix).to(self.x_unit).value,(xy1[1]*self.units.ypix).to(self.y_unit).value)
                xy2 = ((xy2[0]*self.units.xpix).to(self.x_unit).value,(xy2[1]*self.units.ypix).to(self.y_unit).value)
                print ("CIRCLE: centr:{0}, rad:{1}".format(xy1,rad))
                self.patch = currentAxis.add_patch(Circle(xy1, rad, facecolor=[0.5,0.5,0.5,0.5]))
            else :
                if (len (self.d.point_of_interest) == 2):
                    xy2 = self.d.point_of_interest[1]
                    xy3 = self.d.point_of_interest[1]
                    centr = ((xy1[0]*self.units.xpix).to(self.x_unit).value,(xy1[1]*self.units.ypix).to(self.y_unit).value)
                    long = ((xy2[0]*self.units.xpix).to(self.x_unit).value,(xy2[1]*self.units.ypix).to(self.y_unit).value)
                    short = ((xy3[0]*self.units.xpix).to(self.x_unit).value,(xy3[1]*self.units.ypix).to(self.y_unit).value)
                    A = centr
                    B = long
                    C = short
                    AB = (B[0]-A[0],B[1]-A[1])
                    AC = (C[0]-A[0],C[1]-A[1])
                    width = 2
                else :
                    xy2 = self.d.point_of_interest[1]
                    xy3 = self.d.point_of_interest[2]
    
                    centr = ((xy1[0]*self.units.xpix).to(self.x_unit).value,(xy1[1]*self.units.ypix).to(self.y_unit).value)
                    long = ((xy2[0]*self.units.xpix).to(self.x_unit).value,(xy2[1]*self.units.ypix).to(self.y_unit).value)
                    short = ((xy3[0]*self.units.xpix).to(self.x_unit).value,(xy3[1]*self.units.ypix).to(self.y_unit).value)
                    A = centr
                    B = long
                    C = short
                    AB = (B[0]-A[0],B[1]-A[1])
                    AC = (C[0]-A[0],C[1]-A[1])
                    sinc = ((1-((AB[0]*AC[0]+AB[1]*AC[1])/rangexy(A,B)/rangexy(C,A))**2)**0.5)
                    print (sinc)
                    width = rangexy(C,A)*sinc

                print ("ELLIPSE: centr:{0}, width:{1}, height:{2}, angle:{3}".format(centr,rangexy(centr,short)*2,width, atan((long[0]-centr[0])/(long[1]-centr[1]))))
                self.patch = currentAxis.add_patch(Ellipse(centr,rangexy(centr,long)*2,width*2, atan((long[1]-centr[1])/(long[0]-centr[0]))*u.rad.to(u.deg) , facecolor=[0.5,0.5,0.5,0.5]))
        elif (self.d.shape == Shapes.polygon):
            self.patch = currentAxis.add_patch(Polygon(np.array(self.d.point_of_interest), facecolor=[0.5,0.5,0.5,0.5]))
        plt.draw()
            
            
    def axes(self, index=0):
        if index == 0:
            return ( 
                np.arange(0,self.data.axes[0]),
                np.arange(0,self.data.axes[1]),
                np.arange(0,self.data.axes[2])
            )
        else:
            return ( 
                np.arange(0,self.axes[i]),
            )
        
    def removeBaseline(self):
        for (x,y),val in np.ndenumerate(self.d.xyfdata[0]):
            if self.d.xyfdata[:,x,y].mean() != 0: 
                spectra = self.d.xyfdata[:,x,y].copy()
                tmp = np.zeros(spectra.shape)
                for f,i in enumerate(spectra[:-20]):
                    tmp[f] = np.std(spectra[f:f+20])
                    
                hold = 0.05
                print('1:({0},{1}): {2}!'.format(x,y,np.count_nonzero(tmp > hold)/len(tmp)))
                while  np.count_nonzero(tmp > hold)/len(tmp) > 0.1 :
                    hold=2*hold
                    print('2:({0},{1}): {3}|{2}!'.format(x,y,np.count_nonzero(tmp > hold)/len(tmp),hold))
    
                for f,i in enumerate(spectra[0:-20]):
                    if (tmp[f]>hold): spectra[f+10]=spectra[f+9]
                model = baseline(spectra,fitfunc,500)
                self.d.xyfdata[:,x,y]-=model
            else:
                print('({0},{1}): skip!'.format(x,y))

