
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 17:15:38 2016

@author: pete
"""

#a
def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]

def line(p1,p2):
#    k = (p2[1]-p1[1])/(p2[0]-p1[0])
#    b = p1[1] - k*p1[0]
    
    b = (p2[0]-p1[0])/(p2[1]*p1[0]-p2[0]/p1[1])
    a = -(b*p1[0]+1)/p1[0]
    return (a,b)
    
def pv(k, b, dv, flux, sigma, rmax, vrange, points):
    y = lambda x: k*x+b
    r = lambda x,y: (k*x-y+b)/(k**2+1)
    rax = linspace(-rmax,rmax,points)
    vax = linspace(vrange[0],vrange[1],points)
    pv = ma.zeros((points,points))
    p2v = [[[] for j in range(points)] for i in range(points)]
    for (ix,iy),val in ndenumerate(dv):
        vv = dv[ix,iy]
        if dv.mask[ix,iy] or flux.mask[ix,iy]: continue
        rv = r(iy,ix)
        ri = np.abs(rax-rv).argmin()
        vi = np.abs(vax-vv).argmin()
        pv[ri,vi] += flux[ix][iy]
        p2v[ri][vi].append(sigma[ix][iy])
    pv.mask |= pv ==0
    return (pv,p2v,rax,vax)

def maskLine(k,b,delta,dv):
    r = lambda x,y: (k*x-y+b)/(k**2+1)
    for (ix,iy),val in ndenumerate(dv):
        if dv.mask[ix,iy]: continue
        rv = r(iy+0.5,ix+0.5) #center of a pixel
        if abs(rv) > delta: dv.mask[ix,iy] = True
                
                
#def line(p1,p2):
#    k = (p2[0]-p1[0])/(p2[1]-p1[1])
#    b = p1[0] - k*p1[1]
#s    return (k,b)
                
def line(p1,p2):
    k = (p2[1]-p1[1])/(p2[0]-p1[0])
    b = p1[1] - k*p1[0]
    return (k,b)             
#for (ix,iy),val in ndenumerate(pvd[0]):
#    for err in pvd[1][ix][iy]:
#        plt.errorbar(pvd[3][iy], pvd[2][ix], xerr=err/2,color='red')


class Positionerovelociter() :
    '''
    Need to be self.armed = True
    '''
    def __init__(self,rec,cube):
        self.rec= rec
        self.cube = cube
        self.armed = False
        self.points = []
        
        
    def loadImg(self,img):
        self.img = img
        
    def onclick(self,event):
        self.levent = event
        if self.armed:
            ion()
            xc = event.ydata
            yc = event.xdata
            self.points.append((yc,xc))
            self.cpoint = self.ax.plot(yc,xc,'xr')
            draw()

    def onhover(self,event):
        self.levent = event
        if self.armed:
            xc = event.ydata
            yc = event.xdata
            if (len(self.points)%2==1):
                if (hasattr(self,'cline')) and self.cline in self.ax.lines:
                    self.cline.remove()
                xt = 2*self.points[-1][1]-xc
                yt = 2*self.points[-1][0]-yc
                self.cline = self.ax.plot([yc,yt],[xc,xt],'-k')[0]
                draw()

            

    def pvd(self,width,step,number,chans_range=(0,0)):
        from numpy import ma, array, vdot
        from numpy.linalg import norm
        from astropy import units as u
        delt = abs(self.rec[0].header['cdelt1']*u.deg.to(u.arcsec))
        width /= delt
        step /= delt
        number = int(number/delt/step)
        p1 = self.points[-2]
        p2 = self.points[-1]
        n1 = array((p2[0]-p1[0],p2[1]-p1[1]))
        n1 = n1/norm(n1)
        n2 = array((n1[1],-n1[0]))
        point = array(p1)
        if chans_range[1] == 0:
            chans_range = (0,self.cube.shape[0])
        ispec = len(self.cube[chans_range[0]:chans_range[1],0,0])

        pvd = zeros((number*2,ispec))
        for i in range(-number,number):
            co = n1*i*step+n1*step/2
            spec = [zeros(ispec)]
            for (xi,yi),val in ndenumerate(self.cube[0]):
                vec = array((yi,xi))-point-co
                if abs(vdot(vec,n1)) < step and abs(vdot(vec,n2)) < width :
                    spec.append(self.cube[chans_range[0]:chans_range[1],xi,yi])
            spec = array(spec).mean(0)
            pvd[i+number]=spec
        
        self.cr = array(chans_range)
        if self.cr[-1] == -1:
            self.cr[-1]=ispec-1
        self.number = number
        self.step   = step
        self.width  = width
        return pvd
            
            

        
    def target(self,fig=None):
        if fig is None: 
            fig,ax = subplots()
        else :
            ax = fig.gca()
            ax.imshow(self.img)
        self.fig = fig
        self.ax = ax
        self.armed=True
        self.cid = fig.canvas.mpl_connect('button_press_event', self.onclick)
        self.cid = fig.canvas.mpl_connect('motion_notify_event', self.onhover)
        
    def extent(self, rec):
        from numpy.linalg import norm
        from numpy import arange, zeros, array
        from astropy import wcs
        
        vel = []
        ra = []
        dec = []    
        ind = ['1','2','3']
        rec = rec.header
        for c in ind:
            if (rec['CTYPE'+c].startswith('VELO')):
                ind.remove(c)
                vel = [rec['crval'+c],rec['cdelt'+c],rec['crpix'+c]]
                break
        
        c = ind[0]
        ra = [rec['crval'+c],rec['cdelt'+c],rec['crpix'+c]]
        c = ind[1]
        dec = [rec['crval'+c],rec['cdelt'+c],rec['crpix'+c]]

        w = wcs.WCS(rec)

        chan = arange(self.cube.shape[0])
        fr = zeros((chan.shape[0],4))
        fr[:,2] = chan
        fr = w.all_pix2world(fr,0)
        v = fr[:,2]#(x*u.Hz).to(u.km/u.s, equ        p1 = self.points[-2]
        
        p2 = self.points[-1]
        p1 = self.points[-2]
        n1 = array((p2[0]-p1[0],p2[1]-p1[1]))
        n1 /= norm(n1)
        n1*=self.number*self.step
        coord = (n1[0]*ra[1]*u.deg.to(u.arcsec),n1[1]*dec[1]*u.deg.to(u.arcsec))
        leng = lambda c: (c[0]**2+c[1]**2)**0.5
        extent = [v[self.cr[0]],v[self.cr[1]-1],-leng(coord),leng(coord)]
        return extent
    
    def plot(self):
        from numpy import arange, zeros, array
        import matplotlib.patches as patches

        p1 = self.points[-2]
        p2 = self.points[-1]
        n1 = array((p2[0]-p1[0],p2[1]-p1[1]))
        n1 = n1/norm(n1)
        n2 = array((n1[1],-n1[0]))

        tr = p1+n1*self.number-n2*self.width        
        tl = p1+n1*self.number+n2*self.width        
        lr = p1-n1*self.number-n2*self.width        
        ll = p1-n1*self.number+n2*self.width        
        
        poly = array([tr,tl,ll,lr])
        
        self.ax.plot(self.points[-2][0],self.points[-2][1],'xr',markersize=12)
        self.ax.plot(self.points[-1][0],self.points[-1][1],'xr',markersize=12)
        self.ax.add_patch(patches.Polygon(poly,closed=True,alpha=0.5))
        
        
class PolyPositionerovelociter() :
    '''
    Need to be self.armed = True
    '''
    def __init__(self,cube):
        self.cube = cube
        self.armed = False
        self.points = []
        self.number = 0
        
    def loadImg(self,img):
        self.img = img
        
    def onclick(self,event):
        self.levent = event
        if self.armed:
            ion()
            xc = event.ydata
            yc = event.xdata
            if (len(self.points) > 0):
                xt = self.points[-1][1]
                yt = self.points[-1][0]
                self.cline = self.ax.plot([yc,yt],[xc,xt],'-k')[0]
            self.points.append((yc,xc))
            self.cpoint = self.ax.plot(yc,xc,'xr')
            draw()

    def onhover(self,event):
        self.levent = event
        if self.armed:
            xc = event.ydata
            yc = event.xdata
            if (hasattr(self,'cline')) and self.cline in self.ax.lines:
                self.cline.remove()
            xt = self.points[-1][1]
            yt = self.points[-1][0]
            self.cline = self.ax.plot([yc,yt],[xc,xt],'-k')[0]
            draw()

    def pvd(self,width,step,chans_range=(0,0)):
        from numpy import ma, array, vdot, concatenate
        from numpy.linalg import norm
        if chans_range[1] == 0:
            chans_range = (0,self.cube.shape[0])
        ispec = len(self.cube[chans_range[0]:chans_range[1],0,0])
        self.len=0

        apvd = zeros((0,ispec))
        
        for i in range(len(self.points)-1):
            p1 = self.points[i]
            p2 = self.points[i+1]
            
            n1 = array((p2[0]-p1[0],p2[1]-p1[1]))
            number = int(norm(n1)/step)
            n1 = n1/norm(n1)
            n2 = array((n1[1],-n1[0]))
            point = array(p1)
            
            pvd = zeros((number,ispec))
            for i in range(0,number):
                co = n1*i*step+n1*step/2
                spec = [zeros(ispec)]
                for (xi,yi),val in ndenumerate(self.cube[0]):
                    vec = array((yi,xi))-point-co
                    if abs(vdot(vec,n1)) < step and abs(vdot(vec,n2)) < width :
                        spec.append(self.cube[chans_range[0]:chans_range[1],xi,yi])
                spec = array(spec).mean(0)
                pvd[i]=spec
            self.len += norm(array((p2[0]-p1[0],p2[1]-p1[1])))
            apvd = concatenate((apvd,pvd),axis=0)
            
        self.cr = array(chans_range)
        if self.cr[-1] == -1:
            self.cr[-1]=ispec-1
        self.step   = step
        self.width  = width
        return apvd
                
                
    
            
    def target(self,fig=None):
        if fig is None: 
            fig,ax = subplots()
        else :
            ax = fig.gca()
            ax.imshow(self.img)
        self.fig = fig
        self.ax = ax
        self.armed=True
        self.cid = fig.canvas.mpl_connect('button_press_event', self.onclick)
        self.cid = fig.canvas.mpl_connect('motion_notify_event', self.onhover)
        
    def extent(self, rec):
        from numpy.linalg import norm
        from numpy import arange, zeros, array
        from astropy import wcs
        
        vel = []
        ra = []
        dec = []    
        ind = ['1','2','3']
        rec = rec.header
        for c in ind:
            if (rec['CTYPE'+c].startswith('VELO')):
                ind.remove(c)
                vel = [rec['crval'+c],rec['cdelt'+c],rec['crpix'+c]]
                break
        
        c = ind[0]
        ra = [rec['crval'+c],rec['cdelt'+c],rec['crpix'+c]]
        c = ind[1]
        dec = [rec['crval'+c],rec['cdelt'+c],rec['crpix'+c]]

        w = wcs.WCS(rec)

        chan = arange(self.cube.shape[0])
        fr = zeros((chan.shape[0],4))
        fr[:,2] = chan
        fr = w.all_pix2world(fr,0)
        v = fr[:,2]#(x*u.Hz).to(u.km/u.s, equ        p1 = self.points[-2]
        
        coord = self.len*abs(ra[1]*u.deg.to(u.arcsec))
        leng = lambda c: (c[0]**2+c[1]**2)**0.5
        extent = [v[self.cr[0]],v[self.cr[1]-1],0,coord]
        return extent
    
    def plot(self,rec,step,color='r'):
        from numpy import arange, zeros, array
        import matplotlib.patches as patches
        pts = array(self.points)
        self.ax.plot(pts[:,0],pts[:,1],'-',color=color)
        ind = ['1','2','3']
        rec = rec.header
        for c in ind:
            if (rec['CTYPE'+c].startswith('VELO')):
                ind.remove(c)
                break
        
        c = ind[0]
        ra = [rec['crval'+c],rec['cdelt'+c],rec['crpix'+c]]
        p2as = abs(ra[1]*u.deg.to(u.arcsec))
        stepp = step/p2as
        
        leng=0
        last = self.points[0]
        prestep = 1
        plot(last[0],last[1],'o',color=color)
        for i in range(len(self.points)-1):
            p1 = self.points[i]
            p2 = self.points[i+1]
            n1 = array((p2[0]-p1[0],p2[1]-p1[1]))
            number = int(norm(n1)/step)
            n1 = n1/norm(n1)
            n2 = array((n1[1],-n1[0]))

            leng += norm(array((p2[0]-p1[0],p2[1]-p1[1])))
            if leng > stepp:
                while leng > stepp:
                    last = last+n1*stepp*prestep
                    self.ax.plot((last[0]-n2[0],last[0]+n2[0]),(last[1]-n2[1],last[1]+n2[1]),'-',color=color)
                    prestep = 1
                    leng-=stepp
#                leng-=stepp
            prestep=(stepp-leng)/stepp

            last = p2

    

    
    
    

