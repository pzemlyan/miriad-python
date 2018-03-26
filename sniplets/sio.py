p1 = Parameters()
p1.add('sigma1',400,min=2e2,  max=1e4)
p1.add('dv1',9000,   min=3376,max=13252,vary=True)
p1.add('I1',1,   min=0,  max=10)

g1 = lambda p,x : p['I1'].value*exp(-(x-p['dv1'].value)**2/2/(p['sigma1'].value)**2)

def res1(p,x,y,c):
    return (y - g1(p,x))**2

resc = zeros(cube2.shape)
ix = arange(cube.shape[1])
iy = arange(cube.shape[2])
parnh = [['' for z in range( 0,cube.shape[2])] for y in range(0,cube.shape[1])]
figure(1)
import time
for xi in ix:
    for yi in iy:
        spec = cube2[:,xi,yi]
        if (spec.max()>spec[80:].std()*factor):# and parnh[xi][yi] == '':
            p = p1
            mini = lmfit.Minimizer(res1,p,fcn_args=(velx,spec,[]))
            par2hned = mini.minimize(method='nelder')

            parnh[xi][yi] = minimize(res1,par2hned.params,args=(velx,spec,[]),)
            resc[:,xi,yi] = parnh[xi][yi].residual




import pywt
from statsmodels.robust import stand_mad
wings = zeros(cube.shape)
for (x,y),val in ndenumerate(cube[0]):
    spectra2 = cube[:,x,y]
    noisy_coefs = pywt.wavedec(spectra2, 'db4', level=11, mode='per')
    sigma = stand_mad(noisy_coefs[-1])
    uthresh = sigma*np.sqrt(2*np.log(len(spectra2)))
    denoised = noisy_coefs[:]
    
    denoised[:] = (pywt.thresholding.soft(i, value=uthresh) for i in denoised[:])
    
    signal = pywt.waverec(denoised, 'db4', mode='per')
    wings[:,x,y] = signal[:-1]
    
    
def getImages(par):
    dv1 =ma.zeros((nh3bg.shape[0],nh3bg.shape[1]))
    sigma1 =ma.zeros((nh3bg.shape[0],nh3bg.shape[1]))
    I1 =ma.zeros((nh3bg.shape[0],nh3bg.shape[1]))
    
    dv1.mask=False
    sigma1.mask=False
    I1.mask=False
    
    xi,yi=0,0
    for xi in ix:
        for yi in iy:
            p = par[xi][yi]
            if (p == '') or (p == 'none'):
                sigma1.mask[xi,yi]|= True
                dv1.mask[xi,yi]|= True
                I1.mask[xi,yi]|= True
    
    
            else :
#                tau[xi,yi,0]=g1.tau(x,p.params).max()
                sigma1[xi,yi] = p.params['sigma1']
                dv1[xi,yi] = p.params['dv1']
                I1[xi,yi] = p.params['I1']
                
    return (sigma1,dv1,I1)
