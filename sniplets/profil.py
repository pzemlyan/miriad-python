def profil(hco, dvr,vlsr,deltas,dist):
    from scipy.constants import pi
    from numpy import zeros
    from math import sqrt,exp
    
    
    
    c=2.997925e5
    idlev = 512
    v0 = zeros((idlev,))
    t0 = zeros((idlev,))
    iline = 1
    dvi = dvr/idlev
    vi = vlsr-dvr/2
    r0 = hco.scp.r0
    for iv in range(idlev):
        v0[iv]=vi
        vi+=dvi
    
    rim = deltas/3600/180*pi*dist*3.08e18
    isp = hco.scp.isp
    rl = hco.scp.rl
    im = 0
    for k in range(1,isp):
        if (rim > rl[k-1] and rim < rl[k]):
            im=k
            
    isp2 = (isp+1-im)*2-1
    rk = sqrt(rl[isp-1]**2-rim**2)
    
    ab = zeros((isp,))
    sab = 0.
    
    for k in range(im,isp):
        ab[k]=abs(sqrt(rl[k]**2-rim**2))-sab
        sab+=ab[k]
    
    hck = hco.const.hck
    tbg = hco.const.tbg
    tir = hco.ir.tir
    fir = hco.ir.fir
    freq = hco.lin.freq
    sc = hco.scp.sc
    tc = hco.scp.tc
    drl =hco.scp.dr


    vsys0 = hco.vsys.vsys0
    rvsys = hco.vsys.rvsys
    iout = hco.vsys.iout
    isys = hco.vsys.isys
    istep = hco.vsys.istep
    
    for iv in range(idlev):
        vi=v0[iv]
        ti = hck*freq[iline]/c*10.*(1-vi/c)
        ubg=1./(exp(ti/tbg)-1.)
        uir=1./(exp(ti/tir)-1.)
        u0=fir*uir+(1-fir)*ubg
        
        t0r=0
        tau=0
        sab=0
        
        for il in range(isp2):
            i = abs(isp-im-il)+im-1 ###+1?
            sl = sc[iline,i]
            abi = ab[i]
            iistep = int(istep)
            if i == im:
                abi*=2
                iistep*=2
            
            dr = abi/iistep
            rs=sab+dr
            for ist in range(iistep):
                drc=rk-rs
                rc=sqrt(drc**2+rim**2)
                cont=-dr/rc
                vs=vsys0*(rc/r0)**rvsys*cont
                if i == 1: vs=0
                if iout!=0 and i==isp and isys == 0: vs=0
                
                tl=tc[iline,i]*hco.prof1(vi-vs-vlsr,i)*dr/drl[i]
                if tl>80: dext=0.
                if abs(tl)<80  and abs(tl) > 1e-4: dext=exp(-tl)
                if abs(tl)<=1e-4: dext = 1.-tl
                if tl <= -80: raise ValueError('level temperature is too low!')
                t0r=t0r+ti*sl*(1.-dext)*exp(-tau)
                tau+=tl
                rs+=dr
                
            sab+=abi
            print('tl:{0}'.format(tl))
        t0[iv]=t0r-ti*u0*(1-exp(-tau))
#        taut(iv)=tau
    return v0,t0

                