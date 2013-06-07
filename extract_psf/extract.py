import numpy as np

import fortran.extrpsf as extrpsf

# turn output into astropy table instead?


class PSF:
    """Point-spread function for extraction"""

    def __init__(self, form=0, guesses=[0., 4., (2.5,True)]):
        self.form = form
        self.setguess(guesses)
        self.par = []
        self.err = []
        self.chi2 = []

    def setguess(self, guesses):
        assert len(guesses) == (3 if self.form == 0 or self.form == -1
                                else 2+(self.form-1)*3)
        self.npoldisp = np.zeros(len(guesses),dtype=np.int)
        self.guess = np.array([], dtype=[('par',np.float), ('fix',np.bool)])
        for i,g in enumerate(guesses):
            try:  # Element(s) given in full?
                npg = np.array(g, dtype=self.guess.dtype)
                # Catch something like ((0.,(1.,True)),...)
                # where the tuple would be taken to be a Boolean
                # (strictly, this should be given as [[0.,(1.,True)],...])
                if isinstance(g[1], tuple):
                    raise TypeError
            except (TypeError, ValueError):
                try:  # Perhaps a number with default choice of not fixed?
                    npg = np.array((g,False), dtype=self.guess.dtype)
                except (TypeError, ValueError):  # Now it must be a list
                    npg = np.array([], dtype=self.guess.dtype)
                    for gg in g:
                        try:  # Element is given in full?
                            npgg = np.array(gg, dtype=self.guess.dtype)
                        except (TypeError, ValueError):  # Just a number?
                            npgg = np.array((gg,False), dtype=self.guess.dtype)
                        npg = np.hstack([npg,npgg])
            self.npoldisp[i] = npg.size-1
            self.guess = np.hstack([self.guess,npg])

    def evalpar(self, y=0., useguess=False):
        if useguess:
            pars = self.guess['par']
        else:
            pars = self.par
        evalpar = np.zeros(len(self.npoldisp))
        k = 0
        for i in range(len(self.npoldisp)):
            for j in range(self.npoldisp[i]+1):
                evalpar[i] += pars[k]*pow(y, j)
                k += 1
        return evalpar

    def eval(self, x, y=0., useguess=False, normalise=True):
        """Evaluate PSF using fitted parameters"""
        # ignore derivatives from extrpsf.psf (all zero if Normalize=True)
        val,_ = extrpsf.psfmany(ipsf=self.form, x=x,
                                par=self.evalpar(y, useguess),
                                normalise=normalise)
        return val

    def offsets(self):
        if len(self.par) == len(self.guess):
            return self.par[0]
        else:
            return np.vstack([self.par[0],
                              self.par[0]+self.par[len(self.guess):]])


def ccd_noise(d, ron, fferr):
    return np.sqrt(np.abs(d) + ron*ron + (fferr*d)**2)


def extract(d, tracepos, psf,
            skypol=2, fitrange=None, clip=[5.,2.],
            e=None, ron=4., fferr=0.03,
            spatial='y', itesttype=103, ibadlimit=3, squeeze_dims=True):
    """Extract spectrum by fitting sky and PSF to stellar traces"""
    assert d.ndim <= 3  # should be at most stacks of images
    # estimate of uncertainties (if not passed on)
    if e is None:
        e = ccd_noise(d, ron, fferr)
    # number of images in stack (treated as orders for echelle extraction)
    norder = 1 if d.ndim < 3 else d.shape[0]
    if d.ndim == 1:
        ndisp, nspat = 1, d.shape[0]
        istep_d = np.array([1,1])
    else:
        if spatial == 'y':
            nspat, ndisp = d.shape[-2:]
            istep_d = np.array([ndisp,1])
        else:
            ndisp, nspat = d.shape[-2:]
            istep_d = np.array([1,nspat])

    # Define the "trace" of the star, offset by a chip size for each "order"
    # and defined for every pixel in the dispersion direction
    # (used to trace curved orders in echelle spectra)
    xord = np.arange(norder, dtype=float)*nspat
    xord = np.repeat(xord[np.newaxis,:], ndisp, axis=0)
    # Slit can in principle be tilted (but not used here)
    tilt1 = np.zeros_like(xord)
    if tracepos is None:
        # No stars to extract, define trace as middle
        nstar = 0
        tracepos = np.array([nspat/2])
    else:
        nstar = len(tracepos)
        # Get guesses for offset, FWHM, etc. from psf
        # as well as whether or not they should be fit for
        psfguess = psf.guess['par']
        ifit = np.array(np.logical_not(psf.guess['fix']), dtype=int)
        npoldisp = psf.npoldisp
        if nstar > 1:
            psfguess = np.hstack([psfguess, tracepos[1:]-tracepos[0]])
            ifit = np.hstack([ifit, np.ones(nstar-1,dtype=int)])
            npoldisp = np.hstack([npoldisp, np.zeros(nstar-1,dtype=int)])
    xord += tracepos[0] + 1  # fortran is 1-based
    # Define physical size of the chip
    # For multiple, use whole stack with different orders as a hack.
    kbeg = 1
    kend = nspat*norder
    lbeg = np.ones(norder, dtype=int)
    lend = lbeg-1+ndisp
    if fitrange is None:
        ipbeg = 1-int(np.round(xord[0,0]))
        ipend = ipbeg-1+nspat
    else:
        ipbeg = int(np.round(fitrange[0]+np.min(tracepos-tracepos[0])))
        ipend = int(np.round(fitrange[1]+np.max(tracepos-tracepos[0])))

    # Number of output pixels
    # and offsets for 1 dispersion pixel, 1 star, 1 order
    nout = ndisp*norder*max(nstar,1)
    istep_out = np.array([1, ndisp, ndisp*norder])

    # Fortran routine wants arrays 1-dimensional
    d_orig_shape = d.shape
    d.shape = e.shape = d.size,

    if nstar > 0:
        mindof = 1
        nitermax = 20
        (psf.par, psf.err, psf.chi2,
         out, eout, back, chi2, test,
         ntdiscard, ntbadl, ntbadh, nproblems) = extrpsf.extract_psf(
             d, istep_in=istep_d, ein=e, istart_ein=0, istep_ein=istep_d,
             ipbeg=ipbeg, ipend=ipend,
             kbeg=kbeg, kend=kend, lbeg=lbeg, lend=lend,
             xord=xord, tilt1=tilt1, tilt2=tilt1, tilt_flag=False, npoltilt=0,
             mindof=mindof, nitermax=nitermax, clip1=clip[0], clip2=clip[1],
             ipsf=psf.form, psfguess=psfguess, ifix=ifit,
             npoldisp=npoldisp, nshape=len(psf.npoldisp),
             nstar=nstar, npolsky=skypol,
             nout=nout, nback=nout, nchi2=nout, istep_out=istep_out,
             ntest=d.size, istart_test=0, istep_test=istep_d,
             o_flag=True, e_flag=True, b_flag=True, c_flag=True,
             p_flag=False, s_flag=True,
             itesttype=itesttype, rnull=0., ibadlimit=ibadlimit)
        # Restore proper shapes of input arrays
        d.shape = e.shape = test.shape = d_orig_shape
        # Set output array shape, removing unnecessary dimensions if wanted
        out.shape = eout.shape = back.shape = chi2.shape = [
            n for n in (nstar, norder, ndisp) if n > 1 or not squeeze_dims]

        return (out, eout, back, chi2, test,
                ntdiscard, ntbadl, ntbadh, nproblems)

    else:
        (_,_,_,_,_, back, chi2, test,
         _, ntbadl, ntbadh, _) = extrpsf.extract_psf(
             d, istep_in=istep_d, ein=e, istart_ein=0, istep_ein=istep_d,
             ipbeg=ipbeg, ipend=ipend,
             kbeg=kbeg, kend=kend, lbeg=lbeg, lend=lend,
             xord=xord, tilt1=tilt1, tilt2=tilt1, tilt_flag=False, npoltilt=0,
             mindof=0, nitermax=0, clip1=clip, clip2=0.,
             ipsf=0, psfguess=np.array([0.]), ifix=np.array([0]),
             npoldisp=np.array([0]), nshape=1,
             nstar=0, npolsky=skypol,
             nout=1, nback=nout, nchi2=nout, istep_out=istep_out,
             ntest=d.size, istart_test=0, istep_test=istep_d,
             o_flag=False, e_flag=False, b_flag=True, c_flag=True,
             p_flag=False, s_flag=False,
             itesttype=itesttype, rnull=0., ibadlimit=2)
        # Restore proper shapes of input arrays
        d.shape = e.shape = test.shape = d_orig_shape
        # Set output array shape, removing unnecessary dimensions
        back.shape = chi2.shape = [n for n in (norder, ndisp)
                                   if n > 1 or not squeeze_dims]

        return (back, chi2, test, ntbadl, ntbadh)


def fitsky(d, skypol=2, clip=5., itesttype=101, **kwargs):
    return extract(d, skypol=skypol, clip=clip, itesttype=itesttype,
                   tracepos=None, psf=None, fitrange=None, **kwargs)
