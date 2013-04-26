module fastexp
  ! turns out not to help much
  implicit none
  double precision :: xlim=log(1.e-10)
contains
  double precision function myexp(x)
    double precision, intent(in) :: x
    if(x.lt.xlim)then
       myexp = 0.d0
    else
       myexp = exp(x)
    endif
  end function myexp
end module fastexp

module psfpol
  integer :: ipsf, nshape
  integer, allocatable :: npoldisp(:)
end module psfpol

module frame
  ! numbers to help find a pixel in a frame (orientation dependent, etc.)
  ! kbeg,end = start,end pix spatial dir.
  ! step_in_space is number of pixels between two successive pixels along slit,
  ! step_in_disp is number of pixels between two successive pix along dispersion
  integer :: kbeg,kend, &
       istep_in_space,istep_in_disp, &
       istart_ein,istep_ein_space,istep_ein_disp, &
       istart_test,istep_test_space,istep_test_disp, &
       istart_prob,istep_prob_space,istep_prob_disp, &
       istep_out_disp,istep_out_order,istep_out_star
end module frame

module log
  integer :: ibadlimit
end module log

module slit
  integer :: ipbeg,nwidth
end module slit

module profile
  integer :: nstar,npolsky
end module profile

module shape_fdp
  double precision, allocatable :: fdp(:)
  save
end module shape_fdp

subroutine Extract_PSF(in,nin,istep_in, &
     ein,nein,istart_ein,istep_ein, &
     ipbeg,ipend,norder,kbeg,kend,lbeg,lend, &
     xord,tilt1,tilt2,tilt_flag,npoltilt, &
     ndisp,mindof,nitermax,clip1,clip2, &
     ipsf,psfguess,ifix,npar,npoldisp,nsppar, &
     nshape,nstar,npolsky, &
     psfpar,psferr,psfchi2, &
     out,eout,nout,back,nback,chi2,nchi2,istep_out, &
     test,ntest,istart_test,istep_test, &
     o_flag,e_flag,b_flag,c_flag,p_flag,s_flag,itesttype, &
     rnull,ntdiscard,ntbadl,ntbadh,nproblems,ibadlimit)
  !
  !--------- extract a long-slit or echelle spectrum by PSF fitting -----------
  !------------------- orders are extracted individually ----------------------
  !
  !--- passed on parameters:
  ! in(nin):              array with input data. All positions are relative to
  !                       first element of this array
  ! nin:                  total number of pixels in input frame
  ! istep_in(2):          step sizes in pixels for input frame. First element
  !                       should hold the offset for making one step in the 
  !                       spatial direction, second element for dispersion 
  !                       direction. Should be 1,nx for spatial direction along
  !                       the X-axis, or nx,1 for along the Y-axis
  ! ein(nein):            array with errors on input data
  ! nein:                 number of pixels in error frame
  ! istart_ein:           offset to get from spatial pixel 1 on the input frame
  !                       to pixel 1 on the error frame
  ! istep_ein(2):         step sizes for the error frame (see istep_in)
  ! ipbeg:                width down the slit for echelle mode
  ! ipend:                 ..    up   ..  ..   ..    ..    ..
  ! norder:               number of orders to extract (unity for long-slit)
  ! kbeg:                 first useful pixel in spatial direction on input frame
  ! kend:                 last    ..    ..   ..    ..      ..     ..  ..    ..
  ! lbeg(norder):         first useful pixel in dispersion dir. for each order
  ! lend(norder):         last    ..    ..   ..     ..      ..  ..   ..   ..
  ! xord(ndisp,norder)    central pos. of each order (for echelle extraction)
  ! tilt1(ndisp,norder)   first-order tilt of each order (for echelle extr.)
  ! tilt2(ndisp,norder)   second-order tilt of each order (for echelle extr.)
  ! tilt_flag             whether tilt is used or not
  ! npoltilt              order of tilt
  ! ndisp                 number of pixels along the dispersion direction
  ! mindof:               minimum # degrees of freedom that is still acceptable
  ! nitermax:             maximum number of iterations
  ! clip1:                number of deviations above which to reject pixels
  ! clip2:                  ..   ..     ..      ..    ..   ..   ..   profiles
  !                       for inclusion in the profile fit
  ! ipsf                  PSF type: 0=moffat, 1=gaussian, 2.. = multiple gauss.
  ! psfguess(npar)        Initial guesses for parameters
  ! ifix(npar)            Whether parameters should be held fixed (0) or not (1)
  ! npar                  number of parameters that define PSF
  ! npoldisp(npar)        polynomial degree with which parameters vary
  ! nsppar                # parameters for individual profile = nshape+nstar-1
  ! nshape                number of parameters defining shape of PSF
  ! nstar                 number of stars to fit (can be 0 to just fit sky)
  ! npolsky               polynomial degree to use for sky fit
  !--- passed back
  ! psfpar(npar,norder)   fitted parameters
  ! psferr(npar,norder)   corresponding errors
  ! psfchi2(norder)       chi2 for each order
  ! out(nout)             array to hold extracted orders (if o_flag is .true.)
  ! eout(nout)             ..   ..  ..  errors on same  (if e_flag is .true.)
  ! nout                  number of pixels in output and output error frames
  ! back(nback)           array to hold background values at extracted positions
  ! nback                 number of pixels in background frame
  ! chi2(nchi2)           array to hold chi2 values for extracted profiles
  ! nchi2                 number of pixels in chi2 frame
  ! istep_out             step sizes for disp. dir. and 1 order (like istep_in)
  ! test(ntest):          array to hold output test results (if itesttype.ne.0)
  ! ntest:                number of pixels in test frame
  ! istart_test:          offset for test frame (see istart_ein)
  ! istep_test(2):        step sizes for the test frame (see istep_in)
  ! o_flag                .true. if output fluxes are wanted
  ! e_flag                .true. if output errors are wanted
  ! b_flag                .true. if output background is wanted
  ! c_flag                .true. if output chi2 is wanted
  ! p_flag                .true. if PSF profile is given rather than TBD
  ! s_flag                .true. if PSF should be fit even w/ given PSF profile
  !                            (uses input PSF prof. as guesses for each order)
  ! itesttype:            indicate kind of test data that are wanted: 
  !                    0: none (no test frame needed at all)
  !                 1-99: predicted fluxes for star 1,..,n, skyconst,lin,quad...
  !                  101: predicted total fluxes
  !                  102: differences of obs-total flux
  !                  103: differences in sigma
  !                +1000: same, but with rejected pixels * factor +/-1.e15, 
  !                       with +/- for too high/too low for predicted fluxes
  ! rnull                 value used to indicate NULL result
  ! ntdiscard             number of profiles rejected for inclusion in the fit
  ! ntbadl                number of pixels rejected because they were too low
  ! ntbadh                  ..   ..   ..      ..      ..     ..   ..  ..  high
  ! nproblems             number of orders for which there were problems
  ! ibadlimit             level ranging from 0..7 to indicate what to log
  !                       (using call writeout(string)). 0 logs all, 7 nothing
  !-----------------------------------------------------------------------------
  !
  use psfpol, only : s_ipsf=>ipsf, s_nshape=>nshape, s_npoldisp=>npoldisp
  use frame, only: s_kbeg=>kbeg,s_kend=>kend, &
       istep_in_space,istep_in_disp, &
       s_istart_ein=>istart_ein,istep_ein_space,istep_ein_disp, &
       s_istart_test=>istart_test,istep_test_space,istep_test_disp, &
       istep_out_disp,istep_out_order,istep_out_star
  use log, only: addlog=>ibadlimit
  use slit, only: s_ipbeg=>ipbeg,nwidth
  use profile, only: s_nstar=>nstar,s_npolsky=>npolsky
  
  implicit none
  !---- passed-on integer parameters, and corresponding common blocks
  integer, intent(in) :: ibadlimit       ! what-to-log badness limit
  ! ... slit parameters
  integer, intent(in) :: ipbeg,ipend ! width of slit down, up
  ! ... frame parameters, indicating how to access frames
  integer, intent(in) :: kbeg,kend,istep_in(2),istart_ein,istep_ein(2), &
       istart_test,istep_test(2),istep_out(3)
  ! minimum D.O.F, max # iter., bad profile, pixel low/high, problem counters
  integer, intent(in) :: mindof,nitermax
  integer, intent(out) :: ntdiscard,ntbadl,ntbadh,nproblems
  ! sizes
  integer, intent(in) :: nin,nein,nout,nback,nchi2,ntest    ! #pix in images
  integer, intent(in) :: ndisp,norder           ! #pix in disp.dir, # orders
  ! order definition (given)
  integer, intent(in) :: lbeg(norder),lend(norder)  ! start,end pixels
  integer, intent(in) :: npoltilt                   ! degree of tilt (if used)
  real, intent(in) :: xord(ndisp,norder)            ! order position
  real, intent(in) :: tilt1(ndisp,norder),tilt2(ndisp,norder) ! first/second degree tilt
  logical, intent(in) :: tilt_flag                 ! slit is tilted or aligned
  real, intent(in) :: clip1,clip2                  ! pixel,profile clip level
  ! null value
  real, intent(in) :: rnull
  ! input images
  real, intent(in) :: in(nin),ein(nein)          ! input data and errors
  ! output images
  real, intent(out) :: out(nout),eout(nout)      ! output fluxes and errors
  real, intent(out) :: back(nback),chi2(nchi2)   ! output backgr. (sky), chi2
  real, intent(out) :: test(ntest)               ! output test image
  ! PSF definition (to be fitted)
  integer, intent(in) :: npar                      ! # of paramaters
  integer, intent(in) :: nsppar                    ! # of pars for single prof.
  integer, intent(in) :: ifix(npar)                ! fix flags
  double precision, intent(out) :: psfpar(npar,norder),psferr(npar,norder)
  real, intent(out) :: psfchi2(norder)             ! chi2 of fit for each order
  double precision, intent(in) :: psfguess(npar)   ! initial guesses
  integer, intent(in) :: ipsf,nshape               ! PSF type, #shape par
  integer, intent(in) :: npoldisp(nsppar)          ! pol.deg. for each shape p.
  integer, intent(in) :: nstar,npolsky             ! # stars, pol.deg. backgr.
  ! flags indicating what to store (output, error, sky, chi2, slit & test)
  logical, intent(in) :: o_flag,e_flag,b_flag,c_flag,p_flag,s_flag
  integer, intent(in) :: itesttype
  !
  !---- locally used variables
  character :: output*80
  integer :: ndiscard,nbadl,nbadh   ! counters for bad profs, low/high pix.
  ! ... temporary variables
  integer :: psfndof                   ! degrees of freedom
  integer :: iorder,iout,istatus
  !
  !---- initialise module variables
  !
  s_nstar = nstar
  s_npolsky = npolsky
  s_ipsf = ipsf
  s_nshape = nshape
  if(allocated(s_npoldisp))deallocate(s_npoldisp)
  allocate(s_npoldisp(nsppar))
  s_npoldisp = npoldisp
  ! ... what to log limit
  addlog = ibadlimit
  ! ... slit width
  s_ipbeg  = ipbeg
  nwidth = ipend-ipbeg+1
  ! ... parameters for how to access frames
  s_kbeg = kbeg
  s_kend = kend
  istep_in_space = istep_in(1)
  istep_in_disp = istep_in(2)
  s_istart_ein = istart_ein
  istep_ein_space = istep_ein(1)
  istep_ein_disp = istep_ein(2)
  s_istart_test = istart_test
  istep_test_space = istep_test(1)
  istep_test_disp = istep_test(2)
  istep_out_disp = istep_out(1)
  istep_out_order = istep_out(2)
  istep_out_star = istep_out(3)
  !---- initialise counters
  ntbadl = 0
  ntbadh = 0
  ntdiscard = 0
  nproblems = 0
  !
  !---- extract each order individually
  !
  orders: do iorder = 1,norder
     if(nstar.eq.0)goto 10    ! if no star being fit, just go fit sky
     if(p_flag.and..not.s_flag)goto 10   ! if profile known & fixed, go extract
     if(.not.p_flag)then  ! w/o input profile, copy guesses for parameters
        psfpar(1:npar,iorder)  = psfguess(1:npar)
     endif
     !---- fit parameters describing profile
     call Fit_Profiles(in,nin,ein,nein,xord(1,iorder), &
          tilt1(1,iorder),tilt2(1,iorder),tilt_flag,npoltilt, &
          ndisp,lbeg(iorder),lend(iorder),ifix,npar, &
          clip2,mindof,nitermax, &
          psfpar(1,iorder),psferr(1,iorder), &
          psfchi2(iorder),psfndof,ndiscard,istatus)
     if(istatus.ne.0)then
        nproblems = nproblems+1
        if(ibadlimit.le.5)then
           write(output,93)iorder
93         format('Order ',i2,' extraction FAILED.')
           call writeout(output)
        endif
        cycle orders
     endif
     ntdiscard = ntdiscard+ndiscard
     !
     !---- use fit to extract amplitudes of all stars, as well as sky
10   iout = (iorder-1)*istep_out(2)+1
     call Extract_Order(in,nin,ein,nein,xord(1,iorder), &
          tilt1(1,iorder),tilt2(1,iorder),tilt_flag,npoltilt, &
          ndisp,lbeg(iorder),lend(iorder), &
          psfpar(1,iorder),npar,clip1, &
          iout,out,eout,nout,back,nback,chi2,nchi2, &
          test,ntest,nbadl,nbadh, &
          o_flag,e_flag,b_flag,c_flag,itesttype,rnull)
     !---- report progress
     if(ibadlimit.le.5)then
        if(p_flag.or.nstar.eq.0)then
           write(output,91)iorder,nbadh,nbadl
91         format('Order ',i2,' extracted. ', &
                i6,'/',i6,' pixels too high/low.')
        else
           write(output,92)iorder,ndiscard,nbadh,nbadl
92         format('Order ',i2,' extracted. ', &
                i4,' prof.s discarded; ', &
                i6,'/',i6,' pixels too high/low.')
        endif
        call writeout(output)
     endif
     !---- adjust bad pixel counters
     ntbadl    = ntbadl+nbadl
     ntbadh    = ntbadh+nbadh
  enddo orders
end subroutine Extract_PSF
!
subroutine Fit_Profiles(in,nin,ein,nein,xord, &
     tilt1,tilt2,tilt_flag,npoltilt, &
     ndisp,lbeg,lend, &
     ifix,npar,chi2clip,mindof,nitermax, &
     psfpar,psferr,psfchi2, &
     psfndof,ndiscard,istatus)
  !
  !---- Fit the PSF for a given order
  !
  ! Important parameters:
  ! ifix(npar): which parameters to hold fixed (0=fixed)
  ! chi2clip: if a profile has an individual reduced chi2 more than 
  !           chi2clip times the total reduced chi2, it is rejected
  ! mindof:   minimum degrees of freedom needed to make fit acceptable
  !
  ! istatus:  0=OK
  !          -1=no convergence
  !          -2=too few degrees of freedom
  !
  use log
  implicit none
  ! input images
  integer, intent(in) :: nin,nein
  real, intent(in) :: in(nin),ein(nein)          ! input data and errors
  ! order definition (given)
  integer, intent(in) :: ndisp                   ! #pix in dispersion dir.
  integer, intent(in) :: lbeg,lend               ! start,end pix in each order
  real, intent(in) :: xord(ndisp),tilt1(ndisp),tilt2(ndisp) ! slit center,tilt (disppos)
  logical, intent(in) :: tilt_flag               ! slit tilted or aligned
  integer, intent(in) :: npoltilt                ! degree of tilt
  ! PSF definition (given)
  integer, intent(in) :: npar,ifix(npar)         ! #par, fix flags (0=fix)
  real, intent(in) :: chi2clip                   ! clip factor for profile chi2
  integer, intent(in) :: mindof                  ! min. dof for fit to be OK
  integer, intent(in) :: nitermax                ! max. # of iterations
  ! PSF definition (to be fitted)
  double precision, intent(inout) :: psfpar(npar),psferr(npar) ! pars,errs
  real, intent(out) :: psfchi2                   ! red. chi2
  integer, intent(out) :: psfndof,ndiscard       ! d.o.f., #discarded profiles
  integer, intent(out) :: istatus                ! !=0 -> something wrong
  ! temporary arrays for fitting routine
  double precision :: alpha(npar,npar),covar(npar,npar)
  double precision :: alambda,alambda_old,chi2,chi2old
  ! external procedure to fit 1 profile
  external :: do1profile
  ! limiting chi2 for single profile (calculated here; passed on to do1profile)
  double precision :: chi2lim
  common/profclip/chi2lim
! loop indices, etc.
  integer :: niter,ip,nparprint
  character :: fmt*80,output*320
  !
  !... initialize for iteration
  !
  istatus = 0
  alambda = -1.d0
  niter = 0
  chi2 = 99.d0
  if(ibadlimit.le.2)then
     call writeout('Iter.    Red.Chi2/dof     Lambda   par.s')
  endif
  !--- iteration loop
  fit: do niter=1,nitermax
     chi2old = chi2
     alambda_old = alambda
     chi2lim = chi2old*dble(chi2clip)
     call mrqmin(in,nin,ein,nein,xord,tilt1,tilt2,tilt_flag,npoltilt, &
          ndisp,lbeg,lend, &
          psfpar,ifix,npar, &
          covar,alpha,npar,chi2,psfndof,ndiscard, &
          do1profile,alambda)
     if(alambda.eq.0.d0)exit fit  ! alambda=0 means this was final run
     if(ibadlimit.le.2)then
        nparprint = min(npar,(len(output)-(4+12+1+8+4+10))/14)
        write(fmt,940)nparprint
        write(output,fmt)niter,chi2,psfndof,alambda, &
             (psfpar(ip),ip=1,nparprint)
940     format('(i4,f12.8,''/'',i8,4x,e10.2,',i2,'e14.6)')
        call writeout(output(1:4+12+1+8+4+10+nparprint*14))
     endif
     !... check degrees of freedom
     if(psfndof.lt.mindof)then
        call writeout('Too few data points')
        istatus = -2
        return
     endif
     !... if converged, set alambda=zero for final run (Num.Rec. for details)
     if((chi2old-chi2.lt.1.d-4*chi2.or.chi2.lt.1.d-4) &
          .and.alambda.le.alambda_old &
          .and.alambda.lt.1.d-4)then
        alambda = 0.d0
     endif
  enddo fit
  if(niter.gt.nitermax)then
     call writeout('No convergence')
     istatus = -1
     return
  endif
  ! ... set output chi2
  psfchi2 = sngl(chi2)
  forall(ip=1:npar) psferr(ip) = sqrt(covar(ip,ip))
  if(ibadlimit.le.5)then                   ! write out results
     write(output,950)psfchi2,psfndof
950  format('Fit results; chi2/dof=',F9.4,'/',i8)
     call writeout(output(:80))
     do ip = 1,npar
        write(output,955)ip,psfpar(ip),psferr(ip)
955     format('par(',i2,')=',f9.4,'+/-',f7.4)
        call writeout(output(:80))
     enddo
  endif
  return
end subroutine Fit_Profiles
!
subroutine Extract_Order(in,nin,ein,nein,xord, &
     tilt1,tilt2,tilt_flag,npoltilt, &
     ndisp,lbeg,lend,psfpar,npar,clip,iout1, &
     out,eout,nout,back,nback,chi2,nchi2, &
     test,ntest,ntbadl,ntbadh, &
     o_flag,e_flag,b_flag,c_flag,itesttype,rnull)
  !
  !--- measure fluxes and sky with given profile and offsets
  !--- fill the appropriate output and test frames
  !
  use psfpol
  use frame
  use slit
  use profile ! nstar, npolsky profile definition numbers 
  !
  implicit none
  ! input images (given)
  integer, intent(in) :: nin,nein
  real, intent(in) :: in(nin),ein(nein)      ! input data and errors
  ! order definition (given)
  integer, intent(in) :: ndisp               ! # pixels in disp. dir.
  integer, intent(in) :: lbeg,lend           ! start,end pix in each order
  real, intent(in) :: xord(ndisp),tilt1(ndisp),tilt2(ndisp) ! slit center/tilt (disppos,order)
  logical, intent(in) :: tilt_flag           ! slit tilted or aligned
  integer, intent(in) :: npoltilt            ! degree of tilt
  ! PSF definition (as found by Fit_Order)
  integer, intent(in) :: npar                ! # of PSF paramaters
  double precision, intent(in) :: psfpar(npar) ! PSF parameters
  ! clipping level
  real, intent(in) :: clip                   ! clipping level for bad pixels
  ! output images
  integer, intent(in) :: iout1               ! first pix for current order
  integer, intent(in) :: nout,nback,nchi2
  real, intent(out) :: out(nout),eout(nout)    ! output data and errors
  real, intent(out) :: back(nback),chi2(nchi2) ! output background, chi2
  integer, intent(in) :: ntest
  real, intent(out) :: test(ntest)           ! output difference/test im.
  integer, intent(out) :: ntbadl,ntbadh      ! #of low,high bad pixels
  logical, intent(in) :: o_flag,e_flag,b_flag,c_flag  ! which outputs to store
  integer, intent(in) :: itesttype           ! what test to make
  real, intent(in) :: rnull                  ! value indicating undefined
  ! variables and arrays used in calc1profile
  integer :: istatus,ndof,nbadl,nbadh
  double precision :: dchi2
  double precision, allocatable :: sppar(:)
  ! amplitudes of PSFs (and sky)
  double precision, allocatable :: amp(:),covar(:,:)
  ! loop indices, etc.
  integer :: idisp,iout,k,nsppar,nlsq
  !     
  nsppar = nshape+nstar-1          ! total # pars for single prof.
  allocate(sppar(nsppar))
  nlsq = nstar+npolsky+1               ! total number of amplitudes
  if(tilt_flag)nlsq = nlsq*2           ! for tilted slits, also derivatives
  allocate(amp(nlsq),covar(nlsq,nlsq))
  ntbadl = 0                       ! zero bad pixel counters
  ntbadh = 0
  do idisp = lbeg,lend             ! loop over all pixels in order
     if(nstar.gt.0)then
!--- find shape parameters at this position from polynomial expressions
        call getshape(idisp,ndisp,psfpar,npar,npoldisp,nsppar,sppar)
     endif
!--- use parameters to extract signals and subtract signal from input frame
     call calc1profile(in,nin,ein,nein,xord, &
          tilt1,tilt2,tilt_flag,npoltilt, &
          ndisp,idisp,lbeg,lend,nstar,npolsky, &
          ipsf,nshape,sppar,nsppar,clip, &
          amp,covar,nlsq,dchi2,ndof,nbadl,nbadh, &
          test,ntest,itesttype,istatus)
     ntbadl = ntbadl+nbadl
     ntbadh = ntbadh+nbadh
     if(istatus.eq.0)then
        !... output pixel and error
        !*** note: for tilted slits, derivatives are returned as well, but these
        !***       are not used here.  they are in amp(nstar+npolsky+1+1...)
        iout = iout1+(idisp-1)*istep_out_disp
        if(c_flag)chi2(iout) = sngl(dchi2)
        do k = 1,nstar
           if(o_flag)out(iout)  = sngl(amp(k))
           if(e_flag)eout(iout) = sngl(sqrt(covar(k,k)))
           iout = iout+istep_out_star
        enddo
        iout = iout1+(idisp-1)*istep_out_disp
        if(b_flag)then
           if(npolsky.eq.0)then
              back(iout) = sngl(amp(nstar+1))
           else if(npolsky.gt.0)then
              !*** TODO: add position dependence of sky for npolsky>0
              do k = 1,nstar
                 back(iout) = sngl(amp(nstar+1))
                 iout = iout+istep_out_star
              enddo
           endif
        endif
     else  ! if extraction failed, set output undefined
        iout = iout1+(idisp-1)*istep_out_disp
        if(c_flag)chi2(iout) = rnull
        do k = 1,nstar
           if(o_flag)out(iout) = rnull
           if(e_flag)eout(iout) = rnull
           iout = iout+istep_out_star
        enddo
        iout = iout1+(idisp-1)*istep_out_disp
        if(b_flag)then
           if(npolsky.eq.0)then
              back(iout) = rnull
           else
              do k = 1,nstar
                 back(iout) = rnull
                 iout = iout+istep_out_star
              enddo
           endif
        endif
     endif
  enddo
end subroutine Extract_Order
!
subroutine do1profile(in,nin,ein,nein,xord, &
     tilt1,tilt2,tilt_flag,npoltilt, &
     ndisp,idisp,lbeg,lend,par,npar, &
     chi2,ndof,beta,alpha,nab,istatus)
  !
  !--- front end for get1profile.  Evaluates model parameters at current pos.
  !--- and gets chi2 and the relevant derivatives
  !
  ! ... assumption is that all parameters required to describe a single profile
  ! ... can depend on dispersion position in a polynomial fashion, i.e.,
  ! ... par_i =  a_i,0 + a_i,1 * x + a_i,2 * x^2 ...
  !
  use psfpol
  use profile ! nstar, npolsky profile definition numbers 
  use shape_fdp
  use log
  !
  implicit none
  !--- INPUT
  integer, intent(in) :: nin,nein            ! # of pixels in input, error arr.
  real, intent(in) :: in(nin),ein(nein)      ! input, error arrays
  integer, intent(in) :: ndisp               ! total # of dispersion pixels
  real, intent(in) :: xord(ndisp),tilt1(ndisp),tilt2(ndisp) ! centre of slit, tilt
  logical, intent(in) :: tilt_flag           ! slit aligned or tilted
  integer, intent(in) :: npoltilt            ! degree of tilt
  integer, intent(in) :: idisp               ! dispersion pixel #
  integer, intent(in) :: lbeg,lend           ! start,end pix in each order
  !--- MODEL PARAMETERS
  integer, intent(in) :: npar                ! total # par
  double precision, intent(in) :: par(npar)  ! model parameters, including disp. dep.
  !--- OUTPUT
  double precision, intent(out) :: chi2      ! chi2
  integer, intent(out) :: ndof               ! # degrees of freedom
  integer, intent(in) :: nab                 ! derivative array dimensions
  double precision, intent(out) :: beta(nab) ! -0.5 d chi2/d par_i
  double precision, intent(out) :: alpha(nab,nab)  ! 0.5 d2 chi2/d par_i d par_j
  integer, intent(out) :: istatus            ! output status; 0=OK
  !
  double precision, allocatable :: sppar(:),tbeta(:),talpha(:,:)
  double precision :: chi2lim               ! limiting reduced chi2
  common/profclip/chi2lim
  !
  double precision :: f1
  !
  integer :: nsppar,is,iy,ip,is2,iy2,ip2
  !--- get model parameters at current profile from polynomial fits
  if(nstar.eq.0)stop'nstar=0 in do1profile. Should not happen!'
  nsppar = nshape+nstar-1   ! # model parameters for 1 profile
  allocate(sppar(nsppar),tbeta(nsppar),talpha(nsppar,nsppar))
  ! evaluate model pars at idisp
  call getshape(idisp,ndisp,par,npar,npoldisp,nsppar,sppar)
  !      if(istat.eq.-10)then
  !         print*,'idisp=',idisp
  !         print*,'Shape=',(sppar(is),is=1,nsppar)
  !         print*,'xord=',xord(idisp)
  !         if(tilt_flag)print*,'tilt1=',tilt1(idisp)
  !         if(tilt_flag.and.npoltilt.eq.2)print*,'tilt2=',tilt2(idisp)
  !      endif
  !--- use these to fit the profile
  if(npar.gt.20.and.ibadlimit.le.2)print*,'do1profile',idisp
  call get1profile(in,nin,ein,nein,xord, &
       tilt1,tilt2,tilt_flag,npoltilt, &
       ndisp,idisp,lbeg,lend,nstar,npolsky, &
       ipsf,nshape,sppar,nsppar, &
       chi2,ndof,tbeta,talpha,nsppar,istatus)
  !--- do something with bad profiles
  !*** KLUDGE !!!!  Just throw them out.
  if(chi2/dble(ndof).gt.chi2lim)istatus = -9
  !--- maybe should use weighting scheme???      chi2 = min(chi2,dble(ndof*10))
  if(istatus.ne.0)return
  !      if(istat.eq.-10)then
  !         print*,'do1profile: tbeta=',(tbeta(ip),ip=1,nsppar)
  !         print*,'do1profile: talpha=', &
  !             ((talpha(ip2,ip),ip2=1,nsppar),ip=1,nsppar)
  !      endif
  !--- evaluate derivatives for pars including the polynomial dependence on 
  !--- dispersion position
  ip = 0
  do is = 1,nsppar              ! for each individual model parameter
     do iy = 1,npoldisp(is)+1   ! loop over all its polynomial ingredients
        ip = ip+1 
        f1 = fdp(iy)
        beta(ip) = tbeta(is)*f1 ! and set beta with the appropriate factor
        !            if(istat.eq.-10)print*,is,tbeta(is),iy,fdp(iy),ip,beta(ip)
        ip2 = 0
        do is2 = 1,nsppar       ! same for alpha
           do iy2 = 1,npoldisp(is2)+1
              ip2 = ip2+1
              alpha(ip2,ip) = talpha(is2,is)*f1*fdp(iy2)
           enddo
        enddo
     enddo
  enddo
  !      if(istat.eq.-10)then
  !         print*,'do1profile: beta=',(beta(ip),ip=1,npar)
  !         print*,'do1profile: alpha=', &
  !             ((alpha(ip2,ip),ip2=1,npar),ip=1,npar)
  !      endif
end subroutine do1profile
!
subroutine getshape(idisp,ndisp,par,npar,npoldisp,nsppar,sppar)
  !
  ! Given polynomial functions defining the parameters as fu. of position
  ! along the dispersion axis, calculate the parameters at a given position
  !
  use shape_fdp
  implicit none
  integer, intent(in) :: idisp,ndisp        ! disp. position, # disp. pixels
  integer, intent(in) :: npar               ! # of fit parameters 
  double precision, intent(in) :: par(npar) ! input parameters
  integer, intent(in) :: nsppar             ! # of parameters for single profile
  integer, intent(in) :: npoldisp(nsppar)   ! polyn. deg. of single prof. pars
  double precision, intent(out) :: sppar(nsppar)  ! output single profile pars
  ! ... local variables
  integer :: imiddisp=-1,npoldispmax,is,iy,ip
  double precision :: dispfact,tsppar
  !
  save  imiddisp,dispfact,npoldispmax
  !  data imiddisp/-1/
  !---- on first call, set some useful numbers that do not change
  if(imiddisp.lt.0)then
     imiddisp = ndisp/2
     dispfact = 2.d0/dble(ndisp)
     npoldispmax = maxval(npoldisp)
     if(allocated(fdp))deallocate(fdp)  ! not guaranteed with f2py
     allocate(fdp(npoldispmax))
  endif
  !---- calculate required polynomial factors in dispersion position
  fdp(1) = 1.d0
  if(npoldispmax.gt.0)then
     fdp(2) = dble(idisp-imiddisp)*dispfact
     do iy = 3,npoldispmax+1
        fdp(iy) = fdp(iy-1)*fdp(2)
     enddo
  endif
  !---- calculate the shape parameters at the current position
  !---- using the polynomial coefficients
  ip = 0
  do is = 1,nsppar
     tsppar = 0.d0
     do iy = 1,npoldisp(is)+1
        ip = ip+1
        tsppar = tsppar+fdp(iy)*par(ip)
     enddo
     sppar(is) = tsppar
  enddo
end subroutine getshape
!
subroutine get1profile(in,nin,ein,nein,xord, &
     tilt1,tilt2,tilt_flag,npoltilt, &
     ndisp,idisp,lbeg,lend, &
     nstar,npolsky,ipsf,nshape,par,npar, &
     chi2,ndof,beta,alpha,nab,istatus)
  !
  !---- fit a slit profile with a nstar PSFs plus a polynomial for the sky
  !---- and return derivatives and cross-derivatives relative to the fit pars
  !
  ! The model is written as a linear function, with a number of amplitudes
  ! The amplitudes are found by finding the minimum
  ! 
  !   chi2 = Sum (y-STARS-SKY)**2/sigma**2, 
  !        where STARS = Sum a_star,i * PSF(x-x_i)
  !         and  SKY  = Sum a_sky,i * x**i
  !
  ! The derivatives of chi2 with respect to a_star and a_sky,i are linear in
  ! the amplitudes, and thus the solution can be found by linear least squares.
  ! Cholesky decomposition is used for speed.
  ! 
  ! The routine returns chi2 as well as the derivatives of chi2 to the fit pars
  ! The amplitudes themselves are not of interest here, as this routine is used
  ! only to find an optimum PSF.
  !
  ! The parameters are ordered as follows:
  ! par(1..nshape) = PSF parameters
  ! par(nshape+1..nshape+nstar-1) = relative positions of stars 2 and higher
  !
  ! npar:    total number of fit paramers:     npar = nshape+nstar-1
  ! nsky:    total number of sky parameters:   nsky = npolsky+1
  ! nlsq:    total number of `amplitudes' fit: nlsq = nstar+npolsky+1
  ! 
  use frame
  use slit    !... slit defs: where it starts rel. to the center and size
  !
  implicit none
  !--- INPUT
  integer, intent(in) :: nin,nein            ! # of pixels in input, error arr.
  real, intent(in) :: in(nin),ein(nein)      ! input, error arrays
  integer, intent(in) :: ndisp               ! total # of dispersion pixels
  real, intent(in) :: xord(ndisp),tilt1(ndisp),tilt2(ndisp) ! slit centre, tilt
  logical, intent(in) :: tilt_flag
  integer, intent(in) :: npoltilt
  integer, intent(in) :: idisp               ! dispersion pixel #
  integer, intent(in) :: lbeg,lend           ! start,end pix in each order
  !--- MODEL PARAMETERS
  integer, intent(in) :: nstar,npolsky       ! #stars in profile,pol.deg.sky
  integer, intent(in) :: ipsf,nshape         ! PSF type, #of shape pars
  integer, intent(in) :: npar                ! total # par = nshape+nstar-1
  double precision, intent(in) :: par(npar)  ! input parameters
  !--- OUTPUT
  double precision, intent(out) :: chi2      ! chi2
  integer, intent(out) :: ndof               ! # degrees of freedom
  integer, intent(in) :: nab                 ! derivative array dimensions
  double precision, intent(out) :: beta(nab) ! -0.5 d chi2/d par_i
  double precision, intent(out) :: alpha(nab,nab) ! 0.5 d2 chi2/d par_i d par_j
  integer, intent(out) :: istatus            ! output status; 0=OK
  !                                            -2 if not completely on chip
  real, allocatable :: px(:),py(:),pe(:),pd(:)
  double precision :: spar(nshape)
  double precision, allocatable :: p(:), dp(:,:)
  double precision, allocatable :: amp(:), dampdpar(:,:)
  double precision :: swyy
  double precision, allocatable :: swyf(:),swff(:,:),swffd(:), &
       swydf(:,:),swfdf(:,:,:),swdfdf(:,:,:,:),swydfc(:)
  !
  double precision :: x
  integer :: ip1,ip2,ix,i,j,k,jp,kp,nlsq,jstar,kstar
  double precision :: w,y,tbeta,talpha
  logical :: average_flag
  !
  !      print*,'get1profile ',idisp,(par(j),j=1,npar)
  !      print*,'ipsf,nshape,nstar=',ipsf,nshape,nstar
  if(nstar.eq.0)stop'nstar=0 in get1profile. Should not happen!'
  istatus = 0
  allocate(px(nwidth),py(nwidth),pe(nwidth),pd(nwidth))
  !
  !... get spatial positions and check that profile is complete
  !*** need to check whether average_flag should be .true.!!!
  !*** (.false. could give some problem with ip1,ip2=1,nwidth check...)
  average_flag = .true.
  call Access_Profile(0,xord,tilt1,tilt2,tilt_flag,npoltilt, &
       average_flag,idisp,lbeg,lend,ip1,ip2,pd,px)
  if(ip1.gt.1.or.ip2.lt.nwidth)then
     istatus = -1
     return
  endif
  !... get errors and check whether there are known bad pixels
  call Access_Profile(2,xord,tilt1,tilt2,tilt_flag,npoltilt, &
       average_flag,idisp,lbeg,lend,ip1,ip2,ein,pe)
  if(any(pe.le.0.))then
     istatus = minloc(pe,dim=1)
     return
  endif
  !... get actual fluxes
  call Access_Profile(1,xord,tilt1,tilt2,tilt_flag,npoltilt, &
       average_flag,idisp,lbeg,lend,ip1,ip2,in,py)
  !
  nlsq = nstar+npolsky+1                    ! total number of amplitudes
  allocate(p(nlsq),dp(nshape,nstar))
  allocate(amp(nlsq),dampdpar(nlsq,npar))
  allocate(swyf(nlsq),swff(nlsq,nlsq),swffd(nlsq),&
       swydf(nshape,nstar),swfdf(nshape,nstar,nlsq), &
       swdfdf(nshape,nstar,nshape,nstar), &
       swydfc(nlsq))
  !--- initialise sums to 0
  swyy = 0.d0 ; swyf = 0.d0 ; swff = 0.d0 
  swydf = 0.d0 ; swfdf = 0.d0 ; swdfdf = 0.d0
  p = 0.d0; dp = 0.d0
  !--- further initialisation
  p(nstar+1) = 1.d0                         ! const.part of sky indep. of x
  !                                         ! other parts zero'd above
  !--- get sums needed for linear least squares solution
  ndof = -nlsq+nwidth                     ! degrees of freedom (not used)
  slitloop: do ix = 1,nwidth              ! consider all pixels along slit
     x = dble(px(ix))                     ! pos. along slit rel. to centre
     y = dble(py(ix))                     ! data
     w = 1.d0/dble(pe(ix))**2             ! weight = 1/sigma**2
     swyy = swyy+w*y*y                    ! add data to y^2 sum
     !... get profile for first star; returns PSF in p(1), derivatives in dp
     !... no normalisation, since that gets taken care of by amplitudes
     !... (which are not interesting here; just trying to get a good PSF)
     call psf(ipsf,x,nshape,par,p(1),dp(:,1),.false.)
     if(nstar.gt.1)then
        !... get profiles for other stars; parameters the same except for pos.
        spar(2:nshape) = par(2:nshape)
        do k = 2,nstar
           !... add offset to position
           spar(1) = par(1)+par(nshape+k-1)
           call psf(ipsf,x,nshape,spar,p(k),dp(:,k),.false.)
        enddo
     endif
     !... constant sky has p=1,dp=0 (set above)
     !... further functions have x**(k-1) and dp=0 (zeroing was done above)
     if(npolsky.gt.0)then
        p(nstar+2) = x*2.d0/dble(nwidth)
        do k = nstar+3,nlsq
           p(k) = p(k-1)*p(nstar+2)
        enddo
     endif
     !--- add to sums, looping over all star probabilities and all sky pars
     ! swyf(k) = sum_i w y f_k; k=1,nlsq
     swyf = swyf+w*y*p
     !... swff(j,k) = sum_i w f_j f_k (symmetric, fill only on/below diag.)
     forall(k=1:nlsq) swff(1:k,k) = swff(1:k,k)+w*p(k)*p(1:k)
     ! swydf(i,j) = sum_i w y d f_j/d p_i; i=1,nshape; j=1,nstar
     swydf = swydf+w*y*dp
     ! swfdf(i,j,k) = sum_i w f_k  d f_j/d p_i; i=1,nshape; j=1,nstar; k=1,nlsq
     forall(k=1:nlsq) swfdf(:,:,k) = swfdf(:,:,k)+w*p(k)*dp
     ! swdfdf(i1,j,l,k) = sum_i w d f_j/d p_i d f_k / d p_l
     ! i,l=1,nshape; j,k=1,nstar (sky derivatives are always zero)
     ! (fill in on/below diag only, since this is slowest step for large nstar)
     forall(k=1:nstar, j=1:nshape) &
          swdfdf(:,1:k,j,k) = swdfdf(:,1:k,j,k)+w*dp(j,k)*dp(:,1:k)
  enddo slitloop
  !---- now Cholesky decompose the matrix
  !---- (upper triangle gets filled, and diagonal stored in swffd)
  call choldc(swff,nlsq,nlsq,swffd,istatus)
  if(istatus.ne.0)then
     print*,'Choldc failed in get1profile!!! Should not happen'
     stop
  endif
  !---- get solution for amplitudes
  call cholsl(swff,nlsq,nlsq,swffd,swyf,amp)
  !---- the amplitudes are not of interest here, but we
  !---- use solution to get derivatives of amplitudes wrt parameters
  !---- by the apropriate matrix multiplications
  do kp = 1,npar
     !---- set up new vectors to solve
     !     sum w y d f_k/dp_i 
     !     - a_j sum_k sum_i w f_j d f_k/d p_i 
     !     - a_j sum_k sum_i w f_k d f_j/d p_i 
     !     + sum_i d a_j/d p_i f_k f_j =0
     !     here the last term was also used to infer the amplitudes,
     !     so just needs new input vector
     if(kp.le.nshape)then
        !--- for shape parameters, all stars contribute
        forall(k=1:nstar) &
             swydfc(k) = swydf(kp,k)-sum(amp(:)*swfdf(kp,k,:))
        swydfc(nstar+1:) = 0.d0      ! sky terms of above all zero
        forall(k=1:nlsq) &
             swydfc(k) = swydfc(k)-sum(amp(:nstar)*swfdf(kp,:,k))
     else
        ! for offsets, only the relevant star has non-zero derivative
        swydfc = 0.d0
        k = kp-nshape+1   ! offset to star #2 is the par following nshape
        swydfc(k) = swydf(1,k) &   ! offset deriv. same as pos. deriv.
             -sum(amp(:)*swfdf(1,k,:)) &
             -sum(amp(k)*swfdf(1,:,k))
     endif
     call cholsl(swff,nlsq,nlsq,swffd,swydfc,dampdpar(:,kp))
  enddo
  !---- fill in above the diagonal from symmetry for use below
  ! ... (this removes the part set by the Cholesky decomposition)
  do k = 2,nlsq
     forall(j=1:k-1) swff(k,j) = swff(j,k)
  enddo
  do k = 2,nstar
     forall(j=1:k-1, i=1:nshape) swdfdf(:,k,i,j) = swdfdf(i,j,:,k)
  enddo
  !---- now calculate chi2
  ! chi2 = Sum (y-STARS-SKY)**2/sigma**2, 
  !     where STARS = Sum a_star,i * PSF(x-x_i)
  !       and  SKY  = Sum a_sky,i * x**i
  ! write as chi2 = Sum(y- Sum(A_k*F_k))**2/sigma**2
  ! ...  [Y**2] chi2 = sum(y**2/sigma**2) + ...
  chi2 = swyy
  ! ...  [2YAF]       +sum(-2*y*A_k*F_k/sigma**2) )
  chi2 = chi2-2.d0*sum(amp(:)*swyf(:))
  ! ! ...  [AFjAFk]   +sum( sum( A_k*F_k*A_j*F_j ))
  do k = 1,nlsq
     chi2 = chi2+amp(k)*sum(amp(:)*swff(:,k))
  enddo
  !---- calculate beta == -0.5 * d chi2 / d par
  ! ...  write as d chi2/d par = d [ Sum(y- Sum(A_k*F_k))**2/sigma**2 ]/ d par
  ! ... Here, d y/d par = 0; Have d A_k/d par, d F_k/d par. 
  ! ... Hence beta = -0.5 * sum( -2 y [d A_k F_k/d par] /sigma**2) 
  ! ...              -0.5 * sum( [d(A_k F_k)**2/d par] /sigma**2
  !                              +sum( [d(A_k F_k A_j F_j)/d par]/sigma**2 ))
  do kp = 1,npar
     ! ...           -0.5 * sum( -2 y [d A_k/d par] F_k /sigma**2) 
     tbeta = sum(dampdpar(:,kp)*swyf)
     do j = 1,nlsq ! ... symmetry is used in the sums
        tbeta = tbeta-sum(dampdpar(:,kp)*amp(j)*swff(j,:))
     enddo
     if(kp.le.nshape)then
        ! ...        -0.5 * sum( -2 y A_k [d F_k/d par] /sigma**2)
        tbeta = tbeta+sum(amp(:nstar)*swydf(kp,:nstar))
        ! ...        -0.5 * sum(  2 [d A_k/d par] A_k F_k**2 
        ! ...                     2 A_k**2 [d F_k/d par] F_k
        ! ...                     sum ( [d A_k/d par] A_j F_k F_j 
        ! ...                          +A_k d [A_j/d par] F_k F_j
        ! ...                          +A_k A_j d [F_k/d par] F_j
        ! ...                          +A_k A_j F_k d [F_j/d par] )/sigma**2 )
        do j = 1,nlsq ! ... symmetry is used in the sums
           tbeta = tbeta-sum(amp(:nstar)*amp(j)*swfdf(kp,:nstar,j))
        enddo
     else
        kstar = kp-nshape+1
        tbeta = tbeta+amp(kstar)*swydf(1,kstar) &
                -sum(amp(kstar)*amp(:)*swfdf(1,kstar,:))
     endif
     beta(kp) = tbeta
  enddo
  !---- and finally alpha == 0.5 * d2 chi2/ d par1 dpar2
  ! ... write as d2 [ Sum(y- Sum(A_k*F_k))**2/sigma**2 ]/ d par1 d par2
  ! ... Here, d y/d par = 0, and second derivatives unknown, 
  ! ... so all terms w/ y disappear
  ! ... Do have d A_k/d par, d F_k/d par.  Hence,
  ! ... alpha = 0.5* sum( [d(A_k F_k)/d par1] sum( [d(A_j F_j)/dpar2])/sigma**2
  do kp = 1,npar
     do jp = 1,kp    ! calculate on and below diag. and fill in from symmetry
        talpha = 0.d0
        if(kp.le.nshape)then
           do j = 1,nlsq
              talpha = talpha &
                   +dampdpar(j,jp)*(sum(dampdpar(:,kp)*swff(:,j)) &
                                   +sum(amp(:nstar)*swfdf(kp,:nstar,j)))
           enddo
           if(jp.le.nshape)then
              do j = 1,nstar
                 talpha = talpha &
                      +amp(j)*(sum(dampdpar(:,kp)*swfdf(jp,j,:)) &
                              +sum(amp(:nstar)*swdfdf(kp,:,jp,j)))
              enddo
           else
              jstar = jp-nshape+1
              talpha = talpha &
                      +amp(j)*(sum(dampdpar(:,kp)*swfdf(1,jstar,:)) &
                              +sum(amp(:nstar)*swdfdf(kp,:,1,jstar)))
           endif
        else
           kstar = kp-nshape+1
           talpha = talpha &
                   +dampdpar(j,jp)*(sum(dampdpar(:,kp)*swff(:,j)) &
                                   +amp(kstar)*swfdf(1,kstar,j))
           if(jp.le.nshape)then
              do j = 1,nstar
                 talpha = talpha &
                      +amp(j)*(sum(dampdpar(:,kp)*swfdf(jp,j,:)) &
                              +amp(kstar)*swdfdf(1,kstar,jp,j))
              enddo
           else
              jstar = jp-nshape+1
              talpha = talpha &
                      +amp(jstar)*(sum(dampdpar(:,kp)*swfdf(1,jstar,:)) &
                              +amp(kstar)*swdfdf(1,jstar,1,kstar))
           endif
        endif
        alpha(jp,kp) = talpha
     enddo
  enddo
  do kp = 2,npar  ! complete matrix from symmetry
     forall(jp=1:kp-1) alpha(kp,jp) = alpha(jp,kp)
  enddo
end subroutine get1profile
!
subroutine calc1profile(in,nin,ein,nein,xord, &
     tilt1,tilt2,tilt_flag,npoltilt, &
     ndisp,idisp,lbeg,lend, &
     nstar,npolsky,ipsf,nshape,par,npar,clip, &
     amp,covar,nlsq,chi2,ndof,nbadl,nbadh, &
     test,ntest,itesttype,istatus)
  !
  !---- fit a slit profile with a nstar PSFs plus a polynomial for the sky
  !---- and return chi2, amplitudes, and associated errors
  !
  ! The amplitudes are found by finding the minimum
  ! 
  !   chi2 = Sum (y-STARS-SKY)**2/sigma**2, 
  !        where STARS = Sum a_star,i * PSF(x-x_i)
  !          and  SKY  = Sum a_sky,i * x**i
  !
  ! The derivatives of chi2 with respect to a_star and a_sky,i are linear in
  ! the amplitudes, and thus the solution can be found by linear least squares. 
  ! Cholesky decomposition is used for speed.
  !
  ! itesttype:      indicate kind of test data that are wanted: 
  !                   0: none (no test frame needed at all)
  !                1-99: predicted fluxes for star 1,..,n, skyconst,lin,quad...
  !                 101: predicted total fluxes
  !                 102: differences of obs-total flux
  !                 103: differences in sigma
  !               +1000: same, but with rejected pixels * factor +/-1.e15, 
  !                      with +/- for too high/too low for predicted fluxes
  ! istatus:         0: OK; 
  !                 -1: too few good points
  !                 +1: too few points on chip
  ! 
  use frame
  use log
  use slit
  implicit none
  !---- INPUT
  integer, intent(in) :: nin,nein           ! # of pixels in in,ein
  real, intent(in) :: in(nin),ein(nein)     ! input, error arrays
  integer, intent(in) :: ndisp              ! # of pixels in xord
  real, intent(in) :: xord(ndisp),tilt1(ndisp),tilt2(ndisp) ! centre of slit position, tilt
  logical, intent(in) :: tilt_flag          ! slit tilted or aligned w/ axis
  integer, intent(in) :: npoltilt           ! degree of tilt model
  integer, intent(in) :: idisp              ! dispersion pixel #
  integer, intent(in) :: lbeg,lend          ! start,end pix in each order
  !---- MODEL PARAMETERS
  integer, intent(in) :: nstar,npolsky      ! #stars in profile,pol.deg.sky
  integer, intent(in) :: ipsf,nshape        !  PSF type, #of shape pars
  integer, intent(in) :: npar               ! # psf parameters = nshape+nstar-1
  double precision, intent(inout) :: par(npar) ! PSF par.s & offsets between stars
  real, intent(in) :: clip                  ! clipping level for bad pixels
  !---- OUTPUT
  integer, intent(in) :: nlsq             ! # lsq par array size
  double precision, intent(out) :: amp(nlsq) ! amplitudes
  double precision, intent(out) :: covar(nlsq,nlsq) ! covariances
  double precision, intent(out) :: chi2     ! chi2 of fit
  integer, intent(out) :: ndof              ! degrees of freedom
  integer, intent(out) :: nbadl,nbadh       ! #bad low,high pixels
  integer, intent(in) :: ntest              ! # of pixels in test
  real, intent(out) :: test(ntest)          ! output test image
  integer, intent(in) :: itesttype          ! type of test data required
  integer, intent(out) :: istatus           ! output status
  !---- local variables
  real, allocatable :: px(:),py(:),pe(:),pd(:),pt(:)
  double precision, allocatable :: w(:),y(:)
  double precision, allocatable :: p(:,:),dp(:)
  double precision, allocatable :: swyf(:), swff(:,:), swffd(:)
  !
  double precision :: x,dd,ddmin,ddmax,par1,spar(nshape)
  double precision :: dqmax,dqmin,dq,dqworst,dt
  integer :: ix,k,ip1,ip2,ix1,ix2
  integer :: itype,imax,imin,iworst,nlsqc
  character :: output*80
  !
  logical :: markbad
  !
  istatus = 0
  nbadl = 0
  nbadh = 0
  if(tilt_flag)then
     ddmin = 1.d0                      ! need to check range in
     ddmax = -1.d0                     ! offsets is sufficient for derivs
  endif
  allocate(px(nwidth),py(nwidth),pe(nwidth),pd(nwidth),pt(nwidth))
  allocate(w(nwidth),y(nwidth))
  allocate(p(nlsq,nwidth),dp(nshape+nstar-1))
  allocate(swyf(nlsq),swff(nlsq,nlsq),swffd(nlsq))
  !
  ! get spatial positions, data, errors for this profile
  call Access_Profile(0,xord,tilt1,tilt2,tilt_flag,npoltilt,.false., &
       idisp,lbeg,lend,ip1,ip2,pd,px)
  call Access_Profile(1,xord,tilt1,tilt2,tilt_flag,npoltilt,.false., &
       idisp,lbeg,lend,ip1,ip2,in,py)
  call Access_Profile(2,xord,tilt1,tilt2,tilt_flag,npoltilt,.false., &
       idisp,lbeg,lend,ip1,ip2,ein,pe)
  !---- see whether there are sufficient data points in the profile
  ! ... want at least half the slit, and more data than amplitudes
  if(count(pe(ip1:ip2).gt.0).le.min(nlsq,nwidth/2).and.nstar.gt.0)then
     if(ibadlimit.le.0)then
        write(output,902)idisp
902     format('Insufficient points at row',i7,'; skipping.')
        call writeout(output)
     endif
     istatus = -1
     return
  endif
  !---- what to put in testing data
  markbad = (itesttype.ge.1000)
  if(markbad)then
     itype = itesttype-1000
  else
     itype = itesttype
  endif
  !
  par1 = par(1)                        ! Star1 pos.rel.to slit centre
  !---- initialise all sums to 0
  swyf = 0.d0; swff = 0.d0
  !---- get input data and PSF probabilities for the stars and sky all along slit
  slit_in: do ix = ip1,ip2                      ! loop over slit
     w(ix) = dble(pe(ix))              ! get error
     if(w(ix).le.0.d0)then             ! if a bad pixel
        w(ix) = 0.d0                   ! weight=0 (neg.weights used to subtr.)
        cycle slit_in
     endif
     ! good pixel, add it
     x = dble(px(ix))               ! spatial pos. rel. to slit centre
     y(ix) = dble(py(ix))           ! data
     w(ix) = 1.d0/w(ix)**2          ! weight = 1/sigma**2
     if(nstar.gt.0)then
        !... get profile for 1st star; returns PSF in p(1) (derivs not used)
        !... normalise=.true., since we want amplitudes in counts
        call psf(ipsf,x,nshape,par,p(1,ix),dp,.true.)
        if(nstar.gt.1)then
           ! ... get profiles for other stars
           spar(2:nshape) = par(2:nshape)
           do k = 2,nstar
              ! ... add offset to position (other pars don't change)
              spar(1) = par1+par(nshape+k-1)
              call psf(ipsf,x,nshape,spar,p(k,ix),dp,.true.)
           enddo
        endif
     endif
     ! ... constant sky has p=1 (indep. of x), further funcs have x**(k-1)
     p(nstar+1,ix) = 1.d0
     if(npolsky.gt.0)then
        do k = nstar+2,nstar+npolsky+1
           p(k,ix) = p(k-1,ix)*x*2.d0/dble(nwidth)
        enddo
     endif
     if(tilt_flag)then
        dd = dble(pd(ix))             ! offset in dispersion direction
        ddmin = min(ddmin,dd)
        ddmax = max(ddmax,dd)
        p(nlsq/2+1:,ix) = p(:nlsq/2,ix)*dd
     endif
  enddo slit_in
  nlsqc = nlsq
  if(tilt_flag)then                      ! check range in disp. offsets
     if(ddmax-ddmin.lt.0.1d0)then        ! is large enough to make fitting
        forall(k=nlsq/2+1:nlsq)          ! sensible; if not, amp.s only
           amp(k) = 0.d0
           covar(k,k) = 0.d0
        end forall
        nlsqc = nlsq/2
     endif
  endif
  ix1 = ip1                              ! set loop beginning and end to
  ix2 = ip2                              ! go over whole profile
  !---- return point here for subtraction of bad points
  !---- these will have negative w(ix) and ix1=ix2
  clean: do
     do ix = ix1,ix2
        if(w(ix).ne.0.d0)then
           !---- add to sums, looping over all star probs and all sky pars
           swyf = swyf+w(ix)*y(ix)*p(:,ix)
           ! ... swff matrix is symmetric; fill only on and below diagonal
           forall(k=1:nlsqc) swff(1:k,k) = swff(1:k,k)+w(ix)*p(k,ix)*p(1:k,ix)
        endif
     enddo
     !---- Cholesky decompose the matrix (upper triangle gets filled)
     call choldc(swff,nlsqc,nlsq,swffd,istatus)
     !---- if it failed, it is likely because there were no pixels left
     !---- underneath the star, so just give up on the pixel
     if(istatus.ne.0)then
        if(ibadlimit.le.0)then
           write(output,900)idisp
900        format('Choldc failed at row',i7,'; giving up')
           call writeout(output)
        endif
        istatus = -1
        return
     endif
     !---- get solution for amplitudes
     call cholsl(swff,nlsqc,nlsq,swffd,swyf,amp)
     !---- check all valid pixels and find the one with the largest difference
     !---- with respect to what is expected
     dqmax = 0.d0
     dqmin = 0.d0
     do ix = ip1,ip2
        if(w(ix).gt.0.)then               !  for all good pixels, ...
           dq = y(ix)-sum(amp(:nlsqc)*p(:nlsqc,ix)) ! data - stars - sky
           dq = dq*sqrt(w(ix))            ! /sqrt(variance) -> dev. in sigma's
           if(dq.gt.dqmax)then            ! if exceeding most positive so far
              dqmax = dq                  ! use this one as most positive
              imax  = ix
           else if(dq.lt.dqmin)then       ! if exceeding most negative so far
              dqmin = dq                  ! use this one as most negative
              imin  = ix
           endif
        endif
     enddo
     !---- check whether largest difference is small enough; if so, we're done
     if(max(dqmax,-dqmin).lt.clip)exit clean
     ! if not, we'll kick it out, but should have enough points left
     if(ip2-ip1+1-nbadl-nbadh.le.nlsqc)then
        if(ibadlimit.le.0)then
           write(output,901)idisp
901        format('Too few points left at row',i7,'; giving up')
           call writeout(output)
        endif
        istatus = -1
        return
     endif
     ! ... check if deviation is positive or negative; kick out very bad 
     ! ... positive deviations first, since these are expected to occur, 
     ! ... while negative ones usually just result from them
     if(dqmax.gt.-dqmin.or.dqmax.gt.clip*clip)then
        iworst  = imax
        dqworst = dqmax
        nbadh = nbadh+1
     else
        iworst  = imin
        dqworst = dqmin
        nbadl = nbadl+1
     endif
     if(ibadlimit.le.0)then
        write(output,903)idisp,iworst,y(iworst),dqworst
903     format('Bad pixel at row',i7,';ip,y,sig:',i4,f10.2,f8.2)
        call writeout(output)
     endif
     ! ... reverse sign on weight, and go subtract this pixel from the sums
     w(iworst) = -w(iworst)
     ix1 = iworst
     ix2 = iworst
  enddo clean
  !---- calculate data for test image and chi2
  chi2 = 0.d0
  ndof = -nlsqc                       ! will add one for every good pixel
  slit_out: do ix = ip1,ip2           ! loop over slit
     ! ... for all pix, calculate difference between input and predicted value
     dq = y(ix)-sum(amp(1:nlsqc)*p(1:nlsqc,ix)) ! data - stars - sky
     if(w(ix).gt.0.d0)then            ! for good pixels
        ndof = ndof+1                 ! increase d.o.f.
        chi2 = chi2+w(ix)*dq**2       ! and update chi2
     endif
     ! ... if test data wanted...
     select case(itype)
     case(1:99)  ! individual components
        if(markbad)then   ! markbad -> predicted flux
           pt(ix) = sngl(amp(itype)*p(itype,ix))
        else              ! .not.markbad -> probability
           pt(ix) = sngl(p(itype,ix))
        endif
        cycle slit_out
     case(101)            ! predicted total flux
        dt = -dq+y(ix)
     case(102)            ! difference betw. obs & total flux
        dt = dq
     case(103)            ! difference in sigma's
        dt = dq*sqrt(abs(w(ix)))
     case default
        stop 'calc1profile: unknown test type'
     end select
     ! ... mark as bad if required
     if(markbad)then
        if(w(ix).lt.0.d0)then    ! for pixels found to be bad
           if(itype.eq.101)then
              dt = dt*sign(1.d15,dq) ! for pred. value, sign fr. diff.
           else
              dt = dt*1.d15          ! for differences, just multiply
           endif
        else if(w(ix).eq.0.d0)then ! pixels originally marked bad
           dt = -1.d30             ! set differently
        endif
     endif
     pt(ix) = sngl(dt)                ! fill test frame
  enddo slit_out
  call Access_Profile(4,xord,tilt1,tilt2,tilt_flag,npoltilt,.false., &
       idisp,lbeg,lend,ip1,ip2,test,pt)
  !---- calculate covariances (using swyf as temporary array for unit vectors)
  swyf(1:nlsqc) = 0.d0
  do k = 1,nlsqc
     swyf(k) = 1.d0
     call cholsl(swff,nlsqc,nlsq,swffd,swyf,covar(1,k))
     swyf(k) = 0.d0
  enddo
end subroutine calc1profile
!
subroutine psfmany(ipsf,n,x,np,par,y,deriv,normalise)
  implicit none
  integer, intent(in) :: ipsf,n,np
  double precision, intent(in) :: x(n),par(np)
  double precision, intent(out) :: y(n),deriv(np,n)
  logical :: normalise
  integer :: i
  do i = 1,n
     call psf(ipsf,x(i),np,par,y(i),deriv(:,i),normalise)
  enddo
end subroutine psfmany

!--- fitting routine written similarly to what is used in applic/fit
subroutine psf(ipsf,x,np,par,y,deriv,normalise)
  ! --- Evaluate the psf of type ipsf with with psf parameters par(np) 
  !     at given x.  Return psf value in y, and derivatives in deriv(np).
  !     Always: par(1)=offset, par(2)=FWHM
  !--- Only if normalise=.true., the PSF used will be normalised to unity
  !    This is needed in extraction only; in this case, no derivatives 
  !    will be returned
  !--- PSF TYPE
  !    ipsf=0: moffat
  !    ipsf=1: Gaussian
  !    ipsf>1: multiple gaussians (usually ipsf gaussians; but inferred from np)
  !    ipsf=-1: Gaussian integrated over pixel for testing
  !             (little effect: is just Gaussian box-car smoothed by 1 pix)
  implicit none
  integer, intent(in) :: ipsf,np
  double precision, intent(in) :: x,par(np)
  double precision, intent(out) :: y,deriv(np)
  logical :: normalise
  if(ipsf.eq.0)then
     call psfmoffat(x,np,par,y,deriv,normalise)
  else if(ipsf.eq.1)then
     call psfgauss(x,np,par,y,deriv)
  else if(ipsf.eq.-1)then
     call psfintgauss(x,np,par,y,deriv)
  else
     call psfmultigauss(x,np,par,y,deriv,normalise)
  endif
end subroutine psf
!
subroutine psfmoffat(x,np,par,y,deriv,normalise)
  !
  !---- Moffat function; par(1)=b=offset; par(2)=c=FWHM; par(3)=d=exponent
  !                        1
  ! ... general: y = --------------
  !                  /    x-b 2 \ d
  !                  | 1+(---)  |
  !                  \     c    /
  !                                                     1/p3
  ! ... to get c from FWHM, divide by          sqrt( 4(2    -1) )
  !
  !                                                     Gamma(d-1/2)
  ! ... to give unit area, divide whole by   c sqrt(pi) ------------
  !                                                       Gamma(d)
  !
  ! y = (1+4(2^(1/p3)-1)*((x-p1)/p2)^2)^(-p3)  = exp(-p3*ln(w2))
  !        -----wn------  -----w----
  !      ---------------w2------------
  ! 
  ! d y/d p1 = y/(...) * -p3 * wn * 2 * (x-p1)/p2 * -1/p2 = y/w2 *2*p3*wn*w/p2
  !
  ! d y/d p2 = y/(...) * -p3 * wn * 2 * (x-p1)/p2 * -(x-p1)/p2**2
  !          = y/w2 *2*p3*wn*w/p2 *w = d y/d p1 * w
  ! d y/d p3 = y*(-ln(w2) -p3/(...) 4*w*w *exp(ln2/p3) * -ln2/p3/p3)
  !          = y*(-lnw2 +4*w*w*exp(ln2/p3)*ln2/p3/w2)
  !
  implicit none
  integer, intent(in) :: np
  double precision, intent(in) :: x,par(np)
  double precision, intent(out) :: y,deriv(np)
  logical, intent(in) :: normalise
  logical :: newp3
  double precision :: p3,ln2,pi,ln2byp3,wn,wn2,w,w2,lnw2,norm
  parameter(ln2=0.6931471805599d0,pi=3.141592653589793d0)
  double precision, external :: gammln
  save p3,ln2byp3,wn2,wn,norm,newp3
  data p3/-1.d10/
  !
  if(par(3).ne.p3)then
     newp3 = .true.
     p3 = par(3)
     ln2byp3 = ln2/p3
     wn2 = exp(ln2byp3)
     wn = 4.d0*(wn2-1.d0)
     wn2 = wn2*ln2byp3
  endif
  w = (x-par(1))/par(2)
  w2 = 1.d0+wn*w*w
  lnw2 = log(w2)
  y = exp(-p3*lnw2)
  if(normalise)then
     if(newp3)then
        norm = sqrt(pi/wn)*exp(gammln(p3-0.5d0)-gammln(p3))
        newp3 = .false.
     endif
     y = y/norm/par(2)
     return
  endif
  deriv(1) = y*p3/w2*2.d0*wn*w/par(2)
  deriv(2) = deriv(1)*w
  deriv(3) = y*(4.d0*wn2*w*w/w2-lnw2)
  return
end subroutine psfmoffat
!
subroutine psfgauss(x,np,par,y,deriv)
  !--- normalised Gaussian; par(1) = offset; par(2) = FWHM
  ! y = exp(-ln2*(2(x-p1)/p2)**2))/sqrt(pi/4ln2)/p2
  !                ----w----
  ! d y/d p1 = y * -ln2 * 2 * 2(x-p1)/p2 * -2/p2 = y * 8ln2 (x-p1)/p2 /p2
  !                                                         ----w----
  ! d y/d p2 = -(y/p2) + y * -ln2 * 2 * 2(x-p1)/p2 * -2(x-p1)/p2**2
  !          = -(y/p2) + y * 8ln2 * (x-p1)/p2 /p2 * (x-p1)/p2
  !                      -------d y/d p1---------   ----w----
  implicit none
  integer, intent(in) :: np
  double precision, intent(in) :: x,par(np)
  double precision, intent(out) :: y,deriv(np)
  double precision :: ln2,sqrtpibyln16,ln16,ln256,w
  parameter(ln2=0.6931471805599d0,ln16=4.d0*ln2,ln256=8.d0*ln2, &
       sqrtpibyln16=1.06446701943d0)
  !
  w = (x-par(1))/par(2)
  y = exp(-ln16*w*w)/par(2)/sqrtpibyln16
  deriv(1) = y*ln256*w/par(2)
  deriv(2) = -y/par(2)+deriv(1)*w
  return
end subroutine psfgauss
!
subroutine psfmultigauss(x,np,par,y,deriv,normalise)
  !--- multiple gaussians; par(1)=offset1, par(2)=fwhm1
  !       par(3)=rel-amp2, par(4)=rel-offset2, par(5)=rel-fwhm2,
  !       par(6)=rel-amp3, etc.
  !       [multiplicative]     [additive]      [multiplicative]
  !       e
  !
  ! y =        exp(-ln2*(2(x-p1)/p2)**2)
  !    +par(3)*exp(-ln2*(2(x-p1-par(4))/p2/par(5))**2)
  !    +par(6)...
  !
  implicit none
  integer, intent(in) :: np
  double precision, intent(in) :: x,par(np)
  double precision, intent(out) :: y,deriv(np)
  logical, intent(in) :: normalise
  double precision :: ln2,sqrtpibyln16,ln16,ln256,w,yt,dt,p1,p2,norm
  parameter(ln2=0.6931471805599d0,ln16=4.d0*ln2,ln256=8.d0*ln2, &
       sqrtpibyln16=1.06446701943d0)
  integer :: ip
  !
  p1 = par(1)
  p2 = par(2)
  w = (x-p1)/p2
  y = exp(-ln16*w*w)
  deriv(1) = y*ln256*w/p2
  deriv(2) = deriv(1)*w
  do ip = 5,np,3
     w = (x-p1-par(ip-1))/(p2*par(ip))
     deriv(ip-2) = exp(-ln16*w*w)
     yt = par(ip-2)*deriv(ip-2)
     y = y+yt
     dt = yt*ln256*w
     deriv(ip-1) = dt/(p2*par(ip))
     deriv(ip)   = dt*w/par(ip)
     deriv(1) = deriv(1)+deriv(ip-1)
     deriv(2) = deriv(2)+dt*w/p2
  enddo
  if(normalise)then
     norm = 1.d0+sum(par(3:np:3)*par(5:np:3))
     y = y/norm/par(2)/sqrtpibyln16
  endif
end subroutine psfmultigauss
!
subroutine psfintgauss(x,np,par,y,deriv)
  !--- normalised Gaussian integrated over some width (in pixels); 
  ! par(1) = offset; par(2) = FWHM; par(3) = full width to integrate over
  ! normalised Gaussian as in psfgauss
  ! y = int_(x-p3/2)^(x+p3/2) dx exp(-(sqrt4ln2(x-p1)/p2)**2))
  !                                         
  ! y = int_(x-p3/2)^(x+p3/2) exp(-g{x}**2))/sqrt(pi/4ln2)/p2 /p3
  !
  ! using that int_x exp(-x**2) = sqrt(pi)/2 erf(x)
  ! y = [ erf(g{x+p3/2})-erf(g{x-p3/2}) ] *  sqrt(pi/16ln2) * p2
  !                                         /sqrt(pi/4ln2)  / p2 /p3
  !   = [ erf(g{x+p3/2})-erf(g{x-p3/2}) ] / 2 / p3
  !
  ! with g{x'} = sqrt4ln2*(x'-p1)/p2
  ! d erf(x)/ dx = 2/sqrt(pi) * exp(-x**2)
  ! d g/d p1 = -sqrt4ln2/p2
  ! d g/d p2 = -sqrt4ln2*(x'-p1)/p2/p2 = -g{x'}/p2
  ! d g/d x' =  sqrt4ln2/p2
  ! d y/d p1 = [exp(-g(x+p3/2)**2) -exp(-g(x-p3/2)**2)]* -sqrt4ln2/p2 
  !                                           *2/sqrt(pi)/ 2p3 
  !          = [exp(-g(x+p3/2)**2) -exp(-g(x-p3/2)**2)]* -sqrt4ln2 
  !                                                  / (p2*p3*sqrtpi) 
  ! d y/d p2 = -(y/p2) + [exp(-g(x+p3/2)**2) * g(x+p3/2)
  !                      -exp(-g(x-p3/2)**2) * g(x-p3/2)] * -1/p2 /(sqrtpi*p3)
  ! d y/d p3 = -(y/p3) + [exp(-g(x+p3/2)**2) * 1/2
  !                      -exp(-g(x-p3/2)**2) * -1/2] *sqrt4ln2/p2 /(sqrtpi*p3) 
  implicit none
  integer, intent(in) :: np
  double precision, intent(in) :: x,par(np)
  double precision, intent(out) :: y,deriv(np)
  double precision :: gm,gaussm,gp,gaussp,derivnorm
  double precision, parameter :: sqrt4ln2=1.66510922232d0, &
       sqrtpi=1.77245385091d0
  !
  gm = sqrt4ln2*(x-0.5d0*par(3)-par(1))/par(2)
  gaussm = exp(-gm**2)
  gp = sqrt4ln2*(x+0.5d0*par(3)-par(1))/par(2)
  gaussp = exp(-gp**2)
  y = (erf(gp)-erf(gm))/(2.d0*par(3))
  derivnorm = par(2)*par(3)*sqrtpi
  deriv(1) = -(gaussp-gaussm)*sqrt4ln2/derivnorm
  deriv(2) = -y/par(2)-(gaussp*gp-gaussm*gm)/derivnorm
  deriv(3) = -y/par(3)+0.5d0*(gaussp+gaussm)*sqrt4ln2/derivnorm
  return
end subroutine psfintgauss

!  (C) Copr. 1986-92 Numerical Recipes Software :)z%+.
SUBROUTINE mrqmin(in,nin,ein,nein,xord, &
     tilt1,tilt2,tilt_flag,npoltilt, &
     ndisp,lbeg,lend, &
     a,ia,ma,covar,alpha,nca,chisq,ndof,ndiscard,funcs,alamda)
  !... adapted so that it can pass on image data and positions
  INTEGER, intent(in) :: nin,nein,ma,nca,ndisp,lbeg,lend,ia(ma)
  INTEGER,intent(out)  :: ndof,ndiscard
  REAL, intent(in) :: in(nin),ein(nein),xord(ndisp),tilt1(ndisp),tilt2(ndisp)
  INTEGER, intent(in) :: npoltilt
  LOGICAL, intent(in) :: tilt_flag
  DOUBLE PRECISION, intent(inout) :: alamda,chisq,a(ma),alpha(nca,nca),covar(nca,nca)
  !U    USES covsrt,gaussj,mrqcof
  INTEGER :: j,l,mfit,ondof
  DOUBLE PRECISION :: ochisq
  double precision, allocatable :: atry(:),beta(:),da(:)
  EXTERNAL funcs
  SAVE ochisq,atry,beta,da,mfit,ondof
  if(.not.allocated(atry))allocate(atry(ma),beta(ma),da(ma))
  if(alamda.lt.0.d0)then
     mfit=0
     do j=1,ma
        if (ia(j).ne.0) mfit=mfit+1
     enddo
     alamda=1.d0
     call mrqcof(in,nin,ein,nein,xord,tilt1,tilt2,tilt_flag,npoltilt, &
          ndisp,lbeg,lend, &
          a,ia,ma,mfit,alpha,beta,nca,chisq,ndof,ndiscard, &
          funcs)
     ochisq=chisq
     atry(1:ma)=a(1:ma)
  endif
  do j=1,mfit
     covar(j,1:mfit)=alpha(j,1:mfit)
     covar(j,j)=alpha(j,j)*(1.+alamda)
     da(j)=beta(j)
  enddo
  call gaussj(covar,mfit,nca,da,1,1)
  if(alamda.eq.0.d0)then
     call covsrt(covar,nca,ma,ia,mfit)
     call covsrt(alpha,nca,ma,ia,mfit)
     return
  endif
  j=0
  do l=1,ma
     if(ia(l).ne.0) then
        j=j+1
        atry(l)=a(l)+da(j)
     endif
  enddo
  call mrqcof(in,nin,ein,nein,xord,tilt1,tilt2,tilt_flag,npoltilt, &
       ndisp,lbeg,lend, &
       atry,ia,ma,mfit,covar,da,nca,chisq,ndof,ndiscard, &
       funcs)
  if(chisq.le.ochisq.or.ndof.gt.ondof)then
     alamda=0.1d0*alamda
     ochisq=chisq
     do j=1,mfit
        alpha(j,1:mfit)=covar(j,1:mfit)
        beta(j)=da(j)
     enddo
     a(1:ma)=atry(1:ma)
  else
     alamda=10.d0*alamda
     chisq=ochisq
     ndof=ondof
  endif
  return
END SUBROUTINE mrqmin
!  (C) Copr. 1986-92 Numerical Recipes Software :)z%+.
SUBROUTINE mrqcof(in,nin,ein,nein,xord, &
     tilt1,tilt2,tilt_flag,npoltilt,ndisp, &
     lbeg,lend, &
     a,ia,ma,mfit,alpha,beta,nalp,chisq,ndof,ndiscard,funcs)
  !... adapted so that it can pass on image data and positions
  INTEGER, intent(in) :: nin,nein,ndisp,lbeg,lend,ma,mfit,nalp,ia(ma)
  INTEGER, intent(out) :: ndof,ndiscard
  REAL, intent(in) :: in(nin),ein(nein),xord(ndisp),tilt1(ndisp),tilt2(ndisp)
  LOGICAL, intent(in) :: tilt_flag
  INTEGER, intent(in) :: npoltilt
  DOUBLE PRECISION, intent(out) :: chisq,a(ma),alpha(nalp,nalp),beta(ma)
  EXTERNAL funcs
  !INTEGER :: nparmax
  !PARAMETER (nparmax=128)
  INTEGER :: i,j,k,l,m,istatus,tndof
  DOUBLE PRECISION :: tchisq,tbeta(ma),talpha(ma,ma)
  do j=1,mfit
     alpha(j,1:j)=0.d0
     beta(j)=0.d0
  enddo
  chisq=0.d0
  ndof=0
  ndiscard=0
  profile: do i=lbeg,lend
     call funcs(in,nin,ein,nein,xord,tilt1,tilt2,tilt_flag,npoltilt, &
          ndisp,i,lbeg,lend, &
          a,ma,tchisq,tndof,tbeta,talpha,ma,istatus)
     if(istatus.ne.0)then   ! don't use if something wrong
        if(istatus.lt.0)ndiscard = ndiscard+1  ! keep track of bad data
        cycle profile
     endif
     j=0
     do l=1,ma
        if(ia(l).ne.0) then
           j=j+1
           k=0
           do m=1,l
              if(ia(m).ne.0) then
                 k=k+1
                 alpha(j,k)=alpha(j,k)+talpha(l,m)
              endif
           enddo
           beta(j)=beta(j)+tbeta(l)
        endif
     enddo
     chisq=chisq+tchisq
     ndof=ndof+tndof
  enddo profile
  do j=2,mfit
     do k=1,j-1
        alpha(k,j)=alpha(j,k)
     enddo
  enddo
  chisq = chisq/dble(ndof)
  return
END SUBROUTINE mrqcof
!  (C) Copr. 1986-92 Numerical Recipes Software :)z%+.
SUBROUTINE gaussj(a,n,np,b,m,mp)
  INTEGER, intent(in) :: m,mp,n,np
  DOUBLE PRECISION, intent(inout) :: a(np,np),b(np,mp)
  INTEGER :: i,icol,irow,j,k,l,ll,indxc(n),indxr(n),ipiv(n)
  DOUBLE PRECISION :: big,dum,pivinv
  ipiv(1:n)=0
  do i=1,n
     big=0.d0
     do j=1,n
        if(ipiv(j).ne.1)then
           do k=1,n
              if (ipiv(k).eq.0) then
                 if (abs(a(j,k)).ge.big)then
                    big=abs(a(j,k))
                    irow=j
                    icol=k
                 endif
              endif
           enddo
        endif
     enddo
     ipiv(icol)=ipiv(icol)+1
     if (irow.ne.icol) then
        do l=1,n
           dum=a(irow,l)
           a(irow,l)=a(icol,l)
           a(icol,l)=dum
        enddo
        do l=1,m
           dum=b(irow,l)
           b(irow,l)=b(icol,l)
           b(icol,l)=dum
        enddo
     endif
     indxr(i)=irow
     indxc(i)=icol
     if (a(icol,icol).eq.0.d0) stop 'singular matrix in gaussj'
     pivinv=1.d0/a(icol,icol)
     a(icol,icol)=1.
     a(icol,1:n)=a(icol,1:n)*pivinv
     b(icol,1:m)=b(icol,1:m)*pivinv
     do ll=1,n
        if(ll.ne.icol)then
           dum=a(ll,icol)
           a(ll,icol)=0.d0
           a(ll,1:n)=a(ll,1:n)-a(icol,1:n)*dum
           b(ll,1:m)=b(ll,1:m)-b(icol,1:m)*dum
        endif
     enddo
  enddo
  do l=n,1,-1
     if(indxr(l).ne.indxc(l))then
        do k=1,n
           dum=a(k,indxr(l))
           a(k,indxr(l))=a(k,indxc(l))
           a(k,indxc(l))=dum
        enddo
     endif
  enddo
  return
END SUBROUTINE gaussj
!  (C) Copr. 1986-92 Numerical Recipes Software :)z%+.
SUBROUTINE covsrt(covar,npc,ma,ia,mfit)
  INTEGER, intent(in) :: ma,mfit,npc,ia(ma)
  DOUBLE PRECISION, intent(inout) :: covar(npc,npc)
  INTEGER :: i,j,k
  DOUBLE PRECISION :: swap
  do i=mfit+1,ma
     covar(i,1:i)=0.d0
     covar(1:i,i)=0.d0
  enddo
  k=mfit
  do j=ma,1,-1
     if(ia(j).ne.0)then
        do i=1,ma
           swap=covar(i,k)
           covar(i,k)=covar(i,j)
           covar(i,j)=swap
        enddo
        do i=1,ma
           swap=covar(k,i)
           covar(k,i)=covar(j,i)
           covar(j,i)=swap
        enddo
        k=k-1
     endif
  enddo
  return
END SUBROUTINE covsrt
!  (C) Copr. 1986-92 Numerical Recipes Software :)z%+.
SUBROUTINE choldc(a,n,np,p,istatus)
  IMPLICIT NONE
  INTEGER, intent(in) :: n,np
  INTEGER, intent(out) :: istatus
  DOUBLE PRECISION, intent(inout) :: a(np,np)
  DOUBLE PRECISION, intent(out) :: p(n)
  INTEGER :: i,j,k
  DOUBLE PRECISION :: sum
  istatus = 0
  do i=1,n
     do j=i,n
        sum=a(i,j)
        do k=i-1,1,-1
           sum=sum-a(i,k)*a(j,k)
        enddo
        if(i.eq.j)then
           if(sum.le.0.d0)then
              istatus = -1
              return
           endif
           p(i)=sqrt(sum)
        else
           a(j,i)=sum/p(i)
        endif
     enddo
  enddo
  return
END SUBROUTINE choldc
!  (C) Copr. 1986-92 Numerical Recipes Software :)z%+.
SUBROUTINE cholsl(a,n,np,p,b,x)
  IMPLICIT NONE
  INTEGER, intent(in) :: n,np
  DOUBLE PRECISION, intent(in) :: a(np,np),b(n),p(n)
  DOUBLE PRECISION, intent(out) :: x(n)
  INTEGER :: i,k
  DOUBLE PRECISION :: sum
  do i=1,n
     sum=b(i)
     do  k=i-1,1,-1
        sum=sum-a(i,k)*x(k)
     enddo
     x(i)=sum/p(i)
  enddo
  do i=n,1,-1
     sum=x(i)
     do k=i+1,n
        sum=sum-a(k,i)*x(k)
     enddo
     x(i)=sum/p(i)
  enddo
  return
END SUBROUTINE cholsl
!  (C) Copr. 1986-92 Numerical Recipes Software :)z%+.
FUNCTION gammln(xx)
  DOUBLE PRECISION :: gammln,xx
  INTEGER :: j
  DOUBLE PRECISION :: ser,stp,tmp,x,y,cof(6)
  SAVE cof,stp
  DATA cof,stp/76.18009172947146d0,-86.50532032941677d0, &
       24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2, &
       -.5395239384953d-5,2.5066282746310005d0/
  x=xx
  y=x
  tmp=x+5.5d0
  tmp=(x+0.5d0)*log(tmp)-tmp
  ser=1.000000000190015d0
  do j=1,6
     y=y+1.d0
     ser=ser+cof(j)/y
  enddo
  gammln=tmp+log(stp*ser/x)
  return
END FUNCTION gammln
!  (C) Copr. 1986-92 Numerical Recipes Software :)z%+.
!
subroutine Access_Profile(directive,xord, &
     tilt1,tilt2,tilt_flag,npoltilt, &
     average_flag,iout,lbeg,lend,ip1,ip2,a,p)
  !
  !----------------------------------------------------------------------------
  !
  !  Get a profile from the input frame, checking whether it is completely on
  !  the chip
  !
  ! ... parameters
  !
  ! directive [i] what to do? 
  !               0=setup, 1=get data, 2=get variance, 3=get prob., 4=write test
  ! xord()  [r]   Spatial pos. of the order in pixel units on the input frame
  !               (for echelle only)
  ! tilt1()  [r]  tilt, first degree (pixel per pixel)
  ! tilt2()  [r]  tilt, second degree (pixel per pixel^2)
  ! tilt_flag [l] whether slit is tilted or not
  ! npoltilt [i]  polynomial degree of tilt
  ! average_flag [l] .true. -> interpolate between pixels to get flux/variance
  !                            at requested xord
  !                  .false. ->get pixels w/i 0.5 pix of requested xord
  ! iout    [i]   pixel number in the dispersion direction
  ! lbeg,lend [i] start,end of order in dispersion direction
  ! ip1     [i]   start point in profile (1      unless partly off chip)
  ! ip2     [i]   end    ..   ..    ..   (nwidth   ..     ..    ..  .. )
  ! a()     [r]   image array holding data to be accessed (except for setup)
  ! p()     [r]   profile to be read or written
  !
  !----------------------------------------------------------------------------
  !
  use frame
  use slit
  implicit none
  !---- passed parameters
  integer, intent(in) :: directive,iout,lbeg,lend,npoltilt
  integer, intent(inout) :: ip1,ip2
  real, intent(in) :: xord(1),tilt1(1),tilt2(1)
  !*** intent can be in or out; really should just output indices!
  real :: a(1),p(1)
  logical, intent(in) :: tilt_flag,average_flag
  !---- other variables that are used
  integer :: i0,ip,ia0,ia,istep_space,istep_disp,it1,it2,id1,idx,ipx
  real :: x0,x,d0,d,ddc,f2,slope,a1,a2
  !*** temporary while fiddling with this!
  integer :: isample
  !
  isample = 3
  if(isample.gt.2)then
     i0 = nint(xord(iout))+ipbeg-1
     x0 = float(i0)-xord(iout)
  else
     i0 = kbeg-1
     x0 = float(iout)
  endif
  tilt: if(.not.tilt_flag)then
     !
     !**** access data for a properly aligned slit
     !
     select case(directive)
     case(0)
        !---- directive=0: setup and slit x positions
        if(isample.gt.2)then
           !---- only use the part on the chip
           ip1 = max(1,kbeg-i0)
           ip2 = min(nwidth,kend-i0)
        else
           ip1 = 1
           ip2 = nwidth
        endif
        !---- calculate spatial positions for each pixel along slit
        do ip = ip1,ip2
           p(ip) = x0+float(ip)
        enddo
        return
     case(4)   !---- put test data
        ia = (i0-1+istart_test+ip1)*istep_test_space+&
             (iout-1)*istep_test_disp+1
        istep_space = istep_test_space
        do ip = ip1,ip2
           a(ia) = p(ip)
           ia = ia+istep_space
        enddo
        return
     case(1)   !---- get input data
        ia = (i0-1+ip1)*istep_in_space+(iout-1)*istep_in_disp+1
        istep_space = istep_in_space
     case(2)   !---- get variances
        ia = (i0-1+istart_ein+ip1)*istep_ein_space+(iout-1)*istep_ein_disp+1
        istep_space = istep_ein_space
     case(3)   !---- get probabilities
        ia = (i0-1+istart_prob+ip1)*istep_prob_space+ &
             (iout-1)*istep_prob_disp+1
        istep_space = istep_prob_space
     case default
        stop'Access_Profile: unknown directive'
     end select
     do ip = ip1,ip2
        p(ip) = a(ia)
        ia = ia+istep_space
     enddo
     return
  else
     !
     !**** access data for a tilted slit
     !**** interpolate linearly to find flux at requested dispersion position
     !
     average: if(average_flag)then
        if(isample.gt.2)then
           d0 = x0*tilt1(iout)
           if(npoltilt.eq.2)d0 = d0+x0**2*tilt2(iout)
        else
           d0 = float(ipbeg-1)*tilt1(iout)
           if(npoltilt.eq.2)d0 = d0+float((ipbeg-1)**2)*tilt2(iout)
        endif
        avtilt: select case(directive)
        case(0) !---- setup and slit x positions
           if(isample.gt.2)then
              !---- find part of the slit on chip in the spatial dir.
              ip1 = max(1,kbeg-i0)
              ip2 = min(nwidth,kend-i0)
           else
              ip1 = 1
              ip2 = nwidth
           endif
           !---- slit may run off chip in dispersion dir.
           !*** NOTE: NOT updated for 2nd degree tilt!!
           ddc = (nint(xord(iout))-xord(iout))*tilt1(iout)
           it1 = nint((ddc+float(iout-lbeg))/abs(tilt1(iout))-0.5)
           it2 = nint((ddc+float(lend-1-iout))/abs(tilt1(iout))-0.5)
           if(tilt1(iout).gt.0)then
              ip1 = max(ip1,-ipbeg+1-it1)
              ip2 = min(ip2,-ipbeg+1+it2)
           else
              ip1 = max(ip1,-ipbeg+1-it2)
              ip2 = min(ip2,-ipbeg+1+it1)
           endif
           !---- calculate spatial positions for each pixel along slit
           do ip = ip1,ip2
              p(ip) = x0+float(ip)
           enddo
           return
        case(1)  !---- get input data
           ia0 = (i0-1)*istep_in_space+(iout-1)*istep_in_disp+1
           do ip = ip1,ip2
              d = d0+float(ip)*tilt1(iout)
              if(npoltilt.eq.2)d = d+float(ip*ip)*tilt2(iout)
              id1 = nint(d-0.5)
              if(npoltilt.eq.1.and. &
                   (iout+id1.lt.lbeg.or.iout+id1.ge.lend))goto 199
              f2 = d-float(id1)
              ia = ia0+ip*istep_in_space+id1*istep_in_disp
              a1 = a(ia)
              a2 = a(ia+istep_in_disp)
              if(a1.lt.-1.e10.or.a2.lt.-1.e10)then
                 p(ip) = -1.e30
              else
                 p(ip) = a1*(1.-f2)+a2*f2
              endif
           enddo
           return
        case(2)  !---- get variances
           ia0 = (i0-1+istart_ein)*istep_ein_space+(iout-1)*istep_ein_disp+1
           do ip = ip1,ip2
              d = d0+float(ip)*tilt1(iout)
              if(npoltilt.eq.2)d = d+float(ip*ip)*tilt2(iout)
              id1 = nint(d-0.5)
              if(npoltilt.eq.1.and. &
                   (iout+id1.lt.lbeg.or.iout+id1.ge.lend))goto 199
              f2 = d-float(id1)
              ia = ia0+ip*istep_ein_space+id1*istep_ein_disp
              a1 = a(ia)
              a2 = a(ia+istep_in_disp)
              if(a1.gt.1.e10)then
                 !            print*,'a1,a2:',iout,xord(iout),ip,id1,a1,a2
                 p(ip) = max(a2,1.e10)
              else if(a2.gt.1.e10)then
                 !            print*,'a1,a2:',iout,xord(iout),ip,id1,a1,a2
                 p(ip) = a1
              else
                 p(ip) = a1*(1.-f2)**2+a2*f2**2
              endif
           enddo
           return
        case(3)  !---- get input probabilities
           stop'Reading input prob. for interp. tilted slit makes no sense'
        case(4)
           stop'Writing test for interpolated tilted slit makes no sense.'
        case default
           stop'Access_Profile: unknown directive'
        end select avtilt
199     if(npoltilt.eq.1)then
           print*,'x0,d0,tilt,ip1,ip2,ip=',x0,d0,tilt1(iout),ip1,ip2,ip
        else
           print*,'x0,d0,tilt1,tilt2,ip1,ip2,ip=', &
                x0,d0,tilt1(iout),tilt2(iout),ip1,ip2,ip
        endif
        print*,'iout,d,id1,lbeg,lend=',iout,d,id1,lbeg,lend
        stop'Tilted slit/average: Out of bounds! Should not happen!'
     else
        !
        !**** access data for a tilted slit
        !**** use all pixels w/i 0.5 pix of dispersion position
        !
        tilted: select case(directive)
        case(0) !---- setup and slit x positions
           if(isample.gt.2)then
              !---- only use the part on the chip
              ip1 = max(1,kbeg-i0)
              ip2 = min(nwidth,kend-i0)
              it1 = max(lbeg,iout-1)
              it2 = min(lend,iout+1)
              slope = (xord(it2)-xord(it1))/float(it2-it1)
           else
              ip1 = 1
              ip2 = nwidth
              slope = 0.
           endif
           !---- slit may run off chip in dispersion dir.
           ! NOTE: not corrected for 2nd degree!
           ddc = (nint(xord(iout))-xord(iout))*tilt1(iout)
           it1 = nint((ddc+float(iout-lbeg))/abs(tilt1(iout)))
           it2 = nint((ddc+float(lend-1-iout))/abs(tilt1(iout)))
           if(tilt1(iout).gt.0)then  ! slanted to right
              ip1 = max(ip1,-ipbeg+1-it1)  ! bottom left can be cut off
              ip2 = min(ip2,-ipbeg+1+it2)
           else
              ip1 = max(ip1,-ipbeg+1-it2)
              ip2 = min(ip2,-ipbeg+1+it1)
           endif
           !---- going to a next column may cause one to go off chip
           !---- correct for that
           idx = nint(float(ip1-1+ipbeg)*tilt1(iout))
           ipx = nint(xord(iout+idx))-nint(xord(iout))
           if(ipx.lt.0)ip1 = max(ip1,kbeg-i0-ipx)
           idx = nint(float(ip2-1+ipbeg)*tilt1(iout))
           ipx = nint(xord(iout+idx))-nint(xord(iout))
           if(ipx.gt.0)ip2 = min(ip2,kend-i0-ipx)
           do ip = ip1,ip2
              d = float(ip-1+ipbeg)*tilt1(iout)
              if(npoltilt.eq.2)d = d+float((ip-1+ipbeg)**2)*tilt2(iout)
              idx = nint(d)
              ipx = nint(xord(iout+idx))-nint(xord(iout))
              x = x0+float(ip+ipx)  !-slope*(float(idx)-d)
              d = x*tilt1(iout)
              if(npoltilt.eq.2)d = d+x**2*tilt2(iout)
              a(ip) = float(idx)-d
              !***CHECK***
              p(ip) = x-slope*a(ip)
              !         print*,iout,ip,idx,ipx,p(ip),a(ip)
           enddo
           return
        case(4)  !---- put test values
           ia0 = (i0-1+istart_test)*istep_test_space &
                +(iout-1)*istep_test_disp+1
           istep_space = istep_test_space
           istep_disp = istep_test_disp
           do ip = ip1,ip2
              d = float(ip-1+ipbeg)*tilt1(iout)
              if(npoltilt.eq.2)d = d+float((ip-1+ipbeg)**2)*tilt2(iout)
              idx = nint(d)
              if(iout+idx.lt.lbeg.or.iout+idx.gt.lend)goto 299
              ipx = nint(xord(iout+idx))-nint(xord(iout))
              if(i0+ip+ipx.lt.kbeg.or.i0+ip+ipx.gt.kend)goto 299
              ia = ia0+(ip+ipx)*istep_space+idx*istep_disp
              a(ia) = p(ip)
           enddo
           return
        case(1)  !---- get input data
           ia0 = (i0-1)*istep_in_space+(iout-1)*istep_in_disp+1
           istep_space = istep_in_space
           istep_disp = istep_in_disp
        case(2)  !---- get variances
           ia0 = (i0-1+istart_ein)*istep_ein_space+(iout-1)*istep_ein_disp+1
           istep_space = istep_ein_space
           istep_disp = istep_ein_disp
        case(3)  !---- get input probabilities
           ia0 = (i0-1+istart_prob)*istep_prob_space &
                +(iout-1)*istep_prob_disp+1
           istep_space = istep_prob_space
           istep_disp = istep_prob_disp
        case default
           stop'Access_Profile: unknown directive'
        end select tilted
        do ip = ip1,ip2
           d = float(ip-1+ipbeg)*tilt1(iout)
           if(npoltilt.eq.2)d = d+float((ip-1+ipbeg)**2)*tilt2(iout)
           idx = nint(d)
           if(iout+idx.lt.lbeg.or.iout+idx.gt.lend)goto 299
           ipx = nint(xord(iout+idx))-nint(xord(iout))
           if(i0+ip+ipx.lt.kbeg.or.i0+ip+ipx.gt.kend)goto 299
           ia = ia0+(ip+ipx)*istep_space+idx*istep_disp
           p(ip) = a(ia)
        enddo
        return
299     print*,x0,ip1,ip2,ip
        print*,'D: iout,idx,lbeg,lend=',iout,idx,lbeg,lend
        print*,'X:i0+ip,ipx,kbeg,kend=',i0+ip,ipx,kbeg,kend
        stop'Tilted slit/close-by: Out of bounds! Should not happen!'
     endif average
  endif tilt
end subroutine Access_Profile
subroutine writeout(out)
  character :: out*(*)
  print*,out
end
