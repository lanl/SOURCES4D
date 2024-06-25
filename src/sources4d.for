c***********************************************************************
c                                                                      *
c   SOURCES: Version 4D       (Last Edited: 04/30/24 JAF)              *
c                                                                      *
c (c) 2024. Triad National Security, LLC. All rights reserved.         *
c This program was produced under U.S. Government contract             *
c 89233218CNA000001 for Los Alamos National Laboratory (LANL), which   *
c is operated by Triad National Security, LLC for the U.S. Department  *
c of Energy/National Nuclear Security Administration. All rights in    *
c the program are reserved by Triad National Security, LLC, and the    *
c U.S. Department of Energy/National Nuclear Security Administration.  *
c The Government is granted for itself and others acting on its        *
c behalf a nonexclusive, paid-up, irrevocable worldwide license in     *
c this material to reproduce, prepare. derivative works, distribute    *
c copies to the public, perform publicly and display publicly, and to  *
c permit others to do so.                                              *
c***********************************************************************
c                                                                      *
c SOURCES was written by:                                              *
c                                                                      *
c   W. B. Wilson   (LANL, T-16)                                        *
c   R. T. Perry    (LANL, ESH-12)                                      *
c   E. F. Shores   (LANL, ESH-12)                                      *
c   W. S. Charlton (University of Texas at Austin)                     *
c   T. A. Parish   (Texas A&M University)                              *
c   G. P. Estes    (LANL, X-5)                                         *
c   T. H. Brown    (LANL, ESH-12)                                      *
c   J. E. Stewart  (LANL, NIS-5)                                       *
c                                                                      *
c Nuclear Data Contributors:                                           *
c                                                                      *
c   T. R. England  (LANL, T-16)                                        *
c   E. D. Arthur   (LANL, D-DO)                                        *
c   D. G. Madland  (LANL, T-16)                                        *
c                                                                      *
c The last edit of SOURCES4C was 03/27/02 by E. F. Shores.             *
c                                                                      *
c Mods for SOURCES4D by J. A. Favorite (LANL, XCP-3), 2015, 2017.      *
c Disable "stop" on some errors.                                       *
c Print more digits.                                                   *
c Compute first derivatives of (alpha,n) sources w.r.t. densities.     *
c Increase max number of elements in subroutine homog.                 *
c Correct units in format 2230, subroutine homog.                      *
c Error checking for nng and nag.                                      *
c                                                                      *
c Mods for SOURCES4D by J. A. Favorite (LANL, XCP-3), 2018.            *
c Compute second derivatives of (alpha,n) sources w.r.t. densities.    *
c Compute first derivatives of (alpha,n) sources w.r.t. (alpha,n)      *
c  cross sections from tape3 and stopping power data from tape2.       *
c Print more digits.                                                   *
c Correct units in format 310, subroutine homog.                       *
c Initialize dummy1, dummy2, dummy3.                                   *
c Write summary spectrum info to tape8 instead of tape7.               *
c Replace "character*n variable" with "character variable*n".          *
c                                                                      *
c Mods for SOURCES4D by J. A. Favorite (LANL, XCP-3), 2019.            *
c More digits to pi in subroutine homog. Corrects negative spont. fiss.*
c  sources (fix discovered by Ruixian Fang, U. South Carolina).        *
c Likewise, replaced 2.71828 with rnapier=2.7182818284590d+0           *
c                                                                      *
c compile with                                                         *
c   ifort -o ../bin/sources4d -r8 -i8 -traceback -check bounds -static -diag-disable 8290 -diag-disable 8291 sources4d.for
c***********************************************************************
c***********************************************************************
c                                                                      *
c Overview:                                                            *
c                                                                      *
c   The SOURCES code calculates (alpha,n), spontaneous fission,        *
c and delayed neutron source magnitudes and spectra using available    *
c library data. Data sources are identified on leading comment         *
c records in the four input data library files.                        *
c                                                                      *
c   The (alpha,n) neutron sources are calculated with 4th.-order fits  *
c to stopping cross sections and linear-linear cross section interpo-  *
c lation points. Spectra are calculated with the approximation of      *
c Whitmore [see Phys. Rev. 78, 799 (1950)], assuming isotropic neutron *
c emission in the center-of-mass system. Branchings to product levels  *
c were produced with model code calculations.                          *
c https://doi.org/10.1103/PhysRev.78.799                               *
c                                                                      *
c   Spontaneous fission sources are calculated with s.f. branching and *
c s.f. nu-bar values gathered from available sources. Spectra are      *
c calculated with watt spectrum parameters.                            *
c                                                                      *
c   Delayed neutron sources and spectra are calculated with pn (i.e.,  *
c branching) and independent d.n. spectra of England, et al.           *
c (see LANL Report LA-UR-82-841).                                      *
c                                                                      *
c   First and second derivatives of (alpha,n) sources in a homogeneous *
c material can be calculated from outputs in the pdata file. SOURCES   *
c does not know the total material atom density, so user post-         *
c processing is needed.                                                *
c                                                                      *
c***********************************************************************

c***********************************************************************
c Description of File Structure
c***********************************************************************
c tape1 = user input.
c tape2 = stopping cross section expansion coefficients library.
c tape3 = targets (alpha,n) cross section library.
c tape4 = targets (alpha,n) product level branching library.
c tape5 = sources decay data library.
c tape6 = neutron source magnitudes output file.
c tape7 = absolute neutron spectra output file.
c tape8 = normalized neutron spectra output file.
c tape9 = neutron source spectra by product level output file.
c outp  = summary of input and output.  
c outp2 = supplementary output.  
c JAF for derivatives (3 lines)
c pdata = data needed to compute derivatives of the (alpha,n) source
c w.r.t densities.
c sdata = derivatives of the (alpha,n) source w.r.t nuclear data.

c***********************************************************************
c Definition of Variables (List may be incomplete)
c***********************************************************************
c
c a          = watt spec. parameter: watt(e)=w*exp(-e/a)*sinh(sqrt(be))
c adef       = default value of watt spectrum parameter a
c alam       = source nuclide decay constant (/sec).
c alph       = mass of alpha particle (4).
c aneut      = mass of neutron (1).
c apro       = mass of product nuclide.
c aps        = source rate of alphas at eal(l) from source nuclide k
c aq(k)      = density (atoms/cc) of source nuclide k
c at         = input fraction of target atoms that are nuclide i
c atar       = mass of target nuclide.
c azm(j)     = fraction of target material atoms that are element j
c b          = watt spec. parameter: watt(e)=w*exp(-e/a)*sinh(sqrt(be))
c barnu      = source nuclide spontaneous fission nu-bar value
c bdef       = default value of watt spectrum parameter b
c beamn      = neutrons/sec per microamp of beam on target
c bfsf       = source nuclide decay branching fraction for s.f.
c bx         = branching fraction of alphas at ea reacting with target i
c              and producing product level il
c c(jt)      = solid stopping cross section coef. for terms jt=1,5;
c              gas stopping cross section coef. for terms jt=6,10.
c c1         =          erf() arguments in
c c2         =          watt spectrum integral
c c3         =          calculation.
c c4         =
c c1s        =          exp() arguments in
c c2s        =          watt spectrum integral
c c3s        =          calculation.
c c4s        =
c cx(m)      = (alpha,n) cross section (mb.) at alpha energy pt. m
c czm(jt,j)  = scx(j) coef. jt of elemental constituent
c dcx        = ln of element j stopping cross section. (I think it is
c the stopping cross section, not the ln of it--JAF)
c dcxe(j,m)  = stopping power for element j, energy m
c de         = that part (mev) of dele in neutron energy group n
c dea        = alpha energy grid width (mev) used for beam or source i
c dele       = width (mev) of neutron energy range possible from alphas
c              reacting at ea with target i populating product nuclide l
c e(ip)      = energy at point ip of library cross-section evaluation
c e90        = minimum alpha energy (mev) at which neutrons may be produ
c              at 90 degrees in the center-of-mass system in reactions
c              with target i populating product level il.
c ea         = midpoint of alpha energy group
c eal(l)     = energy (mev) of source alpha l
c eamax      = maximum alpha energy (mev) of source nuclide k
c eamin      = 0.001 mev or, if greater, target (alpha,n) threshold
c ebeam      = alpha beam energy (mev), used only if nq=0
c ee(m)      = lower energy bound (mev) of alpha energy group m,
c              groups numbered in ascending alpha energy.
c eft        = collection of erf() terms used in watt spec. integral
c el(il)     = excitation energy (mev) of product nuclide level il
c elg        = ln. of alpha energy
c en(n)      = upper energy bound (mev) of neutron energy group n,
c              groups numbered in descending neutron energy
c enmax      = maximum energy of neutrons produced by alphas at ea
c              reacting with target i and populating product level il.
c enmin      = minimum energy (mev) of neutrons produced by alphas at ea
c              reacting with target i and populating level il.
c enrsep     = cross section linear interpolation intercept
c ep(ip)     = energy value ip at which branching fractions
c              f(il,ip) are evaluated for target i
c epsel      = electronic component of stopping power
c epsnu      = nuclear component of stopping power
c erg        = >1, output spectra in ascending order (e.g. for MCNP input)
c              <1, output spectra in descending order 
c f(il,ip)   = fraction of (alpha,n) reactions with target i at energy
c              e(ip) resulting in the production of product level il.
c fact       = factor used in the calculation of p(m) to account
c              for units of ee,cx,and scx.
c fal(l)     = fraction of source nuclide decays by any mode resulting
c              in the emission of an alpha at eal(l).
c fdng(nd)   = fractional delayed neutron spectrum contribution in
c              input 100kev-wide bin structure.
c fwatt      = fraction of source k s.f. watt spectrum included in
c              neutron energy group structure
c gtmg       = sum of (alpha,n) + s.f. + d.n. multigroup neutron
c              spectra values for all sources and targets
c gtmga      = sum of (alpha,n) multigroup neutron spectra values for al
c              sources and targets
c gttqan     = total (alpha,n) neutron source from all sources and targe
c gtsan(n)   = group n spectrum value of total (alpha,n) spectrum
c i          = target nuclide index
c icon       = 0, indicates last library title card
c            > 0, indicates more library title card(s) to follow.
c id         = 1, calculate source magnitudes only
c            = 2, calculate source spectra also
c idd        = 1. homogeneous problem
c            = 2. interface problem
c            = 3. beam problem
c            = 4. three-region interface problem
c idnnq      = number of source nuclides that are d.n. precursors.
c idq        = source nuclide ident from library(lsq+10*laq+10000*lzq)
c idt        = target nuclide ident from input(lst+10*lat+10000*lzt)
c il         = product nuclide level index(1=ground,2=1st excited,etc.)
c ip         = index for input cross section or level branching data
c isfnq      = number of source nuclides decaying by s.f.
c isg        = 0, use solid stopping cross sections,
c            = 1, use gas stopping cross sections if available,
c              otherwise use solid stopping cross sections.
c iz         = element z value
c j          = elemental stopping cross section constituent index.
c jdt        = target ident read from library
c jl         = number of product nuclide levels
c jn         = temp element index used in ordering
c jp         = number of product level branching data points
c jps        = number of points in library (alpha,n) cross-section
c              evaluation for nuclide i.
c jq(k)      = source nuclide k ident read from user input;
c            = state+10*mass+10000*z
c jsm(iz)    = chemical symbol of element with z=iz
c jt         = scx coefficient index
c ajwd1      =       table i column heading words for incident
c ajwd2      =       alpha beam or decay alpha source
c jzm(j)     = z value of stopping cross section elemental constituent
c k          = source nuclide index
c kn         = temp source nuclide index used in ordering
c l          = source nuclide alpha spectrum index
c laq        = source nuclide a value
c lat        = target nuclide a value
c lp         = index of level branching data frame containing ea
c lsq        = source nuclide state (1=gnd.,1=1st.iso.,etc.)
c lst        = target nuclide state (0=gnd.,1=1st.iso.,etc.)
c lzq        = source nuclide z value
c lzt        = target nuclide z value
c m          = alpha energy group index
c mm         = alpha energy group containing energy of source nuclide k
c              alpha eal(l).
c amq        = source nuclide isomer tag (blank or m)
c amt        = target nuclide isomer tag (blank or m)
c n          = neutron energy group index
c nag        = number of alpha energy groups
c nal        = number of alpha energies in source nuclide k library spec
c nd         = delayed neutron energy structure group index.
c ndn        = number of 100kev-wide delayed neutron energy groups
c              required to describe the delayed neutron spectrum
c              associated with precursor k.
c nh         = nng/2
c nng        = number of neutron spectrum energy groups
c nq         = number of source nuclides to be evaluated;if=0,source is
c              incident alpha particle beam at ebeam.
c nt         = number of target nuclides specified in user input;
c              if=0, neutrons from s.f. only.
c nz         = number of stopping cross section elemental constituents
c p(m)       = probability of an alpha at ee(m) reacting with target i
c              and populating product level il rather than slowing
c              below eamin.
c JAF 
c pi         = 3.1415926535898d+0
c pval       = interpolated value of p at eal(l) or at ebeam.
c q          = target nuclide i (alpha,n) reaction q value (mev).
c qan        = neutron source rate (#/cc/sec) due to source nuclide k
c              alphas at eal(l) reacting at some energy above eamin with
c              target i and populating product level il.
c qlev       = target nuclide i (alpha,n) reaction q value (mev) for
c              production of product nuclide level il.
c qsf        = neutron source k rate due to s.f. of source k
c r(m)       = ratio of (alpha,n) cross section to alpha stopping cross
c              at ee(m); this ratio is the integrand of p at ee(m).
c JAF 
c rnapier    = napier's constant; 2.7182818284590d+0
c rr(m)      = calculated fraction of target i product level il
c              reactions of source k alphas l occuring in alpha energy g
c s(n)       = accumulated group n neutron spectrum contribution from
c              source k alpha l on target i
c sa         = sqrt(a), used in watt spectrum integration
c san(n)     = accumulated neutron group n (alpha,n) spec. contribution
c              for source k alphas on target i.
c sbaso4     = sqrt(b * a**2 / 4)
c se1        = sqrt of maximum energy of neutron structure
c seh        = sqrt of neutron group upper energy bound
c sel        = sqrt of neutron group lower energy bound
c sen1       = sqrt of minimum energy of neutron structure
c senmax     = sqrt of maximum neutron energy from reactions at ea
c senmin     = sqrt of minimum neutron energy from reactions at ea
c scx(m)     = stopping cross section (ev/10**15 atoms/cm**2) at ee(m).
c sdn(n)     = group n contribution to delayed neutron spectrum of
c              nuclide k.
c sl(n,il)   = accumulated group n neutron spectrum from production of
c              target i product nuclide level il by source k alphas.
c slope      = cross section linear interpolation slope
c smga       = sum of source k on target i (alpha,n) multigroup
c              neutron spectrum values
c smgs       = check sum of source k multigroup neutron spectrum values
c spba       = sqrt(pi * b * a)
c ssf(n)     = group n source k s.f. neutron spectrum contribution.
c sbtqan     = total (alpha,n) source for target i / source k
c term1      =       combinations of mass and energy
c term2      =       quantities used in the calculation of
c term3      =       values of enmin, and enmax.
c thre       = threshold energy (mev) of (alpha,n) reactions on target i
c              populating product nuclide energy level il.
c tmga       = sum of target i (alpha,n) multigroup neutron spectra
c              values for all sources
c tmgs       = check sum of multigroup spectra values for all sources
c totqan     = accumulated total of qan values all sources on target i.
c totqsf     = accumulated total of qsf values.
c ts(n)      = total (alpha,n)+ s.f. neutron spectrum in group n
c tsan(n)    = total group n target i (alpha,n) neutron spectrum
c tsdn(n)    = group n contribution to total delayed neutron spectrum
c tssf(n)    = total s.f. neutron spectrum in group n
c w          = watt spec. parameter: watt(e)=w*exp(-e/a)*sinh(sqrt(be))
c x(ip)      = cross section value at library evaluation point ip.
c **********************************************************************

c=======================================================================
c  SOURCES Main Program (6/97)
c=======================================================================

      program source
      implicit none

c----------------------
c  Common Data Storage
c----------------------
      character title*7,jsm(105)*2
      real amass
      integer erg, id, idd, ios, nerr
      dimension title(11), amass(105)
      common /names/ jsm
      common /masses/ amass

c-----------------
c  File Structure
c-----------------
      nerr=0
      open(unit=1,file='tape1',status='old',iostat=ios)
      if(ios.ne.0)then
        write(*,'("error. tape1 (input file) not found.")')
        nerr=nerr+1
      end if
      open(unit=2,file='tape2',status='old',iostat=ios)
      if(ios.ne.0)then
        write(*,'("error. tape2 (data file) not found.")')
        nerr=nerr+1
      end if
      open(unit=3,file='tape3',status='old',iostat=ios)
      if(ios.ne.0)then
        write(*,'("error. tape3 (data file) not found.")')
        nerr=nerr+1
      end if
      open(unit=4,file='tape4',status='old',iostat=ios)
      if(ios.ne.0)then
        write(*,'("error. tape4 (data file) not found.")')
        nerr=nerr+1
      end if
      open(unit=5,file='tape5',status='old',iostat=ios)
      if(ios.ne.0)then
        write(*,'("error. tape5 (data file) not found.")')
        nerr=nerr+1
      end if
      if(nerr.gt.0)then
        close(1)
        close(2)
        close(3)
        close(4)
        close(5)
        stop
      end if
      open(unit=6,file='tape6',status='unknown')
      open(unit=7,file='tape7',status='unknown')
      open(unit=8,file='tape8',status='unknown')
      open(unit=10,file='tape9',status='unknown')
      open(unit=11,file='outp',status='unknown')
      open(unit=12,file='outp2',status='unknown')
      
c------------------
c  Read From Tape1
c------------------
      read (1,17,end=99) title
      read(1,*)idd,id,erg
      if (idd.eq.1) then
       call homog(title,idd,id,erg)
      elseif (idd.eq.2) then
       call interf(title,idd,id,erg)
      elseif (idd.eq.3) then
       call homog(title,idd,id,erg)
      elseif (idd.eq.4) then
       call three(title,idd,id,erg)
      else
       stop 'Error in idd input!'
      endif
   99 continue

c-------------------------
c  Close All Output Files
c-------------------------
      close(unit=1)
      close(unit=2)
      close(unit=3)
      close(unit=4)
      close(unit=5)
      close(unit=6)
      if (id.eq.1) then
       close(unit=7, status='delete')
       close(unit=8, status='delete')
       close(unit=10, status='delete')
       close(unit=12, status='delete')
      elseif (idd.eq.2) then
       close(unit=10, status='delete')
      else
       close(unit=7)
       close(unit=8)
       close(unit=10)
       close(unit=12)
      endif
c JAF for derivatives
      if(idd.eq.1.or.idd.eq.3)then
        close(unit=13)
        close(unit=14)
      else
        close(unit=13, status='delete')
        close(unit=14, status='delete')
      end if

c--------------------
c  Format Statement
c--------------------
   17 format(11a7)

c------------------------------
c  End of Execution of SOURCES
c------------------------------
      Stop 'Normal Execution'
      End



c=======================================================================
c  Homogeneous Problem Subroutine (6/97)
c=======================================================================

      subroutine homog(title,idd,id,erg)



c---------
c Storage
c---------
      character title*7,jsm(105)*2,amt*1,amq*1,ajwd1*1,ajwd2*4,ajwd3*4
      integer erg
c JAF introduce parameters
      integer maxel,maxlv,maxbr,maxnq,maxnag,maxnng,maxnd,maxxs,maxnal
      parameter (maxel=100,maxlv=20,maxbr=200,maxnq=300,maxnag=2001,
     1 maxnng=750,maxnd=1000,maxxs=1100,maxnal=30)
c JAF replaced some 20's in this list with maxel, some with maxlv. some might be wrong.
c JAF replaced both 200's in this list with maxbr.
c JAF replaced all 300's in this list with maxnq. some might be wrong.
c JAF replaced all 4001's in this list with maxnag.
c JAF replaced all 750's in this list with maxnng.
c JAF replaced the 1000 in this list with maxnd.
c JAF replaced both 1100's in this list with maxxs.
c JAF replaced both 30's in this list with maxnal.
      dimension title(11), jzm(maxel), azm(maxel), czm(9,maxel), c(14),
     1 jq(maxnq),
     1 aq(maxnq),el(maxlv), ep(maxbr),f(maxlv,maxbr), eal(maxnal),
     2 fal(maxnal), e(maxxs), x(maxxs),
     2 ee(maxnag), scx(maxnag), en(maxnng+1), san(maxnng), ssf(maxnng),
     3 ts(maxnng), cx(maxnag), r(maxnag), rr(maxnag), p(maxnag),
     3 tsan(maxnng),
     4 fdng(maxnd), tsdn(maxnng), sdn(maxnng), gtsan(maxnng), s(maxnng),
     4 tssf(maxnng),
     5 totlev(maxlv),sl(maxnng,maxlv),etopl(maxlv),ebotl(maxlv)
c JAF for derivatives
c san_der dimensions from jzm(maxel),jq(maxnq). will go to nz,nq
c iz_tar dimensions from jzm(maxel). will go to nz
c iz_src dimensions from jq(maxnq). will go to nq
c dcxe dimensions from jzm(maxel),scx(maxnag). will go to nz,nagp1
c r_drv_n dimensions from r(maxnag). will go to nagp1
c p_drv_n dimensions from jzm(maxel),p(maxnag). will go to nz,nagp1
c pval_drv_n dimensions from jzm(maxel). will go to nz
c rk_drv_n dimensions from jzm(maxel). will go to nz
c rk_drv_n_tot dimensions from jzm(maxel). will go to nz
c fp_tot dimensions from jzm(maxel). will go to nz
c r_drv2_n dimensions from r(maxnag). will go to nagp1
c p_drv2_n dimensions from jzm(maxel),p(maxnag). will go to nz,nagp1
c pval_drv2_n dimensions from jzm(maxel). will go to nz
c rk_drv2_n dimensions from jzm(maxel). will go to nz
c dqdd(maxnag); will go to nagp1
c dqdx1(maxnag); will go to nagp1
c depsdb(6,maxel,maxnag); will go to 6,nz,nagp1
c depsdc(9,maxel,maxnag); will go to 9,nz,nagp1
c dq3c(maxnng,maxnag); will go to nng,nagp1
c dq3s(maxnng,maxnag,maxel); will go to nng,nagp1,nz
      dimension san_der(maxel,maxnq),iz_tar(maxel),iz_src(maxnq),
     1 dcxe(maxel,maxnag),niso(2),r_drv_n(maxnag),
     2 p_drv_n(maxel,maxnag),pval_drv_n(maxel),rk_drv_n(maxel),
c    3 rk_drv_n_tot(maxel),fp_tot(2,maxel),
     4 r_drv2_n(maxnag),p_drv2_n(maxel,maxel,maxnag),
     5 pval_drv2_n(maxel,maxel),rk_drv2_n(maxel,maxel),
     6 dqdd(maxnag),dqdx1(maxnag),depsdb(6,maxel,maxnag),
     7 depsdc(9,maxel,maxnag),dqdb(6),dqdc(9,maxel),
     8 dq3c(maxnng,maxnag),dq3s(maxnng,maxnag,maxel)
c JAF end mod
      common /names/ jsm
      common /masses/ amass(105)
c     dimension frclev(maxel) 
c JAF for derivatives
c l_sdata = .true./.false. calculate and write sdata/don't
c l_2nd = .true./.false. calculate and write 2nd derivatives/don't
      logical l_sdata,l_2nd
      l_sdata=.false.
      l_2nd=.true.
      open(unit=13,file='pdata',status='unknown')
      if(l_sdata)open(unit=14,file='sdata',status='unknown')
c JAF end mod


c-------------------
c Format Statements
c-------------------
   10 format(i3,11a7)
   15 format(3x,11a7)
   17 format(11a7)
   20 format(2i3)
   30 format(i8,e12.5)
   40 format(i3,8e8.3,f10.3,/,3x,5f8.3)
   50 format(' failed to find stopping cross section coefficients on',
     1' tape 2 for iz=',i3,' stop')
   60 format(6e12.5)
   70 format(/,1x,'Neutron Multigroup Structure (MeV)')
c JAF for derivatives: 2 lines, more digits. (was 10.3)
c also in 81 remove '.', and change ' MeV.' to ' MeV'
   80 format(1p8e15.7)
   81 format(35x,'Total (all groups): ',1pE15.7,' neutrons',a1,'sec-',
     1a4,a4,/,31x,'Average Neutron Energy: ',1pE15.7,' MeV')
   83 format(' Title    gttqan  ebaran   totqsf  ebarsf   totqdn  ',
     1 'ebardn   qtotal  ebrall')
   84 format(1x,a7,4(1pe10.3,7x))
   85 format(1x,a7,4(1pe10.3,0pf7.3))
   90 format (i3,5e12.5)
   95 format(2i6)
   96 format(/,26x,'Neutron Source Magnitudes',/,26x,25(1h~),/)
   97 format(/,23x,'Absolute Neutron Source Spectra',/,23x,31(1h~),/)
   98 format(/,22x,'Normalized Neutron Source Spectra',/,22x,
     133(1h~),/)
   99 format(/,16x,'Neutron Source Spectra by Nuclide Energy Level',
     1/,16x,46(1h~),/)
  100 format(//,1h1,30x,'Table I',/,31x,7(1h=),//,12x,'(alpha,n) ',
     1'Neutron Production by Target Per Source Alpha',///,17x,
     2'target    alpha   alpha   alphas/sec     p(e)     neuts/sec',
     3/,7x,'target  atom frac.  source  energy  ',a1,a4,a4,2x,
     4' neut/alpha  ',a1,a4,a4,/,1h+,6x,6(1h_),2x,10(1h_),2x,6(1h_),2x,
c JAF more digits
c    5 6(1h_),3(12h  __________))
     5 6(1h_),2(12h  __________),14h  ____________)
  102 format(1h+,77x,37hfractional ground/excited level split,/,
     1 1h+,77x,49(1h_))
  105 format(i8,i4)
  110 format (i8,2i4)
  120 format (8e10.3)
  130 format (52h failed to find alpha-n cross section on tape 3 for ,i8
c JAF no need to stop
c    1 ,5h stop)
     1         )
c JAF end mod
  140 format (i8,2i3)
  150 format (66h failed to find comp.nucl.neut.emiss.lev. branching on 
c JAF no need to stop
c    1tape 4 for ,i8,5h stop)
     1tape 4 for ,i8        )
c JAF end mod
  160 format (56h failed to find alpha/s.f./dn source data on tape 5 for
c JAF no need to stop
c    1 ,i8,5h stop)
     1 ,i8        )
c JAF end mod
  170 format (17h alpha energy of ,1pe12.5,62h exceeds current 6.5 mev.
     1 limit of most cross sections, stop. )
  180 format (7h target,i8,40h cross-section data not available at ee(
     1 ,i4,3h)= ,1pe12.5,6h stop.)
  190 format(7x,a2,i3,a1,1pe12.4,3x,5hbeam ,0pf8.3,1p3e12.4)
c JAF more digits
c 200 format(33x,f8.3,1p3e12.4)
  200 format(33x,f8.3,1p2e12.4,e14.6)
  205 format(1h+,(76x,7f7.5,/,1x))
c JAF more digits
c 210 format(7x,a2,i3,a1,1pe12.4,2x,a2,i3,a1,0pf8.3,1p3e12.4)
  210 format(7x,a2,i3,a1,1pe12.4,2x,a2,i3,a1,0pf8.3,1p2e12.4,e14.6)
  220 format (14h alpha energy ,1pe12.5,45h exceeds range of level branc
     1hing data, stop.)
  230 format (//,' (a,n) neutrons/sec-microamp ',f6.3,' mev a on ',a2,
     1i3,a1,' in target')
  240 format (//,' (a,n) neutrons/sec-cc from ',1pe12.5,' at/cc ',a2,
     1i3,a1,' alphas on ',a2,i3,a1,' in target')
c JAF more digits
c 250 format (1h+,66x,10(1h_),/,55x,'Total:',4x,1pe12.4,/)
c 251 format (1h+,66x,10(1h_),/,41x,'Total (this target):',
c    14x,1pe12.4,/)
c 252 format (1h+,66x,10(1h_),/,41x,'Total (all targets):',
c    14x,1pe12.4,/)
  250 format (1h+,66x,12(1h_),/,55x,'Total:',4x,1pe14.6,/)
  251 format (1h+,66x,12(1h_),/,41x,'Total (this target):',
     14x,1pe14.6,/)
  252 format (1h+,66x,12(1h_),/,41x,'Total (all targets):',
     14x,1pe14.6,/)
  255 format(1h )
  260 format (///,' Total (alpha,n) neutron spectrum this target')
  261 format(//,' Neutron spectrum from ',a2,i3,a1,' alphas on ',a2,i3,
     1a1,' via the ',f6.2,'-MeV product level')
  262 format(//,' Neutron spectrum from ',F5.2,' MeV alphas on ',a2,i3,
     1a1,' via the ',f6.2,'-MeV product level')
  270 format(///,' Grand total (alpha,n) neutron spectrum, all targets,'
     1,1x,'all sources')
  280 format (////,1h2,35x,'Table II',/,36x,8(1h=),//,21x,
     1'Spontaneous Fission Neutron Production',///,8x,2(8h source ),
     2'atoms  dk constant  sf decay    nu    neutrons',/,8x,'nuclide
     3    per cm**3    (/second)   branching   bar   sec/cm**3',
     4/,1h+,7x,7(1h_),2x,12(1h_),2x,11(1h_),2x,9(1h_),2x,5(1h_),2x,
c JAF more digits
c    59(1h_))
c 290 format (9x,a2,i3,a1,1p2e13.4,e12.3,0pf7.3,1pe11.3)
     512(1h_))
  290 format (9x,a2,i3,a1,1p2e13.4,e12.3,0pf7.3,1pe14.6)
  300 format (39h no watt spectrum parameters found for ,a2,i3,a1,22h, d
     1efault values used:,/,3x,3ha= ,1pe12.5,5h, b= ,e12.5)
c JAF correct units; S.F. was neutrons/cc
  310 format (//,' S.F. neutrons/sec-cc from ',1pe12.5,' at/cc ',
     1 a2,i3,a1)
c JAF more digits
c 320 format (1h+,61x,9(1h_),/,52x,'Total:',4x,1pe9.3)
  320 format (1h+,61x,12(1h_),/,52x,'Total:',4x,1pe12.6)
  330 format (///,' Total S.F. neutron spectrum')
  340 format (////,1h3,35x,'Table III',/,36x,9(1h=),//,22x,
     1'Delayed Neutron Production',///,8x,2(8h source ),
     2'atoms  dk constant  branching neutrons'/,/,8x,
     3'nuclide   per cm**3     (/second)     (pn)sec/cm**3',
     4/,1h+,7x,7(1h_),2x,12(1h_),2x,11(1h_),2x,9(1h_),2x,9(1h_))
  350 format (9x,a2,i3,a1,1p2e13.4,e12.3,e11.3)
  360 format (//,' Delayed neutrons/sec-cc from',1pe12.5,' at/cc ',
     1a2,i3,a1)
  370 format (1h+,54x,9(1h_),/,53x,1pe11.3)
  380 format (///,' Total delayed neutron spectrum')
  390 format (////,' Total Neutron Spectrum')

 1900 format('SOURCES 4D Calculation',/,'<<<<<<<<<<<>>>>>>>>>>>',//)
 2000 format(/,22x,'Summary of Input')
 2010 format(22x,'================')
 2020 format(' ')
 2021 format(/,80(1h-))
 2030 format('Title: ', 11a7)
 2040 format('Homogeneous problem input (idd =',i2,')')
 2041 format('Ascending energy structure for output (erg =',i3,')')
 2042 format('Descending energy structure for output (erg =',i3,')')
 2043 format('Interface problem input (idd =',i2,')')
 2047 format('Beam problem input (idd =',i2,')')
 2048 format('Magnitudes and spectra computed (id =',i2,')')
 2049 format('Magnitudes only computed (id =',i2,')')
 2050 format('Number of elemental constituents:',i3)
 2060 format('Solid stopping cross-sections used (isg =',i2,')')
 2070 format('Gas stopping cross-sections used (isg =',i2,')')
 2080 format('Elemental Constituents:')
 2090 format(' ')
 2100 format('  Z-value  Atom Fraction')
 2110 format('  -------  -------------')
 2120 format(4x,i3,5x,f12.10)
 2130 format('Number of neutron spectrum energy groups:', i4)
 2140 format('Maximum neutron energy is ',1pE10.3,' MeV.')
 2150 format('Minimum neutron energy is ',1pE10.3,' MeV.')
 2160 format('Energy Group Structure:')
 2170 format('  Group  Upper-Bound  Lower-Bound')
 2180 format('  -----  -----------  -----------')
 2190 format(2x,i4,4x,1pE10.3,3x,1pE10.3)
 2193 format(/,'Warning: First Group Upper-Bound Does Not Coincide ',
     1'With Maximum Neutron Energy.')
 2195 format(/,'Warning: Energy Group Structure Not Entered in ',
     1'Decreasing Order.',/,5x,'Group ',i3,' Upper-Bound is Greater ',
     2'Than Group ',i3,' Upper-Bound.') 
 2200 format('Number of source nuclides to be evaluated:',i3)
 2210 format('Alpha beam energy is ',1pE10.3,' Mev.')
 2220 format('Source Nuclides:')
c JAF correct units; was g/cc
 2230 format('    ZAID     Atom Density (atoms/cc)')
 2240 format('    ----     -----------------------')
 2250 format(3x,i6,4x,1pE10.3)
 2260 format('Number of target nuclides to be used:',i3)
 2270 format(i5,' Alpha energy groups used.')
 2280 format('Target Nuclides:')
 2290 format('    ZAID     Atom Fraction')
 2300 format('    ----     -------------')
 2310 format(3x,i6,4x,1pE10.3)

 3000 format(//,22x,'Summary of Output')
 3010 format(22x,'=================')
 3020 format('Total (alpha,n) neutron source from all sources and ',
     1'targets: ',1pE10.3,' n/sec-cm^3.')
 3021 format('Total (alpha,n) neutron source from all sources and ',
     1'targets: ',1pE10.3,' n/sec-microamp.')
 3022 format(/,'Total spontaneous fission neutron source from all ',
     1'sources and targets: ',1pE10.3,' n/sec-cm^3.')
 3024 format(/,'Total delayed neutron source from all sources and ',
     1'targets: ',1pE10.3,' n/sec-cm^3.')
 3026 format(/,'Total neutron source from all sources and ',
     1'targets: ',1pE10.3,' n/sec-cm^3.')
 3030 format(/,'Average (alpha,n) neutron energy: ',1pE10.3,' MeV.')
 3032 format(/,'Average spontaneous fission neutron energy: ',
     11pE12.5,' MeV.')
 3034 format(/,'Average delayed neutron energy: ',1pE10.3,' MeV.')
 3036 format(/,'Average neutron energy from all sources: ',1pE10.3,
     1' MeV.')
 3050 format(/,'Normalized Neutron Energy Spectrum by Energy Group ',
     1'for All Sources:',/'(Note: descending group structure is',
     2' independent of erg record.)',//,10x,'Group',5x,'Contribution',
     3/,10x,'-----',5x,'------------')
 3051 format(11x,i3,7x,1pE10.3)
 3500 format(/,'Portion of Total Neutron Source Rate Accounted for ',
     1'in the Total Energy Spectrum: ',F5.1,'%.')
 3501 format(/,'Portion of (a,n) Neutron Source Rate Accounted for ',
     1'in the (a,n) Energy Spectrum: ',F5.1,'%.')
 3502 format(/,'Portion of Spontaneous Fission Neutron Source Rate ',
     1'Accounted For in the Spontaneous Fission Energy Spectrum: ',
     2F5.1,'%.')
 3503 format(/,'Portion of Delayed Neutron Neutron Source Rate ',
     1'Accounted for in the Delayed Neutron Energy Spectrum: ',
     2F5.1,'%.')
 3510 format(/,'Warning: Energy Group Structure May Be Poorly Chosen.')
 3600 format(' Title: ', 11a7) 
 3601 format(' Neutron Source Spectra in Columnar Format.'/
     1' Examine other tape files for additional information.'/
     2' Recall the Energy values are bin bounds.'/)  
 3605 format(/,'            Normalized Totals (neutrons/sec-basis)')
 3610 format('            ==========================================')
 3611 format(/,'            Absolute Totals (neutrons/sec-basis)')
c JAF more digits
c3615 format(' E (MeV)    (a,n)      (sf)       (dn)       Total')
c3620 format(1x,9('-'),2x,9('-'),2x,9('-'),2x,9('-'),2x,9('-'))
c3630 format(1p1e10.3,1p4e11.3)
c3635 format(' Total    ',1p4e11.3)
 3615 format(' E (MeV)       (a,n)         (sf)          (dn)',
     1 '          Total')
 3620 format(1x,12('-'),2x,12('-'),2x,12('-'),2x,12('-'),2x,12('-'))
 3630 format(1p1e13.6,1p4e14.6)
 3635 format(' Total       ',1p4e14.6)
c***********************************************************************

c-------------
c  Parameters
c-------------
      rewind 1
c JAF d0 added to reals here
c JAF pi changed from 3.14159 to 3.1415926535898d+0
c JAF added rnapier and changed e from to 2.71828 to 2.7182818284590d+0
      aneut=1.d0
      alph=4.d0
      zalp=2.d0
      pi=3.1415926535898d+0
      rnapier=2.7182818284590d+0
      adef=0.82d0
      bdef=4.6d0
c JAF line label 400 removed
      totqsf=0.d0
      isfnq=0
      idnnq=0
      gttqan=0.d0

c-----------------------------------
c  Read Input Parameters from Tape1
c-----------------------------------
      read (1,17,end=1440) title
      read(1,*)idd,id,erg
      write(6,1900)
      write(7,1900)
      write(8,1900)
      write(10,1900)
      write(11,1900)
      write(6,96)
      write(6,2030)title
      write(11,2000)
      write(11,2010)
      write(11,2020)
      write(11,2030)title
      if (idd.eq.1) then
        write(11,2040)idd
      elseif (idd.eq.3) then
        write(11,2047)idd
      endif
c JAF for derivatives (12 lines)
      if (idd.eq.1 .or. idd.eq.3) then
        write(13,1900)
        write(13,2020)
        write(13,2030)title
        write(13,2040)idd
        write(13,2020)
        if(l_sdata)then
          write(14,1900)
          write(14,2020)
          write(14,2030)title
          write(14,2040)idd
          write(14,2020)
        endif
      endif
      if (id.eq.1) then
        write(11,2049)id
      elseif (id.eq.2) then
        write(11,2048)id
        if (erg.ge.1) then
          write(11,2041)erg
        elseif (erg.le.-1) then
          write(11,2042)erg
        else
          stop 'Error in erg input!'
        endif       
      else
        stop 'Error in id input!'
      endif
  410 if(id.eq.2) then
        write(7,97)
        write(7,2030)title
        write(8,98)
        write(8,2030) title
        write(10,99)
        write(10,2030)title
      endif

      read (1,*) nz,isg
c JAF for derivatives (2 lines)
      write(13,'(/,"number of stopping elems:",i6)')nz
      if(l_sdata)write(14,'(/,"number of stopping elems:",i6)')nz
      write(11,2050)nz
c JAF error check (18 lines)
      if(nz.gt.maxel)then
        write(*,'("error. too many elements."/,
     1   "nz=",i0," maxel=",i0)')nz,maxel
        close(unit=1)
        close(unit=2)
        close(unit=3)
        close(unit=4)
        close(unit=5)
        close(unit=6, status='delete')
        close(unit=7, status='delete')
        close(unit=8, status='delete')
        close(unit=10, status='delete')
        close(unit=11)
        close(unit=12, status='delete')
        close(unit=13, status='delete')
        if(l_sdata)close(unit=14, status='delete')
        stop
      end if
      if (isg.eq.0) then
        write(11,2060)isg
      elseif (isg.eq.1) then
        write(11,2070)isg
      else
        stop 'Error in isg input!'
      endif
      if (nz.le.0) go to 525

c---------------------------------
c read input material constituents
c---------------------------------
      write(11,2090)
      write(11,2080)
      write(11,2090)
      write(11,2100)
      write(11,2110)
      do 420 j=1,nz
        read (1,*) jzm(j),azm(j)
        write(11,2120) jzm(j),azm(j)
  420 continue
      write(11,2090)

c----------------------------------------
c order the nz material constituents by z
c----------------------------------------
      do 440 j=1,nz
        isw=0
        do 430 jn=2,nz
          if (jzm(jn).gt.jzm(jn-1)) go to 430
          jtem=jzm(jn-1)
          atem=azm(jn-1)
          jzm(jn-1)=jzm(jn)
          azm(jn-1)=azm(jn)
          jzm(jn)=jtem
          azm(jn)=atem
          isw=1
  430   continue
        if (isw.eq.0) go to 450
  440 continue

c----------------------------------------------------
c read tape 2 for stopping cross-section coefficients
c----------------------------------------------------
  450 rewind 2
  460 read (2,10) icont
      if (icont.ne.0) go to 460
      do 520 j=1,nz
  470   read (2,40,end=510) iz,c
  480   if (jzm(j)-iz) 510,490,470
  490   continue
        do 500 jt=1,9
  500   czm(jt,j)=c(jt)
        if (isg.eq.1.and.c(10).gt.0.) then
          do 505 jt=1,5
  505     czm(jt,j)=c(jt+9)
        endif
        go to 520
  510   write (6,50) jzm(j)
        stop 'Element not found in tape2'
  520 continue
c JAF DEBUG (5 lines)
c     if(l_sdata)then
c       write(14,'("czm(3,15)",1pe16.8)')czm(3,15)
c       czm(3,15)=czm(3,15)*(1.d0+0.010d0)
c       write(14,'("czm(3,15)",1pe16.8)')czm(3,15)
c     end if
  525 if (id.eq.1) go to 580

c------------------------------------------------------------
c construct neutron group structure desired from tape 1 input
c neutron groups ordered in decreasing energy
c------------------------------------------------------------
      read (1,*) nng,enmax,enmin
      if (nng.gt.0) then
        write(11,2130)nng
      else
        idummy=-nng
        write(11,2130)idummy
      endif
      write(11,2140)enmax
      write(11,2150)enmin
      write(11,2090)
      write(11,2160)
      write(11,2090)
      write(11,2170)
      write(11,2180)
c JAF error check (18 lines)
      if(abs(nng).gt.maxnng)then
        write(*,'("error. too many neutron energy groups."/,
     1   "nng=",i0," maxnng=",i0)')nng,maxnng
        close(unit=1)
        close(unit=2)
        close(unit=3)
        close(unit=4)
        close(unit=5)
        close(unit=6, status='delete')
        close(unit=7, status='delete')
        close(unit=8, status='delete')
        close(unit=10, status='delete')
        close(unit=11)
        close(unit=12, status='delete')
        close(unit=13, status='delete')
        if(l_sdata)close(unit=14, status='delete')
        stop
      end if
      nngp1=nng+1
      if (nng.gt.0) go to 540
      nng=-nng
      nngp1=nng+1
      en(nngp1)=enmin
      read (1,*) (en(n),n=1,nng)
      do 527 n=1,nng
        write(11,2190) n,en(n),en(n+1)
  527 continue
      if (en(1).ne.enmax) write(11,2193)
      do 528 n=2,nng
        nm1=n-1
        if (en(n).gt.en(n-1)) then
          write(11,2195)n,nm1
          stop 'Energy Structure Incorrectly Entered!'
        endif
  528 continue
      if (en(1).gt.en(2)) go to 565
      nh=nngp1/2
      do 530 n=1,nh
        etem=en(nngp1-n+1)
        en(nngp1-n+1)=en(n)
  530 en(n)=etem
      go to 565
  540 fnng=nng
      den=(enmax-enmin)/fnng
      en(nngp1)=enmin
      en(1)=enmax
      do 550 n=2,nng
        nnm1=n-1
        fnm1=n-1
        en(n)=en(1)-fnm1*den
        write(11,2190)nnm1,en(nnm1),en(n)
  550 continue
      write(11,2190)nnm1,en(nng),en(nngp1)
  565 if (erg.ge.1) then
        write(7,70)
        write(7,80) (en(n),n=nngp1,1,-1)
        write(7,2021)
        write(8,70)
        write(8,80) (en(n),n=nngp1,1,-1)
        write(8,2021)
        write(10,70)
        write(10,80) (en(n),n=nngp1,1,-1)
        write(10,2021)      
      elseif (erg.le.-1) then        
  560   write(7,70)
        write(7,80) (en(n),n=1,nngp1)
        write(7,2021)
        write(8,70)
        write(8,80) (en(n),n=1,nngp1)
        write(8,2021)
        write(10,70)
        write(10,80) (en(n),n=1,nngp1)
        write(10,2021)
      endif
      
c----------------------------
c zero total spectrum storage
c----------------------------
      do 570 n=1,nng
        ts(n)=0.
        gtsan(n)=0.
  570 tssf(n)=0.
      etpall=0.
      ebtall=0.
      etopan=0.
      ebotan=0.
      etopsf=0.
      ebotsf=0.
      etopdn=0.
      ebotdn=0.
      gttqan=0.
      totqsf=0.
      totqdn=0.
      qtotal=0.
      ebarsf=0.
      ebardn=0.
      ebrall=0.
c JAF for derivatives (7 lines)
      san_der(1:maxel,1:maxnq)=0.
      iz_tar(1:maxel)=0
      iz_src(1:maxnq)=0
c     rk_drv_n_tot(1:maxel)=0.
c     fp_tot(1:2,1:maxel)=0.
      depsdb(1:6,1:maxel,1:maxnag)=0.
      depsdc(1:9,1:maxel,1:maxnag)=0.
c
c---------------------------------------------
c read alpha/s.f./d.n. sources from user input
c---------------------------------------------
  580 if (idd.eq.1) then
        read(1,*)nq
      elseif (idd.eq.3) then
        read(1,*)ebeam
        nq=0
        nal=0
      endif
      write(11,2090)
      if (idd.eq.1) then
        write(11,2200)nq
c JAF for derivatives (2 lines)
        write(13,'("number of sources input: ",i6)')nq
        if(l_sdata)write(14,'("number of sources input: ",i6)')nq
      else
        write(11,2210)ebeam
        go to 620
      endif
      if (nq.eq.0) then
        stop 'Error in nq input!'
      endif
      write(11,2090)
      write(11,2220)
      write(11,2090)
      write(11,2230)
      write(11,2240)
      do 590 k=1,nq
        read (1,*) jq(k),aq(k)
        write(11,2250) jq(k),aq(k)
  590 continue

c-----------------------------------------------
c order sources by z-a-state id = s+10*a+10000*z
c-----------------------------------------------
      if (nq.eq.1) go to 620
      do 610 k=1,nq
        isw=0
        do 600 kn=2,nq
          if (jq(kn).gt.jq(kn-1)) go to 600
          jtem=jq(kn-1)
          atem=aq(kn-1)
          jq(kn-1)=jq(kn)
          aq(kn-1)=aq(kn)
          jq(kn)=jtem
          aq(kn)=atem
          isw=1
  600   continue
        if (isw.eq.0) go to 620
  610 continue

c----------------------------------------
c read target information from user input
c----------------------------------------
  620 if (nz.eq.0) go to 1140
      read (1,*) nt,nag
      write(11,2090)
      write(11,2260)nt
      write(11,2270)nag
c JAF error check (18 lines)
      if(nag.gt.maxnag-1)then
        write(*,'("error. too many alpha energy groups."/,
     1   "nag=",i0," maxnag=",i0)')nag,maxnag-1
        close(unit=1)
        close(unit=2)
        close(unit=3)
        close(unit=4)
        close(unit=5)
        close(unit=6, status='delete')
        close(unit=7, status='delete')
        close(unit=8, status='delete')
        close(unit=10, status='delete')
        close(unit=11)
        close(unit=12, status='delete')
        close(unit=13, status='delete')
        if(l_sdata)close(unit=14, status='delete')
        stop
      end if
c JAF for derivatives (6 lines)
      write(13,'("number of targets input: ",i6)')nt
      if(idd.eq.1)then
        write(13,'("in this file, i is targets, j is stopping ",
     1   "elements, k is sources, l is alpha levels")')
      else if(idd.eq.3)then
        write(13,'("in this file, i is targets, j is stopping ",
     1   "elements")')
      end if
      if(l_sdata)then
        write(14,'("number of targets input: ",i6)')nt
        if(idd.eq.1)then
          write(14,'("in this file, i is targets, j is stopping ",
     1     "elements, k is sources, l is alpha levels")')
        else if(idd.eq.3)then
          write(14,'("in this file, i is targets, j is stopping ",
     1     "elements")')
        end if
      end if
      write(11,2090)
      write(11,2280)
      write(11,2090)
      write(11,2290)
      write(11,2300)
c JAF move this block to here from after "if (id.eq.2) rewind 4"
      if (idd.eq.3) then
        ajwd1='/'
        ajwd2='micr'
        ajwd3='oamp'
      else
        ajwd1='/'
        ajwd2='cm**'
        ajwd3='3   '
      endif
      if (nt.eq.0) go to 1140

c----------------------------------------------------
c beginning of big loop over target nuclides
c         for nt=0 calculate spontaneous fission only
c----------------------------------------------------
      rewind 3
      if (id.eq.2) rewind 4
      write (6,100) ajwd1,ajwd2,ajwd3,ajwd1,ajwd2,ajwd3

c--------------------------
c loop on target nuclides i
c--------------------------
  630 read (3,10) icont
      if (icont.ne.0) go to 630
  640 if (id.eq.2) read (4,10) icont
      if (icont.ne.0) go to 640
      do 1120 i=1,nt
      etpant=0.
      ebtant=0.
      totqan=0.
      if(id.eq.1) go to 660
      do 650 n=1,nng
  650 tsan(n)=0.

c------------------------------------------------------
c find alpha-n cross sections for this target on tape 3
c------------------------------------------------------
  660 read (1,*) idt,at
      write(11,2310)idt,at
  665 read (3,105,end=680) jdt,jps
  670 read (3,120) (e(ip),x(ip),ip=1,jps)
      if (idt-jdt) 680,690,665
  680 write (6,130) idt
c JAF no need to stop
c     stop 'Target nuclide not found in tape3'
      write(*,'("Target nuclide not found in tape3")')
      cycle
  690 lzt=idt/10000
c JAF DEBUG
c     if(l_sdata)then
c       write(14,2310)idt,at
c       write(14,'("table (alpha,n) cross sections")')
c       do m=1,jps
c         write(14,'(i6,1p2e16.8)')m,e(m),x(m)
c       end do ! m
c       write(14,'("x(0007)",1pe16.8)')x(0007)
c       x(0007)=x(0007)+10.0d0
c       write(14,'("x(0007)",1pe16.8)')x(0007)
c     end if
c JAF end mod
      lat=idt/10-lzt*1000
      atar=lat
      lst=idt-lzt*10000-lat*10
      amt=' '
      if (lst.ne.0) amt='m'
      apro=atar+3.
      do 700 ip=2,jps
        eamin=e(ip-1)
        if (x(ip).ne.0.) go to 710
  700 continue
  710 if (eamin.lt.0.001) eamin=0.001
c JAF set jl and jp
c     if (id.eq.1) go to 760
      if (id.eq.1)then
        jl=0
        jp=0
        go to 760
      end if

c----------------------------------------------------------------
c find product nuclide level branchings for this target on tape 4
c f(il,ip)   = fraction of (alpha,n) reactions with target i at energy
c              e(ip) resulting in the production of product level il.
c jl         = number of product nuclide levels
c jp         = number of product level branching data points
c----------------------------------------------------------------
  720 read (4,140,end=750) jdt,jl,jp
  730 read (4,120) q,(el(il),il=1,jl)
      do 740 ip=1,jp
  740 read (4,120) ep(ip),(f(il,ip),il=1,jl)
      if (idt-jdt) 750,760,720
  750 write (6,150) idt
c JAF no need to stop
c     stop 'Target nuclide not found in tape4'
      write(*,'("Target nuclide not found in tape4")')
      cycle
  760 rewind 5
      if (idd.eq.3) then 
        sbtqan=0.
        if(id.ne.1) then
          do 765 il=1,jl
  765     totlev(il)=0.
        endif
        nq=1
        go to 771
      endif
  770 read (5,10) icont
      if (icont.ne.0) go to 770
  771 continue
c--------------------------------------------
c loop on source nuclides k for this target i
c contribution to (alpha,n) neutrons
c--------------------------------------------
      do 1100 k=1,nq
      if (idd.eq.3) goto 850
      etop=0.
      ebot=0.
      do 775 il=1,jl
        etopl(il)=0.
        ebotl(il)=0.
        do 774 n=1,nng
  774   sl(n,il)=0.
  775 totlev(il)=0.
      sbtqan=0.
      if(id.eq.1) go to 790
      do 780 n=1,nng
  780 san(n)=0.
      if(nq.le.0) go to 850

c-------------------------------------------
c read source nuclide decay data from tape 5
c-------------------------------------------
  790 read (5,110,end=830) idq,nal,ndn
  800 read (5,60) alam,bfsf,barnu,a,b,bfdn
      if (nal.eq.0) go to 810
      read (5,120) (eal(l),fal(l),l=1,nal)
  810 if (ndn.eq.0) go to 820
      read (5,120) (fdng(nd),nd=1,ndn)
  820 if (idq-jq(k)) 790,840,830
  830 write (6,160) jq(k)
c JAF no need to stop
c     stop 'Source nuclide not found in tape5'
      write(*,'("Source nuclide not found in tape5")')
      cycle
  840 if (i.eq.1.and.bfsf.gt.0.) isfnq=isfnq+1
      if (i.eq.1.and.bfdn.gt.0.) idnnq=idnnq+1
      if (nal.eq.0) go to 1100
      if(eal(nal).gt.e(1)) go to 845
      write(6,842) idq,idt
  842 format(i8,13h alphas below,i8,20h (alpha,n) threshold)
      go to 1100
  845 lzq=idq/10000
      laq=idq/10-lzq*1000
      lsq=idq-lzq*10000-laq*10
      amq=' '
      if (lsq.ne.0) amq='m'

c---------------------------------------------------------
c set alpha energy grid ee(m) for this target/source combo
c---------------------------------------------------------
      eamax=eal(nal)
  850 continue
      if (idd.eq.3) eamax=ebeam
      if (eamax.le.eamin) then
        qan=0.
        mm=0
c JAF this is a bug. nagp1 is not yet defined. add this line:
        nagp1=nag+1
        p(nagp1)=0.
        pval=0.
        go to 945
      endif
      if (eamax.le.6.5) go to 860
      write (6,170) eamax
      stop 'Maximum alpha energy above 6.5 MeV'
  860 fnag=nag
      nagp1=nag+1
      dea=(eamax-eamin)/fnag
c JAF for derivatives
      if(idd.eq.3)then
        idq=0     ! otherwise undefined
        alam=0.d0 ! set in case a derivative tries to use it
      end if
      iz_tar(i)=idt
      iz_src(k)=idq
      write(13,'(/,"target nuclide",i6,i12)')i,idt
      write(13,'("atom frac.",2x,1pe15.7," (target)")')at
      if(idd.eq.1)then
        write(13,'("source nuclide",i6,i12)')k,idq
c       write(13,'("number of alphas",i4)')nal
        write(13,'("atom density",1pe15.7," (source)")')aq(k)*1.d-24
        write(13,'("lambda",6x,1pe15.7," (source)")')alam
      end if
c     write(13,'("eamax,eamin",1x,1p2e15.7)')eamax,eamin
c     write(13,'("delta-e",5x,1pe15.7," [(eamax-eamin)/nag]")')dea
      if(l_sdata)then
        write(14,'(/,"target nuclide",i6,i12)')i,idt
        write(14,'("atom frac.",2x,1pe15.7," (target)")')at
        if(idd.eq.1)then
          write(14,'("source nuclide",i6,i12)')k,idq
          write(14,'("atom density",1pe15.7," (source)")')aq(k)*1.d-24
          write(14,'("lambda",6x,1pe15.7," (source)")')alam
          write(14,'(3x,"l",2x,"erg_alpha",6x,"fal")')
          do l=1,nal
            write(14,'(i4,1p2e15.7)')l,eal(l),fal(l)
          end do ! l
        end if
      end if
c JAF end mod
      ee(1)=eamin
      do 870 m=2,nag
        fm=m-1
  870 ee(m)=ee(1)+dea*fm
      ee(nagp1)=eamax

c----------------------------------------------------
c for each energy ee(m) calc. and store cross section
c----------------------------------------------------
      ip=1
      do 900 m=1,nagp1
  880   if (ee(m).ge.e(ip).and.ee(m).le.e(ip+1)) go to 890
        ip=ip+1
        if (ip.lt.jps) go to 880
        write (6,180) idt,m,ee(m)
        stop 'Error in energy structure!' 
  890   slope=(x(ip+1)-x(ip))/(e(ip+1)-e(ip))
        enrsep=(x(ip)*e(ip+1)-x(ip+1)*e(ip))/(e(ip+1)-e(ip))
  900 cx(m)=enrsep+ee(m)*slope

c-------------------------------------------------------------
c for each energy ee(m) calc. and store stopping cross section
c-------------------------------------------------------------
      tmp1=0.
      tmp2=0.
      do 930 m=1,nagp1
        scx(m)=0.
        do 920 j=1,nz
          zmat=jzm(j)
          amat=amass(jzm(j))
c-----  nuclear stopping---
          term=zalp**0.66667+zmat**0.66667
          bot=zalp*zmat*(alph+amat)*sqrt(term)
          b1=32530.
          b2=1.593
          b3=1.7
          b4=6.8
          b5=3.4
          b6=0.47
c JAF DEBUG (5 lines)
c         if(l_sdata)then
c           if(m*j.eq.1)write(14,'("b6=",1pe16.8)')b6
c           b6=b6*(1.d0-0.03d0)
c           if(m*j.eq.1)write(14,'("b6=",1pe16.8)')b6
c         end if
          rep=b1*amat*ee(m)/bot
          if(rep.lt.0.001) then
            rep12=sqrt(rep)
            epsnu=b2*rep12
            depsdb(2,j,m)=rep12
            depsdb(1,j,m)=depsdb(2,j,m)*0.5d0*b2/b1
            go to 905
          endif
          if(rep.lt.10.) then
            rep12=sqrt(rep)
            rep15=rep**1.5
            bot=1.+b4*rep+b5*rep15
            al=alog(rep+rnapier)
            epsnu=b3*rep12*al/bot
            depsdb(1,j,m)=b3*rep15/(b1*bot)*(0.5d0*al/rep
     1       +1.d0/(rep+rnapier)-(b4+1.5d0*b5*rep12)*al/bot)
            depsdb(3,j,m)=rep12*al/bot
            depsdb(4,j,m)=-b3*rep15*al/bot**2
            depsdb(5,j,m)=depsdb(4,j,m)*rep12
            go to 905
          endif
          al=alog(b6*rep)
          epsnu=al/(2.*rep)
          depsdb(1,j,m)=(1.d0-al)/(2.*b1*rep)
          depsdb(6,j,m)=0.5d0/(b6*rep)
c-----  electronic stopping---
  905     if(ee(m).gt.30.) go to 906
          slow=czm(1,j)*(1000.*ee(m))**czm(2,j)
          shigh=(czm(3,j)/ee(m))*alog(1.+czm(4,j)/ee(m)+czm(5,j)*ee(m))
          epsel=slow*shigh/(slow+shigh)
          depsdc(1,j,m)=epsel**2/(slow*czm(1,j))
          depsdc(2,j,m)=epsel**2/slow*alog(1000.*ee(m))
          depsdc(3,j,m)=epsel**2/(shigh*czm(3,j))
          depsdc(5,j,m)=(epsel/shigh)**2*czm(3,j)
     1     /(1.+czm(4,j)/ee(m)+czm(5,j)*ee(m))
          depsdc(4,j,m)=depsdc(5,j,m)/ee(m)**2
          go to 907
  906     eil=alog(1/ee(m))
          arg=czm(6,j)
          arg=arg+eil*(czm(7,j)+czm(8,j)*eil+czm(9,j)*eil*eil)
          epsel=exp(arg)
          depsdc(6,j,m)=epsel
          depsdc(7,j,m)=eil*depsdc(6,j,m)
          depsdc(8,j,m)=eil*depsdc(7,j,m)
          depsdc(9,j,m)=eil*depsdc(8,j,m)
  907     continue
          dcx=epsnu+epsel
c JAF for derivatives
c         if(m.eq.1.and.j.eq.1)then
c           write(13,'(/,3x,"m indexes alpha groups, ",
c    1       "j indexes stopping elements",/,"   m   j jzm",
c    2       2x,"at.frac.",7x,"eps_j")')
c         end if
c         write(13,'(3i4,1p20e15.7)')m,j,jzm(j),azm(j),dcx
c JAF DEBUG (5 lines)
c         if(l_sdata.and.j.eq. 1.and.m.eq.101)then
c           write(14,'("dcx",2i6,1p2e16.8)')j,m,dcx,rep
c           dcx=dcx*(1.d0-0.010d0)
c           write(14,'("dcx",2i6,1p2e16.8)')j,m,dcx,rep
c         end if
          if(l_sdata)then
            if(m*j.eq.1)write(14,'("stopping power",/,
     1       4x,"j",4x,"g",2x,"ee",14x,"nuclear",9x,"nuc.frac.",
     2       7x,"electronic",6x,"elec.frac.",6x,"total")')
            write(14,'(2i5,1p6e16.8)')j,m,ee(m),epsnu,epsnu/dcx,
     1       epsel,epsel/dcx,dcx
          end if
          dcxe(j,m)=dcx ! stopping power, element j, energy m
          if(j.eq.1)then
            tmp1=max(tmp1,epsnu/dcx)
          else if(j.eq.2)then
            tmp2=max(tmp2,epsnu/dcx)
          end if
c JAF end mod
  920   scx(m)=scx(m)+azm(j)*dcx ! Eq. (10) -- JAF
  930 continue
c JAF DEBUG (1 line)
c     if(l_sdata)write(14,'(1p2e16.8)')tmp1,tmp2

c-----------------------------------------------------------
c for each energy ee(m) calc. ratio r(m), then integral p(m)
c-----------------------------------------------------------
      r(1)=cx(1)/scx(1)
      p(1)=0.
      fact=1.0d-06*at
c JAF DEBUG (5 lines)
c     if(l_sdata)then
c       write(14,'("cx(77-79)",1p3e16.8)')cx(77),cx(78),cx(79)
c       cx(79)=cx(79)*(1.d0+0.500d0)
c       write(14,'("cx(77-79)",1p3e16.8)')cx(77),cx(78),cx(79)
c     end if
      do 940 m=2,nagp1
        r(m)=cx(m)/scx(m) ! sigma/epsilon, Eq. (14) -- JAF
  940 p(m)=p(m-1)+fact*(r(m-1)+r(m))*dea/2. ! Eq. (14) -- JAF
c JAF for derivatives
c     write(13,'(/,3x,"fact=",1pe15.7)')fact
c     write(13,'(/,3x,"eamax=",1pe15.7)')eamax
c     if(l_sdata)then
c       write(14,'(/,3x,"i indexes target nuclides, ",
c    1   "k indexes source nuclides, m indexes alpha groups",
c    2   /,"   i   k   m  energy",
c    3   9x,"cross-sec.",5x,"stop.-power",4x,"p_i(energy)")')
c       do m=1,nagp1
c         write(14,'(3i4,1p20e15.7)')i,k,m,ee(m),cx(m),scx(m),p(m)
c       end do ! m
c     end if
c
c derivatives of source rate density w.r.t. (alpha,n) cross sections
c and stopping power data. first deal with the interpolation of pval
c that happens below (label 970).
      dqdd(1:nagp1)=0.d0
      do m=1,nagp1
        h=1.d0
        if(m.eq.1.or.m.eq.nagp1)then
          h=0.5d0
        end if
        do l=1,nal
          if(m.gt.1)then
c m is immediately above l
            if(eal(l).gt.ee(m-1).and.eal(l).lt.ee(m))then
              dqdd(m)=dqdd(m)+0.5d0*fal(l)*(eal(l)-ee(m-1))/dea
            end if
          end if             ! this "end if/if" looks like it should be
          if(m.lt.nagp1)then ! be an "else if" but it should not
c m is immediately below l
            if(eal(l).gt.ee(m).and.eal(l).lt.ee(m+1))then
              dqdd(m)=dqdd(m)+0.5d0*fal(l)*(1.d0+(eal(l)-ee(m))/dea)
c m is anywhere else below l
            else if(eal(l).ge.ee(m+1))then
              dqdd(m)=dqdd(m)+h*fal(l)
            end if
          else if(m.eq.nagp1.and.eal(l).eq.ee(m))then
            dqdd(m)=dqdd(m)+h*fal(l)
          end if
        end do ! l
      end do ! m
c derivative of the source rate density w.r.t. (alpha,n) cross sections
c in the data table. dqdx1 is the derivative of the source rate density
c w.r.t. the (alpha,n) cross sections on the alpha energy grid.
      if(l_sdata)then
        write(14,'(2x,"derivative of source rate density ",
     1   "w.r.t. (alpha,n) cross section comp. grid values")')
        write(14,'(4x,"g  ee",14x,"cx",14x,"dqdsig")')
        do m=1,nagp1
          dqdx1(m)=dqdd(m)*aq(k)*alam*fact/scx(m)*dea
          write(14,'(i5,1p3e16.8)')m,ee(m),cx(m),dqdx1(m)
        end do ! m
        write(14,'(2x,"derivative of source rate density ",
     1   "w.r.t. (alpha,n) cross section table values")')
        write(14,'(3x,"jps=",i6)')jps
        write(14,'(4x,"g  energy",10x,"xsec",12x,"dqdx")')
c identify all points on alpha energy grid affected by each point on
c data table. here dq2 is the desired derivative.
        m2=1
        do ip=1,jps
          dq2=0.d0
          m1=m2
c         if(ip.eq.1)then
c           write(14,'(i5,1pe14.6,i5,e14.6)')ip,e(ip),ip+1,e(ip+1)
c         else if(ip.eq.jps)then
c           write(14,'(i5,1pe14.6,i5,e14.6)')ip-1,e(ip-1),ip,e(ip)
c         else
c           write(14,'(i5,1pe14.6,i5,e14.6,i5,e14.6)')ip-1,e(ip-1),
c    1       ip,e(ip),ip+1,e(ip+1)
c         end if
          do m=m1,nagp1
            if(ip.eq.1)then
              if(ee(m).lt.e(ip+1))then
c               write(14,'("ee(m)",2i5,1pe14.6," included")')ip,m,ee(m)
                if(m2.eq.m1.and.ee(m).gt.e(ip))m2=m
              else
c               write(14,'("ee(m)",2i5,1pe14.6," exit")')ip,m,ee(m)
                exit
              end if
            else if(ip.eq.jps)then
              if(ee(m).gt.e(ip-1))then
c               write(14,'("ee(m)",2i5,1pe14.6," included")')ip,m,ee(m)
                continue
              else
c               write(14,'("ee(m)",2i5,1pe14.6," exit")')ip,m,ee(m)
                exit
              end if
            else
              if(ee(m).gt.e(ip-1).and.ee(m).lt.e(ip+1))then
c               write(14,'("ee(m)",2i5,1pe14.6," included")')ip,m,ee(m)
                if(m2.eq.m1.and.ee(m).gt.e(ip))m2=m
              else if(ee(m).gt.e(ip+1))then
c               write(14,'("ee(m)",2i5,1pe14.6," exit")')ip,m,ee(m)
                exit
              else
c               write(14,'("ee(m)",2i5,1pe14.6," cycle")')ip,m,ee(m)
                cycle
              end if
            end if
            if(ee(m).lt.e(ip))then
              dq2=dq2+dqdx1(m)*(ee(m)-e(ip-1))/(e(ip)-e(ip-1))
            else if(ee(m).gt.e(ip))then
              dq2=dq2+dqdx1(m)*(e(ip+1)-ee(m))/(e(ip+1)-e(ip))
            else if(ee(m).eq.e(ip))then
              dq2=dq2+dqdx1(m)
            end if
          end do ! m
          write(14,'(i5,1p3e16.8)')ip,e(ip),x(ip),dq2
        end do ! ip
c
c derivative of source rate density w.r.t. stopping power data.
c fact is 1E-6 times atom fraction of target nuclide.
c aq(k) is the atom density of source nuclide.
c here dq2 is the derivative of the source rate density w.r.t.
c the stopping power of element j at alpha energy m.
        dqdb(1:6)=0.
        dqdc(1:9,1:nz)=0.
        write(14,'(4x,"j",4x,"g",2x,"ee",14x,"eps",13x,"dqde")')
        do j=1,nz
          do m=1,nagp1
            dq2=-dqdd(m)*aq(k)*azm(j)*alam*fact*dea*r(m)/scx(m)
            write(14,'(2i5,1p20e16.8)')j,m,ee(m),dcxe(j,m),dq2
            dqdb(1:6)=dqdb(1:6)+dq2*depsdb(1:6,j,m)
            dqdc(1:9,j)=dqdc(1:9,j)+dq2*depsdc(1:9,j,m)
          end do ! m
        end do ! j
        write(14,'(2x,"derivative of (alpha,n) source rate density ",
     1   "w.r.t. nuclear stopping power constants")')
        write(14,'(4x,"b",2x,"dqdb")')
        do m=1,6
          write(14,'(i5,1pe16.8)')m,dqdb(m)
        end do ! j
        write(14,'(2x,"derivative of (alpha,n) source rate density ",
     1   "w.r.t. electronic stopping power data")')
        write(14,'(4x,"c",4x,"j",2x,"dqdc")')
        do m=1,9
          do j=1,nz
            write(14,'(2i5,1pe16.8)')m,j,dqdc(m,j)
          end do ! j
        end do ! m
c       write(14,'(3x,"c",4x,"j",4x,"g",3x,"depsdc")')
c       do m=1,9
c         do j=1,nz
c           do l=1,nagp1
c             write(14,'(3i5,1pe16.8)')m,j,l,depsdc(m,j,l)
c           end do ! l
c         end do ! j
c       end do ! m
      end if ! l_sdata
c
      do j=1,nz
        r_drv_n(1)=r(1)*dcxe(j,1)/scx(1)
        p_drv_n(j,1)=0.
c fact is 1E-6*at = 1E-6*Ni/N where i is target
        do m=2,nagp1
          r_drv_n(m)=r(m)*dcxe(j,m)/scx(m)
          p_drv_n(j,m)=p_drv_n(j,m-1)+fact*(r_drv_n(m-1)+r_drv_n(m))
     1     *dea/2.
        end do ! m
c second derivatives
        if(l_2nd)then
          do j2=1,nz
            r_drv2_n(1)=r_drv_n(1)*dcxe(j2,1)/scx(1)
            p_drv2_n(j,j2,1)=0.
            do m=2,nagp1
              r_drv2_n(m)=r_drv_n(m)*dcxe(j2,m)/scx(m)
              p_drv2_n(j,j2,m)=p_drv2_n(j,j2,m-1)+fact
     1         *(r_drv2_n(m-1)+r_drv2_n(m))*dea/2.
            end do ! m
          end do ! j2
        end if
      end do ! j
c JAF end mod
  945 if (idd.ne.3) go to 950

c---------------------------------------------------------------------
c if nq=0 calculate neutrons from alpha beam: 1 microamp=3.1209+12 a/s
c---------------------------------------------------------------------
      aps=3.1209d+12
      beamn=aps*p(nagp1)
      pval=p(nagp1)
      write (6,190) jsm(lzt),lat,amt,at,ebeam,aps,pval,beamn
      nal=1
      eal(1)=ebeam
      qan=beamn
      sbtqan=sbtqan+qan ! Eq. (45), summing over targets -- JAF
      totqan=totqan+qan
c JAF for derivatives
c the derivative of p(nagp1) w.r.t Nj is p_drv_n(j,nagp1)
      write(13,'(/,3x,"(alpha,n) source rate density for this ",
     1 "target",1pe15.7)')beamn
      if(l_sdata)then
        write(14,'(/,3x,"(alpha,n) source rate density for this ",
     1   "target",1pe15.7)')beamn
      end if
      write(13,'(/,3x,"j   Z  at.frac.",7x,"dp_i/dN_j*microamp")')
      do j=1,nz
        write(13,'(2i4,1p2e15.7)')j,jzm(j),azm(j),
     1   aps*p_drv_n(j,nagp1)*1.d24
      end do ! j
      if(l_2nd)then
        write(13,'(/,2x,"j1  Z1  at.frac.",7x,"j2  Z2  at.frac.",7x,
     1   "d^2p_i/dN_j1 dN_j2*microamp")')
        do j=1,nz
          do j2=1,nz
            write(13,'(2i4,1pe15.7,2i4,2e15.7)')j,jzm(j),azm(j),
     1       j2,jzm(j2),azm(j2),p_drv2_n(j,j2,nagp1)*1.d24
          end do ! j2
        end do ! j
      end if
c JAF end mods
      if (id.eq.1) go to 1101
      mm=nag

c----------------------------------------------
c if nq>0 calculate neutrons from source alphas
c----------------------------------------------
  950 do 1080 l=1,nal
c JAF for derivatives (5 lines)
      if(l.eq.1)then
        fp=0.
        rk_drv_n(1:maxel)=0.
        if(l_2nd)rk_drv2_n(1:maxel,1:maxel)=0.
      end if
      if (idd.eq.3) go to 1000
      do 960 m=1,nag
        mm=m
  960 if (eal(l).ge.ee(m).and.eal(l).le.ee(m+1)) go to 970
      mm=0
      pval=0.
      go to 975
  970 pval=p(mm)+(eal(l)-ee(mm))*(p(mm+1)-p(mm))/(ee(mm+1)-ee(mm))
c JAF for derivatives
      do j=1,nz
        pval_drv_n(j)=p_drv_n(j,mm)+(eal(l)-ee(mm))*
     1   (p_drv_n(j,mm+1)-p_drv_n(j,mm))/(ee(mm+1)-ee(mm))
        if(l_2nd)then
          do j2=1,nz
            pval_drv2_n(j,j2)=p_drv2_n(j,j2,mm)+(eal(l)-ee(mm))*
     1       (p_drv2_n(j,j2,mm+1)-p_drv2_n(j,j2,mm))/(ee(mm+1)-ee(mm))
          end do ! j2
        end if
      end do ! j
      fp=fp+fal(l)*pval
c     fp_tot(1,i)=fp_tot(1,i)+fal(l)*pval
c     fp_tot(2,i)=fp_tot(2,i)+fal(l)*pval*alam*aq(k)
      do j=1,nz
        rk_drv_n(j)=rk_drv_n(j)+fal(l)*pval_drv_n(j)
c       rk_drv_n_tot(j)=rk_drv_n_tot(j)+fal(l)*pval_drv_n(j)*alam*aq(k)
        if(l_2nd)then
          do j2=1,nz
            rk_drv2_n(j,j2)=rk_drv2_n(j,j2)+fal(l)*pval_drv2_n(j,j2)
          end do ! j2
        end if
      end do ! j
c k indexes source nuclide, l indexes alpha level
c JAF end mod
  975 aps=aq(k)*alam*fal(l)
c JAF for derivatives
c     if(l.eq.1)then
c       write(13,'(/,3x,"target nuclide",i6,i12)')i,idt
c       write(13,'(3x,"source nuclide",i6,i12)')k,idq
c       write(13,'(3x,"l indexes alpha lines, ",
c    1   "i indexes target nuclides,",
c    2   " j (across) indexes stopping elements",/,
c    3   "   divide dp_i/dN_j by material atom density",/,
c    4   "   l  erg_alpha",6x,"fal",12x,"p_i(erg_alpha)",
c    5   1x,"dp_i/dN_j(erg_alpha)")')
c     end if
c     write(13,'(i4,1p30e15.7)')l,eal(l),fal(l),pval,
c    1 (pval_drv_n(j),j=1,nz)
      if(l.eq.nal)then
c after this, fp is the (a,n) source rate density for this
c source and target, same as sbtqan below.
        fp=fp*aq(k)*alam
        write(13,'(/,3x,"(alpha,n) source rate density for this ",
     1   "source and target",1pe15.7)')fp
        if(l_sdata)then
          write(14,'(/,3x,"(alpha,n) source rate density for this ",
     1     "source and target",1pe15.7)')fp
        end if
        write(13,'(/,3x,"j   Z  at.frac.",7x,"sum_l {fal*dp_i/dN_j}")')
        do j=1,nz
          write(13,'(2i4,1p2e15.7)')j,jzm(j),azm(j),rk_drv_n(j)*1.d24
        end do ! j
        if(l_2nd)then
          write(13,'(/,2x,"j1  Z1  at.frac.",7x,"j2  Z2  at.frac.",7x,
     1     "sum_l {fal*d^2p_i/dN_j1 dN_j2}")')
          do j=1,nz
            do j2=1,nz
              write(13,'(2i4,1pe15.7,2i4,2e15.7)')j,jzm(j),azm(j),
     1         j2,jzm(j2),azm(j2),rk_drv2_n(j,j2)*1.d24
            end do ! j2
          end do ! j
        end if
      end if
c JAF end mod
      qan=aps*pval
      sbtqan=sbtqan+qan
      totqan=totqan+qan
      if (l.eq.1) go to 980
      write (6,200) eal(l),aps,pval,qan
      if(nal.ne.1.and.l.eq.nal) write(6,250) sbtqan
      go to 990
  980 write(6,210)jsm(lzt),lat,amt,at,jsm(lzq),laq,amq,eal(l),aps,
     1 pval,qan
      if(nal.eq.1) write(6,255)
  990 if (id.eq.1.or.mm.eq.0) go to 1080

c-----------------------------------------------------------------
c calculate (alpha,n) neutron spectrum contribution in multigroups
c-----------------------------------------------------------------
c JAF mm is the number of bins, not endpoints
 1000 mmm1=mm-1
      do 1010 m=1,mmm1
c JAF Eq. (32), except does not depend on the product level.
 1010 rr(m)=(p(m+1)-p(m))/pval
      rr(mm)=(pval-p(mm))/pval
      do 1020 n=1,nng
 1020 s(n)=0.
c JAF
      dq3c(1:nng,1:nagp1)=0.d0
      dq3s(1:nng,1:nagp1,1:nz)=0.d0
c fact is 1E-6*at = 1E-6*Ni/N where i is target
c     write(*,'("i,k,l,mm",3i4,i8,1pe17.9)')i,k,l,mm,pval ! DEBUG
c     if(l_sdata)then
c       write(14,'("i,k,l,mm",3i4,i8,1pe17.9)')i,k,l,mm,pval
c       do m=1,mmm1
c         tjaf=0.5d0*dea*fact*(cx(m+1)/scx(m+1)+cx(m)/scx(m))
c         write(14,'(i6,1p2e19.11)')m,rr(m),tjaf/pval
c       end do ! m
c       m=mm
c       if(l.lt.nal)then
c         tjaf=0.5d0*dea*fact*(cx(m+1)/scx(m+1)+cx(m)/scx(m))
c    1     *(eal(l)-ee(mm))/(ee(mm+1)-ee(mm))
c         write(14,'(i6,1p2e19.11,a)')m,rr(m),tjaf/pval," l.lt.nal"
c       else
c         tjaf=0.5d0*dea*fact*(cx(m+1)/scx(m+1)+cx(m)/scx(m))
c         write(14,'(i6,1p2e19.11,a)')m,rr(m),tjaf/pval," l.eq.nal"
c       end if
c     end if
c JAF
      do 1060 il=1,jl
        qlev=q-el(il)
        e90=-qlev*apro/(apro-alph)
        thre=-qlev*(aneut+apro)/(aneut+apro-alph)
        if (qlev.gt.0.) thre=0.
        do 1055 m=1,mm
          ea=(ee(m)+ee(m+1))/2.
          if (m.eq.mm) ea=(eal(l)+ee(mm))/2.
          if (ea.lt.thre) go to 1055
          do 1030 ip=2,jp
            lp=ip
 1030     if (ep(ip-1).le.ea.and.ep(ip).ge.ea) go to 1040
          write (6,220)
          stop 'Error in alpha energy structure'
 1040     lp1=lp-1
c JAF f is read from tape4 at line 740
          bx=f(il,lp1)
c JAF Eq. (33)
          bx=bx+(ea-ep(lp1))*(f(il,lp)-f(il,lp1))/(ep(lp)-ep(lp1))
          term1=sqrt(4.*alph*aneut*ea)/(2.*(aneut+apro))
          term2=alph*aneut*ea/((aneut+apro)*(aneut+apro))
          term3=(apro*ea+apro*qlev-alph*ea)/(aneut+apro)
          senmax=term1+sqrt(term2+term3)
          enmax=senmax*senmax
          if (ea.le.e90) senmin=term1-sqrt(term2+term3)
          if (ea.gt.e90) senmin=-term1+sqrt(term2+term3)
          enmin=senmin*senmin
          ebar=(enmin+enmax)/2.
          val=qan*rr(m)*bx
          etop=etop+val*ebar
          etopl(il)=etopl(il)+val*ebar
          ebotl(il)=ebotl(il)+val
          ebot=ebot+val
          dele=enmax-enmin
          do 1050 n=1,nng
            if (en(n).lt.enmin) go to 1050
            if (en(n+1).gt.enmax) go to 1050
            de=en(n)-en(n+1)
            if (en(n+1).lt.enmin) de=de-(enmin-en(n+1))
            if (en(n).gt.enmax) de=de-(en(n)-enmax)
c JAF Eq. (35). 
c qan=R_k
c rr=H_i,k I think
c bx=S_i,k I think
c           write(6,'("i,k,l,il,lp,m,n",7i4,1p3e14.6)')i,k,l,il,lp,m,n,
c    1       ep(ip),ea,ep(ip-1)
            gpadd=qan*rr(m)*bx*de/dele
            s(n)=s(n)+gpadd
c JAF for derivatives (2 lines)
c           write(13,'("il,m,n,qan,rr(m),bx,de,dele,s(n),gpadd",3i4,
c     1      1p20e15.7)')il,m,n,qan,rr(m),bx,de,dele,s(n),gpadd
            sl(n,il)=sl(n,il)+gpadd
            totlev(il)=totlev(il)+gpadd
c JAF cycle through alpha groups to compute dH/dp, where p
c is (alpha,n) cross section and stopping power.
c fact is 1E-6*at = 1E-6*Ni/N where i is target
c td1c is the derivative of the numerator of H wrt cross section
c td2c is the derivative of the denominator of H, pval wrt cross sec
c td1s is the derivative of the numerator of H wrt stopping power
c td2s is the derivative of the denominator of H, pval wrt stop. pow.
c     rr(m)=(p(m+1)-p(m))/pval
c     rr(mm)=(pval-p(mm))/pval
c correct this if needed.
c           dqdf(1)=dqdf(1)-(ea-ep(lp1))/(ep(lp)-ep(lp1))*rr(m)
c    1       *de/dele*sbtqan
c           dqdf(2)=dqdf(2)+(ea-ep(lp1))/(ep(lp)-ep(lp1))*rr(m)
c    1       *de/dele*sbtqan
            do j=1,nz
              do mp=1,mm+1 ! index boundaries, not bins
                if(mp.eq.m.or.mp.eq.m+1)then
                  if(j.eq.1)td1c=0.5d0*fact/scx(mp)*dea
                  td1s=-0.5d0*azm(j)*fact*dea*r(mp)/scx(mp)
                else
                  if(j.eq.1)td1c=0.d0
                  td1s=0.d0
                end if
                if(j.eq.1)td2c=fact/scx(mp)*dea
                td2s=-azm(j)*fact*dea*r(mp)/scx(mp)
                if(mp.eq.1)then
                  if(j.eq.1)td2c=td2c*0.5d0
                  td2s=td2s*0.5d0
                else if(l.lt.nal.and.mp.eq.mm)then
                  if(j.eq.1)td2c=td2c*0.5d0
                  td2s=td2s*0.5d0
                else if(l.lt.nal.and.mp.eq.mm+1)then
                  if(j.eq.1)td2c=0.d0
                  td2s=0.d0
                else if(l.eq.nal.and.mp.eq.mm+1)then
                  if(j.eq.1)td2c=td2c*0.5d0
                  td2s=td2s*0.5d0
                end if
                if((mp.eq.mm.or.mp.eq.mm+1).and.l.lt.nal)then
                  if(j.eq.1)td2c=td2c+0.5d0*fact/scx(mp)*dea
     1             *(eal(l)-ee(mm))/(ee(mm+1)-ee(mm))
                  td2s=td2s-0.5d0*azm(j)*fact*dea*r(mp)/scx(mp)
     1             *(eal(l)-ee(mm))/(ee(mm+1)-ee(mm))
                  if(m.eq.mm)then
                    if(j.eq.1)td1c=td1c*(eal(l)-ee(mm))
     1               /(ee(mm+1)-ee(mm))
                    td1s=td1s*(eal(l)-ee(mm))/(ee(mm+1)-ee(mm))
                  end if
                end if
                if(j.eq.1)dq3c(n,mp)=dq3c(n,mp)+bx/pval
     1           *(td1c-td2c*rr(m))*de/dele
                dq3s(n,mp,j)=dq3s(n,mp,j)+bx/pval*(td1s-td2s*rr(m))
     1           *de/dele
              end do ! mp
            end do ! j
c JAF
 1050     continue ! JAF n, neutron energy groups
 1055   continue ! JAF m, alpha energy groups
 1060 continue ! JAF il, product nuclide levels
c
      do 1070 n=1,nng
        san(n)=san(n)+s(n)
        gtsan(n)=gtsan(n)+s(n)
        ts(n)=ts(n)+s(n)
 1070 tsan(n)=tsan(n)+s(n)

c1075 continue JAF, unused
 1080 continue ! JAF l, alpha energy level
c JAF
      if(l_sdata)then
        write(14,'(2x,"derivative of source rate density in g w.r.t. ",
     1   "(alpha,n) cross section comp. grid values")')
        write(14,'(4x,"g    m",2x,"dqgdsig")')
        do n=1,nng
          do m=1,nagp1
            dq3c(n,m)=dqdx1(m)*san(n)/sbtqan+sbtqan*dq3c(n,m)
            write(14,'(2i5,1pe19.11)')n,m,dq3c(n,m)
          end do ! m
        end do ! n
c identify all points on alpha energy grid affected by each point on
c data table. here dq2 is the desired derivative.
        write(14,'(2x,"derivative of source rate density in g w.r.t. ",
     1   "(alpha,n) cross section table values")')
        write(14,'(3x,"jps=",i6," nng=",i6)')jps,nng
        write(14,'(4x,"g    m  energy",10x,"xsec",12x,"dqgdx")')
c assume standard output is desired; put in a switch for this.
        do n=1,nng
c make a subroutine for this
          m2=1
          do ip=1,jps
            dq2=0.d0
            m1=m2
            do m=m1,nagp1
              if(ip.eq.1)then
                if(ee(m).lt.e(ip+1))then
                  if(m2.eq.m1.and.ee(m).gt.e(ip))m2=m
                else
                  exit
                end if
              else if(ip.eq.jps)then
                if(ee(m).gt.e(ip-1))then
                  continue
                else
                  exit
                end if
              else
                if(ee(m).gt.e(ip-1).and.ee(m).lt.e(ip+1))then
                  if(m2.eq.m1.and.ee(m).gt.e(ip))m2=m
                else if(ee(m).gt.e(ip+1))then
                  exit
                else
                  cycle
                end if
              end if
              if(ee(m).lt.e(ip))then
                dq2=dq2+dq3c(n,m)*(ee(m)-e(ip-1))/(e(ip)-e(ip-1))
              else if(ee(m).gt.e(ip))then
                dq2=dq2+dq3c(n,m)*(e(ip+1)-ee(m))/(e(ip+1)-e(ip))
              else if(ee(m).eq.e(ip))then
                dq2=dq2+dq3c(n,m)
              end if
            end do ! m
            write(14,'(2i5,1p3e16.8)')n,ip,e(ip),x(ip),dq2
          end do ! ip
        end do ! n
c
        write(14,'(2x,"derivative of source rate density in g w.r.t. ",
     1   "nuclear stopping power constants")')
        write(14,'(4x,"g    b",2x,"dqgdb")')
        do n=1,nng
          do nb=1,6
            dq2=dqdb(nb)*san(n)/sbtqan
            do j=1,nz
              do m=1,nagp1
                dq2=dq2+sbtqan*dq3s(n,m,j)*depsdb(nb,j,m)
              end do ! m
            end do ! j
           write(14,'(2i5,1pe19.11)')n,nb,dq2
          end do ! nb
        end do ! n
c
        write(14,'(2x,"derivative of source rate density in g w.r.t. ",
     1   "electronic stopping power data")')
        write(14,'(4x,"g    c    j",2x,"dqgdc")')
        do n=1,nng
          do nb=1,9
            do j=1,nz
              dq2=dqdc(nb,j)*san(n)/sbtqan
              do m=1,nagp1
                dq2=dq2+sbtqan*dq3s(n,m,j)*depsdc(nb,j,m)
              end do ! m
              write(14,'(3i5,1pe19.11)')n,nb,j,dq2
            end do ! j
          end do ! m
        end do ! n
      end if ! l_sdata
c JAF
c
c----------------------------------------------------
c output this (alpha,n) neutron spectrum contribution
c for this source/target combination
c----------------------------------------------------
      ebarqt=0.
      if(ebot.le.0.) go to 1085
      ebarqt=etop/ebot
      etpant=etpant+etop
      ebtant=ebtant+ebot
 1085 continue
      smga=0.
      do 1090 n=1,nng
 1090 smga=smga+san(n)
c     fracgp=smga/sbtqan
c     do 1095 il=1,jl
c1095 frclev(il)=totlev(il)/sbtqan
c     write(6,205) (frclev(il),il=1,jl)
      if (idd.eq.3) write(7,230)ebeam,jsm(lzt),lat,amt
      if (idd.ne.3) write(7,240)aq(k),jsm(lzq),laq,amq,jsm(lzt),lat,amt
      if (erg.ge.1) then
        write(7,80) (san(n),n=nng,1,-1)
      elseif (erg.le.-1) then
        write(7,80) (san(n),n=1,nng)
      endif
      write(7,81) smga,ajwd1,ajwd2,ajwd3,ebarqt
      smga=0.
      if (idd.eq.3) write(8,230)ebeam,jsm(lzt),lat,amt
      if (idd.ne.3) write(8,240)aq(k),jsm(lzq),laq,amq,jsm(lzt),lat,amt
      do 1097 n=1,nng
        san(n)=san(n)/sbtqan
 1097 smga=smga+san(n)
      if (erg.ge.1) then
        write(8,80) (san(n),n=nng,1,-1)
      elseif (erg.le.-1) then
        write(8,80) (san(n),n=1,nng)
      endif
      write(8,81) smga,ajwd1,ajwd2,ajwd3,ebarqt
      do 1098 il=1,jl
        if(ebotl(il).le.0.) go to 1098
        ebarl=etopl(il)/ebotl(il)
        if (idd.eq.3) then
          write(10,262) ebeam,jsm(lzt),lat,amt,el(il)
        else
          write(10,261) jsm(lzq),laq,amq,jsm(lzt),lat,amt,el(il)
        endif
        write(10,80) (sl(n,il),n=1,nng)
        write(10,81) totlev(il),ajwd1,ajwd2,ajwd3,ebarl
 1098 continue
      if(idd.eq.3) go to 1101
 1100 continue ! JAF k, source nuclide
c
 1101 gttqan=gttqan+totqan

c--------------------------------------
c output total (alpha,n) source for all
c alpha sources on this target nuclide
c--------------------------------------
      if (nq.gt.1) write (6,251) totqan
      if (id.eq.1) go to 1120
      tmga=0.
      do 1110 n=1,nng
 1110 tmga=tmga+tsan(n)
c     fracgp=tmga/totqan
      ebart=etpant/ebtant
      etopan=etopan+etpant
      ebotan=ebotan+ebtant
      if(nq.le.1) go to 1120
      if (erg.ge.1) then
        write(7,260)
        write(7,80) (tsan(n),n=nng,1,-1)
      elseif (erg.le.-1) then      
        write(7,260)      
        write(7,80) (tsan(n),n=1,nng)
      endif
      write(7,81) tmga,ajwd1,ajwd2,ajwd3,ebart
      write(7,2090)
      write(7,2021)
      write(7,2090)
      tmga=0.
      do 1115 n=1,nng
        tsan(n)=tsan(n)/totqan
 1115 tmga=tmga+tsan(n)
      write(8,260)
      if (erg.ge.1) then
        write(8,80) (tsan(n),n=nng,1,-1)
      elseif (erg.le.-1) then
        write(8,80) (tsan(n),n=1,nng)
      endif
      write(8,81) tmga,ajwd1,ajwd2,ajwd3,ebart
      write(8,2090)
      write(8,2021)
      write(8,2090)
      write(10,2090)
      write(10,2021)
      write(10,2090)
 1120 continue ! JAF i, target nuclide
c
c JAF for derivatives
c     write(13,'(/,"totals")')
c     write(13,'(/,3x,"i  izm     sum_k sum_l {fal*p_i}",
c    1 2x,"sum_k sum_l {alam_k*N_k*fal*p_i}")')
c     do i=1,nt
c       write(13,'(i4,i8,1pe15.7,8x,e15.7)')i,iz_tar(i),fp_tot(1,i),
c    1   fp_tot(2,i)
c     end do ! i
c     write(13,'(/,3x,"j jzm  azm",12x,
c    1 "sum_i sum_k sum_l {alam_k*N_k*fal*dp_i/dN_j}")')
c     do j=1,nz
c       write(13,'(2i4,1p2e15.7)')j,jzm(j),azm(j),rk_drv_n_tot(j)
c     end do ! j
      niso(1:2)=0
      do i=1,nt
        if(iz_tar(i).eq.0)cycle
        niso(1)=niso(1)+1
      end do ! i
      do k=1,nq
        if(iz_src(k).eq.0)cycle
        niso(2)=niso(2)+1
      end do ! i
      write(13,'(/,"(alpha,n) targets and sources")')
      write(13,'("number of targets, sources",2i6)')niso(1:2)
      do i=1,nt
        if(iz_tar(i).eq.0)cycle
        write(13,'("target nuclide",i6,i12)')i,iz_tar(i)
      end do ! i
      do k=1,nq  
        if(iz_src(k).eq.0)cycle
        write(13,'("source nuclide",i6,i12)')k,iz_src(k)
      end do ! k
c JAF end mod

c----------------------------------------------------------------------
c output grand total (alpha,n) neutron source: all targets, all sources
c----------------------------------------------------------------------
      if (nt.gt.1) write (6,252) gttqan
      if(id.eq.1) go to 1140
      gtmga=0.
      do 1130 n=1,nng
 1130 gtmga=gtmga+gtsan(n)
      ebaran=etopan/ebotan
c     fracgp=gtmga/gttqan
      write(7,2021)
      write(7,270)
      if (erg.ge.1) then
        write(7,80) (gtsan(n),n=nng,1,-1)
      elseif (erg.le.-1) then
        write(7,80) (gtsan(n),n=1,nng)
      endif      
      write(7,81) gtmga,ajwd1,ajwd2,ajwd3,ebaran
      write(7,2090)
      write(7,2021)
      write(7,2021)
      write(7,2090)
      gtmga=0.
c JAF bug. initialize dummy1.
      dummy1=0.
      do 1135 n=1,nng
        dummy1=dummy1+gtsan(n)
        gtsan(n)=gtsan(n)/gttqan
 1135 gtmga=gtmga+gtsan(n)
      write(8,2021)
      write(8,270)
      if (erg.ge.1) then
        write(8,80) (gtsan(n),n=nng,1,-1)
      elseif (erg.le.-1) then
        write(8,80) (gtsan(n),n=1,nng)
      endif
      write(8,81) gtmga,ajwd1,ajwd2,ajwd3,ebaran
      write(8,2090)
      write(8,2021)
      write(8,2021)
      write(8,2090)
      write(8,2090)

c-----------------------------
c calculate s.f. neutron source
c-----------------------------
 1140 if (idd.eq.3) go to 1400
      if (isfnq.eq.0.and.nt.gt.0) go to 1260
      rewind 5
 1150 read (5,10) icont
      if (icont.ne.0) go to 1150
      write (6,280)
      tmgs=0.
      etopsf=0.
      ebotsf=0.
      do 1250 k=1,nq
 1160   read (5,110,end=1200) idq,nal,ndn
 1170   read (5,60) alam,bfsf,barnu,a,b,bfdn
        if (nal.eq.0) go to 1180
        read (5,120) (eal(l),fal(l),l=1,nal)
 1180   if (ndn.eq.0) go to 1190
        read (5,120) (fdng(nd),nd=1,ndn)
 1190   if (idq-jq(k)) 1160,1210,1200
 1200   write (6,160) jq(k)
c JAF no need to stop
c       stop 'S.F. source nuclide not found on tape5'
        cycle
c JAF bug--weird things happen to s.f. summary if no (alpha,n).
c   possibly fixed with "isfnq=isfnq+1" after 1210. check this.
 1210   continue
c       isfnq=isfnq+1
        qsf=aq(k)*alam*bfsf*barnu
        if(qsf.le.0.) go to 1250
        ebar=0.25*a*a*b + 1.5*a
        etopsf=etopsf+ebar*qsf
        ebotsf=ebotsf+qsf
        totqsf=totqsf+qsf
        lzq=idq/10000
        laq=idq/10-lzq*1000
        lsq=idq-lzq*10000-laq*10
        amq=' '
        if (lsq.ne.0) amq='m'

c---------------------------------------------
c output this s.f. neutron source contribution
c---------------------------------------------
        write (6,290) jsm(lzq),laq,amq,aq(k),alam,bfsf,barnu,qsf
        if (id.eq.1) go to 1250
        do 1220 n=1,nng
 1220   ssf(n)=0.

c----------------------------------------------------------------
c calc. frac. watt spectrum contained in neutron energy structure
c----------------------------------------------------------------
        if (a.gt.0.0.and.b.gt.0.0) go to 1230
        a=adef
        b=bdef
        write (7,300) jsm(lzq),laq,amq,a,b
 1230   se1=sqrt(en(1))
        sen1=sqrt(en(nng+1))
        sa=sqrt(a)
        sbaso4=sqrt(b*a*a/4.)
        c1=(sen1-sbaso4)/sa
        c2=(se1-sbaso4)/sa
        c3=(sen1+sbaso4)/sa
        c4=(se1+sbaso4)/sa
c JAF erf changed to derf
        if (c1.lt.0..and.c2.lt.0.) eft=0.5*(DERF(-c1)-DERF(-c2))
        if (c1.lt.0..and.c2.ge.0.) eft=0.5*(DERF(-c1)+DERF(c2))
        if (c1.ge.0..and.c2.ge.0.) eft=0.5*(DERF(c2)-DERF(c1))
        eft=eft+0.5*(DERF(c4)-DERF(c3))
        c1s=c1*c1
        c2s=c2*c2
        c3s=c3*c3
        c4s=c4*c4
        spba=sqrt(pi*b*a)
c       fwatt=eft+(1./spba)*(exp(-c1s)-exp(-c2s)-exp(-c3s)+exp(-c4s))

c-------------------------------------------
c calculate multigroup s.f. neutron spectrum
c-------------------------------------------
        w=qsf
        smgs=0.
        do 1240 n=1,nng
          seh=sqrt(en(n))
          sel=sqrt(en(n+1))
          c1=(sel-sbaso4)/sa
          c2=(seh-sbaso4)/sa
          c3=(sel+sbaso4)/sa
          c4=(seh+sbaso4)/sa
          c1s=c1*c1
          c2s=c2*c2
          c3s=c3*c3
          c4s=c4*c4
c JAF erf changed to derf
          if (c1.lt.0..and.c2.lt.0.) eft=0.5*(DERF(-c1)-DERF(-c2))
          if (c1.lt.0..and.c2.ge.0.) eft=0.5*(DERF(-c1)+DERF(c2))
          if (c1.ge.0..and.c2.ge.0.) eft=0.5*(DERF(c2)-DERF(c1))
          eft=eft+0.5*(DERF(c4)-DERF(c3))
          ssf(n)=(1./spba)*(exp(-c1s)-exp(-c2s)-exp(-c3s)+exp(-c4s))
          ssf(n)=w*(eft+ssf(n))
          ts(n)=ts(n)+ssf(n)
          smgs=smgs+ssf(n)
 1240   tssf(n)=tssf(n)+ssf(n)
        tmgs=tmgs+smgs

c-----------------------------------------------
c output this s.f. neutron spectrum contribution
c-----------------------------------------------
        write (7,310) aq(k),jsm(lzq),laq,amq
        if (erg.ge.1) then
          write (7,80) (ssf(n),n=nng,1,-1)
        elseif (erg.le.-1) then
          write (7,80) (ssf(n),n=1,nng)
        endif
        write (7,81) smgs,ajwd1,ajwd2,ajwd3,ebar
        smgs=0.
        do 1245 n=1,nng
          ssf(n)=ssf(n)/qsf
 1245   smgs=smgs+ssf(n)
        write (8,310) aq(k),jsm(lzq),laq,amq
        if (erg.ge.1) then
          write (8,80) (ssf(n),n=nng,1,-1)
        elseif (erg.le.-1) then
          write (8,80) (ssf(n),n=1,nng)
        endif
        write (8,81) smgs,ajwd1,ajwd2,ajwd3,ebar
        sfoutp2=smgs      
 1250 continue

c----------------------------------------------
c output total s.f. neutron source and spectrum
c----------------------------------------------
      if(isfnq.gt.1) write(6,320) totqsf
      if(id.eq.1.or.isfnq.lt.2) go to 1260
      ebarsf=etopsf/ebotsf
c     fracgp=tmgs/totqsf
      write(7,2021)
      write(7,330)
      if (erg.ge.1) then
        write(7,80) (tssf(n),n=nng,1,-1)
      elseif (erg.le.-1) then
        write(7,80) (tssf(n),n=1,nng)
      endif
      write(7,81) tmgs,ajwd1,ajwd2,ajwd3,ebarsf
      write(7,2090)
      write(7,2021)
      write(7,2090)
      tmgs=0.
c JAF bug. initialize dummy2.
      dummy2=0.
      do 1255 n=1,nng
        dummy2=dummy2+tssf(n)
        tssf(n)=tssf(n)/totqsf
 1255 tmgs=tmgs+tssf(n)
      write(8,2021)
      write(8,330)
      if (erg.ge.1) then
        write(8,80) (tssf(n),n=nng,1,-1)
      elseif (erg.le.-1) then
        write(8,80) (tssf(n),n=1,nng)
      endif
      write(8,81) tmgs,ajwd1,ajwd2,ajwd3,ebarsf
      write(8,2090)
      write(8,2021)
      write(8,2090)

c----------------------------------
c calculate delayed neutron sources
c----------------------------------
 1260 if (idnnq.eq.0) go to 1400
      rewind 5
      totqdn=0.
      if (id.eq.1) go to 1280
      do 1270 n=1,nng
 1270 tsdn(n)=0.
 1280 read (5,10) icont
      if (icont.ne.0) go to 1280
      write (6,340)
      tmgs=0.
      etopdn=0.
      ebotdn=0.
      do 1390 k=1,nq
 1290   read (5,110,end=1330) idq,nal,ndn
 1300   read (5,60) alam,bfsf,barnu,a,b,bfdn
        if (nal.eq.0) go to 1310
        read (5,120) (eal(l),fal(l),l=1,nal)
 1310   if (ndn.eq.0) go to 1320
        read (5,120) (fdng(nd),nd=1,ndn)
 1320   if (idq-jq(k)) 1290,1340,1330
 1330   write (6,160)
c JAF no need to stop
c       stop 'D.N. source nuclide not found on tape5'
        cycle
 1340   qdn=aq(k)*alam*bfdn
        etpdnq=0.
        ebtdnq=0.
        if(qdn.le.0.) go to 1390
        totqdn=totqdn+qdn
        lzq=idq/10000
        laq=idq/10-lzq*1000
        lsq=idq-laq*10-lzq*10000
        amq=' '
        if (lsq.eq.1) amq='m'
        if (lsq.eq.2) amq='n'

c-----------------------------------------------
c output the delayed neutron source contribution
c-----------------------------------------------
        write (6,350) jsm(lzq),laq,amq,aq(k),alam,bfdn,qdn
        if (id.eq.1) go to 1390

c---------------------------------------------------
c calculate fraction of this d.n. spectrum contained
c in the neutron energy structure.
c---------------------------------------------------
        fdn=1.
        fndn=ndn
        top=fndn*.01
        if (en(nngp1).eq.0..and.en(1).ge.top) go to 1355
        fdn=0.
        do 1350 nd=1,ndn
          fnd=nd
          top=fnd*.01
          bot=top-.01
          ebar=(top+bot)/2.
          val=qdn*fdng(nd)
          etpdnq=etpdnq+val*ebar
          ebtdnq=ebtdnq+val
          if (top.le.en(nngp1)) go to 1350
          if (bot.ge.en(1)) go to 1350
          frac=1.
          botlim=bot
          if (bot.lt.en(nngp1)) botlim=en(nngp1)
          toplim=top
          if (top.gt.en(1))  toplim=en(1)
          frac=(toplim-botlim)/.01
          fdn=fdn+frac*fdng(nd)
 1350   continue

c----------------------------------------------
c calculate multigroup delayed neutron spectrum
c----------------------------------------------
 1355   smgs=0.
        do 1380 n=1,nng
          sdn(n)=0.
          do 1370 nd=1,ndn
            fnd=nd
            top=fnd*.01
            bot=top-.01
            if(bot.ge.en(n)) go to 1375
            if(top.lt.en(n+1)) go to 1370
            frac=1.
            if(bot.ge.en(n+1).and.top.le.en(n)) go to 1360
            if(top.gt.en(n)) top=en(n)
            if(bot.lt.en(n+1)) bot=en(n+1)
            frac=(top-bot)/.01
 1360       sdn(n)=sdn(n)+frac*fdng(nd)*qdn
 1370     continue
 1375     smgs=smgs+sdn(n)
          tsdn(n)=tsdn(n)+sdn(n)
          ts(n)=ts(n)+sdn(n)
 1380   continue
        tmgs=tmgs+smgs

c--------------------------------------------------
c output this delayed neutron spectrum contribution
c--------------------------------------------------
        etopdn=etopdn+etpdnq
        ebotdn=ebotdn+ebtdnq
        ebrdnq=etpdnq/ebtdnq
        write (7,360) aq(k),jsm(lzq),laq,amq
        if (erg.ge.1) then
          write (7,80) (sdn(n),n=nng,1,-1)
        elseif (erg.le.-1) then
          write (7,80) (sdn(n),n=1,nng)
        endif
        write (7,81) smgs,ajwd1,ajwd2,ajwd3,ebrdnq
        smgs=0.
        do 1385 n=1,nng
          sdn(n)=sdn(n)/qdn
 1385   smgs=smgs+sdn(n)
        write (8,360) aq(k),jsm(lzq),laq,amq
        if (erg.ge.1) then
          write (8,80) (sdn(n),n=nng,1,-1)
        elseif (erg.le.-1) then
          write (8,80) (sdn(n),n=1,nng)
        endif
        write (8,81) smgs,ajwd1,ajwd2,ajwd3,ebrdnq
 1390 continue ! k

c-------------------------------------------------
c output total delayed neutron source and spectrum
c-------------------------------------------------
      if(idnnq.gt.1) write(6,370) totqdn
      if(id.eq.1.or.idnnq.lt.2) go to 1400
      ebardn=etopdn/ebotdn
c     fracgp=tmgs/totqdn
      write(7,2021)
      write(7,380)
      if (erg.ge.1) then
        write(7,80) (tsdn(n),n=nng,1,-1)
      elseif (erg.le.-1) then
        write(7,80) (tsdn(n),n=1,nng)
      endif
      write(7,81) tmgs,ajwd1,ajwd2,ajwd3,ebardn
      write(7,2090)
      write(7,2021)
      write(7,2090)
      tmgs=0.
c JAF bug. initialize dummy3.
      dummy3=0.
      do 1395 n=1,nng
        dummy3=dummy3+tsdn(n)
        tsdn(n)=tsdn(n)/totqdn
 1395 tmgs=tmgs+tsdn(n)
      write(8,2021)
      write(8,380)
      if (erg.ge.1) then
        write(8,80) (tsdn(n),n=nng,1,-1)
      elseif (erg.le.-1) then
        write(8,80) (tsdn(n),n=1,nng)
      endif
      write(8,81) tmgs,ajwd1,ajwd2,ajwd3,ebardn
      write(8,2090)
      write(8,2021)
      write(8,2090)
 1400 continue
      if (id.eq.1) go to 1430
      itest=0
      if(nt.gt.0) itest=1
      if(isfnq.gt.0) itest=itest+1
      if(idnnq.gt.0) itest=itest+1
      if(itest.lt.2) then
c JAF remove dummy; it is unused and immediately set to 0 in the
c next section.
c JAF bug. initialize dummy.
c       dummy=0.
        do 1405 n=1,nng
c         dummy=dummy+ts(n)
c JAF bug fix: gttqan is 0 if no (alpha,n) (or maybe the code
c never gets here?)
c         ts(n)=ts(n)/gttqan
          if(gttqan.ne.0.)ts(n)=ts(n)/gttqan
 1405   gtmg=gtmg+ts(n)
c JAF bug fix: why go to 1430?
c       go to 1430
      endif

c------------------------------------------------------------
c output grand total (alpha,n) + s.f. + d.n. neutron spectrum
c------------------------------------------------------------
      gtq=gttqan+totqsf+totqdn
      dummy=0.0
      gtmg=0.
      do 1420 n=1,nng
 1420 gtmg=gtmg+ts(n)
c     fracgp=gtmg/gtq
      etpall=etopan+etopsf+etopdn
      ebtall=ebotan+ebotsf+ebotdn
      ebrall=etpall/ebtall
      write(7,2090)
      write(7,2021)
      write(7,390)
      if (erg.ge.1) then
        write(7,80) (ts(n),n=nng,1,-1)
      elseif (erg.le.-1) then
        write(7,80) (ts(n),n=1,nng)
      endif
      write(7,81) gtmg,ajwd1,ajwd2,ajwd3,ebrall
      gtmg=0.
      do 1425 n=1,nng
        dummy=dummy+ts(n)
        ts(n)=ts(n)/gtq
 1425 gtmg=gtmg+ts(n)
      write(8,2090)
      write(8,2021)
      write(8,390)
      if (erg.ge.1) then
c JAF bug. was 7
        write(8,80) (ts(n),n=nng,1,-1)
      elseif (erg.le.-1) then
c JAF bug. was 7
        write(8,80) (ts(n),n=1,nng)
      endif
      write(8,81) gtmg,ajwd1,ajwd2,ajwd3,ebrall
 1430 continue

c---------------------------------------------------
c spectra for outp2 in ascending or descending order
c---------------------------------------------------
      if(isfnq.lt.2) then
        if (erg.ge.1) then
          write(12,1900)
          write(12,3600)title
          write(12,3601)
          write(12,3605)
          write(12,3610)
          write(12,3615)
          write(12,3620)
          write(12,3630)(en(n),gtsan(n),ssf(n),tsdn(n),ts(n),
     1     n=nngp1,1,-1)
          write(12,3635)gtmga,sfoutp2,totqdn,gtmg
          write(12,3611)
          write(12,3610)
          write(12,3615)
          write(12,3620)
          write(12,3630)(en(n),gtsan(n)*gttqan,ssf(n)*totqsf,
     1    tsdn(n)*totqdn,ts(n)*gtq,n=nngp1,1,-1)  
          write(12,3635)gttqan,totqsf,totqdn,gtq   
        elseif (erg.le.-1) then
          write(12,1900)
          write(12,3600)title
          write(12,3601)
          write(12,3605)
          write(12,3610)
          write(12,3615)
          write(12,3620)
          write(12,3630)(en(n),gtsan(n),ssf(n),tsdn(n),ts(n),n=1,nngp1)
c JAF bug. why would sfoutp2 become tmgs just because group order
c changes? but this works unless nngp1=2?
c         write(12,3635)gtmga,tmgs,totqdn,gtmg
          write(12,3635)gtmga,sfoutp2,totqdn,gtmg
          write(12,3611)
          write(12,3610)
          write(12,3615)
          write(12,3620)
          write(12,3630)(en(n),gtsan(n)*gttqan,ssf(n)*totqsf,
     1    tsdn(n)*totqdn,ts(n)*gtq,n=1,nngp1)
          write(12,3635)gttqan,totqsf,totqdn,gtq
        endif

      elseif(isfnq.ge.2) then
        nngp1=nng+1

        if (erg.ge.1) then
          write(12,1900)
          write(12,3600)title
          write(12,3601)
          write(12,3605)
          write(12,3610)
          write(12,3615)
          write(12,3620)
          write(12,3630)(en(n),gtsan(n),tssf(n),tsdn(n),ts(n),
     1     n=nngp1,1,-1)
          write(12,3635)gtmga,tmgs,totqdn,gtmg
          write(12,3611)
          write(12,3610)
          write(12,3615)
          write(12,3620)
          write(12,3630)(en(n),gtsan(n)*gttqan,tssf(n)*totqsf,
     1    tsdn(n)*totqdn,ts(n)*gtq,n=nngp1,1,-1)  
          write(12,3635)gttqan,totqsf,totqdn,gtq     
        elseif (erg.le.-1) then
          write(12,1900)
          write(12,3600)title
          write(12,3601)
          write(12,3605)
          write(12,3610)
          write(12,3615)
          write(12,3620)
          write(12,3630)(en(n),gtsan(n),tssf(n),tsdn(n),ts(n),n=1,nngp1)
          write(12,3635)gtmga,tmgs,totqdn,gtmg
          write(12,3611)
          write(12,3610)
          write(12,3615)
          write(12,3620)
          write(12,3630)(en(n),gtsan(n)*gttqan,tssf(n)*totqsf,
     1    tsdn(n)*totqdn,ts(n)*gtq,n=1,nngp1)
          write(12,3635)gttqan,totqsf,totqdn,gtq
        endif
      endif ! isfnq
c-----------------
c  Output Summary
c-----------------
 1440 continue
      qtotal=gttqan+totqsf+totqdn
      ebrall=ebaran*gttqan
      if (isfnq.ne.0) ebrall=ebrall+(ebarsf*totqsf)
      if (idnnq.ne.0) ebrall=ebrall+(ebardn*totqdn)
      ebrall=ebrall/qtotal
      if (nq.eq.0) ebrall=ebaran
      write(11,2090)
      write(11,2090)
      write(11,3000)
      write(11,3010)
      write(11,2090)
      if (idd.eq.1) then
        write(11,3020)gttqan
        write(11,3022)totqsf
        write(11,3024)totqdn
        write(11,3026)qtotal
      else
        write(11,3021)gttqan
      endif
      if (id.gt.1) then
        write(11,3030)ebaran
        if (idd.eq.1) then
          write(11,3032)ebarsf
          write(11,3034)ebardn
          write(11,3036)ebrall
        endif
        write(11,3050)
        do 1450 n=1,nng
          write(11,3051)n,ts(n)
1450    continue
        frac=100.0*(dummy/qtotal)
        write(11,3500)frac
        if (frac.le.90.0) write(11,3510)
        if (itest.gt.2) then
          frac=100.0*(dummy1/gttqan)
          write(11,3501)frac
        endif
        if (isfnq.ne.0) then
          frac=100.0*(dummy2/totqsf)
          write(11,3502)frac
          if (frac.le.90.0) write(11,3510)
        endif
        if (idnnq.ne.0) then
          frac=100.0*(dummy3/totqdn)
          write(11,3503)frac
          if (frac.le.90.0) write(11,3510)
        endif
      endif ! id

c-------------------------
c  End Execution of HOMOG
c-------------------------
      return
      end



c=======================================================================
c  Interface Problem Main Subroutine (8/96)
c=======================================================================

      subroutine interf(title,idd,id,erg)

c----------
c  Storage
c----------
      
      character title*7
      integer erg
c JAF introduce parameters
      integer maxel,maxnag
      parameter (maxel=100,maxnag=2001)
c JAF replaced both 4001's in this list with maxnag.
      dimension title(11), ast(maxnag), eal(maxnag)

c------------------
c  Call Statements
c------------------

c JAF pass in maxel,maxnag
      call sainf(maxel,maxnag,nag,title,idd,id,erg,ast,eal)
      call sninf(maxel,maxnag,nag,title,id,erg,ast,eal)

c--------------------------------
c  End of Execution of INTERFACE
c--------------------------------

      return
      end


c=======================================================================
c  Interface Problem Alpha Source Subroutine (8/96)
c=======================================================================

      subroutine sainf(maxel,maxnag,nag,title,idd,id,erg,ast,eal)

c---------
c storage
c---------

      character title*7,jsm(105)*2
      integer erg
c JAF introduce parameters
      integer maxel,maxnag
c JAF replaced all 4001's in this list with maxnag.
c JAF replaced all 20's in this list with maxel.
      dimension title(11), jzm(maxel), azm(maxel), czm(9,maxel), c(14),
     1 jq(300),
     1 aq(300), eala(30), eal(maxnag), ee(maxnag), scx(maxnag),
     2 r(maxnag), 
     3 p(maxnag), fdng(1000), ast(maxnag),fal(30),ea(maxnag),ps(maxnag),
     4 faq(300)
       common /names/ jsm
       common /masses/ amass(105)

c--------------------
c  Format Statements
c--------------------
   10 format (i3,11a7)
   15 format(11a7)
   16 format(11a7,//)
   20 format (2i3,f7.3)
   30 format (i8,e12.5)
   40 format(i3,8e8.3,f10.3,/,3x,5f8.3)
   50 format (69h failed to find stopping cross section coefficients on
     1tape 2 for iz=,i3,5h stop)
   60 format (6e12.5)
   90 format (i3,5e12.5)
   95 format(2i6)
   96 format(/,26x,'Neutron Source Magnitudes',/,26x,25(1h~),/)
   97 format(/,23x,'Absolute Neutron Source Spectra',/,23x,31(1h~),/)
   98 format(/,22x,'Normalized Neutron Source Spectra',/,22x,
     133(1h~),/)
   99 format(/,16x,'Neutron Source Spectra by Nuclide Energy Level',
     1/,16x,46(1h~),/)
  110 format (i8,2i4)
  120 format (8e10.3)
  130 format (1x,'SOURCE DOES NOT HAVE ALPHA EMITTER')
  131 format(' WARNING-AN ALPHA ENERGY THAT EXCEEDED ',f7.3,' MEV WAS',
     1' SET TO:',f7.3,' MeV')
  132 format(' NUMBER OF STOPPING ELEMENTS ON SOURCE SIDE:',i3,/,
     1' SOLID OR GAS STOPPING POWER TRIGER:',i3,/,
     2' MAXIMUM ENERGY FOR ALPHA SPECTRA:',f7.3,/)
  133 format(' ISOTOPE IDENTIFICATION:',i6,/
     1' ISOTOPE FRACTION:',e15.6,/)
  134 format(' NUMBER OF SOURCES:',i3)
  135 format(/,' SOURCE IDENTIFICATION: ',i6,/
     1' SOURCE ISOTOPIC FRACTION:',e15.6)
  136  format(//,' ALPHA SPECTRA ENERGY BOUNDARIES',/)
  137  format(5e14.6)
  138 format(' MAXIMUM ALPHA DECAY ENERGIES')
  160 format (56h failed to find alpha/s.f./dn source data on tape 5 for
     1 ,i8,5h stop)

 1900 format('SOURCES 4D Calculation',/,'<<<<<<<<<<<>>>>>>>>>>>',//)
 2000 format(/,22x,'Summary of Input')
 2010 format(22x,'================')
 2020 format(' ')
 2021 format(/,80(1h-))
 2030 format('Title: ', 11a7,/)
 2041 format('Ascending energy structure for output (erg =',i3,').')
 2042 format('Descending energy structure for output (erg =',i3,').')
 2043 format('Interface problem input (idd =',i2,')')
 2048 format('Magnitudes and spectra computed (id =',i2,').')
 2049 format('Magnitudes only computed (id =',i2,').')
 2050 format('Number of elemental constituents on source side:',i3)
 2060 format('Solid stopping cross-sections used (isg =',i2,') on',
     1 1x,'source side.')
 2070 format('Gas stopping cross-sections used (isg =',i2,')',
     1 1x,'on sources side.')
 2072 format('Maximum energy for alpha spectra: ',1pE10.3,' MeV.')
 2074 format('Minimum energy for alpha spectra: ',1pE10.3,' MeV.')
 2080 format('Elemental Constituents on Source Side:')
 2090 format(' ')
 2100 format('  Z-value  Atom Fraction')
 2110 format('  -------  -------------')
 2120 format(4x,i3,5x,f12.10)
 2200 format('Number of source nuclides to be evaluated:',i3)
 2220 format('Source Nulcides:')
 2230 format('    ZAID     Atom Fraction')
 2240 format('    ----     -------------')
 2250 format(3x,i6,4x,1pE10.3)
 2270 format(i5,' alpha energy groups used at interface.')

c-----------------
c  Set Parameters
c-----------------
      rewind 2
      rewind 3
      rewind 4
      rewind 5
c JAF aneut and pi (already commented) deleted
c JAF added rnapier; "d0" added to reals
      alph=4.d0
      zalp=2.d0
      rnapier=2.7182818284590d+0
c JAF remove--unused
c 400 continue

c------------------
c  Read from Tape1
c------------------
      write(6,1900)
      write(7,1900)
      write(8,1900)
      write(10,1900)
      write(11,1900)
      write(11,2000)
      write(11,2010)
      write(11,2020)
      write(11,2030)title
      write(11,2043)idd
      if (id.eq.1) then
       write(11,2049)id
      elseif (id.eq.2) then
       write(11,2048)id
       write(7,97)
       write(8,98)
       if (erg.ge.1) then
         write(11,2041)erg
         elseif (erg.le.-1) then
         write(11,2042)erg
         else
         stop 'Error in erg input!'
         endif       
      else
       stop 'Error in id input!'
      endif
      write(6,96)
      write(6,2030)title
      write(7,2030)title
      write(8,2030)title
      read (1,*) nz,isg,em,eamin
      write(11,2050)nz
c JAF error check (16 lines)
      if(nz.gt.maxel)then
        write(*,'("error. too many elements."/,
     1   "nz=",i0," maxel=",i0)')nz,maxel
        close(unit=1)
        close(unit=2)
        close(unit=3)
        close(unit=4)
        close(unit=5)
        close(unit=6, status='delete')
        close(unit=7, status='delete')
        close(unit=8, status='delete')
        close(unit=10, status='delete')
        close(unit=11)
        close(unit=12, status='delete')
        stop
      end if
      if (isg.eq.0) then
       write(11,2060)isg
      elseif (isg.eq.1) then
       write(11,2070)isg
      else
       stop 'Error in isg input!'
      endif
      write(11,2072)em
      write(11,2074)eamin
      write(11,2090)
      write(11,2080)
      write(11,2090)
      write(11,2100)
      write(11,2110)
      if (eamin.lt.0.001) eamin=0.001

c---------------------------------
c read input material constituents
c---------------------------------
      do 420 j=1,nz
       read (1,*) jzm(j),azm(j)
       write(11,2120)jzm(j),azm(j)
  420 continue

c----------------------------------------
c order the nz material constituents by z
c----------------------------------------
      do 440 j=1,nz
       isw=0
       do 430 jn=2,nz
        if (jzm(jn).gt.jzm(jn-1)) go to 430
        jtem=jzm(jn-1)
        atem=azm(jn-1)
        jzm(jn-1)=jzm(jn)
        azm(jn-1)=azm(jn)
        jzm(jn)=jtem
        azm(jn)=atem
        isw=1
  430  continue
       if (isw.eq.0) go to 450
  440 continue

c----------------------------------------------------
c read tape 2 for stopping cross-section coefficients
c----------------------------------------------------
  450 rewind 2
  460 read (2,10) icont
      if (icont.ne.0) go to 460
      do 520 j=1,nz
  470 read (2,40,end=510) iz,c
  480 if (jzm(j)-iz) 510,490,470
  490 continue
      do 500 jt=1,9
  500 czm(jt,j)=c(jt)
      if (isg.eq.1.and.c(10).gt.0.) then
      do 505 jt=1,5
  505 czm(jt,j)=c(jt+9)
      endif
      go to 520
  510 write (6,50) jzm(j)
      stop 1
  520 continue
      read (1,*)nag
      write(11,2090)
      write(11,2270)nag
c JAF error check (16 lines)
      if(nag.gt.maxnag-1)then
        write(*,'("error. too many alpha energy groups."/,
     1   "nag=",i0," maxnag=",i0)')nag,maxnag-1
        close(unit=1)
        close(unit=2)
        close(unit=3)
        close(unit=4)
        close(unit=5)
        close(unit=6, status='delete')
        close(unit=7, status='delete')
        close(unit=8, status='delete')
        close(unit=10, status='delete')
        close(unit=11)
        close(unit=12, status='delete')
        stop
      end if
      do 948 m=1,nag
  948 ast(m)=0.0

c---------------------------------------------------------
c set alpha energy grid ee(m) for this target/source combo
c---------------------------------------------------------
      eamax=em
  860 fnag=nag
      nagp1=nag+1
      dea=(eamax-eamin)/fnag
      ee(1)=eamin
      do 870 m=2,nag
       fm=m-1
  870 ee(m)=ee(1)+dea*fm
      ee(nagp1)=eamax
      do 871 m=1,nag
  871 ea(m)=(ee(m+1) + ee(m))/2.

c-------------------------------------------------------------
c for each energy ee(m) calc. and store stopping cross section
c-------------------------------------------------------------
      do 930 m=1,nagp1
       scx(m)=0.
       do 920 j=1,nz
        zmat=jzm(j)
        amat=amass(jzm(j))
c-----  nuclear stopping---
        term=zalp**0.66667 + zmat**0.66667
        sterm=sqrt(term)
        bot=zalp*zmat*(alph+amat)*sterm
        rep=32530.*amat*ee(m)/bot
        if(rep.lt.0.001) then
         dcx=1.593*sqrt(rep)
         go to 905
         endif
        if(rep.lt.10.) then
         bot=1.+6.8*rep+3.4*rep**1.5
         dcx=1.7*sqrt(rep)*alog(rep+rnapier)/bot
         go to 905
        endif
        dcx=alog(0.47*rep)/(2.*rep)
c-----  electronic stopping---
  905   if(ee(m).gt.30.) go to 906
        slow=czm(1,j)*(1000.*ee(m))**czm(2,j)
        shigh=(czm(3,j)/ee(m))
        shigh=shigh*alog(1.+czm(4,j)/ee(m)+czm(5,j)*ee(m))
        dcx=dcx+slow*shigh/(slow+shigh)
        go to 907
  906   eil=alog(1/ee(m))
        arg=czm(6,j)+eil*(czm(7,j)+czm(8,j)*eil+czm(9,j)*eil*eil)
        dcx=dcx+exp(arg)
  907  continue
  920  scx(m)=scx(m)+azm(j)*dcx
  930 continue

c-------------------------------------------
c read alpha sources from user input
c---------------------------------------------
  580 read (1,*) nq
      write(11,2200)nq
      write(11,2090)
      write(11,2220)
      write(11,2090)
      write(11,2230)
      write(11,2240)
      do 590 k=1,nq
       read (1,*) jq(k),faq(k)
       write(11,2250)jq(k),faq(k)
  590 continue

c-----------------------------------------------
c order sources by z-a-state id = s+10*a+10000*z
c-----------------------------------------------
      if (nq.eq.1) go to 620
      do 610 k=1,nq
       isw=0
       do 600 kn=2,nq
        if (jq(kn).gt.jq(kn-1)) go to 600
        jtem=jq(kn-1)
        atem=aq(kn-1)
        jq(kn-1)=jq(kn)
        aq(kn-1)=aq(kn)
        jq(kn)=jtem
        aq(kn)=atem
        isw=1
  600  continue
       if (isw.eq.0) go to 620
  610 continue

c----------------------------------------
c read target information from user input
c----------------------------------------
  620 continue

c-------------------------------------------
c read source nuclide decay data from tape 5
c-------------------------------------------
  770 read (5,10) icont
      if (icont.ne.0) go to 770
      do 1100 k=1,nq
  790 read (5,110,end=830) idq,nal,ndn
  800 read (5,60) alam,bfsf,barnu,a,b,bfdn
      if (nal.eq.0) go to 810
      read (5,120) (eala(l),fal(l),l=1,nal)
  810 if (ndn.eq.0) go to 820
      read (5,120) (fdng(nd),nd=1,ndn)
  820 if (idq-jq(k)) 790,840,830
  830 write (6,160) jq(k)
      stop 4
  840 continue
      if(nal.eq.0)write(6,130)
      if(nal.eq.0)go to 841
  841 continue

c---------------------------------------------------------
c for each energy group cal integral 1/scx
c---------------------------------------------------------
      r(1)=1./scx(1)
      fact=alam*faq(k)*2.5e+20
      do 940 m=1,nag
      ps(m)= 0.0
      r(m+1)= 1.0/scx(m+1)
      p(m)=fact*(r(m+1)+r(m))*dea/2.
  940 continue
c---------------------------------------------------------
      do 946 mm = 1,nal
       if(eala(mm).gt.em)write(6,131)em,em
       if(eala(mm).gt.6.5)eala(mm)=6.5
      do 945  m = 1,nag
      f=1.
      if(eala(mm).ge.ee(m+1))go to 944
      if(eala(mm).gt.ee(m))f=(eala(mm)-ee(m))/dea
      if(eala(mm).lt.ee(m))f=0.0
  944 ps(m) = ps(m) +fal(mm)*p(m)*f
  945 continue
  946 continue
      do 947 m=1,nag
      ast(m) = ast(m) + ps(m)
  947 continue
c----------------------------------------------------------
 1100 continue
      knt=nag
      do 2 i=1,nag
       if(ast(knt).gt.0.0)go to 3
       knt=knt-1
 2    continue
 3    nag=knt
      if(nag.eq.0)write(6,130)
      if(nag.eq.0)stop 20
      do 1 i=1,nag
       eal(i)=ea(i)
 1    continue

c-------------------------
c  End Execution of SAINF
c-------------------------
      return
 1440 continue
      stop
      end



c=======================================================================
c  Interface Problem Neutron Source Subroutine (8/96)
c=======================================================================

      subroutine sninf(maxel,maxnag,nal,title,id,erg,ast,eal)

c---------
c storage
c---------
      character title*7,title2*7,jsm(105)*2,amt*1
      integer erg
c JAF introduce parameters
      integer maxel,maxnag
c JAF replaced all 4001's in this list with maxnag.
c JAF replaced some 20's in this list with maxel.
      dimension title(11),title2(11),jzm(maxel),azm(maxel),czm(9,maxel),
     1 c(14),
     1 el(20), ep(200), f(20,200), e(1100), x(1100),
     2 ee(maxnag), scx(maxnag), en(751), san(750), eal(maxnag),
     3 ts(750), cx(maxnag), r(maxnag), rr(maxnag), p(maxnag), tsan(750),
     4 gtsan(750), s(750), totlev(20), ast(maxnag),
     5 sl(750,20), etopl(20), ebotl(20)
      common /names/ jsm
      common /masses/ amass(105)
c     dimension frclev(20)

c--------------------
c formats statements
c--------------------
   10 format (11a7)
   15 format(//,11a7,/)
   20 format (3i3)
   30 format (i8,e12.5)
   40 format(i3,8e8.3,f10.3,/,3x,5f8.3)
   50 format (69h failed to find stopping cross section coefficients on
     1tape 2 for iz=,i3,5h stop)
   60 format (6e12.5)
   70 format(/,1x,'Neutron Multigroup Structure (MeV)')
   80 format (1p8e10.3)
   81 format(31x,'Total (all groups): ',1pE10.3,' neutrons/sec-cm^2.',
     1/,27x,'Average Neutron Energy: ',1pE10.3,' MeV.')
   90 format (i3,5e12.5)
   95 format(2i6)
  100 format(//,1h1,30x,'Table I',/,31x,7(1h=),//,12x,'(alpha,n) ',
     1'Neutron Production by Target and Alpha Energy',///,17x,
     2'target            alpha   alphas/sec     p(e)     neuts/sec',
     3/,7x,'target  atom frac.          energy     /cm^2',4x,
     4'neut/alpha',4x,'/cm^2',/,1h+,6x,6(1h_),2x,10(1h_),2x,6(1h_),2x,
     5 6(1h_),3(12h  __________))
c 102 format(1h+,77x,37hfractional ground/excited level split,/,
c    1 1h+,77x,49(1h_))
  105 format(i8,i4)
  120 format (8e10.3)
  130 format (52h failed to find alpha-n cross section on tape 3 for ,i8
     1 ,5h stop)
  140 format (i8,2i3)
  150 format (66h failed to find comp.nucl.neut.emiss.lev. branching on 
     1tape 4 for ,i8,5h stop)
  170 format (/,17h alpha energy of ,1pe12.5,62h exceeds current 6.5 MeV
     1 limit of most cross sections, stop.  )
  180 format (7h target,i8,40h cross-section data not available at ee(
     1 ,i4,3h)= ,1pe12.5,6h stop.)
  200 format(33x,f8.3,1p3e12.4)
c 205 format(1h+,(76x,7f7.5,/,1x))
  210 format(7x,a2,i3,a1,1pe12.4,8x,         0pf8.3,1p3e12.4)
  220 format (14h alpha energy ,1pe12.5,45h exceeds range of level branc
     1hing data, stop.)
  240 format(/,' (a,n) neutrons from alphas on ',a2,i3,a1,' in target.')
  250 format(1h+,66x,10(1h_),/,55x,'Total:',4x,1pe12.4,/)
  251 format(1h+,66x,10(1h_),/,41x,'Total (this target):',
     14x,1pe12.4,/)
  252 format (1h+,66x,10(1h_),/,41x,'Total (all targets):',
     14x,1pe12.4,/)
  255 format(1h )
  260 format (///,' Total (alpha,n) neutron spectrum this target')
  270 format(///,' Grand total (alpha,n) neutron spectrum, all targets,'
     1 ,1x,'all sources')
  390 format (////,' Total Neutron Spectrum')
 1900 format('SOURCES 4D Calculation',/,'<<<<<<<<<<<>>>>>>>>>>>',//)
 2021 format(/,80(1h-))
 2030 format('Target title: ', 11a7)
 2050 format('Number of elemental constituents on target side:',i3)
 2060 format('Solid stopping cross-sections used (isg =',i2,') on',
     1 1x,'target side.')
 2070 format('Gas stopping cross-sections used (isg =',i2,')',
     1 1x,'on target side.')
 2082 format('Elemental Constituents on Target Side:')
 2090 format(' ')
 2100 format('  Z-value  Atom Fraction')
 2110 format('  -------  -------------')
 2120 format(4x,i3,5x,f12.10)
 2130 format('Number of neutron spectrum energy groups:', i4)
 2140 format('Maximum neutron energy is ',1pE12.5,' MeV.')
 2150 format('Minimum neutron energy is ',1pE12.5,' MeV.')
 2160 format('Energy Group Structure:')
 2170 format('  Group  Upper-Bound  Lower-Bound')
 2180 format('  -----  -----------  -----------')
 2190 format(2x,i4,4x,1pE12.5,3x,1pE12.5)
 2193 format(/,'Warning: First Group Upper-Bound Does Not Coincide ',
     1'With Maximum Neutron Energy.')
 2195 format(/,'Warning: Energy Group Structure Not Entered in ',
     1'Decreasing Order.',/,5x,'Group ',i3,' Upper-Bound is Greater ',
     2'Than Group ',i3,' Upper-Bound.') 
 2260 format(/,'Number of target nuclides to be used:',i3)
 2270 format(i5,' alpha energy groups used in target calculation.')
 2280 format('Target Nuclides:')
 2290 format('    ZAID     Atom Fraction')
 2300 format('    ----     -------------')
 2310 format(3x,i6,4x,1pE12.5)

 3000 format(//,22x,'Summary of Output')
 3010 format(22x,'=================')
 3020 format('Total (alpha,n) neutron source from all sources and ',
     1'targets: ',1pE12.5,' n/sec-cm^2.')
 3030 format(/,'Average (alpha,n) neutron energy: ',1pE12.5,' MeV.')
 3050 format(/,'Normalized Neutron Energy Spectrum by Energy Group ',
     1'for All Sources:',/'(Note: group structure is independent of',
     2' erg record!)',//,10x,'Group',5x,'Contribution',/,10x,
     3'-----',5x,'------------')
 3051 format(11x,i3,7x,1pE12.5)
 3500 format(/,'Portion of Total Neutron Source Rate Accounted for ',
     1'in the Energy Spectrum: ',F5.1,'%.')
 3510 format(/,'Warning: Energy Group Structure May Be Poorly Chosen.')
 3600 format(' Title: ', 11a7)
 3601 format(' Neutron Source Spectra in Columnar Format.'/
     1' Examine other tape files for additional information.'/
     2' Recall the Energy values are bin bounds.'/) 
 3605 format(/,'           Normalized Total (neutrons/sec-cm^2)')
 3610 format('           ========================================')
 3611 format(/,'            Absolute Total (neutrons/sec-cm^2)')
 3615 format(' Energy     (a,n)')
 3620 format(1x,9('-'),2x,9('-'))
 3625 format(1p5e10.3)
 3630 format(1p1e10.3,1p1e11.3)
 3635 format(' Total    ',1p1e11.3) 
c-----------------
c  Set Parameters
c-----------------
      rewind 2
      rewind 3
      rewind 4
      rewind 5
c JAF added rnapier; "d0" added to reals
      aneut=1.d0
      alph=4.d0
      zalp=2.d0
      rnapier=2.7182818284590d+0
  400 gttqan=0.

c------------------
c  Read from Tape1
c------------------
      read (1,10)title2
      write(11,2090)
      write(11,2030)title2
      read (1,*) nz,isg
      write(11,2050)nz
c JAF error check (16 lines)
      if(nz.gt.maxel)then
        write(*,'("error. too many elements."/,
     1   "nz=",i0," maxel=",i0)')nz,maxel
        close(unit=1)
        close(unit=2)
        close(unit=3)
        close(unit=4)
        close(unit=5)
        close(unit=6, status='delete')
        close(unit=7, status='delete')
        close(unit=8, status='delete')
        close(unit=10, status='delete')
        close(unit=11)
        close(unit=12, status='delete')
        stop
      end if
      if (isg.eq.0) then
      write(11,2060)isg
      elseif (isg.eq.1) then
       write(11,2070)isg
      else
       stop 'Error in isg input!'
      endif
      write(11,2090)
      write(11,2082)
      write(11,2090)
      write(11,2100)
      write(11,2110)
  410 continue

c---------------------------------
c read input material constituents
c---------------------------------
      do 420 j=1,nz
       read (1,*) jzm(j),azm(j)
  420 write(11,2120)jzm(j),azm(j)

c----------------------------------------
c order the nz material constituents by z
c----------------------------------------
      do 440 j=1,nz
       isw=0
       do 430 jn=2,nz
        if (jzm(jn).gt.jzm(jn-1)) go to 430
        jtem=jzm(jn-1)
        atem=azm(jn-1)
        jzm(jn-1)=jzm(jn)
        azm(jn-1)=azm(jn)
        jzm(jn)=jtem
        azm(jn)=atem
        isw=1
  430  continue
       if (isw.eq.0) go to 450
  440 continue

c----------------------------------------------------
c read tape 2 for stopping cross-section coefficients
c----------------------------------------------------
  450 rewind 2
  460 read (2,20) icont
      if (icont.ne.0) go to 460
      do 520 j=1,nz
  470 read (2,40,end=510) iz,c
  480 if (jzm(j)-iz) 510,490,470
  490 continue
      do 500 jt=1,9
  500 czm(jt,j)=c(jt)
      if (isg.eq.1.and.c(10).gt.0.) then
      do 505 jt=1,5
  505 czm(jt,j)=c(jt+9)
      endif
      go to 520
  510 write (6,50) jzm(j)
      stop 1
  520 continue
  525 if (id.eq.1) go to 580

c------------------------------------------------------------
c construct neutron group structure desired from tape 1 input
c neutron groups ordered in decreasing energy
c------------------------------------------------------------
      read (1,*) nng,enmax,enmin
      write(11,2090)
      write(11,2130)nng
      write(11,2140)enmax
      write(11,2150)enmin
      write(11,2090)
      write(11,2160)
      write(11,2090)
      write(11,2170)
      write(11,2180)
      nngp1=nng+1
      if (nng.gt.0) go to 540
      nng=-nng
      nngp1=nng+1
      read (1,*) (en(n),n=1,nng)
      en(nngp1)=enmin
      if (en(1).gt.en(2)) go to 560
      nh=nngp1/2
      do 530 n=1,nh
       etem=en(nngp1-n+1)
       en(nngp1-n+1)=en(n)
  530 en(n)=etem
      go to 560
  540 en(nngp1)=enmin
      fnng=nng
      den=(enmax-enmin)/fnng
      en(1)=enmax
      do 550 n=2,nng
       fnm1=n-1
  550 en(n)=en(1)-fnm1*den
  560 do 562 n=1,nng
       write(11,2190)n,en(n),en(n+1)
  562 continue
      if (en(1).ne.enmax) write(11,2193)
        do 528 n=2,nng
         nm1=n-1
         if (en(n).gt.en(n-1)) then
          write(11,2195)n,nm1
          stop 'Energy Structure Incorrectly Entered!'
         endif
  528 continue
      write (7,70)
      if (erg.ge.1) then
      write (7,80) (en(n),n=nngp1,1,-1)
      elseif (erg.le.-1) then
      write (7,80) (en(n),n=1,nngp1)
      endif
      write (8,70)
      if (erg.ge.1) then
      write (8,80) (en(n),n=nngp1,1,-1)
      elseif (erg.le.-1) then
      write (8,80) (en(n),n=1,nngp1)
      endif

c----------------------------
c zero total spectrum storage
c----------------------------
      do 570 n=1,nng
       ts(n)=0.
       gtsan(n)=0.
  570 etopan=0.
      ebotan=0.
      gttqan=0.

c----------------------------------------------------
c    alpha sources treated as if one isotope ie nq=1
c----------------------------------------------------
  580 continue
      nq=1

c----------------------------------------
c read target information from user input
c----------------------------------------
  620 continue
      read (1,*) nt, nag
      write(11,2260)nt
      write(11,2270)nag
      write(11,2090)
      write(11,2280)
      write(11,2090)
      write(11,2290)
      write(11,2300)
c JAF error check (16 lines)
      if(nag.gt.maxnag-1)then
        write(*,'("error. too many alpha energy groups."/,
     1   "nag=",i0," maxnag=",i0)')nag,maxnag-1
        close(unit=1)
        close(unit=2)
        close(unit=3)
        close(unit=4)
        close(unit=5)
        close(unit=6, status='delete')
        close(unit=7, status='delete')
        close(unit=8, status='delete')
        close(unit=10, status='delete')
        close(unit=11)
        close(unit=12, status='delete')
        stop
      end if
      if (nt.eq.0) go to 1440

c----------------------------------------------------
c beginning of big loop over target nuclides
c         for nt=0 calculate spontaneous fission only
c----------------------------------------------------
      rewind 3
      if (id.eq.2) rewind 4
      write (6,100)
c     if(id.eq.2) write(6,102)

c--------------------------
c loop on target nuclides i
c--------------------------
  630 read (3,20) icont
      if (icont.ne.0) go to 630
  640 if (id.eq.2) read (4,20) icont
      if (icont.ne.0) go to 640
      do 1120 i=1,nt
      etpant=0.
      ebtant=0.
      totqan=0.
      if(id.eq.1) go to 660
      do 650 n=1,nng
  650 tsan(n)=0.

c------------------------------------------------------
c find alpha-n cross sections for this target on tape 3
c------------------------------------------------------
  660 read (1,*) idt,at
      write(11,2310)idt,at
  665 read (3,105,end=680) jdt,jps
  670 read (3,120) (e(ip),x(ip),ip=1,jps)
      if (idt-jdt) 680,690,665
  680 write (6,130) idt
      stop 2
  690 lzt=idt/10000
      lat=idt/10-lzt*1000
      atar=lat
      lst=idt-lzt*10000-lat*10
      amt=' '
      if (lst.ne.0) amt='m'
      apro=atar+3.
      do 700 ip=2,jps
      eamin=e(ip-1)
      if (x(ip).ne.0.) go to 710
  700 continue
  710 if (eamin.lt.0.001) eamin=0.001
      if (id.eq.1) go to 760

c----------------------------------------------------------------
c find product nuclide level branchings for this target on tape 4
c----------------------------------------------------------------
  720 read (4,140,end=750) jdt,jl,jp
  730 read (4,120) q,(el(il),il=1,jl)
      do 740 ip=1,jp
  740 read (4,120) ep(ip),(f(il,ip),il=1,jl)
      if (idt-jdt) 750,760,720
  750 write (6,150) idt
      stop 3
  760 rewind 5
      if (nq.ne.0) go to 770
      sbtqan=0.
  770 read (5,20) icont
      if (icont.ne.0) go to 770

c--------------------------------------------
c loop on source nuclides k for this target i
c contribution to (alpha,n) neutrons
c--------------------------------------------
      do 1100 k=1,nq
      etop=0.
      ebot=0.
      do 775 il=1,jl
      etopl(il)=0.
      ebotl(il)=0.
      do 774 n=1,nng
  774 sl(n,il)=0.
  775 totlev(il)=0.
      sbtqan=0.
      if(id.eq.1) go to 790
      do 780 n=1,nng
  780 san(n)=0.
  790 continue

c---------------------------------------------------------
c set alpha energy grid ee(m) for this target/source combo
c---------------------------------------------------------
      eamax=eal(nal)
  850 continue
      if(eamax.le.eamin) then
      qan=0.
      mm=0
      p(nagp1)=0.
      pval=0.
      go to 945
      endif
      if (eamax.le.6.5) go to 860
      write (6,170) eamax
      stop 5
  860 fnag=nag
      nagp1=nag+1
      dea=(eamax-eamin)/fnag
      ee(1)=eamin
      do 870 m=2,nag
      fm=m-1
  870 ee(m)=ee(1)+dea*fm
      ee(nagp1)=eamax

c----------------------------------------------------
c for each energy ee(m) calc. and store cross section
c----------------------------------------------------
      ip=1
      do 900 m=1,nagp1
  880 if (ee(m).ge.e(ip).and.ee(m).le.e(ip+1)) go to 890
      ip=ip+1
      if (ip.lt.jps) go to 880
      write (6,180) idt,m,ee(m)
      stop 6
  890 slope=(x(ip+1)-x(ip))/(e(ip+1)-e(ip))
      enrsep=(x(ip)*e(ip+1)-x(ip+1)*e(ip))/(e(ip+1)-e(ip))
  900 cx(m)=enrsep+ee(m)*slope

c-------------------------------------------------------------
c for each energy ee(m) calc. and store stopping cross section
c-------------------------------------------------------------
      do 930 m=1,nagp1
      scx(m)=0.
      do 920 j=1,nz
      zmat=jzm(j)
      amat=amass(jzm(j))
c-----nuclear stopping---
      term=zalp**0.66667 + zmat**0.66667
      sterm=sqrt(term)
      bot=zalp*zmat*(alph+amat)*sterm
      rep=32530.*amat*ee(m)/bot
      if(rep.lt.0.001) then
      dcx=1.593*sqrt(rep)
      go to 905
      endif
      if(rep.lt.10.) then
      bot=1.+6.8*rep+3.4*rep**1.5
      dcx=1.7*sqrt(rep)*alog(rep+rnapier)/bot
      go to 905
      endif
      dcx=alog(0.47*rep)/(2.*rep)
c-----electronic stopping---
  905 if(ee(m).gt.30.) go to 906
      slow=czm(1,j)*(1000.*ee(m))**czm(2,j)
      shigh=(czm(3,j)/ee(m))*alog(1.+czm(4,j)/ee(m) + czm(5,j)*ee(m))
      dcx=dcx+slow*shigh/(slow+shigh)
      go to 907
  906 eil=alog(1/ee(m))
      arg=czm(6,j)+czm(7,j)*eil+czm(8,j)*eil*eil+czm(9,j)*eil*eil*eil
      dcx=dcx+exp(arg)
  907 continue
  920 scx(m)=scx(m)+azm(j)*dcx
  930 continue

c-----------------------------------------------------------
c for each energy ee(m) calc. ratio r(m), then integral p(m)
c-----------------------------------------------------------
      r(1)=cx(1)/scx(1)
      p(1)=0.
      fact=1.0e-06*at
      do 940 m=2,nagp1
      r(m)=cx(m)/scx(m)
  940 p(m)=p(m-1)+fact*(r(m-1)+r(m))*dea/2.
  945 continue

c----------------------------------------------
c  calculate neutrons from alphas
c----------------------------------------------
  950 do 1080 l=1,nal
      if (nq.eq.0) go to 1000
      do 960 m=1,nag
      mm=m
  960 if (eal(l).ge.ee(m).and.eal(l).le.ee(m+1)) go to 970
      mm=0
      pval=0.
      go to 975
  970 pval=p(mm)+(eal(l)-ee(mm))*(p(mm+1)-p(mm))/(ee(mm+1)-ee(mm))
  975 aps=ast(l)
      qan=aps*pval
      sbtqan=sbtqan+qan
      totqan=totqan+qan
      if (l.eq.1) go to 980
      write (6,200) eal(l),aps,pval,qan
      if(nal.ne.1.and.l.eq.nal) write(6,250) sbtqan
      go to 990
  980 write(6,210)jsm(lzt),lat,amt,at,eal(l),aps,pval,qan
      if(nal.eq.1) write(6,255)
  990 if (id.eq.1.or.mm.eq.0) go to 1080

c-----------------------------------------------------------------
c calculate (alpha,n) neutron spectrum contribution in multigroups
c-----------------------------------------------------------------
 1000 mmm1=mm-1
      do 1010 m=1,mmm1
 1010 rr(m)=(p(m+1)-p(m))/pval
      rr(mm)=(pval-p(mm))/pval
      do 1020 n=1,nng
 1020 s(n)=0.
      do 1060 il=1,jl
      qlev=q-el(il)
      e90=-qlev*apro/(apro-alph)
      thre=-qlev*(aneut+apro)/(aneut+apro-alph)
      if (qlev.gt.0.) thre=0.
      do 1060 m=1,mm
      ea=(ee(m)+ee(m+1))/2.
      if (m.eq.mm) ea=(eal(l)+ee(mm))/2.
      if (ea.lt.thre) go to 1055
      do 1030 ip=2,jp
      lp=ip
 1030 if (ep(ip-1).le.ea.and.ep(ip).ge.ea) go to 1040
      write (6,220)ea
      stop 7
 1040 lp1=lp-1
      bx=f(il,lp1)+(ea-ep(lp1))*(f(il,lp)-f(il,lp1))/(ep(lp)-ep(lp1))
      term1=sqrt(4.*alph*aneut*ea)/(2.*(aneut+apro))
      term2=alph*aneut*ea/((aneut+apro)*(aneut+apro))
      term3=(apro*ea+apro*qlev-alph*ea)/(aneut+apro)
      senmax=term1+sqrt(term2+term3)
      enmax=senmax*senmax
      if (ea.le.e90) senmin=term1-sqrt(term2+term3)
      if (ea.gt.e90) senmin=-term1+sqrt(term2+term3)
      enmin=senmin*senmin
      ebar=(enmin+enmax)/2.
      val=qan*rr(m)*bx
      etop=etop+val*ebar
      etopl(il)=etopl(il)+val*ebar
      ebotl(il)=ebotl(il)+val
      ebot=ebot+val
      dele=enmax-enmin
      do 1050 n=1,nng
      if (en(n).lt.enmin) go to 1050
      if (en(n+1).gt.enmax) go to 1050
      de=en(n)-en(n+1)
      if (en(n+1).lt.enmin) de=de-(enmin-en(n+1))
      if (en(n).gt.enmax) de=de-(en(n)-enmax)
      gpadd=qan*rr(m)*bx*de/dele
      s(n)=s(n)+gpadd
      sl(n,il)=sl(n,il)+gpadd
      totlev(il)=totlev(il)+gpadd
 1050 continue
 1055 continue
 1060 continue
      do 1070 n=1,nng
      san(n)=san(n)+s(n)
      gtsan(n)=gtsan(n)+s(n)
      ts(n)=ts(n)+s(n)
 1070 tsan(n)=tsan(n)+s(n)
 1075 continue
 1080 continue
      if(id.eq.1) go to 1100

c----------------------------------------------------
c output this (alpha,n) neutron spectrum contribution
c for this source/target combination
c----------------------------------------------------
      ebarqt=0.
      if(ebot.le.0.) go to 1085
      ebarqt=etop/ebot
      etpant=etpant+etop
      ebtant=ebtant+ebot
 1085 continue
      smga=0.
      do 1090 n=1,nng
 1090 smga=smga+san(n)
c     fracgp=smga/sbtqan
c     do 1095 il=1,jl
c1095 frclev(il)=totlev(il)/sbtqan
c     write(6,205) (frclev(il),il=1,jl)
      if (nq.ne.0) write (7,240) jsm(lzt),lat,amt
      if (erg.ge.1) then
      write (7,80) (san(n),n=nng,1,-1)
      elseif (erg.le.-1) then
      write (7,80) (san(n),n=1,nng)
      else
      stop 'Error in erg input!'
      endif
      write (7,81) smga,ebarqt
      smga=0.
      if (nq.ne.0) write (8,240) jsm(lzt),lat,amt
      do 1097 n=1,nng
      san(n)=san(n)/sbtqan
 1097 smga=smga+san(n)
      if (erg.ge.1) then
      write (8,80) (san(n),n=nng,1,-1)
      elseif (erg.le.-1) then
      write (8,80) (san(n),n=1,nng)
      endif
      write (8,81) smga, ebarqt
      do 1098 il=1,jl
      if(ebotl(il).le.0.) go to 1098
c     ebarl=etopl(il)/ebotl(il)
 1098 continue
 1100 continue
 1101 gttqan=gttqan+totqan

c--------------------------------------
c output total (alpha,n) source for all
c alpha sources on this target nuclide
c--------------------------------------
      if(id.eq.1) go to 1120
      tmga=0.
      do 1110 n=1,nng
 1110 tmga=tmga+tsan(n)
c     fracgp=tmga/totqan
      ebart=etpant/ebtant
      etopan=etopan+etpant
      ebotan=ebotan+ebtant
      if(nq.le.1) go to 1120
      write(7,260)
      if (erg.ge.1) then
      write(7,80) (tsan(n),n=nng,1,-1)
      elseif (erg.le.-1) then
      write(7,80) (tsan(n),n=1,nng)
      endif
      write(7,81) tmga,ebart
      write(7,2090)
      write(7,2021)
      write(7,2090)
      tmga=0.
      do 1115 n=1,nng
       tsan(n)=tsan(n)/totqan
 1115 tmga=tmga+tsan(n)
      write(8,260)
      if (erg.ge.1) then
      write(8,80) (tsan(n),n=nng,1,-1)
      elseif (erg.le.-1) then
      write(8,80) (tsan(n),n=1,nng)
      endif
      write(8,81) tmga,ebart
      write(8,2090)
      write(8,2021)
      write(8,2090)
      write(10,2090)
      write(10,2021)
      write(10,2090)
 1120 continue

c----------------------------------------------------------------------
c output grand total (alpha,n) neutron source: all targets, all sources
c----------------------------------------------------------------------
      if (nt.gt.1) write (6,252) gttqan
      if(id.eq.1) go to 1440
      gtmga=0.
      do 1130 n=1,nng
 1130 gtmga=gtmga+gtsan(n)
      ebaran=etopan/ebotan
c     fracgp=gtmga/gttqan
      write(7,2090)
      write(7,2021)
      write(7,390)
      if (erg.ge.1) then
      write (7,80) (gtsan(n),n=nng,1,-1)
      elseif (erg.le.-1) then
      write (7,80) (gtsan(n),n=1,nng)
      endif
      write (7,81) gtmga,ebaran
      gtmga=0.
        dummy=0.0
      do 1135 n=1,nng
         dummy=dummy+gtsan(n)
       gtsan(n)=gtsan(n)/gttqan
 1135 gtmga=gtmga+gtsan(n)
      write(8,2090)
      write(8,2021)
      write(8,390)
      if (erg.ge.1) then
      write(8,80) (gtsan(n),n=nng,1,-1)
      elseif (erg.le.-1) then
      write(8,80) (gtsan(n),n=1,nng)
      endif
      write(8,81) gtmga,ebaran

c---------------------------------------------------
c spectra for outp2 in ascending or descending order
c---------------------------------------------------
      if (erg.ge.1) then
       write(12,1900)
       write(12,3600)title
       write(12,3601)
       write(12,3605)
       write(12,3610)
       write(12,3615)
       write(12,3620)
       write(12,3630)(en(n),gtsan(n),n=nng,1,-1)
       write(12,3635)gtmga
       write(12,3611)
       write(12,3610)
       write(12,3615)
       write(12,3620)
       write(12,3630)(en(n),gtsan(n)*gttqan,n=nng,1,-1)  
       write(12,3635)gttqan    
      elseif (erg.le.-1) then
       write(12,1900)
       write(12,3600)title
       write(12,3601)
       write(12,3605)
       write(12,3610)
       write(12,3615)
       write(12,3620)
       write(12,3630)(en(n),gtsan(n),n=1,nng)
       write(12,3635)gtmga
       write(12,3611)
       write(12,3610)
       write(12,3615)
       write(12,3620)
       write(12,3630)(en(n),gtsan(n)*gttqan ,n=1,nng)
       write(12,3635)gttqan 
      endif

c-----------------
c  Output Summary
c-----------------
 1440 continue
      write(11,2090)
      write(11,2090)
      write(11,3000)
      write(11,3010)
      write(11,2090)
      write(11,3020)gttqan
      if (id.gt.1) then
       write(11,3030)ebaran
         write(11,3050)
         do 1450 n=1,nng
           write(11,3051)n,gtsan(n)
 1450  continue
       frac=100.0*(1.0-(abs(dummy-gttqan)/gttqan))
         write(11,3500)frac
         if (frac.le.85.0) write(11,3510)
      endif

c-------------------------
c  End Execution of SNINF
c-------------------------
      return
      end


c=======================================================================
c  Three Region Problem Main Subroutine (1/99)
c=======================================================================

      subroutine three(title,idd,id,erg)

c----------
c  Storage
c----------
      
      character title*7
      integer erg
c JAF introduce parameters
      integer maxnag
      parameter (maxnag=2001)
c JAF declare ee,cg,eal
      dimension ee(maxnag), cg(maxnag), eal(maxnag)
      dimension title(11)

c------------------
c  Call Statements
c------------------

c JAF pass in maxnag,ee,cg,eal
      call input(maxnag,title,idd,id,erg,ee,cg,eal)
      call processor(maxnag,title,idd,id,erg,ee,cg,eal)

c----------------------------
c  End of Execution of THREE
c----------------------------

      return
      end


c=======================================================================
c  Three Region Problem Input Subroutine (1/99)
c=======================================================================

      subroutine input(maxnag,title,idd,id,erg,ee,cg,eal)

c----------
c  Storage
c----------
      character title*7, title1*7, title2*7, title3*7, jsm(105)*2,
     1 amt(20)*1, amtc(20)*1
      integer erg
c JAF introduce parameters
      integer maxnag
c JAF replaced all 4001's in this list with maxnag.
      dimension title(11), title1(11), title2(11), title3(11), 
     1 ee(maxnag), ea(maxnag), en(751), cg(maxnag), jza(20), aza(20),
     2 jq(300), faq(300), c(14), cza(9,20), fal(30,20), eala(30,20), 
     3 fdng(1000), jzb(20), azb(20), czb(9,20), idb(20), atb(20),
     4 jps(20), e(1100,20), x(1100,20), lzt(20), lat(20), lst(20),
     5 apro(20), jl(20), jp(20), q(20), el(20,20), ep(200,20),
     6 f(20,200,20), jzc(20), azc(20), czc(9,20), idc(20), atc(20),
     7 jpsc(20), ec(1100,20), xc(1100,20), lztc(20), latc(20), 
     8 lstc(20), aproc(20), jlc(20), jpc(20), qc(20), elc(20,20), 
     9 epc(200,20), fc(20,200,20), nal(20), alam(20), eal(maxnag)

       common /names/ jsm
       common /masses/ amass(105)
       common /vars/ jza, aza, jq, faq, cza, fal, eala, fdng, jzb, 
     1 azb, czb, idb, atb, jps, e, x, lzt, lat, lst, apro, jl, jp, 
     2 q, el, ep, f, jzc, azc, czc, idc, atc, jpsc, ec, xc, lztc, 
     3 latc, lstc, aproc, jlc, jpc, qc, elc, epc, fc, nza, isga, nq,
     5 nzb, isgb, anumb, thick, ntb, nzc, isgc, ntc, amt, amtc, nal,
     6 alam, title1, title2, title3
c JAF remove ee, cg, eal from common
c      common /grids/ ee, en, cg, nag, eamax, eamin, nng, 
c    1 enmax, enmin, ncg, dea, eal
       common /grids/ en, nag, eamax, eamin, nng, 
     1 enmax, enmin, ncg, dea

c--------------------
c  Format Statements
c--------------------
   10 format (i3,11a7)
   15 format(11a7)
   16 format(11a7,//)
   20 format (2i3,f7.3)
   30 format (i8,e12.5)
   40 format(i3,8e8.3,f10.3,/,3x,5f8.3)
   50 format (69h failed to find stopping cross section coefficients on
     1tape 2 for iz=,i3,5h stop)
   60 format (6e12.5)
   70 format(/,1x,'Neutron Multigroup Structure (MeV)')
   80 format (1p8e10.3)
   90 format (i3,5e12.5)
   95 format(2i6)
   96 format(/,26x,'Neutron Source Magnitudes',/,26x,25(1h~),/)
   97 format(/,23x,'Absolute Neutron Source Spectra',/,23x,31(1h~),/)
   98 format(/,22x,'Normalized Neutron Source Spectra',/,22x,
     133(1h~),/)
   99 format(/,16x,'Neutron Source Spectra by Nuclide Energy Level',
     1/,16x,46(1h~),/)
  100 format(//,1h1,30x,'Table I',/,31x,7(1h=),//,12x,'(alpha,n) ',
     1'Neutron Production by Target and Alpha Energy',///,17x,
     2'target            alpha   alphas/sec     p(e)     neuts/sec',
     3/,7x,'target  atom frac.          energy     /cm^2',4x,
     4'neut/alpha',4x,'/cm^2',/,1h+,6x,6(1h_),2x,10(1h_),2x,6(1h_),2x,
     5 6(1h_),3(12h  __________))
  105 format(i8,i4)
  110 format (i8,2i4)
  120 format (8e10.3)
  130 format (1x,'SOURCE DOES NOT HAVE ALPHA EMITTER ', i6)
  131 format(' WARNING-AN ALPHA ENERGY THAT EXCEEDED ',f7.3,' MEV WAS',
     1' SET TO:',f7.3,' MeV')
  132 format(' NUMBER OF STOPPING ELEMENTS ON SOURCE SIDE:',i3,/,
     1' SOLID OR GAS STOPPING POWER TRIGER:',i3,/,
     2' MAXIMUM ENERGY FOR ALPHA SPECTRA:',f7.3,/)
  133 format(' ISOTOPE IDENTIFICATION:',i6,/
     1' ISOTOPE FRACTION:',e15.6,/)
  134 format(' NUMBER OF SOURCES:',i3)
  135 format(/,' SOURCE IDENTIFICATION: ',i6,/
     1' SOURCE ISOTOPIC FRACTION:',e15.6)
  136  format(//,' ALPHA SPECTRA ENERGY BOUNDARIES',/)
  137  format(5e14.6)
  138 format(' MAXIMUM ALPHA DECAY ENERGIES')
  140 format (i8,2i3)
  150 format (66h failed to find comp.nucl.neut.emiss.lev. branching on 
     1tape 4 for ,i8,5h stop)
  160 format (56h failed to find alpha/s.f./dn source data on tape 5 for
     1 ,i8,5h stop)

 1900 format('SOURCES 4D Calculation',/,'<<<<<<<<<<<>>>>>>>>>>>',//)
 2000 format(/,22x,'Summary of Input')
 2010 format(22x,'================')
 2020 format(' ')
 2021 format(/,80(1h-))
 2030 format('Title: ', 11a7,/)
 2043 format('Three region problem input (idd =',i2,')')
 2041 format('Ascending energy structure for output (erg =',i3,')')
 2042 format('Descending energy structure for output (erg =',i3,')')
 2048 format('Magnitudes and spectra computed.')
 2049 format('Magnitudes only computed.')
 2072 format('Maximum energy for alpha spectra: ',1pE10.3,' MeV.')
 2074 format('Minimum energy for alpha spectra: ',1pE10.3,' MeV.')
 2080 format('Elemental Constituents in Region A:')
 2081 format('Elemental Constituents in Region B:')
 2082 format('Elemental Constituents in Region C:')
 2090 format(' ')
 2100 format('  Z-value  Atom Fraction')
 2110 format('  -------  -------------')
 2120 format(4x,i3,5x,f12.10)
 2130 format('Number of neutron spectrum energy groups:', i4)
 2140 format('Maximum neutron energy is ',1pE12.5,' MeV.')
 2150 format('Minimum neutron energy is ',1pE12.5,' MeV.')
 2160 format('Neutron Energy Group Structure:')
 2170 format('  Group    Upper-Bound    Lower-Bound')
 2180 format('  -----    -----------    -----------')
 2190 format(2x,i4,4x,1pE12.5,3x,1pE12.5)
 2193 format(/,'Warning: First Group Upper-Bound Does Not Coincide ',
     1'With Maximum Neutron Energy.')
 2195 format(/,'Warning: Energy Group Structure Not Entered in ',
     1'Decreasing Order.',/,5x,'Group ',i3,' Upper-Bound is Greater ',
     2'Than Group ',i3,' Upper-Bound.') 
 2200 format('Number of source nuclides to be evaluated:',i3)
 2220 format('Source Nuclides in Region A:')
 2230 format('    ZAID     Atom Fraction')
 2240 format('    ----     -------------')
 2250 format(3x,i6,4x,1pE10.3)
 2260 format(/,'Number of target nuclides in region B:',i3)
 2261 format(/,'Number of target nuclides in region C:',i3)
 2270 format(i5,' alpha energy groups used at each interface.')
 2280 format('Target Nuclides in Region B:')
 2281 format('Target Nuclides in Region C:')

 3030 format(/,'Region A Title: ', 11a7,/)
 3040 format(/,'Region B Title: ', 11a7,/)
 3050 format(/,'Region C Title: ', 11a7,/)
 3130 format('Number of angular groups:', i4)
 3160 format('Angular Group Structure:',/)
 3170 format('  Group    Upper-Bound    Lower-Bound')
 3180 format('  -----    -----------    -----------')
 3190 format(2x,i4,4x,1pE12.5,3x,1pE12.5)
 3200 format('Number of elemental constituents in region A:',i3)
 3201 format('Number of elemental constituents in region B:',i3)
 3202 format('Number of elemental constituents in region C:',i3)
 3210 format('Solid stopping cross-sections used (isga=',i2,') in',
     1 1x,'region A.')
 3211 format('Solid stopping cross-sections used (isgb=',i2,') in',
     1 1x,'region B.')
 3212 format('Solid stopping cross-sections used (isgc=',i2,') in',
     1 1x,'region C.')
 3220 format('Gas stopping cross-sections used (isga=',i2,')',
     1 1x,'in region A.')
 3221 format('Gas stopping cross-sections used (isgb=',i2,')',
     1 1x,'in region B.')
 3222 format('Gas stopping cross-sections used (isgc=',i2,')',
     1 1x,'in region C.')
 3231 format('Material B atom density: ',1pE12.5,' atoms/b-cm.')
 3240 format('Interface region thickness: ',1pE12.5,' cm.')

c-----------------
c  Set parameters
c-----------------
c JAF replaced pi=4.0*atan(1.0)
      pi=3.1415926535898d+0

c----------------
c  Write headers
c----------------
      write(6,1900)
      write(7,1900)
      write(8,1900)
      write(10,1900)
      write(11,1900)
      write(11,2000)
      write(11,2010)
      write(11,2020)
      write(11,2030)title
      write(11,2043)idd
      if (id.eq.1) then
       write(11,2049)
      elseif (id.eq.2) then
       write(11,2048)
       write(7,97)
       write(8,98)
      else
       stop 'Error in id input!'
      endif
      write(6,96)
      write(6,2030)title
      write(7,2030)title
      write(8,2030)title

c--------------------------------------
c  Create alpha energy group structure
c--------------------------------------
      read (1,*)nag,eamax,eamin
      write(11,2090)
      write(11,2270)nag
      write(11,2072)eamax
      write(11,2074)eamin
c JAF error check (16 lines)
      if(nag.gt.maxnag-1)then
        write(*,'("error. too many alpha energy groups."/,
     1   "nag=",i0," maxnag=",i0)')nag,maxnag-1
        close(unit=1)
        close(unit=2)
        close(unit=3)
        close(unit=4)
        close(unit=5)
        close(unit=6, status='delete')
        close(unit=7, status='delete')
        close(unit=8, status='delete')
        close(unit=10, status='delete')
        close(unit=11)
        close(unit=12, status='delete')
        stop
      end if
      if (eamin.lt.0.001) eamin=0.001
  200 fnag=nag
      nagp1=nag+1
      dea=(eamax-eamin)/fnag
      ee(1)=eamin
      do 210 m=2,nag
       fm=m-1
  210 ee(m)=ee(1)+dea*fm
      ee(nagp1)=eamax
      do 211 m=1,nag
  211 ea(m)=(ee(m+1) + ee(m))/2.

c----------------------------------------
c  Create neutron energy group structure
c----------------------------------------
      read (1,*) nng,enmax,enmin
      write(11,2090)
      write(11,2130)nng
      write(11,2140)enmax
      write(11,2150)enmin
      write(11,2090)
      write(11,2160)
      write(11,2090)
      write(11,2170)
      write(11,2180)
      nngp1=nng+1
      if (nng.gt.0) go to 240
      nng=-nng
      nngp1=nng+1
      read (1,*) (en(n),n=1,nng)
      en(nngp1)=enmin
      if (en(1).gt.en(2)) go to 260
      nh=nngp1/2
      do 231 n=1,nh
       etem=en(nngp1-n+1)
       en(nngp1-n+1)=en(n)
  231 en(n)=etem
      go to 260
  240 en(nngp1)=enmin
      fnng=nng
      den=(enmax-enmin)/fnng
      en(1)=enmax
      do 250 n=2,nng
       fnm1=n-1
  250 en(n)=en(1)-fnm1*den
  260 do 262 n=1,nng
       write(11,2190)n,en(n),en(n+1)
  262 continue
      if (en(1).ne.enmax) write(11,2193)
      do 268 n=2,nng
       nm1=n-1
       if (en(n).gt.en(n-1)) then
        write(11,2195)n,nm1
        stop 'Energy Structure Incorrectly Entered!'
       endif
  268 continue
      write (7,70)
        if (erg.ge.1) then
      write(7,80) (en(n),n=nngp1,1,-1)
      elseif (erg.le.-1) then
      write(7,80) (en(n),n=1,nngp1)
      endif
      if (erg.ge.1) then
      write(8,80) (en(n),n=nngp1,1,-1)
      elseif (erg.le.-1) then
      write(8,80) (en(n),n=1,nngp1)
      endif

c---------------------------------
c  Create angular group structure
c---------------------------------
      read (1,*)ncg
      write(11,2090)
      write(11,3130)ncg
      write(11,2090)
      write(11,3160)
      write(11,3170)
      write(11,3180)
      fncg=ncg
      cg(1)=0
      do 270 i=1,ncg
        cg(i+1)=cg(i)+(pi/2.0/fncg)
        write(11,3190)i,cg(i+1),cg(i)
  270 continue

c-------------------------------------------
c  Read elemental constituents for region A
c-------------------------------------------
      read(1,15)title1
      write(11,2090)
      write(11,3030)title1
      read(1,*)nza,isga
      write(11,3200)nza
c JAF error check (16 lines)
      if(nza.gt.20)then
        write(*,'("error. too many elements."/,
     1   "nza=",i0," max=",i0)')nza,20
        close(unit=1)
        close(unit=2)
        close(unit=3)
        close(unit=4)
        close(unit=5)
        close(unit=6, status='delete')
        close(unit=7, status='delete')
        close(unit=8, status='delete')
        close(unit=10, status='delete')
        close(unit=11)
        close(unit=12, status='delete')
        stop
      end if
      if (isga.eq.0) then
       write(11,3210)isga
      elseif (isga.eq.1) then
       write(11,3220)isga
      else
       stop 'Error in isga input!'
      endif
      write(11,2080)
      write(11,2090)
      write(11,2100)
      write(11,2110)
      do 280 j=1,nza
       read (1,*) jza(j),aza(j)
       write(11,2120)jza(j),aza(j)
  280 continue

c----------------------------------------
c order the nz material constituents by z
c----------------------------------------
      do 300 j=1,nza
       isw=0
       do 290 jn=2,nza
        if (jza(jn).gt.jza(jn-1)) go to 290
        jtem=jza(jn-1)
        atem=aza(jn-1)
        jza(jn-1)=jza(jn)
        aza(jn-1)=aza(jn)
        jza(jn)=jtem
        aza(jn)=atem
        isw=1
  290  continue
       if (isw.eq.0) go to 310
  300 continue

c----------------------------------------------------
c read tape 2 for stopping cross-section coefficients
c----------------------------------------------------
  310 rewind 2
  320 read (2,10) icont
      if (icont.ne.0) go to 320
      do 390 j=1,nza
  330  read (2,40,end=380) iz,c
  340  if (jza(j)-iz) 380,350,330
  350  continue
       do 360 jt=1,9
  360  cza(jt,j)=c(jt)
       if (isga.eq.1.and.c(10).gt.0.) then
       do 370 jt=1,5
  370  cza(jt,j)=c(jt+9)
       endif
       go to 390
  380  write (6,50) jza(j)
       stop 1
  390 continue

c-------------------------------------------
c read alpha sources from user input
c---------------------------------------------
      read (1,*) nq
      write(11,2090)
      write(11,2200)nq
      write(11,2090)
      write(11,2220)
      write(11,2090)
      write(11,2230)
      write(11,2240)
      do 400 k=1,nq
       read (1,*) jq(k),faq(k)
       write(11,2250)jq(k),faq(k)
  400 continue

c-----------------------------------------------
c order sources by z-a-state id = s+10*a+10000*z
c-----------------------------------------------
      if (nq.eq.1) go to 430
      do 420 k=1,nq
       isw=0
       do 410 kn=2,nq
        if (jq(kn).gt.jq(kn-1)) go to 410
        jtem=jq(kn-1)
        atem=faq(kn-1)
        jq(kn-1)=jq(kn)
        faq(kn-1)=faq(kn)
        jq(kn)=jtem
        faq(kn)=atem
        isw=1
  410  continue
       if (isw.eq.0) go to 430
  420 continue

c-------------------------------------------
c read source nuclide decay data from tape 5
c-------------------------------------------
  430 continue
  440 read (5,10) icont
      if (icont.ne.0) go to 440
      do 510 k=1,nq
  450  read (5,110,end=480) idq,nal(k),ndn
       read (5,60) alam(k),bfsf,barnu,a,b,bfdn
       if (nal(k).eq.0) go to 460
       read (5,120) (eala(l,k),fal(l,k),l=1,nal(k))
  460  if (ndn.eq.0) go to 470
       read (5,120) (fdng(nd),nd=1,ndn)
  470  if (idq-jq(k)) 450,490,480
  480  write (6,160) jq(k)
       stop 'Source nuclide not in library'
  490  continue
       if(nal(k).eq.0)write(6,130)
       if(nal(k).eq.0)go to 500
  500  continue
  510 continue

c------------------------------------------
c  Read material constituents for region B
c------------------------------------------
      read (1,15)title2
      write(11,2090)
      write(11,3040)title2
      read (1,*) nzb,isgb,anumb,thick
      write(11,3201)nzb
c JAF error check (16 lines)
      if(nzb.gt.20)then
        write(*,'("error. too many elements."/,
     1   "nzb=",i0," max=",i0)')nzb,20
        close(unit=1)
        close(unit=2)
        close(unit=3)
        close(unit=4)
        close(unit=5)
        close(unit=6, status='delete')
        close(unit=7, status='delete')
        close(unit=8, status='delete')
        close(unit=10, status='delete')
        close(unit=11)
        close(unit=12, status='delete')
        stop
      end if
      if (isgb.eq.0) then
      write(11,3211)isgb
      elseif (isgb.eq.1) then
       write(11,3221)isgb
      else
       stop 'Error in isgb input!'
      endif
      write(11,3231)anumb
      write(11,3240)thick
      write(11,2090)
      write(11,2081)
      write(11,2090)
      write(11,2100)
      write(11,2110)

c---------------------------------
c read input material constituents
c---------------------------------
      do 520 j=1,nzb
       read (1,*) jzb(j),azb(j)
  520 write(11,2120)jzb(j),azb(j)

c----------------------------------------
c order the nz material constituents by z
c----------------------------------------
      do 540 j=1,nzb
       isw=0
       do 530 jn=2,nzb
        if (jzb(jn).gt.jzb(jn-1)) go to 530
        jtem=jzb(jn-1)
        atem=azb(jn-1)
        jzb(jn-1)=jzb(jn)
        azb(jn-1)=azb(jn)
        jzb(jn)=jtem
        azb(jn)=atem
        isw=1
  530  continue
       if (isw.eq.0) go to 550
  540 continue

c----------------------------------------------------
c read tape 2 for stopping cross-section coefficients
c----------------------------------------------------
  550 rewind 2
  560 read (2,20) icont
      if (icont.ne.0) go to 560
      do 620 j=1,nzb
  570  read (2,40,end=610) iz,c
  580  if (jzb(j)-iz) 610,590,570
  590  continue
       do 600 jt=1,9
  600  czb(jt,j)=c(jt)
       if (isgb.eq.1.and.c(10).gt.0.) then
        do 605 jt=1,5
  605   czb(jt,j)=c(jt+9)
       endif
       go to 620
  610  write (6,50) jzb(j)
       stop 1
  620 continue
  625 if (id.eq.1) go to 720

c----------------------------------------
c read target information from user input
c----------------------------------------
  720 continue
      read (1,*) ntb
      write(11,2260)ntb
      write(11,2090)
      write(11,2280)
      write(11,2090)
      write(11,2230)
      write(11,2240)

c--------------------------
c loop on target nuclides i
c--------------------------
      rewind 3
      if (id.eq.2) rewind 4
      write (6,100)
  730 read (3,20) icont
      if (icont.ne.0) go to 730
  740 if (id.eq.2) read (4,20) icont
      if (icont.ne.0) go to 740
      do 880 i=1,ntb

c------------------------------------------------------
c find alpha-n cross sections for this target on tape 3
c------------------------------------------------------
  760  read (1,*) idb(i),atb(i)
       write(11,2250)idb(i),atb(i)
  765  read (3,105,end=780) jdt,jps(i)
  770  read (3,120) (e(ip,i),x(ip,i),ip=1,jps(i))
       if (idb(i)-jdt) 780,790,765
  780  write (6,130) idb(i)
       stop 2
  790  lzt(i)=idb(i)/10000
       lat(i)=idb(i)/10-lzt(i)*1000
       atar=lat(i)
       lst(i)=idb(i)-lzt(i)*10000-lat(i)*10
       amt(i)=' '
       if (lst(i).ne.0) amt(i)='m'
       apro(i)=atar+3.
c       do 800 ip=2,jps(i)
c        eamin=e(ip-1,i)
c        if (x(ip,i).ne.0.) go to 810
c  800  continue
  810  if (eamin.lt.0.001) eamin=0.001
       if (id.eq.1) go to 860

c----------------------------------------------------------------
c find product nuclide level branchings for this target on tape 4
c----------------------------------------------------------------
  820  read (4,140,end=850) jdt,jl(i),jp(i)
  830  read (4,120) q(i),(el(il,i),il=1,jl(i))
       do 840 ip=1,jp(i)
  840  read (4,120) ep(ip,i),(f(il,ip,i),il=1,jl(i))
       if (idb(i)-jdt) 850,860,820
  850  write (6,150) idb(i)
       stop 3
  860  rewind 5
       if (nq.ne.0) go to 870
       sbtqan=0.
  870  read (5,20) icont
       if (icont.ne.0) go to 870
  880 continue

c------------------------------------------
c  Read material constituents for region C
c------------------------------------------
      read (1,15)title3
      write(11,2090)
      write(11,3050)title3
      read (1,*) nzc,isgc
      write(11,3202)nzc
c JAF error check (16 lines)
      if(nzc.gt.20)then
        write(*,'("error. too many elements."/,
     1   "nzc=",i0," max=",i0)')nzc,20
        close(unit=1)
        close(unit=2)
        close(unit=3)
        close(unit=4)
        close(unit=5)
        close(unit=6, status='delete')
        close(unit=7, status='delete')
        close(unit=8, status='delete')
        close(unit=10, status='delete')
        close(unit=11)
        close(unit=12, status='delete')
        stop
      end if
      if (isgc.eq.0) then
      write(11,3212)isgc
      elseif (isgc.eq.1) then
       write(11,3222)isgc
      else
       stop 'Error in isgc input!'
      endif
      write(11,2090)
      write(11,2082)
      write(11,2090)
      write(11,2100)
      write(11,2110)

c---------------------------------
c read input material constituents
c---------------------------------
      do 1520 j=1,nzc
       read (1,*) jzc(j),azc(j)
 1520 write(11,2120)jzc(j),azc(j)

c----------------------------------------
c order the nz material constituents by z
c----------------------------------------
      do 1540 j=1,nzc
       isw=0
       do 1530 jn=2,nzc
        if (jzc(jn).gt.jzc(jn-1)) go to 1530
        jtem=jzc(jn-1)
        atem=azc(jn-1)
        jzc(jn-1)=jzc(jn)
        azc(jn-1)=azc(jn)
        jzc(jn)=jtem
        azc(jn)=atem
        isw=1
 1530  continue
      if (isw.eq.0) go to 1550
 1540 continue

c----------------------------------------------------
c read tape 2 for stopping cross-section coefficients
c----------------------------------------------------
 1550 rewind 2
 1560 read (2,20) icont
      if (icont.ne.0) go to 1560
      do 1620 j=1,nzc
 1570  read (2,40,end=1610) iz,c
 1580  if (jzc(j)-iz) 1610,1590,1570
 1590  continue
       do 1600 jt=1,9
 1600  czc(jt,j)=c(jt)
       if (isgc.eq.1.and.c(10).gt.0.) then
        do 1605 jt=1,5
 1605   czc(jt,j)=c(jt+9)
       endif
       go to 1620
 1610  write (6,50) jzc(j)
       stop 1
 1620 continue
 1625 if (id.eq.1) go to 1720

c----------------------------------------
c read target information from user input
c----------------------------------------
 1720 continue
      read (1,*) ntc
      write(11,2261)ntc
      write(11,2090)
      write(11,2281)
      write(11,2090)
      write(11,2230)
      write(11,2240)

c--------------------------
c loop on target nuclides i
c--------------------------
      rewind 3
      if (id.eq.2) rewind 4
 1730 read (3,20) icont
      if (icont.ne.0) go to 1730
 1740 if (id.eq.2) read (4,20) icont
      if (icont.ne.0) go to 1740
      do 1880 i=1,ntc

c------------------------------------------------------
c find alpha-n cross sections for this target on tape 3
c------------------------------------------------------
 1760  read (1,*) idc(i),atc(i)
       write(11,2250)idc(i),atc(i)
 1765  read (3,105,end=1780) jdt,jpsc(i)
 1770  read (3,120) (ec(ip,i),xc(ip,i),ip=1,jpsc(i))
       if (idc(i)-jdt) 1780,1790,1765
 1780  write (6,130) idc(i)
       stop 2
 1790  lztc(i)=idc(i)/10000
       latc(i)=idc(i)/10-lztc(i)*1000
       atar=latc(i)
       lstc(i)=idc(i)-lztc(i)*10000-latc(i)*10
       amtc(i)=' '
       if (lstc(i).ne.0) amtc(i)='m'
       aproc(i)=atar+3.
c       do 1800 ip=2,jpsc(i)
c        eamin=ec(ip-1,i)
c        if (xc(ip,i).ne.0.) go to 1810
c 1800  continue
 1810  if (eamin.lt.0.001) eamin=0.001
       if (id.eq.1) go to 1860

c----------------------------------------------------------------
c find product nuclide level branchings for this target on tape 4
c----------------------------------------------------------------
 1820  read (4,140,end=1850) jdt,jlc(i),jpc(i)
 1830  read (4,120) qc(i),(elc(il,i),il=1,jlc(i))
       do 1840 ip=1,jpc(i)
 1840  read (4,120) epc(ip,i),(fc(il,ip,i),il=1,jlc(i))
       if (idc(i)-jdt) 1850,1860,1820
 1850  write (6,150) idc(i)
       stop 3
 1860  rewind 5
       if (nq.ne.0) go to 1870
       sbtqan=0.
 1870  read (5,20) icont
       if (icont.ne.0) go to 1870
 1880 continue

c-------------------------
c  End Execution of INPUT
c-------------------------
      return
      end


c=======================================================================
c  Three Region Problem Processor Subroutine (1/99)
c=======================================================================

      subroutine processor(maxnag,title,idd,id,erg,ee,cg,eal)

c----------
c  Storage
c----------
      character title*7, title1*7, title2*7, title3*7, jsm(105)*2,
     1 amt(20)*1, amtc(20)*1
      integer erg
c JAF introduce parameters
      integer maxnag
c JAF replaced all 4001's in this list with maxnag.
      dimension title(11), r(maxnag), ps(maxnag), p(maxnag), nal(20),
     1 ee(maxnag), en(751), cg(maxnag), jza(20), aza(20),
     2 jq(300), faq(300), cza(9,20), fal(30,20), eala(30,20), 
     3 fdng(1000), jzb(20), azb(20), czb(9,20), idb(20), atb(20),
     4 jps(20), e(1100,20), x(1100,20), lzt(20), lat(20), lst(20),
     5 apro(20), jl(20), jp(20), q(20), el(20,20), ep(200,20),
     6 f(20,200,20), jzc(20), azc(20), czc(9,20), idc(20), atc(20),
     7 jpsc(20), ec(1100,20), xc(1100,20), lztc(20), latc(20), 
     8 lstc(20), aproc(20), jlc(20), jpc(20), qc(20), elc(20,20), 
     9 epc(200,20), fc(20,200,20), scxa(maxnag), scxb(maxnag), 
     1 astab(maxnag), eal(maxnag), alam(20), D(maxnag), 
     2 itrans(maxnag,maxnag), astbc(maxnag), title1(11), title2(11), 
     3 title3(11)
      real nstabb(750), nstbcb(750), nstbcc(750), nst(750), 
     1 nstnorm(750)

       common /names/ jsm
       common /masses/ amass(105)
       common /vars/ jza, aza, jq, faq, cza, fal, eala, fdng, jzb, 
     1 azb, czb, idb, atb, jps, e, x, lzt, lat, lst, apro, jl, jp, 
     2 q, el, ep, f, jzc, azc, czc, idc, atc, jpsc, ec, xc, lztc, 
     3 latc, lstc, aproc, jlc, jpc, qc, elc, epc, fc, nza, isga, nq,
     5 nzb, isgb, anumb, thick, ntb, nzc, isgc, ntc, amt, amtc, nal,
     6 alam, title1, title2, title3
c JAF remove ee, cg, eal from common
c      common /grids/ ee, en, cg, nag, eamax, eamin, nng, 
c    1 enmax, enmin, ncg, dea, eal
       common /grids/ en, nag, eamax, eamin, nng, 
     1 enmax, enmin, ncg, dea

c--------------------
c  Format Statements
c--------------------
   16 format(//,10(8h========),//,'Title: ',11a7)
   80 format (1p8e10.3)
   81 format(31x,'Total (all groups): ',1pE10.3,' neutrons/sec-cm^2.',
     1/,27x,'Average Neutron Energy: ',1pE10.3,' MeV.')
  130 format (1x,'SOURCE DOES NOT HAVE ALPHA EMITTER ', i6)
  131 format(' WARNING-AN ALPHA ENERGY THAT EXCEEDED ',f7.3,' MEV WAS',
     1' SET TO:',f7.3,' MeV')
  390 format (////,' Total Neutron Spectrum')
 1900 format('SOURCES 4D Calculation',/,'<<<<<<<<<<<>>>>>>>>>>>',//)
 2021 format(/,80(1h-))
 3000 format(//,22x,'Summary of Output')
 3010 format(22x,'=================')
 3020 format('Total (alpha,n) neutron source from all sources and ',
     1'targets: ',1pE12.5,' n/sec-cm^2.')
 3030 format(/,'Average (alpha,n) neutron energy: ',1pE12.5,' MeV.')
 3050 format(/,'Normalized Neutron Energy Spectrum by Energy Group ',
     1'for All Sources:',/'(Note: group structure is independent of',
     2' erg record!)',//,10x,'Group',6x,'Contribution',/,10x,
     3'-----',6x,'------------')
 3051 format(11x,i3,7x,1pE12.5)
 3500 format(/,'Portion of Total Neutron Source Rate Accounted for ',
     1'in the Energy Spectrum: ',F5.1,'%.')
 3510 format(/,'Warning: Energy Group Structure May Be Poorly Chosen.')
 3600 format(' Title: ', 11a7)
 3601 format(' Neutron Source Spectra in Columnar Format.'/
     1' Examine other tape files for additional information.'/
     2' Recall the Energy values are bin bounds.'/) 
 3605 format(/,'           Normalized Total (neutrons/sec-cm^2)')
 3610 format('           ========================================')
 3611 format(/,'            Absolute Total (neutrons/sec-cm^2)')
 3615 format(' Energy     (a,n)')
 3620 format(1x,9('-'),2x,9('-'))
 3625 format(1p5e10.3)
 3630 format(1p1e10.3,1p1e11.3)
 3635 format(' Total    ',1p1e11.3) 
 
c-------------------------------------------------------------------
c   Calculate stopping-power cross sections for regions A, B, and C
c-------------------------------------------------------------------
c JAF pass in maxnag
      call stopxs(maxnag, nag, nza, jza, aza, ee, cza, scxa)
      call stopxs(maxnag, nag, nzb, jzb, azb, ee, czb, scxb)
     
c-------------------------------------------------------
c   Calculate alpha source term at the interface for ab
c-------------------------------------------------------
      do 230 m=1,nag
       astab(m)=0.0
  230 continue

      fnag=nag
      do 250 k=1,nq
       r(1)=1./scxa(1)
       fact=alam(k)*faq(k)*2.5e+20
       do 240 m=1,nag
        ps(m)= 0.0
        r(m+1)= 1.0/scxa(m+1)
        p(m)=fact*(r(m+1)+r(m))*dea/2.
  240  continue
       do 246 mm = 1,nal(k)
        if(eala(mm,k).gt.eamax)write(6,131)eamax,eamax
        if(eala(mm,k).gt.6.5)eala(mm,k)=6.5
        do 245 m = 1,nag
         fm=1.
         if(eala(mm,k).ge.ee(m+1))go to 244
         if(eala(mm,k).gt.ee(m))fm=(eala(mm,k)-ee(m))/dea
         if(eala(mm,k).lt.ee(m))fm=0.0
  244    ps(m) = ps(m)+fal(mm,k)*p(m)*fm
  245   continue
  246  continue
       do 247 m=1,nag
        astab(m) = astab(m) + ps(m)
  247  continue
  250 continue

      knt=nag
      do 260 i=1,nag
       if(astab(knt).gt.0.0)go to 270
       knt=knt-1
  260 continue
  270 naga=knt
      if(naga.eq.0)write(6,130)
      if(naga.eq.0)stop 'All astab are zero'

c--------------------------------------
c   Calculate the ab transition energy
c--------------------------------------
      do 305 i=1,ncg
       if (cos(cg(i)).le.0.0) then
        stop 'Error in angular structure'
       endif
       D(i)=thick/cos(cg(i))
       do 304 ig=1,nag
        itrans(i,ig)=ig
        sum1=0.001*(dea/2.)*((1./scxb(ig))+(1./scxb(ig+1)))/anumb
        if (sum1.gt.D(i)) goto 304
  302   sum1=0
        itrans(i,ig)=itrans(i,ig)+1
        if (itrans(i,ig).gt.nag) goto 304
        do 303 m=ig,itrans(i,ig)
         sum1=sum1+0.001*(dea/2.)*((1./scxb(m))+(1./scxb(m+1)))/anumb
  303   continue
        if (sum1.lt.D(i)) goto 302
  304   continue
  305 continue

c-------------------------------------------------------
c   Calculate alpha source term at the interface for bc
c-------------------------------------------------------
      do 330 m=1,nag
       astbc(m)=0.0
  330 continue

      fnag=nag
      do 350 k=1,nq
       fact=alam(k)*faq(k)*1.25e+20
       do 348 i=1,ncg
        r(1)=1./scxa(1)
        do 340 m=1,nag
         ps(m)= 0.0
         r(m+1)= 1.0/scxa(m+1)
         p(m)=fact*(r(m+1)+r(m))*dea/2.
  340   continue
        do 346 mm = 1,nal(k)
         if(eala(mm,k).gt.eamax)write(6,131)eamax,eamax
         if(eala(mm,k).gt.6.5)eala(mm,k)=6.5
         do 345 m = 1,nag
          fm=1.
          if(eala(mm,k).ge.ee(m+1))go to 344
          if(eala(mm,k).gt.ee(m))fm=(eala(mm,k)-ee(m))/dea
          if(eala(mm,k).lt.ee(m))fm=0.0
  344     pm=1.
          if(eala(mm,k).lt.ee(itrans(i,m)))pm=0.0
          ps(m) = ps(m)+fal(mm,k)*p(m)*fm*pm
  345    continue
  346   continue
        do 347 m=1,nag
         astbc(m)=astbc(m)+ps(m)*(cos(2.*cg(i))-cos(2.*cg(i+1)))
  347   continue
  348  continue
  350 continue

      knt=nag
      do 360 i=1,nag
       if(astbc(knt).gt.0.0)go to 370
       knt=knt-1
  360 continue
  370 nagb=knt
      if(nagb.eq.0)write(6,130)
      if(nagb.eq.0)stop 'All astbc are zero'
      do 380 i=1,nag
       eal(i)=(ee(i+1) + ee(i))/2.
  380 continue

c--------------------------------------------------------------
c   Calculate neutron source from ab interface due to region B
c--------------------------------------------------------------
      if (ntb.ne.0) then
       title1(1)='Alphas '
       title1(2)='at inte' 
       title1(3)='rface a'
       title1(4)='b using'
       title1(5)=' region'
       title1(6)=' B mate'
       title1(7)='rials f'
       title1(8)='or neut'
       title1(9)='ron pro'
       title1(10)='duction'
c JAF pass in maxnag
       call neutron(title1,ntb,idb,atb,jp,ep,f,lzt,lat,lst,amt,el,
     1 apro,q,e,x,nzb,jzb,azb,czb,jps,jl,astab,nstabb,naga,erg,maxnag)

c--------------------------------------------------------------
c   Calculate neutron source from bc interface due to region B
c--------------------------------------------------------------
       title1(1)='Alphas '
       title1(2)='at inte' 
       title1(3)='rface b'
       title1(4)='c using'
       title1(5)=' region'
       title1(6)=' B mate'
       title1(7)='rials f'
       title1(8)='or neut'
       title1(9)='ron pro'
       title1(10)='duction'
c JAF pass in maxnag
       call neutron(title1,ntb,idb,atb,jp,ep,f,lzt,lat,lst,amt,el,
     1 apro,q,e,x,nzb,jzb,azb,czb,jps,jl,astbc,nstbcb,nagb,erg,maxnag)
      endif

c--------------------------------------------------------------
c   Calculate neutron source from bc interface due to region C
c--------------------------------------------------------------
      title1(1)='Alphas '
      title1(2)='at inte' 
      title1(3)='rface b'
      title1(4)='c using'
      title1(5)=' region'
      title1(6)=' C mate'
      title1(7)='rials f'
      title1(8)='or neut'
      title1(9)='ron pro'
      title1(10)='duction'
c JAF pass in maxnag
      call neutron(title1,ntc,idc,atc,jpc,epc,fc,lztc,latc,lstc,amtc,
     1 elc,aproc,qc,ec,xc,nzc,jzc,azc,czc,jpsc,jlc,astbc,nstbcc,nagb,
     2 erg,maxnag)

c----------------------------------
c   Calculate total neutron source
c----------------------------------
      tot=0
      title1(1)='Total n'
      title1(2)='eutron ' 
      title1(3)='product'
      title1(4)='ion fro'
      title1(5)='m all i'
      title1(6)='nterfac'
      title1(7)='es     '
      title1(8)='       '
      title1(9)='       '
      title1(10)='       '
      write(7,16)title1
      write(8,16)title1
      ebaran=0.0
      do 500 m=1,nng
       nst(m)=nstbcc(m)+nstabb(m)-nstbcb(m)
       tot=tot+nst(m)
       ebaran=ebaran+((en(m)+en(m+1))*nst(m)/2.0)
  500 continue
      ebaran=ebaran/tot
      write(7,*)
      write(7,2021)
      write(7,390)
      if (erg.ge.1) then
      write(7,80) (nst(n),n=nng,1,-1)
      elseif (erg.le.-1) then
      write(7,80) (nst(n),n=1,nng)
      else
      stop 'Error in erg input!'
      endif
      write(7,81) tot,ebaran
      gtmga=0.0
      do 600 n=1,nng
       nstnorm(n)=nst(n)/tot
       gtmga=gtmga+nstnorm(n)
  600 continue
      write(8,*)
      write(8,2021)
      write(8,390)
      if (erg.ge.1) then
      write(8,80) (nstnorm(n),n=nng,1,-1)
      elseif (erg.le.-1) then
      write(8,80) (nstnorm(n),n=1,nng)
      endif
      write(8,81) gtmga,ebaran

c---------------------------------------------------
c spectra for outp2 in ascending or descending order
c---------------------------------------------------
      if (erg.ge.1) then
       write(12,1900)
       write(12,3600)title
       write(12,3601)
       write(12,3605)
       write(12,3610)
       write(12,3615)
       write(12,3620)
       write(12,3630)(en(n),nst(n)/tot,n=nng,1,-1)
       write(12,3635)gtmga
       write(12,3611)
       write(12,3610)
       write(12,3615)
       write(12,3620)
       write(12,3630)(en(n),nst(n),n=nng,1,-1)  
       write(12,3635)tot    
      elseif (erg.le.-1) then
       write(12,1900)
       write(12,3600)title
       write(12,3601)
       write(12,3605)
       write(12,3610)
       write(12,3615)
       write(12,3620)
       write(12,3630)(en(n),nst(n)/tot,n=1,nng)
       write(12,3635)gtmga
       write(12,3611)
       write(12,3610)
       write(12,3615)
       write(12,3620)
       write(12,3630)(en(n),nst(n),n=1,nng)
       write(12,3635)tot 
      endif


c------------------
c   Output Summary
c------------------
      write(11,*)
      write(11,*)
      write(11,3000)
      write(11,3010)
      write(11,*)
      write(11,3020)tot
      if (id.gt.1) then
       write(11,3030)ebaran
       write(11,3050)
       do 1450 n=1,nng
        write(11,3051)n,nstnorm(n)
 1450  continue
      endif

c------------------------------------
c   Terminate execution of PROCESSOR
c------------------------------------
      return
      end


c=======================================================================
c  Three Region Problem Stopping Cross Section Calculator (1/99)
c=======================================================================

      subroutine stopxs(maxnag,nag,nz,jzm,azm,ee,czm,scx)

c----------
c  Storage
c----------
c JAF introduce parameters
      integer maxnag
      character jsm(105)*2
c JAF replaced both 4001's in this list with maxnag.
      dimension ee(maxnag), jzm(20), azm(20), czm(9,20),  scx(maxnag)

      common /names/ jsm
      common /masses/ amass(105)

c-----------------
c  Set parameters
c-----------------
c JAF added rnapier; "d0" added to reals
      zalp=2.d0
      alph=4.d0
      rnapier=2.7182818284590d+0

c---------------------------------------------------------------
c  For each energy ee(m) calc. and store stopping cross section
c---------------------------------------------------------------
      do 930 m=1,nag+1
       scx(m)=0.
       do 920 j=1,nz
        zmat=jzm(j)
        amat=amass(jzm(j))
c-----  nuclear stopping---
        term=zalp**0.66667 + zmat**0.66667
        sterm=sqrt(term)
        bot=zalp*zmat*(alph+amat)*sterm
        rep=32530.*amat*ee(m)/bot
        if(rep.lt.0.001) then
         dcx=1.593*sqrt(rep)
         go to 905
        endif
        if(rep.lt.10.) then
         bot=1.+6.8*rep+3.4*rep**1.5
         dcx=1.7*sqrt(rep)*alog(rep+rnapier)/bot
         go to 905
        endif
        dcx=alog(0.47*rep)/(2.*rep)
c-----  electronic stopping---
  905   if(ee(m).gt.30.) go to 906
        slow=czm(1,j)*(1000.*ee(m))**czm(2,j)
        shigh=(czm(3,j)/ee(m))
        shigh=shigh*alog(1.+czm(4,j)/ee(m)+czm(5,j)*ee(m))
        dcx=dcx+slow*shigh/(slow+shigh)
        go to 907
  906   eil=alog(1/ee(m))
        arg=czm(6,j)+eil*(czm(7,j)+czm(8,j)*eil+czm(9,j)*eil*eil)
        dcx=dcx+exp(arg)
  907  continue
  920  scx(m)=scx(m)+azm(j)*dcx
  930 continue

c--------------------------------
c  Terminate execution of STOPXS
c--------------------------------
      return
      end



c=======================================================================
c  Three Region Problem Neutron Source Subroutine (1/99)
c=======================================================================

      subroutine neutron(title,nt,idt,at,jp,ep,f,lzt,lat,lst,amt,el,
     1 apro,q,e,x,nz,jz,az,cz,jps,jl,ast,gtsan,nal,erg,maxnag)

c----------
c  Storage
c----------
      character title*7,jsm(105)*2,amt(20)*1
      integer erg
c JAF introduce parameters
      integer maxnag
c JAF replaced all 4001's in this list with maxnag.
      dimension ee(maxnag), en(751), cg(maxnag), e(1100,20), el(20,20),
     1 x(1100,20), idt(20), at(20), scx(maxnag), jps(20), jl(20), 
     2 cx(maxnag), eal(maxnag), apro(20), q(20), jp(20), ep(200,20),
     3 f(20,200,20), lat(20), ast(maxnag), san(750), ts(750), lzt(20),
     4 r(maxnag), rr(maxnag), p(maxnag), tsan(750), gtsan(750), s(750), 
     5 totlev(20), sl(750,20), etopl(20), ebotl(20), jz(20), az(20),
     6 cz(9,20), gtnorm(750), title(11)

       common /names/ jsm
       common /masses/ amass(105)
c JAF remove ee, cg, eal from common
c      common /grids/ ee, en, cg, nag, eamax, eamin, nng, 
c    1 enmax, enmin, ncg, dea, eal
       common /grids/ en, nag, eamax, eamin, nng, 
     1 enmax, enmin, ncg, dea

c--------------------
c formats statements
c--------------------
   10 format (11a7)
   15 format(//,11a7,/)
   16 format(//,10(8h========),//,'Title: ',11a7)
   20 format (3i3)
   30 format (i8,e12.5)
   40 format(i3,8e8.3,f10.3,/,3x,5f8.3)
   50 format (69h failed to find stopping cross section coefficients on
     1tape 2 for iz=,i3,5h stop)
   60 format (6e12.5)
   70 format(/,1x,'Neutron Multigroup Structure (MeV)')
   79 format(/,1x,'Neutron Spectrum (neuts/cm^2-sec)')
   80 format (1p8e10.3)
   81 format(31x,'Total (all groups): ',1pE10.3,' neutrons/sec-cm^2.',
     1/,27x,'Average Neutron Energy: ',1pE10.3,' MeV.')
   90 format (i3,5e12.5)
   95 format(2i6)
  100 format(//,1h1,30x,'Table I',/,31x,7(1h=),//,12x,'(alpha,n) ',
     1'Neutron Production by Target and Alpha Energy',///,17x,
     2'target            alpha   alphas/sec     p(e)     neuts/sec',
     3/,7x,'target  atom frac.          energy     /cm^2',4x,
     4'neut/alpha',4x,'/cm^2',/,1h+,6x,6(1h_),2x,10(1h_),2x,6(1h_),2x,
     5 6(1h_),3(12h  __________))
c 102 format(1h+,77x,37hfractional ground/excited level split,/,
c    1 1h+,77x,49(1h_))
  105 format(i8,i4)
  120 format (8e10.3)
  130 format (52h failed to find alpha-n cross section on tape 3 for ,i8
     1 ,5h stop)
  140 format (i8,2i3)
  150 format (66h failed to find comp.nucl.neut.emiss.lev. branching on 
     1tape 4 for ,i8,5h stop)
  170 format (/,17h alpha energy of ,1pe12.5,62h exceeds current 6.5 MeV
     1 limit of most cross sections, stop.  )
  180 format (7h target,i8,40h cross-section data not available at ee(
     1 ,i4,3h)= ,1pe12.5,6h stop.)
  200 format(33x,f8.3,1p3e12.4)
c 205 format(1h+,(76x,7f7.5,/,1x))
  210 format(7x,a2,i3,a1,1pe12.4,8x,         0pf8.3,1p3e12.4)
  220 format (14h alpha energy ,1pe12.5,45h exceeds range of level branc
     1hing data, stop.)
  240 format(/,' (a,n) neutrons from alphas on ',a2,i3,a1,' in target.')
  250 format(1h+,66x,10(1h_),/,55x,'Total:',4x,1pe12.4,/)
  251 format(1h+,66x,10(1h_),/,41x,'Total (this target):',
     14x,1pe12.4,/)
  252 format (1h+,66x,10(1h_),/,41x,'Total (all targets):',
     14x,1pe12.4,/)
  255 format(1h )
  260 format (///,' Total (alpha,n) neutron spectrum this target')
  270 format(///,' Grand total (alpha,n) neutron spectrum, all targets,'
     1 ,1x,'all sources')
  390 format (////,' Total Neutron Spectrum')
 2021 format(/,80(1h-))
 2030 format('Target title: ', 11a7)
 2050 format('Number of elemental constituents on target side:',i3)
 2060 format('Solid stopping cross-sections used (isg= ',i2,') on',
     1 1x,'target side.')
 2070 format('Gas stopping cross-sections used (isg= ',i2,')',
     1 1x,'on target side.')
 2082 format('Elemental Constituents on Target Side:')
 2090 format(' ')
 2100 format('  Z-value  Atom Fraction')
 2110 format('  -------  -------------')
 2120 format(4x,i3,5x,f12.10)
 2130 format('Number of neutron spectrum energy groups:', i4)
 2140 format('Maximum neutron energy is ',1pE12.5,' MeV.')
 2150 format('Minimum neutron energy is ',1pE12.5,' MeV.')
 2160 format('Energy Group Structure:')
 2170 format('  Group    Upper-Bound    Lower-Bound')
 2180 format('  -----    -----------    -----------')
 2190 format(2x,i4,4x,1pE12.5,3x,1pE12.5)
 2193 format(/,'Warning: First Group Upper-Bound Does Not Coincide ',
     1'With Maximum Neutron Energy.')
 2195 format(/,'Warning: Energy Group Structure Not Entered in ',
     1'Decreasing Order.',/,5x,'Group ',i3,' Upper-Bound is Greater ',
     2'Than Group ',i3,' Upper-Bound.') 
 2260 format(/,'Number of target nuclides to be used:',i3)
 2270 format(i5,' alpha energy groups used in target calculation.')
 2280 format('Target Nuclides:')
 2290 format('    ZAID     Atom Fraction')
 2300 format('    ----     -------------')
 2310 format(3x,i6,4x,1pE12.5)

 3000 format(//,22x,'Summary of Output')
 3010 format(22x,'=================')
 3020 format('Total (alpha,n) neutron source from all sources and ',
     1'targets: ',1pE12.5,' n/sec-cm^2.')
 3030 format(/,'Average (alpha,n) neutron energy: ',1pE12.5,' MeV.')
 3050 format(/,'Normalized Neutron Energy Spectrum by Energy Group ',
     1'for All Sources:',/'(Note: group structure is independent of',
     2' erg record!)',//,10x,'Group',6x,'Contribution',/,10x,
     3'-----',6x,'------------')
 3051 format(11x,i3,7x,1pE12.5)
 3500 format(/,'Portion of Total Neutron Source Rate Accounted for ',
     1'in the Energy Spectrum: ',F5.1,'%.')
 3510 format(/,'Warning: Energy Group Structure May Be Poorly Chosen.')

c---------------
c write headers
c---------------
      write(6,16)title
      write(7,16)title
      write(8,16)title

c--------------------------
c loop on target nuclides i
c--------------------------
      alph=4.
      aneut=1.
      do 1120 i=1,nt
       etpant=0.
       ebtant=0.
       totqan=0.
       if(id.eq.1) go to 790
       do 650 n=1,nng
        tsan(n)=0.
  650  continue
       do 700 ip=2,jps(i)
        eamin=e(ip-1,i)
        if (x(ip,i).ne.0.) go to 710
  700  continue
  710  if (eamin.lt.0.001) eamin=0.001
       etop=0.
       ebot=0.
       do 775 il=1,jl(i)
        etopl(il)=0.
        ebotl(il)=0.
        do 774 n=1,nng
         sl(n,il)=0.
  774   continue
        totlev(il)=0.
  775  continue
       sbtqan=0.
       if(id.eq.1) go to 790
       do 780 n=1,nng
        san(n)=0.
  780  continue
  790  continue

c---------------------------------------------------------
c set alpha energy grid ee(m) for this target/source combo
c---------------------------------------------------------
      eamax=eal(nag)
  850 continue
      if(eamax.le.eamin) then
      qan=0.
      mm=0
      p(nag+1)=0.
      pval=0.
      go to 945
      endif
      if (eamax.le.6.5) go to 860
      write (6,170) eamax
      stop 5
  860 fnag=nag
      nagp1=nag+1
      dea=(eamax-eamin)/fnag
      ee(1)=eamin
      do 870 m=2,nag
       fm=m-1
  870 ee(m)=ee(1)+dea*fm
      ee(nagp1)=eamax

c----------------------------------------------------
c for each energy ee(m) calc. and store cross section
c----------------------------------------------------
       ip=1
       do 900 m=1,nag+1
  880   if (ee(m).ge.e(ip,i).and.ee(m).le.e(ip+1,i)) go to 890
        ip=ip+1
        if (ip.lt.jps(i)) go to 880
        write (6,180) idt(i),m,ee(m)
        stop 6
  890   slope=(x(ip+1,i)-x(ip,i))/(e(ip+1,i)-e(ip,i))
        enrsep=(x(ip,i)*e(ip+1,i)-x(ip+1,i)*e(ip,i))/
     1   (e(ip+1,i)-e(ip,i))
        cx(m)=enrsep+ee(m)*slope
  900  continue

c-------------------------------------------------------------
c for each energy ee(m) calc. and store stopping cross section
c-------------------------------------------------------------
c JAF pass in maxnag
      call stopxs(maxnag, nag, nz, jz, az, ee, cz, scx)

c-----------------------------------------------------------
c for each energy ee(m) calc. ratio r(m), then integral p(m)
c-----------------------------------------------------------
       r(1)=cx(1)/scx(1)
       p(1)=0.
       fact=1.0e-06*at(i)
       do 940 m=2,nag+1
        r(m)=cx(m)/scx(m)
        p(m)=p(m-1)+fact*(r(m-1)+r(m))*dea/2.
  940  continue
  945  continue

c----------------------------------------------
c  calculate neutrons from alphas
c----------------------------------------------
  950 do 1080 l=1,nal
      do 960 m=1,nag
       mm=m
  960 if (eal(l).ge.ee(m).and.eal(l).le.ee(m+1)) go to 970
      mm=0
      pval=0.
      go to 975
  970 pval=p(mm)+(eal(l)-ee(mm))*(p(mm+1)-p(mm))/(ee(mm+1)-ee(mm))
  975 aps=ast(l)
      qan=aps*pval
      sbtqan=sbtqan+qan
      totqan=totqan+qan
      if (l.eq.1) go to 980
      write (6,200) eal(l),aps,pval,qan
      if(nal.ne.1.and.l.eq.nal) write(6,250) sbtqan
      go to 990
  980 write(6,210)jsm(lzt(i)),lat(i),amt(i),at(i),eal(l),aps,pval,qan
      if(nal.eq.1) write(6,255)
  990 if (id.eq.1.or.mm.eq.0) go to 1080

c-----------------------------------------------------------------
c calculate (alpha,n) neutron spectrum contribution in multigroups
c-----------------------------------------------------------------
 1000 mmm1=mm-1
      do 1010 m=1,mmm1
 1010 rr(m)=(p(m+1)-p(m))/pval
      rr(mm)=(pval-p(mm))/pval
      do 1020 n=1,nng
 1020 s(n)=0.
      do 1060 il=1,jl(i)
      qlev=q(i)-el(il,i)
      e90=-qlev*apro(i)/(apro(i)-alph)
      thre=-qlev*(aneut+apro(i))/(aneut+apro(i)-alph)
      if (qlev.gt.0.) thre=0.
      do 1060 m=1,mm
      ea=(ee(m)+ee(m+1))/2.
      if (m.eq.mm) ea=(eal(l)+ee(mm))/2.
      if (ea.lt.thre) go to 1055
      do 1030 ip=2,jp(i)
      lp=ip
 1030 if (ep(ip-1,i).le.ea.and.ep(ip,i).ge.ea) go to 1040
      write (6,220)ea
      stop 7
 1040 lp1=lp-1
      bx=f(il,lp1,i)+(ea-ep(lp1,i))*(f(il,lp,i)-f(il,lp1,i))/
     1 (ep(lp,i)-ep(lp1,i))
      term1=sqrt(4.*alph*aneut*ea)/(2.*(aneut+apro(i)))
      term2=alph*aneut*ea/((aneut+apro(i))*(aneut+apro(i)))
      term3=(apro(i)*ea+apro(i)*qlev-alph*ea)/(aneut+apro(i))
      senmax=term1+sqrt(term2+term3)
      enmax=senmax*senmax
      if (ea.le.e90) senmin=term1-sqrt(term2+term3)
      if (ea.gt.e90) senmin=-term1+sqrt(term2+term3)
      enmin=senmin*senmin
      ebar=(enmin+enmax)/2.
      val=qan*rr(m)*bx
      etop=etop+val*ebar
      etopl(il)=etopl(il)+val*ebar
      ebotl(il)=ebotl(il)+val
      ebot=ebot+val
      dele=enmax-enmin
      do 1050 n=1,nng
      if (en(n).lt.enmin) go to 1050
      if (en(n+1).gt.enmax) go to 1050
      de=en(n)-en(n+1)
      if (en(n+1).lt.enmin) de=de-(enmin-en(n+1))
      if (en(n).gt.enmax) de=de-(en(n)-enmax)
      gpadd=qan*rr(m)*bx*de/dele
      s(n)=s(n)+gpadd
      sl(n,il)=sl(n,il)+gpadd
      totlev(il)=totlev(il)+gpadd
 1050 continue
 1055 continue
 1060 continue
      do 1070 n=1,nng
      san(n)=san(n)+s(n)
      gtsan(n)=gtsan(n)+s(n)
      ts(n)=ts(n)+s(n)
 1070 tsan(n)=tsan(n)+s(n)
 1075 continue
 1080 continue
      if(id.eq.1) go to 1200

c----------------------------------------------------
c output this (alpha,n) neutron spectrum contribution
c for this source/target combination
c----------------------------------------------------
      ebarqt=0.
      if(ebot.le.0.) go to 1085
      ebarqt=etop/ebot
      etpant=etpant+etop
      ebtant=ebtant+ebot
 1085 continue
      smga=0.
      do 1090 n=1,nng
 1090 smga=smga+san(n)
c     fracgp=smga/sbtqan
c     do 1095 il=1,jl(i)
c1095 frclev(il)=totlev(il)/sbtqan
c     write(6,205) (frclev(il),il=1,jl(i))
      if (nq.ne.0) write (7,240) jsm(lzt(i)),lat(i),amt(i)
      write (7,79)
      if (erg.ge.1) then
      write (7,80) (san(n),n=nng,1,-1)
      elseif (erg.le.-1) then
      write (7,80) (san(n),n=1,nng)
      endif
      write (7,81) smga,ebarqt
      smga=0.
      if (nq.ne.0) write (8,240) jsm(lzt(i)),lat(i),amt(i)
      do 1097 n=1,nng
      san(n)=san(n)/sbtqan
 1097 smga=smga+san(n)
      if (erg.ge.1) then
      write (8,80) (san(n),n=nng,1,-1)
      elseif (erg.le.-1) then
      write (8,80) (san(n),n=1,nng)
      endif
      write (8,81) smga, ebarqt
      do 1098 il=1,jl(i)
      if(ebotl(il).le.0.) go to 1098
c     ebarl=etopl(il)/ebotl(il)
 1098 continue
 1101 gttqan=gttqan+totqan

c--------------------------------------
c output total (alpha,n) source for all
c alpha sources on this target nuclide
c--------------------------------------
      if(id.eq.1) go to 1120
      tmga=0.
      do 1110 n=1,nng
 1110 tmga=tmga+tsan(n)
c     fracgp=tmga/totqan
      ebart=etpant/ebtant
      etopan=etopan+etpant
      ebotan=ebotan+ebtant
      if(nq.le.1) go to 1120
      write(7,260)
      if (erg.ge.1) then
      write(7,80) (tsan(n),n=nng,1,-1)
      elseif (erg.le.-1) then
      write(7,80) (tsan(n),n=1,nng)
      endif
      write(7,81) tmga,ebart
      write(7,2090)
      write(7,2021)
      write(7,2090)
      tmga=0.
      do 1115 n=1,nng
       tsan(n)=tsan(n)/totqan
 1115 tmga=tmga+tsan(n)
      write(8,260)
      if (erg.ge.1) then
      write(8,80) (tsan(n),n=nng,1,-1)
      elseif (erg.le.-1) then
      write(8,80) (tsan(n),n=1,nng)
      endif
      write(8,81) tmga,ebart
      write(8,2090)
      write(8,2021)
      write(8,2090)
      write(10,2090)
      write(10,2021)
      write(10,2090)
 1120 continue

c----------------------------------------------------------------------
c output grand total (alpha,n) neutron source: all targets, all sources
c----------------------------------------------------------------------
      if (nt.gt.1) write (6,252) gttqan
      if(id.eq.1) go to 1200
      gtmga=0.
      do 1130 n=1,nng
 1130 gtmga=gtmga+gtsan(n)
      ebaran=etopan/ebotan
c     fracgp=gtmga/gttqan
      write(7,2090)
      write(7,2021)
      write(7,390)
      if (erg.ge.1) then
      write(7,80) (gtsan(n),n=nng,1,-1)
      elseif (erg.le.-1) then
      write(7,80) (gtsan(n),n=1,nng)
      endif
      write(7,81) gtmga,ebaran
      gtmga=0.
      dummy=0.0
      do 1135 n=1,nng
       dummy=dummy+gtsan(n)
       gtnorm(n)=gtsan(n)/gttqan
 1135 gtmga=gtmga+gtnorm(n)
      write(8,2090)
      write(8,2021)
      write(8,390)
      if (erg.ge.1) then
      write(8,80) (gtnorm(n),n=nng,1,-1)
      elseif (erg.le.-1) then
      write(8,80) (gtnorm(n),n=1,nng)
      endif
      write(8,81) gtmga,ebaran
 1200 continue


c------------------------------
c  End of Execution of NEUTRON
c------------------------------
      return
      end


c=======================================================================
c  Block Data Statement (7/97)
c  amass values reflect CRC Handbook (81st Ed. 2000) (EFS 11 Apr 01)
c=======================================================================

      block data initvl
      character jsm(105)*2
      dimension amass(105)
      common /names/ jsm
      common /masses/ amass

      data jsm /' h','he','li','be',' b',' c',' n',' o',' f','ne','na',
     1 'mg','al','si',' p',' s','cl','ar',' k','ca','sc','ti',' v','cr',
     2 'mn','fe','co','ni','cu','zn','ga','ge','as','se','br','kr','rb',
     3 'sr',' y','zr','nb','mo','tc','ru','rh','pd','ag','cd','in','sn',
     4 'sb','te',' i','xe','cs','ba','la','ce','pr','nd','pm','sm','eu',
     5 'gd','tb','dy','ho','er','tm','yb','lu','hf','ta',' w','re','os',
     6 'ir','pt','au','hg','tl','pb','bi','po','at','rn','fr','ra','ac',
     7 'th','pa',' u','np','pu','am','cm','bk','cf','es','fm','md','no',
     8 'lr','rf','db'/

      data amass /1.0079,4.0026,6.941,9.0121,10.811,12.0107,14.0067,
     1 15.9994,18.9984,20.1797,22.9897,24.3050,26.9815,28.0855,30.9737,
     2 32.066,35.4527,39.948,39.0983,40.078,44.9559,47.867,50.9415,
     3 51.9961,54.9380,55.845,58.9332,58.6934,63.546,65.39,69.723,
     4 72.61,74.9216,78.96,79.904,83.80,85.4678,87.62,88.9058,91.224,
     5 92.9063,95.94,98.,101.07,102.9055,106.42,107.8682,112.411,
     6 114.818,118.710,121.760,127.6,126.9044,131.29,132.9054,137.327,
     7 138.9055,140.116,140.9076,144.24,145.,150.36,151.964,157.25,
     8 158.9253,162.5,164.9303,167.26,168.9342,173.04,174.967,178.49,
     9 180.9479,183.84,186.207,190.23,192.217,195.078,196.9665,200.59,
     1 204.3833,207.2,208.9803,209.,217.,222.,223.,226.,227.,232.0381,
     2 231.0358,238.0289,237.,244.,243.,247.,247.,251.,252.,257.,258.,
     3 259.,262.,261.,262./
      end
