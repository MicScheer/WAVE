*CMZ :  4.00/11 17/06/2021  09.33.40  by  Michael Scheer
*CMZ :  3.05/14 27/09/2018  09.23.18  by  Michael Scheer
*CMZ :  3.05/12 29/08/2018  09.36.53  by  Michael Scheer
*CMZ :  3.05/09 01/08/2018  08.45.26  by  Michael Scheer
*CMZ :  3.05/07 20/07/2018  13.05.35  by  Michael Scheer
*-- Author :    Michael Scheer   20/07/2018
*KEEP,BSCALEMOD.
      module mrad_bscalemod

      double precision, save :: bxscale_m=1.0, byscale_m=1.0, bzscale_m=1.0

      double precision, save :: brotmat_m(3,3)=reshape([
     &  1.,0.,0.,0.,1.,0.,0.,0.,1.],[3,3])

      end module mrad_bscalemod
*KEEP,BFOURMOD.
      module mrad_bfourmod

      integer, parameter :: maxfoumagp=2**10,
     &  maxfour_p_m=2**10,
     &  nfourd_p_m=2*maxfour_p_m

      integer, save ::
     &  nfour_m(maxfoumagp)=1,
     &  ifour0_m(maxfoumagp)=0,
     &  nfoumags_m(maxfoumagp)=1,
     &  kbfourini_m=0,
     &  kbfour_m=0

      double precision, save ::
     &  xlenfour_m(maxfoumagp)=0.0d0,
     &  fouentr_m(maxfoumagp)=-9999.0d0,
     &  fouexit_m(maxfoumagp)=-9999.0d0,
     &  xshbfour_m(maxfoumagp)=0.0d0,
     &  zl0four_m(maxfoumagp),zl0four2_m(maxfoumagp),xk0four_m(maxfoumagp),
     &  a0_m(maxfoumagp)=0.0d0,a_m(maxfour_p_m,maxfoumagp)=0.0d0,
     &  xkfour_m(maxfour_p_m,maxfoumagp)=0.0d0,
     &  ykfour_m(maxfour_p_m,maxfoumagp)=0.0d0,
     &  zkfour_m(maxfour_p_m,maxfoumagp)=0.0d0

      character(2048) chfilebfour_m

      end module mrad_bfourmod
*KEEP,BMAPMOD.
      module mrad_bmapmod

      double precision, dimension (:,:), allocatable,save :: bmappe

      double precision,save :: bmxmin=1.0d99,bmxmax=-1.0d99
      double precision,save :: bmymin=1.0d99,bmymax=-1.0d99
      double precision,save :: bmzmin=1.0d99,bmzmax=-1.0d99

      double precision,save :: bmbxmin=1.0d99,bmbxmax=-1.0d99
      double precision,save :: bmbymin=1.0d99,bmbymax=-1.0d99
      double precision,save :: bmbzmin=1.0d99,bmbzmax=-1.0d99

      double precision,save :: scalex=1.0d0, scaley=1.0d0, scalez=1.0d0
      double precision,save :: scalebx=1.0d0, scaleby=1.0d0, scalebz=1.0d0

      double precision, save :: bmapdy,bmapdz

      integer,save :: mbmapini_m=0, intpolbmap=-1, kbmap_m=1, kalloc_bmap_m=0
      integer, save :: ntot,nx,ny,nz,nyz

      end module mrad_bmapmod
*KEEP,PHYCONMOD.
      module phycon

      complex*16, parameter :: zone=(1.0d0,0.0d0), zi=(0.0d0,1.0d0)
      complex*16, save :: vstokes(4,3)

      DOUBLE PRECISION, save ::
     &  CLIGHT1,CGAM1,CQ1,alpha1,
     &  dnull1,done1,sqrttwopi1,
     &  EMASSG1,EMASSE1,ECHARGE1,EMASSKG1,EPS01,ERAD1,
     &  GRARAD1,
     &  HBAR1,HBAREV1,HPLANCK1,
     &  PI1,POL1CON1,POL2CON1,
     &  RADGRA1,rmu01,rmu04pi1,
     &  TWOPI1,HALFPI1,
     &  WTOE1,gaussn1,cK934,
     &  ecdipev1,ecdipkev1,powcon1

      end module phycon
*KEEP,MRADMOD.
      module mrad_mradmod

      integer, parameter :: ndim1elecs_p_m=17, ndim1optics_m=27

      integer, save :: iwritefiles_m=0, iploteps_m=0, iwtra_m=1
      integer, save :: nstepp_m=10000, nphp_m=1000
      integer, save :: minit_m=0, modeobsrndm_m=1, modeobsv_m=1,
     &  nturns_m=1, koptics_m=1, kapert_m=0, nturnstot_m=0, iturnstot_m=0,
     &  n10turns_m=1,kdynaper_m=0,ndynang_m=2, kplstart_m=1, kltetolat_m=0,
     &  ksubco_m=0, ksynpho_m=0, lunout_m, kturn_m=0

      complex*16, dimension (:), save, allocatable ::  aradx_m,arady_m,aradz_m

      double precision, dimension (:,:,:), allocatable :: stoksum_m, ampsum_m
      double precision, dimension (:,:), save, allocatable :: elecs_m
      double precision, dimension (:), allocatable :: powsum_m

      double precision, dimension (:,:,:), save, allocatable ::
     &  traxyz_m

      double precision, dimension (:,:), save, allocatable ::
     &  stokes_m, obsv_m, optics_m, dynr_m, dynaper_m, tref_m

      double precision, dimension (:), save, allocatable :: phener_m

      double precision, save :: betah_m, betav_m, betaph_m, betapv_m,
     &  disph_m, dispv_m, dispph_m, disppv_m,
     &  epsh_m, epsv_m, pathlen_m,circum_m,
     &  tfinvbetah_m(2,2),tfinvbetav_m(2,2),tfbetah_m(2,2),tfbetav_m(2,2),
     &  bmovecut_m, beta_per_h_m, betap_per_h_m, beta_per_v_m, betap_per_v_m,
     &  tuneh_m,tunev_m,xelec_m,yelec_m,zelec_m,vxelec_m,vyelec_m,vzelec_m,
     &  dynrstart_m,dyndrmin_m,dynangi_m,dynange_m,
     &  traref_m(14,2),emithor_m=0.0d0,
     &  xlat_m=0.0d0, ylat_m=0.0d0, zlat_m=0.0d0, dzlat_m=0.0d0,
     &  exlat_m=0.0d0, eylat_m=0.0d0, ezlat_m=0.0d0, ebeam_m, brho_m, emom_m,
     &  aperhor_m,aperver_m, beta_m,v0ele_m,gammai_m,ds_m,
     &  cox_m,coy_m,coz_m,cos_m,coyp_m,cozp_m,cot_m,dering_m,cavigain_m,
     &  dmomtot_m(3),dmom_m(3),quadscal_m,sextscal_m,dipscal_m,octoscal_m

      integer, save :: kelec_m, kthread_m, kstep_m, ksteps_m, mobscen_m,
     &  nphener_m,
     &  nobsv_m, nelec_m, nthphplot_m=1, nthelec_m=1, ngoodelecs_m=0,
     &  nyobsv_m,nzobsv_m,krun_m, kbecho_m=1, nelectra_m,
     &  n1dimtra_m=19, klost_m, kfirstplane_m=1, klastplane_m=0

      integer, save :: kalloc_obsv_m=0, kalloc_phener_m=0,
     &  kalloc_elec_m=0, kalloc_track_m=0, kalloc_arad_m=0,
     &  kalloc_stoksum_m(2)=0, kalloc_powsum_m=0, kalloc_ampsum_m(2)=0,
     &  kalloc_optics_m=0, nmapdat_m=0, kisenall_m=1

      integer :: ignoreaper_m=0

!$OMP THREADPRIVATE(kelec_m,kthread_m,kturn_m,kstep_m,ksteps_m, kapert_m,
!$OMP& kfirstplane_m, klastplane_m,aradx_m,arady_m,aradz_m,stokes_m,kalloc_arad_m,
!$OMP& nphp_m)

      end module mrad_mradmod
*KEEP,MRADLATTICE.
      module mrad_lattice

      integer, parameter :: ntramat_m_p=1000, nmap_m_p=1000, ndimring_m=63,
     &  ndimlat_m=7

      double precision, save ::
     &  tramat_m(7,7,ntramat_m_p),offset_m(7,ntramat_m_p),
     &  tfmat_m(6,6,ntramat_m_p)

      double precision, dimension (:), save, allocatable :: bmaps

      double precision, dimension (:,:,:), save, allocatable :: design_orbit_m
      double precision, dimension (:,:), save, allocatable :: ring_m
      character(128), dimension(:), allocatable :: chring_m

      double precision, save :: delgam_m=0.002, deltf_m=0.00001,
     &  epsbet_m=1.0d-8, btabeps_m=1.0d-4, fringe_m=1000.0d0, cavivolt_m=9999.

      integer, save :: kalloc_lattice_m=0, nplane_m=0, mplane_m=0,
     &  kplane_m=0, kmplane_m=0, kplaneold_m=0, kfouind_m=0, kreftrack_m=0,
     &  koptic_m=0, ktfind_m=0, ntfind_m=0, nosind_m=0, kosind_m=0,
     &  ntrmind_m=0, ktrmind_m=0, nquad_m=0, nmap_m=0, nsext_m=0, nocto_m=0,
     &  nbxyz_m=0,icharge_m=-1, nbend_m=0, ncavi_m=0,
     &  nmapsize_m(3,nmap_m_p),nmappoi_m=0, kelem_m=1, intpolbmaps_m=-1,
     &  koffquad_m=0, koffsext_m=0, koffocto_m=0, nbpm_m, koffcavi_m=0

      integer :: modfringe_m=2

      integer, dimension(:,:), save, allocatable :: lattice_m
      integer, dimension(:), save, allocatable :: kplanidx_m, kbpm_m

!$OMP THREADPRIVATE(nmappoi_m,kfouind_m,kplane_m,kmplane_m,kplaneold_m,kelem_m,kosind_m,ktfind_m,ktrmind_m)

      end module mrad_lattice
*KEEP,MRADERROR.
      module mrad_error

      integer, save :: kwarn_nstep_m=1, debug_m=0

      end module mrad_error
*KEEP,MRADFILES.
      module mrad_files

      integer, save :: kwtrack_m=0, kwtrack1_m=0, kwtrack2_m=0,
     &  luntrack_m=0, luntracks_m=0, lunsyn_m, lunbpm_m,luncav_m

      character(2048), save :: chfiletrack_m="mrad_reference_orbit.dat"
      character(2048), save :: chfiletracks_m="mrad_trajectories.dat"

      end module mrad_files
*KEEP,MRADFIT.
      module mrad_fit

      integer nvar,nparp,npar
      parameter (nparp=100)

      integer, save :: ifitorbit_m=0, nfitpar_m_p=0, nfitcoupmaxelem_m=0,
     &  lunfit_m=0, nfitcoupmaxpar_m=0, nfitelemaxpar_m=0, nfitorbit_m=0,
     &  ifitoptic_m=0, nfitoptic_m, ifitco_m, nfitco_m

      integer, save ::
     &  nfitelem_m=0, nfitpar_m=0, lsave_m=7

      double precision, dimension(:,:), allocatable ,save ::
     &  fitpar_m, bpmfit_m, fitorbit_m

      integer, dimension (:,:), allocatable, save ::
     &  kfitcouplelem_m, nfitcoupelem_m,
     &  kfitcouplpar_m, nfitcouppar_m,
     &  kbpmparidx_m, idxfitelem_m, idxfitpar_m,
     &  idxfitcoupelem_m, idxfitcouppar_m

      character(128), dimension (:), allocatable, save ::
     &  chfitelem_m, chfitbpm_m

      character(10), dimension (:), allocatable, save :: chparam_m

      double precision, save :: chi2, chi2min=1.0d30

      double precision, dimension(:,:), allocatable ,save :: var, varopt_m
      double precision, dimension(:), allocatable ,save ::
     &  paramopt_m

      integer, save ::
     &  minit_fcn,ithread_fcn,ielec_fcn,icharge_fcn,
     &  ivelofield_fcn,itrackback_fcn,nstep_fcn, nthstep_fcn,
     &  ieneloss_fcn,nphener_fcn, maxcalls_fcn, ifitverbose_m, kfit_m

      double precision, save ::
     &  gammai_fcn, dgamtot_fcn, dmomtot_fcn(3),
     &  xelec_fcn, yelec_fcn, zelec_fcn,
     &  vxelec_fcn, vyelec_fcn, vzelec_fcn,
     &  dxelec_fcn, ds_fcn,
     &  xf_fcn, yf_fcn, zf_fcn, efx_fcn,efy_fcn,efz_fcn,
     &  xobsv_fcn, yobsv_fcn, zobsv_fcn,
     &  xexit_fcn, yexit_fcn, zexit_fcn,
     &  vnxex_fcn, vnyex_fcn,vnzex_fcn,sexit_fcn,
     &  texit_fcn,powden_fcn, tolerance_fcn

      character(32) chmethod_fcn

      end module mrad_fit
*KEEP,MRADENE.
      module mrad_ene
      integer ieneloss_m
      end module mrad_ene
