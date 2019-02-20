MODULE asmbkg
   !!======================================================================
   !!                       ***  MODULE asmtrj -> asmbkg  ***
   !! Assimilation trajectory interface: Write to file the background state and the model state trajectory
   !!======================================================================
   !! History :       ! 2007-03 (M. Martin)  Met. Office version
   !!                 ! 2007-04 (A. Weaver)  asm_trj_wri, original code
   !!                 ! 2007-03 (K. Mogensen)  Adapt to NEMOVAR and use IOM instead of IOIPSL
   !!                 ! 2007-04 (A. Weaver)  Name change (formally asmbkg.F90). Distinguish
   !!                                        background states in Jb term and at analysis time.
   !!                                        Include state trajectory routine (currently empty)
   !!                 ! 2007-07 (A. Weaver)  Add tke_rst and flt_rst for case nitbkg=0 
   !!                 ! 2009-03 (F. Vigilant)  Add hmlp (zdfmxl) for no tracer nmldp=2 
   !!                 ! 2009-06 (F. Vigilant) asm_trj_wri: special case when kt=nit000-1
   !!                 ! 2009-07 (F. Vigilant) asm_trj_wri: add computation of eiv at restart
   !!                 ! 2010-01 (A. Vidard) split asm_trj_wri into tam_trj_wri and asm_bkg_wri
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   'key_asminc' : Switch on the assimilation increment interface
   !!----------------------------------------------------------------------
   !!   asm_bkg_wri  : Write out the background state
   !!   asm_trj_wri  : Write out the model state trajectory (used with 4D-Var)
   !!----------------------------------------------------------------------
   USE oce                ! Dynamics and active tracers defined in memory
   USE sbc_oce            ! Ocean surface boundary conditions
   USE zdf_oce            ! Vertical mixing variables
   USE zdfddm             ! Double diffusion mixing parameterization
   USE ldftra_oce         ! Lateral tracer mixing coefficient defined in memory
   USE ldfslp             ! Slopes of neutral surfaces
   USE tradmp             ! Tracer damping
#if defined key_zdftke
   USE zdftke             ! TKE vertical physics
#endif
   USE eosbn2             ! Equation of state (eos_bn2 routine)
   USE zdfmxl             ! Mixed layer depth
   USE dom_oce, ONLY :   ndastp
   USE sol_oce, ONLY :   gcx   ! Solver variables defined in memory
   USE in_out_manager     ! I/O manager
   USE iom                ! I/O module
   USE asmpar             ! Parameters for the assmilation interface
   USE zdfmxl             ! mixed layer depth
#if defined key_traldf_c2d
   USE ldfeiv             ! eddy induced velocity coef.      (ldf_eiv routine)
#endif
#if defined key_lim2
   USE ice_2
#endif
#if defined key_lim3
   USE ice
#endif
   USE asminc, ONLY: ln_avgbkg
   IMPLICIT NONE
   PRIVATE
   
   PUBLIC   asm_bkg_wri   !: Write out the background state

  !! * variables for calculating time means
   REAL(wp),SAVE, ALLOCATABLE,   DIMENSION(:,:,:) ::   tn_tavg  , sn_tavg  
   REAL(wp),SAVE, ALLOCATABLE,   DIMENSION(:,:,:) ::   un_tavg  , vn_tavg
   REAL(wp),SAVE, ALLOCATABLE,   DIMENSION(:,:,:) ::   avt_tavg
#if defined key_zdfgls || key_zdftke
   REAL(wp),SAVE, ALLOCATABLE,   DIMENSION(:,:,:) ::   en_tavg
#endif
   REAL(wp),SAVE, ALLOCATABLE,   DIMENSION(:,:)   ::   sshn_tavg
   REAL(wp),SAVE :: numtimes_tavg     ! No of times to average over

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE asm_bkg_wri( kt )
      !!-----------------------------------------------------------------------
      !!                  ***  ROUTINE asm_bkg_wri ***
      !!
      !! ** Purpose : Write to file the background state for later use in the
      !!              inner loop of data assimilation or for direct initialization
      !!              in the outer loop.
      !!
      !! ** Method  : Write out the background state for use in the Jb term
      !!              in the cost function and for use with direct initialization
      !!              at analysis time.
      !!-----------------------------------------------------------------------
      INTEGER, INTENT( IN ) :: kt               ! Current time-step
      !
      CHARACTER (LEN=50) :: cl_asmbkg
      CHARACTER (LEN=50) :: cl_asmdin
      LOGICAL :: llok          ! Check if file exists
      INTEGER :: inum          ! File unit number
      REAL(wp) :: zdate        ! Date
      INTEGER :: ierror
      !!-----------------------------------------------------------------------

      ! If creating an averaged assim bkg, initialise on first timestep
      IF ( ln_avgbkg .AND. kt == ( nn_it000 - 1) ) THEN
         ! Allocate memory 
         ALLOCATE( tn_tavg(jpi,jpj,jpk), STAT=ierror )
         IF( ierror > 0 ) THEN
            CALL ctl_stop( 'asm_wri_bkg: unable to allocate tn_tavg' )   ;   RETURN
         ENDIF
         tn_tavg(:,:,:)=0
         ALLOCATE( sn_tavg(jpi,jpj,jpk), STAT=ierror )
         IF( ierror > 0 ) THEN
            CALL ctl_stop( 'asm_wri_bkg: unable to allocate sn_tavg' )   ;   RETURN
         ENDIF
         sn_tavg(:,:,:)=0
         ALLOCATE( un_tavg(jpi,jpj,jpk), STAT=ierror )
         IF( ierror > 0 ) THEN
            CALL ctl_stop( 'asm_wri_bkg: unable to allocate un_tavg' )   ;   RETURN
         ENDIF
         un_tavg(:,:,:)=0
         ALLOCATE( vn_tavg(jpi,jpj,jpk), STAT=ierror )
         IF( ierror > 0 ) THEN
            CALL ctl_stop( 'asm_wri_bkg: unable to allocate vn_tavg' )   ;   RETURN
         ENDIF
         vn_tavg(:,:,:)=0
         ALLOCATE( sshn_tavg(jpi,jpj), STAT=ierror )
         IF( ierror > 0 ) THEN
            CALL ctl_stop( 'asm_wri_bkg: unable to allocate sshn_tavg' )   ;   RETURN
         ENDIF
         sshn_tavg(:,:)=0
#if defined key_zdftke
         ALLOCATE( en_tavg(jpi,jpj,jpk), STAT=ierror )
         IF( ierror > 0 ) THEN
            CALL ctl_stop( 'asm_wri_bkg: unable to allocate en_tavg' )   ;   RETURN
         ENDIF
         en_tavg(:,:,:)=0
#endif
         ALLOCATE( avt_tavg(jpi,jpj,jpk), STAT=ierror )
         IF( ierror > 0 ) THEN
            CALL ctl_stop( 'asm_wri_bkg: unable to allocate avt_tavg' )   ;   RETURN
         ENDIF
         avt_tavg(:,:,:)=0
         
         numtimes_tavg = REAL ( nitavgbkg_r -  nn_it000 + 1 )
      ENDIF   

      ! If creating an averaged assim bkg, sum the contribution every timestep
      IF ( ln_avgbkg ) THEN 
         IF (lwp) THEN
              WRITE(numout,*) 'asm_wri_bkg : Summing assim bkg fields at timestep ',kt
              WRITE(numout,*) '~~~~~~~~~~~~ '
         ENDIF

         tn_tavg(:,:,:)        = tn_tavg(:,:,:) + tsn(:,:,:,jp_tem) / numtimes_tavg
         sn_tavg(:,:,:)        = sn_tavg(:,:,:) + tsn(:,:,:,jp_sal) / numtimes_tavg
         sshn_tavg(:,:)        = sshn_tavg(:,:) + sshn (:,:) / numtimes_tavg
         un_tavg(:,:,:)        = un_tavg(:,:,:) + un(:,:,:) / numtimes_tavg
         vn_tavg(:,:,:)        = vn_tavg(:,:,:) + vn(:,:,:) / numtimes_tavg
         avt_tavg(:,:,:)        = avt_tavg(:,:,:) + avt(:,:,:) / numtimes_tavg
#if defined key_zdftke
         en_tavg(:,:,:)       = en_tavg(:,:,:) + en(:,:,:) / numtimes_tavg
#endif
      ENDIF
     

      ! Write out background at time step nitbkg_r or nitavgbkg_r
      IF ( ( .NOT. ln_avgbkg .AND. (kt == nitbkg_r) ) .OR. &
      &          ( ln_avgbkg .AND. (kt == nitavgbkg_r) ) ) THEN
         !
         WRITE(cl_asmbkg, FMT='(A,".nc")' ) TRIM( c_asmbkg )
         cl_asmbkg = TRIM( cl_asmbkg )
         INQUIRE( FILE = cl_asmbkg, EXIST = llok )
         !
         IF( .NOT. llok ) THEN
            IF(lwp) WRITE(numout,*) ' Setting up assimilation background file '// TRIM( c_asmbkg )
            !
            !                                      ! Define the output file        
            CALL iom_open( c_asmbkg, inum, ldwrt = .TRUE., kiolib = jprstlib)
            !
            !
            ! Write the information
            IF ( ln_avgbkg ) THEN
               IF( nitavgbkg_r == nit000 - 1 ) THEN      ! Treat special case when nitavgbkg = 0
                  zdate = REAL( ndastp )
#if defined key_zdftke
                  ! lk_zdftke=T :   Read turbulent kinetic energy ( en )
                  IF(lwp) WRITE(numout,*) ' Reading TKE (en) from restart...'
                  CALL tke_rst( nit000, 'READ' )               ! lk_zdftke=T :   Read turbulent kinetic energy ( en )

#endif
               ELSE
                  zdate = REAL( ndastp )
               ENDIF
               CALL iom_rstput( kt, nitavgbkg_r, inum, 'rdastp' , zdate   )
               CALL iom_rstput( kt, nitavgbkg_r, inum, 'un'     , un_tavg )
               CALL iom_rstput( kt, nitavgbkg_r, inum, 'vn'     , vn_tavg )
               CALL iom_rstput( kt, nitavgbkg_r, inum, 'tn'     , tn_tavg )
               CALL iom_rstput( kt, nitavgbkg_r, inum, 'sn'     , sn_tavg )
               CALL iom_rstput( kt, nitavgbkg_r, inum, 'sshn'   , sshn_tavg)
#if defined key_zdftke
               CALL iom_rstput( kt, nitavgbkg_r, inum, 'en'     , en_tavg )
#endif
               CALL iom_rstput( kt, nitavgbkg_r, inum, 'avt'    , avt_tavg)
               !
            ELSE
               IF( nitbkg_r == nit000 - 1 ) THEN      ! Treat special case when nitbkg = 0
                  zdate = REAL( ndastp )
#if defined key_zdftke
                  ! lk_zdftke=T :   Read turbulent kinetic energy ( en )
                  IF(lwp) WRITE(numout,*) ' Reading TKE (en) from restart...'
                  CALL tke_rst( nit000, 'READ' )               ! lk_zdftke=T :   Read turbulent kinetic energy ( en )

#endif
               ELSE
                  zdate = REAL( ndastp )
               ENDIF
               CALL iom_rstput( kt, nitbkg_r, inum, 'rdastp' , zdate   )
               CALL iom_rstput( kt, nitbkg_r, inum, 'un'     , un                )
               CALL iom_rstput( kt, nitbkg_r, inum, 'vn'     , vn                )
               CALL iom_rstput( kt, nitbkg_r, inum, 'tn'     , tsn(:,:,:,jp_tem) )
               CALL iom_rstput( kt, nitbkg_r, inum, 'sn'     , tsn(:,:,:,jp_sal) )
               CALL iom_rstput( kt, nitbkg_r, inum, 'sshn'   , sshn              )
#if defined key_zdftke
               CALL iom_rstput( kt, nitbkg_r, inum, 'en'     , en                )
#endif
               CALL iom_rstput( kt, nitbkg_r, inum, 'avt'    , avt               )
               !
            ENDIF
            
            CALL iom_close( inum )

         ENDIF
         !
      ENDIF

      !                                !-------------------------------------------
      IF( kt == nitdin_r ) THEN        ! Write out background at time step nitdin_r
         !                             !-----------------------------------========
         !
         WRITE(cl_asmdin, FMT='(A,".nc")' ) TRIM( c_asmdin )
         cl_asmdin = TRIM( cl_asmdin )
         INQUIRE( FILE = cl_asmdin, EXIST = llok )
         !
         IF( .NOT. llok ) THEN
            IF(lwp) WRITE(numout,*) ' Setting up assimilation background file '// TRIM( c_asmdin )
            !
            !                                      ! Define the output file        
            CALL iom_open( c_asmdin, inum, ldwrt = .TRUE., kiolib = jprstlib)
            !
            IF( nitdin_r == nit000 - 1 ) THEN      ! Treat special case when nitbkg = 0

               zdate = REAL( ndastp )
            ELSE
               zdate = REAL( ndastp )
            ENDIF
            !
            !                                      ! Write the information
            CALL iom_rstput( kt, nitdin_r, inum, 'rdastp' , zdate             )
            CALL iom_rstput( kt, nitdin_r, inum, 'un'     , un                )
            CALL iom_rstput( kt, nitdin_r, inum, 'vn'     , vn                )
            CALL iom_rstput( kt, nitdin_r, inum, 'tn'     , tsn(:,:,:,jp_tem) )
            CALL iom_rstput( kt, nitdin_r, inum, 'sn'     , tsn(:,:,:,jp_sal) )
            CALL iom_rstput( kt, nitdin_r, inum, 'sshn'   , sshn              )
#if defined key_lim2 || defined key_lim3
            IF(( nn_ice == 2 ) .OR. ( nn_ice == 3 )) THEN
	       IF(ALLOCATED(frld)) THEN
                  CALL iom_rstput( kt, nitdin_r, inum, 'iceconc', 1.0 - frld(:,:)   )
               ELSE
		  CALL ctl_warn('Ice concentration not written to background as ice variable frld not allocated on this timestep')
	       ENDIF
            ENDIF
#endif
            !
            CALL iom_close( inum )
         ENDIF
         !
      ENDIF
      !                    
   END SUBROUTINE asm_bkg_wri

   !!======================================================================
END MODULE asmbkg
