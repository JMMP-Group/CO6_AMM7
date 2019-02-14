MODULE obs_oper
   !!======================================================================
   !!                       ***  MODULE obs_oper  ***
   !! Observation diagnostics: Observation operators for various observation
   !!                          types
   !!======================================================================

   !!----------------------------------------------------------------------
   !!   obs_pro_opt :    Compute the model counterpart of temperature and
   !!                    salinity observations from profiles
   !!   obs_pro_sco_opt: Compute the model counterpart of temperature and 
   !!                    salinity observations from profiles in generalised 
   !!                    vertical coordinates 
   !!   obs_sla_opt :    Compute the model counterpart of sea level anomaly
   !!                    observations
   !!   obs_sst_opt :    Compute the model counterpart of sea surface temperature
   !!                    observations
   !!   obs_sss_opt :    Compute the model counterpart of sea surface salinity
   !!                    observations
   !!   obs_seaice_opt : Compute the model counterpart of sea ice concentration
   !!                    observations
   !!
   !!   obs_vel_opt :    Compute the model counterpart of zonal and meridional
   !!                    components of velocity from observations.
   !!   obs_logchl_opt : Compute the model counterpart of log10(chlorophyll)
   !!                    observations
   !!   obs_spm_opt :    Compute the model counterpart of spm
   !!                    observations
   !!   obs_fco2_opt :   Compute the model counterpart of fco2
   !!                    observations
   !!   obs_pco2_opt :   Compute the model counterpart of pco2
   !!                    observations
   !!----------------------------------------------------------------------

   !! * Modules used   
   USE par_kind, ONLY : &         ! Precision variables
      & wp
   USE in_out_manager             ! I/O manager
   USE obs_inter_sup              ! Interpolation support
   USE obs_inter_h2d, ONLY : &    ! Horizontal interpolation to the observation pt
      & obs_int_h2d, &
      & obs_int_h2d_init
   USE obs_inter_z1d, ONLY : &    ! Vertical interpolation to the observation pt
      & obs_int_z1d,    &
      & obs_int_z1d_spl
   USE obs_const,  ONLY :     &
      & obfillflt		  ! Fillvalue   
   USE dom_oce,       ONLY : &
      & glamt, glamu, glamv, &
      & gphit, gphiu, gphiv, & 
#if defined key_vvl 
      & gdept_n 
#else 
      & gdept_0 
#endif  
   USE lib_mpp,       ONLY : &
      & ctl_warn, ctl_stop
   USE obs_grid,      ONLY : & 
      & obs_level_search     
      
   IMPLICIT NONE

   !! * Routine accessibility
   PRIVATE

   PUBLIC obs_pro_opt, &  ! Compute the model counterpart of profile observations
      &   obs_pro_sco_opt, &  ! Compute the model counterpart of profile observations 
                              ! in generalised vertical coordinates 
      &   obs_sla_opt, &  ! Compute the model counterpart of SLA observations
      &   obs_sst_opt, &  ! Compute the model counterpart of SST observations
      &   obs_sss_opt, &  ! Compute the model counterpart of SSS observations
      &   obs_seaice_opt, &
      &   obs_vel_opt, &  ! Compute the model counterpart of velocity profile data
      &   obs_logchl_opt, & ! Compute the model counterpart of logchl data
      &   obs_spm_opt, &  ! Compute the model counterpart of spm data
      &   obs_fco2_opt, & ! Compute the model counterpart of fco2 data
      &   obs_pco2_opt    ! Compute the model counterpart of pco2 data

   INTEGER, PARAMETER, PUBLIC :: imaxavtypes = 20 ! Max number of daily avgd obs types

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

   !! * Substitutions 
#  include "domzgr_substitute.h90" 
CONTAINS

   SUBROUTINE obs_pro_opt( prodatqc, kt, kpi, kpj, kpk, kit000, kdaystp, &
      &                    ptn, psn, pgdept, ptmask, k1dint, k2dint, &
      &                    kdailyavtypes )
      !!-----------------------------------------------------------------------
      !!
      !!                     ***  ROUTINE obs_pro_opt  ***
      !!
      !! ** Purpose : Compute the model counterpart of profiles
      !!              data by interpolating from the model grid to the 
      !!              observation point.
      !!
      !! ** Method  : Linearly interpolate to each observation point using 
      !!              the model values at the corners of the surrounding grid box.
      !!
      !!    First, a vertical profile of horizontally interpolated model
      !!    now temperatures is computed at the obs (lon, lat) point.
      !!    Several horizontal interpolation schemes are available:
      !!        - distance-weighted (great circle) (k2dint = 0)
      !!        - distance-weighted (small angle)  (k2dint = 1)
      !!        - bilinear (geographical grid)     (k2dint = 2)
      !!        - bilinear (quadrilateral grid)    (k2dint = 3)
      !!        - polynomial (quadrilateral grid)  (k2dint = 4)
      !!
      !!    Next, the vertical temperature profile is interpolated to the
      !!    data depth points. Two vertical interpolation schemes are
      !!    available:
      !!        - linear       (k1dint = 0)
      !!        - Cubic spline (k1dint = 1)
      !!
      !!    For the cubic spline the 2nd derivative of the interpolating 
      !!    polynomial is computed before entering the vertical interpolation 
      !!    routine.
      !!
      !!    For ENACT moored buoy data (e.g., TAO), the model equivalent is
      !!    a daily mean model temperature field. So, we first compute
      !!    the mean, then interpolate only at the end of the day.
      !!
      !!    Note: the in situ temperature observations must be converted
      !!    to potential temperature (the model variable) prior to
      !!    assimilation. 
      !!??????????????????????????????????????????????????????????????
      !!    INCLUDE POTENTIAL TEMP -> IN SITU TEMP IN OBS OPERATOR???
      !!??????????????????????????????????????????????????????????????
      !!
      !! ** Action  :
      !!
      !! History :
      !!      ! 97-11 (A. Weaver, S. Ricci, N. Daget)
      !!      ! 06-03 (G. Smith) NEMOVAR migration
      !!      ! 06-10 (A. Weaver) Cleanup
      !!      ! 07-01 (K. Mogensen) Merge of temperature and salinity
      !!      ! 07-03 (K. Mogensen) General handling of profiles
      !!-----------------------------------------------------------------------
  
      !! * Modules used
      USE obs_profiles_def ! Definition of storage space for profile obs.

      IMPLICIT NONE

      !! * Arguments
      TYPE(obs_prof), INTENT(INOUT) :: prodatqc  ! Subset of profile data not failing screening
      INTEGER, INTENT(IN) :: kt        ! Time step
      INTEGER, INTENT(IN) :: kpi       ! Model grid parameters
      INTEGER, INTENT(IN) :: kpj
      INTEGER, INTENT(IN) :: kpk
      INTEGER, INTENT(IN) :: kit000    ! Number of the first time step 
                                       !   (kit000-1 = restart time)
      INTEGER, INTENT(IN) :: k1dint    ! Vertical interpolation type (see header)
      INTEGER, INTENT(IN) :: k2dint    ! Horizontal interpolation type (see header)
      INTEGER, INTENT(IN) :: kdaystp   ! Number of time steps per day                    
      REAL(KIND=wp), INTENT(IN), DIMENSION(kpi,kpj,kpk) :: &
         & ptn,    &    ! Model temperature field
         & psn,    &    ! Model salinity field
         & ptmask       ! Land-sea mask
      REAL(KIND=wp), INTENT(IN), DIMENSION(kpk) :: &
         & pgdept       ! Model array of depth levels
      INTEGER, DIMENSION(imaxavtypes), OPTIONAL :: &
         & kdailyavtypes! Types for daily averages
      !! * Local declarations
      INTEGER ::   ji
      INTEGER ::   jj
      INTEGER ::   jk
      INTEGER ::   jobs
      INTEGER ::   inrc
      INTEGER ::   ipro
      INTEGER ::   idayend
      INTEGER ::   ista
      INTEGER ::   iend
      INTEGER ::   iobs
      INTEGER, DIMENSION(imaxavtypes) :: &
         & idailyavtypes
      REAL(KIND=wp) :: zlam
      REAL(KIND=wp) :: zphi
      REAL(KIND=wp) :: zdaystp
      REAL(KIND=wp), DIMENSION(kpk) :: &
         & zobsmask, &
         & zobsk,    &
         & zobs2k
      REAL(KIND=wp), DIMENSION(2,2,kpk) :: &
         & zweig
      REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE :: &
         & zmask, &
         & zintt, &
         & zints, &
         & zinmt, &
         & zinms
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: &
         & zglam, &
         & zgphi
      INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: &
         & igrdi, &
         & igrdj

      !------------------------------------------------------------------------
      ! Local initialization 
      !------------------------------------------------------------------------
      ! ... Record and data counters
      inrc = kt - kit000 + 2
      ipro = prodatqc%npstp(inrc)
 
      ! Daily average types
      IF ( PRESENT(kdailyavtypes) ) THEN
         idailyavtypes(:) = kdailyavtypes(:)
      ELSE
         idailyavtypes(:) = -1
      ENDIF

      ! Initialize daily mean for first timestep
      idayend = MOD( kt - kit000 + 1, kdaystp )

      ! Added kt == 0 test to catch restart case 
      IF ( idayend == 1 .OR. kt == 0) THEN
         IF (lwp) WRITE(numout,*) 'Reset prodatqc%vdmean on time-step: ',kt
         DO jk = 1, jpk
            DO jj = 1, jpj
               DO ji = 1, jpi
                  prodatqc%vdmean(ji,jj,jk,1) = 0.0
                  prodatqc%vdmean(ji,jj,jk,2) = 0.0
               END DO
            END DO
         END DO
      ENDIF

      DO jk = 1, jpk
         DO jj = 1, jpj
            DO ji = 1, jpi
               ! Increment the temperature field for computing daily mean
               prodatqc%vdmean(ji,jj,jk,1) = prodatqc%vdmean(ji,jj,jk,1) &
                  &                        + ptn(ji,jj,jk)
               ! Increment the salinity field for computing daily mean
               prodatqc%vdmean(ji,jj,jk,2) = prodatqc%vdmean(ji,jj,jk,2) &
                  &                        + psn(ji,jj,jk)
            END DO
         END DO
      END DO
   
      ! Compute the daily mean at the end of day
      zdaystp = 1.0 / REAL( kdaystp )
      IF ( idayend == 0 ) THEN
         DO jk = 1, jpk
            DO jj = 1, jpj
               DO ji = 1, jpi
                  prodatqc%vdmean(ji,jj,jk,1) = prodatqc%vdmean(ji,jj,jk,1) &
                     &                        * zdaystp
                  prodatqc%vdmean(ji,jj,jk,2) = prodatqc%vdmean(ji,jj,jk,2) &
                  &                           * zdaystp
               END DO
            END DO
         END DO
      ENDIF

      ! Get the data for interpolation
      ALLOCATE( &
         & igrdi(2,2,ipro),      &
         & igrdj(2,2,ipro),      &
         & zglam(2,2,ipro),      &
         & zgphi(2,2,ipro),      &
         & zmask(2,2,kpk,ipro),  &
         & zintt(2,2,kpk,ipro),  &
         & zints(2,2,kpk,ipro)   &
         & )

      DO jobs = prodatqc%nprofup + 1, prodatqc%nprofup + ipro
         iobs = jobs - prodatqc%nprofup
         igrdi(1,1,iobs) = prodatqc%mi(jobs,1)-1
         igrdj(1,1,iobs) = prodatqc%mj(jobs,1)-1
         igrdi(1,2,iobs) = prodatqc%mi(jobs,1)-1
         igrdj(1,2,iobs) = prodatqc%mj(jobs,1)
         igrdi(2,1,iobs) = prodatqc%mi(jobs,1)
         igrdj(2,1,iobs) = prodatqc%mj(jobs,1)-1
         igrdi(2,2,iobs) = prodatqc%mi(jobs,1)
         igrdj(2,2,iobs) = prodatqc%mj(jobs,1)
      END DO

      CALL obs_int_comm_2d( 2, 2, ipro, igrdi, igrdj, glamt, zglam )
      CALL obs_int_comm_2d( 2, 2, ipro, igrdi, igrdj, gphit, zgphi )
      CALL obs_int_comm_3d( 2, 2, ipro, kpk, igrdi, igrdj, ptmask,zmask )
      CALL obs_int_comm_3d( 2, 2, ipro, kpk, igrdi, igrdj, ptn,   zintt )
      CALL obs_int_comm_3d( 2, 2, ipro, kpk, igrdi, igrdj, psn,   zints )

      ! At the end of the day also get interpolated means
      IF ( idayend == 0 ) THEN

         ALLOCATE( &
            & zinmt(2,2,kpk,ipro),  &
            & zinms(2,2,kpk,ipro)   &
            & )

         CALL obs_int_comm_3d( 2, 2, ipro, kpk, igrdi, igrdj, &
            &                  prodatqc%vdmean(:,:,:,1), zinmt )
         CALL obs_int_comm_3d( 2, 2, ipro, kpk, igrdi, igrdj, &
            &                  prodatqc%vdmean(:,:,:,2), zinms )

      ENDIF

      DO jobs = prodatqc%nprofup + 1, prodatqc%nprofup + ipro

         iobs = jobs - prodatqc%nprofup

         IF ( kt /= prodatqc%mstp(jobs) ) THEN
            
            IF(lwp) THEN
               WRITE(numout,*)
               WRITE(numout,*) ' E R R O R : Observation',              &
                  &            ' time step is not consistent with the', &
                  &            ' model time step'
               WRITE(numout,*) ' ========='
               WRITE(numout,*)
               WRITE(numout,*) ' Record  = ', jobs,                    &
                  &            ' kt      = ', kt,                      &
                  &            ' mstp    = ', prodatqc%mstp(jobs), &
                  &            ' ntyp    = ', prodatqc%ntyp(jobs)
            ENDIF
            CALL ctl_stop( 'obs_pro_opt', 'Inconsistent time' )
         ENDIF
         
         zlam = prodatqc%rlam(jobs)
         zphi = prodatqc%rphi(jobs)
         
         ! Horizontal weights and vertical mask

         IF ( ( prodatqc%npvend(jobs,1) > 0 ) .OR. &
            & ( prodatqc%npvend(jobs,2) > 0 ) ) THEN

            CALL obs_int_h2d_init( kpk, kpk, k2dint, zlam, zphi,     &
               &                   zglam(:,:,iobs), zgphi(:,:,iobs), &
               &                   zmask(:,:,:,iobs), zweig, zobsmask )

         ENDIF

         IF ( prodatqc%npvend(jobs,1) > 0 ) THEN

            zobsk(:) = obfillflt

	    IF ( ANY (idailyavtypes(:) == prodatqc%ntyp(jobs)) ) THEN

               IF ( idayend == 0 )  THEN
                  
                  ! Daily averaged moored buoy (MRB) data
                  
                  CALL obs_int_h2d( kpk, kpk,      &
                     &              zweig, zinmt(:,:,:,iobs), zobsk )
                  
                  
               ELSE
               
                  CALL ctl_stop( ' A nonzero' //     &
                     &           ' number of profile T BUOY data should' // &
                     &           ' only occur at the end of a given day' )

               ENDIF
	       
            ELSE 
               
               ! Point data

               CALL obs_int_h2d( kpk, kpk,      &
                  &              zweig, zintt(:,:,:,iobs), zobsk )

            ENDIF

            !-------------------------------------------------------------
            ! Compute vertical second-derivative of the interpolating 
            ! polynomial at obs points
            !-------------------------------------------------------------
            
            IF ( k1dint == 1 ) THEN
               CALL obs_int_z1d_spl( kpk, zobsk, zobs2k,   &
                  &                  pgdept, zobsmask )
            ENDIF
            
            !-----------------------------------------------------------------
            !  Vertical interpolation to the observation point
            !-----------------------------------------------------------------
            ista = prodatqc%npvsta(jobs,1)
            iend = prodatqc%npvend(jobs,1)
            CALL obs_int_z1d( kpk,                &
               & prodatqc%var(1)%mvk(ista:iend),  &
               & k1dint, iend - ista + 1,         &
               & prodatqc%var(1)%vdep(ista:iend), &
               & zobsk, zobs2k,                   &
               & prodatqc%var(1)%vmod(ista:iend), &
               & pgdept, zobsmask )

         ENDIF

         IF ( prodatqc%npvend(jobs,2) > 0 ) THEN

            zobsk(:) = obfillflt

            IF ( ANY (idailyavtypes(:) == prodatqc%ntyp(jobs)) ) THEN

               IF ( idayend == 0 )  THEN

                  ! Daily averaged moored buoy (MRB) data
                  
                  CALL obs_int_h2d( kpk, kpk,      &
                     &              zweig, zinms(:,:,:,iobs), zobsk )
                  
               ELSE

                  CALL ctl_stop( ' A nonzero' //     &
                     &           ' number of profile S BUOY data should' // &
                     &           ' only occur at the end of a given day' )

               ENDIF

            ELSE
               
               ! Point data

               CALL obs_int_h2d( kpk, kpk,      &
                  &              zweig, zints(:,:,:,iobs), zobsk )

            ENDIF


            !-------------------------------------------------------------
            ! Compute vertical second-derivative of the interpolating 
            ! polynomial at obs points
            !-------------------------------------------------------------
            
            IF ( k1dint == 1 ) THEN
               CALL obs_int_z1d_spl( kpk, zobsk, zobs2k, &
                  &                  pgdept, zobsmask )
            ENDIF
            
            !----------------------------------------------------------------
            !  Vertical interpolation to the observation point
            !----------------------------------------------------------------
            ista = prodatqc%npvsta(jobs,2)
            iend = prodatqc%npvend(jobs,2)
            CALL obs_int_z1d( kpk, &
               & prodatqc%var(2)%mvk(ista:iend),&
               & k1dint, iend - ista + 1, &
               & prodatqc%var(2)%vdep(ista:iend),&
               & zobsk, zobs2k, &
               & prodatqc%var(2)%vmod(ista:iend),&
               & pgdept, zobsmask )

         ENDIF

      END DO
 
      ! Deallocate the data for interpolation
      DEALLOCATE( &
         & igrdi, &
         & igrdj, &
         & zglam, &
         & zgphi, &
         & zmask, &
         & zintt, &
         & zints  &
         & )
      ! At the end of the day also get interpolated means
      IF ( idayend == 0 ) THEN
         DEALLOCATE( &
            & zinmt,  &
            & zinms   &
            & )
      ENDIF

      prodatqc%nprofup = prodatqc%nprofup + ipro 
      
   END SUBROUTINE obs_pro_opt

   SUBROUTINE obs_pro_sco_opt( prodatqc, kt, kpi, kpj, kpk, kit000, kdaystp, & 
      &                    ptn, psn, pgdept, pgdepw, ptmask, k1dint, k2dint, & 
      &                    kdailyavtypes ) 
      !!----------------------------------------------------------------------- 
      !! 
      !!                     ***  ROUTINE obs_pro_opt  *** 
      !! 
      !! ** Purpose : Compute the model counterpart of profiles 
      !!              data by interpolating from the model grid to the  
      !!              observation point. Generalised vertical coordinate version 
      !! 
      !! ** Method  : Linearly interpolate to each observation point using  
      !!              the model values at the corners of the surrounding grid box. 
      !! 
      !!          First, model values on the model grid are interpolated vertically to the 
      !!          Depths of the profile observations.  Two vertical interpolation schemes are 
      !!          available: 
      !!          - linear       (k1dint = 0) 
      !!          - Cubic spline (k1dint = 1)    
      !! 
      !! 
      !!         Secondly the interpolated values are interpolated horizontally to the  
      !!         obs (lon, lat) point. 
      !!         Several horizontal interpolation schemes are available: 
      !!        - distance-weighted (great circle) (k2dint = 0) 
      !!        - distance-weighted (small angle)  (k2dint = 1) 
      !!        - bilinear (geographical grid)     (k2dint = 2) 
      !!        - bilinear (quadrilateral grid)    (k2dint = 3) 
      !!        - polynomial (quadrilateral grid)  (k2dint = 4) 
      !! 
      !!    For the cubic spline the 2nd derivative of the interpolating  
      !!    polynomial is computed before entering the vertical interpolation  
      !!    routine. 
      !! 
      !!    For ENACT moored buoy data (e.g., TAO), the model equivalent is 
      !!    a daily mean model temperature field. So, we first compute 
      !!    the mean, then interpolate only at the end of the day. 
      !! 
      !!    This is the procedure to be used with generalised vertical model  
      !!    coordinates (ie s-coordinates. It is ~4x slower than the equivalent 
      !!    horizontal then vertical interpolation algorithm, but can deal with situations 
      !!    where the model levels are not flat. 
      !!    ONLY PERFORMED if ln_sco=.TRUE.  
      !!       
      !!    Note: the in situ temperature observations must be converted 
      !!    to potential temperature (the model variable) prior to 
      !!    assimilation.  
      !!?????????????????????????????????????????????????????????????? 
      !!    INCLUDE POTENTIAL TEMP -> IN SITU TEMP IN OBS OPERATOR??? 
      !!?????????????????????????????????????????????????????????????? 
      !! 
      !! ** Action  : 
      !! 
      !! History : 
      !!      ! 2014-08 (J. While) Adapted from obs_pro_opt to handel generalised 
      !!                           vertical coordinates
      !!----------------------------------------------------------------------- 
   
      !! * Modules used 
      USE obs_profiles_def   ! Definition of storage space for profile obs. 
      USE dom_oce,  ONLY : & 
#if defined key_vvl  
      &   gdepw_n 
#else 
      &   gdepw_0 
#endif 
       
      IMPLICIT NONE 
 
      !! * Arguments 
      TYPE(obs_prof), INTENT(INOUT) :: prodatqc   ! Subset of profile data not failing screening 
      INTEGER, INTENT(IN) :: kt        ! Time step 
      INTEGER, INTENT(IN) :: kpi       ! Model grid parameters 
      INTEGER, INTENT(IN) :: kpj 
      INTEGER, INTENT(IN) :: kpk 
      INTEGER, INTENT(IN) :: kit000    ! Number of the first time step  
                                       !   (kit000-1 = restart time) 
      INTEGER, INTENT(IN) :: k1dint    ! Vertical interpolation type (see header) 
      INTEGER, INTENT(IN) :: k2dint    ! Horizontal interpolation type (see header) 
      INTEGER, INTENT(IN) :: kdaystp   ! Number of time steps per day                     
      REAL(KIND=wp), INTENT(IN), DIMENSION(kpi,kpj,kpk) :: & 
         & ptn,    &    ! Model temperature field 
         & psn,    &    ! Model salinity field 
         & ptmask       ! Land-sea mask 
      REAL(KIND=wp), INTENT(IN), DIMENSION(kpi,kpj,kpk) :: & 
         & pgdept, &       ! Model array of depth T levels 
         & pgdepw       ! Model array of depth W levels 
      INTEGER, DIMENSION(imaxavtypes), OPTIONAL :: & 
         & kdailyavtypes   ! Types for daily averages 
      
      !! * Local declarations 
      INTEGER ::   ji 
      INTEGER ::   jj 
      INTEGER ::   jk 
      INTEGER ::   iico, ijco 
      INTEGER ::   jobs 
      INTEGER ::   inrc 
      INTEGER ::   ipro 
      INTEGER ::   idayend 
      INTEGER ::   ista 
      INTEGER ::   iend 
      INTEGER ::   iobs 
      INTEGER ::   iin, ijn, ikn, ik   ! looping indices over interpolation nodes 
      INTEGER, DIMENSION(imaxavtypes) :: & 
         & idailyavtypes 
      INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: & 
         & igrdi, & 
         & igrdj 
      INTEGER :: & 
         & inum_obs
      INTEGER, ALLOCATABLE, DIMENSION(:) :: iv_indic    
      REAL(KIND=wp) :: zlam 
      REAL(KIND=wp) :: zphi 
      REAL(KIND=wp) :: zdaystp 
      REAL(KIND=wp), DIMENSION(kpk) :: & 
         & zobsmask, & 
         & zobsk,    & 
         & zobs2k 
      REAL(KIND=wp), DIMENSION(2,2,1) :: & 
         & zweig, & 
         & l_zweig 
      REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE :: & 
         & zmask, & 
         & zintt, & 
         & zints, & 
         & zinmt, & 
         & zgdept,& 
         & zgdepw,& 
         & zinms 
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: & 
         & zglam, & 
         & zgphi    
      REAL(KIND=wp), DIMENSION(1) :: zmsk_1       
      REAL(KIND=wp), DIMENSION(:,:,:), ALLOCATABLE :: interp_corner       
 
      !------------------------------------------------------------------------ 
      ! Local initialization  
      !------------------------------------------------------------------------ 
      ! ... Record and data counters 
      inrc = kt - kit000 + 2 
      ipro = prodatqc%npstp(inrc) 
  
      ! Daily average types 
      IF ( PRESENT(kdailyavtypes) ) THEN 
         idailyavtypes(:) = kdailyavtypes(:) 
      ELSE 
         idailyavtypes(:) = -1 
      ENDIF 
 
      ! Initialize daily mean for first time-step 
      idayend = MOD( kt - kit000 + 1, kdaystp ) 
 
      ! Added kt == 0 test to catch restart case  
      IF ( idayend == 1 .OR. kt == 0) THEN 
          
         IF (lwp) WRITE(numout,*) 'Reset prodatqc%vdmean on time-step: ',kt 
         DO jk = 1, jpk 
            DO jj = 1, jpj 
               DO ji = 1, jpi 
                  prodatqc%vdmean(ji,jj,jk,1) = 0.0 
                  prodatqc%vdmean(ji,jj,jk,2) = 0.0 
               END DO 
            END DO 
         END DO 
       
      ENDIF 
       
      DO jk = 1, jpk 
         DO jj = 1, jpj 
            DO ji = 1, jpi 
               ! Increment the temperature field for computing daily mean 
               prodatqc%vdmean(ji,jj,jk,1) = prodatqc%vdmean(ji,jj,jk,1) & 
               &                        + ptn(ji,jj,jk) 
               ! Increment the salinity field for computing daily mean 
               prodatqc%vdmean(ji,jj,jk,2) = prodatqc%vdmean(ji,jj,jk,2) & 
               &                        + psn(ji,jj,jk) 
            END DO 
         END DO 
      END DO 
    
      ! Compute the daily mean at the end of day 
      zdaystp = 1.0 / REAL( kdaystp ) 
      IF ( idayend == 0 ) THEN 
         DO jk = 1, jpk 
            DO jj = 1, jpj 
               DO ji = 1, jpi 
                  prodatqc%vdmean(ji,jj,jk,1) = prodatqc%vdmean(ji,jj,jk,1) & 
                  &                        * zdaystp 
                  prodatqc%vdmean(ji,jj,jk,2) = prodatqc%vdmean(ji,jj,jk,2) & 
                  &                           * zdaystp 
               END DO 
            END DO 
         END DO 
      ENDIF 
 
      ! Get the data for interpolation 
      ALLOCATE( & 
         & igrdi(2,2,ipro),      & 
         & igrdj(2,2,ipro),      & 
         & zglam(2,2,ipro),      & 
         & zgphi(2,2,ipro),      & 
         & zmask(2,2,kpk,ipro),  & 
         & zintt(2,2,kpk,ipro),  & 
         & zints(2,2,kpk,ipro),  & 
         & zgdept(2,2,kpk,ipro), & 
         & zgdepw(2,2,kpk,ipro)  & 
         & ) 
 
      DO jobs = prodatqc%nprofup + 1, prodatqc%nprofup + ipro 
         iobs = jobs - prodatqc%nprofup 
         igrdi(1,1,iobs) = prodatqc%mi(jobs,1)-1 
         igrdj(1,1,iobs) = prodatqc%mj(jobs,1)-1 
         igrdi(1,2,iobs) = prodatqc%mi(jobs,1)-1 
         igrdj(1,2,iobs) = prodatqc%mj(jobs,1) 
         igrdi(2,1,iobs) = prodatqc%mi(jobs,1) 
         igrdj(2,1,iobs) = prodatqc%mj(jobs,1)-1 
         igrdi(2,2,iobs) = prodatqc%mi(jobs,1) 
         igrdj(2,2,iobs) = prodatqc%mj(jobs,1) 
      END DO 
 
      ! Initiialise depth arrays
      zgdept = 0.0
      zgdepw = 0.0
      
      CALL obs_int_comm_2d( 2, 2, ipro, igrdi, igrdj, glamt, zglam ) 
      CALL obs_int_comm_2d( 2, 2, ipro, igrdi, igrdj, gphit, zgphi ) 
      CALL obs_int_comm_3d( 2, 2, ipro, kpk, igrdi, igrdj, ptmask,zmask ) 
      CALL obs_int_comm_3d( 2, 2, ipro, kpk, igrdi, igrdj, ptn,   zintt ) 
      CALL obs_int_comm_3d( 2, 2, ipro, kpk, igrdi, igrdj, psn,   zints ) 
      CALL obs_int_comm_3d( 2, 2, ipro, kpk, igrdi, igrdj, pgdept(:,:,:), & 
        &                     zgdept ) 
      CALL obs_int_comm_3d( 2, 2, ipro, kpk, igrdi, igrdj, pgdepw(:,:,:), & 
        &                     zgdepw ) 
 
      ! At the end of the day also get interpolated means 
      IF ( idayend == 0 ) THEN 
 
         ALLOCATE( & 
            & zinmt(2,2,kpk,ipro),  & 
            & zinms(2,2,kpk,ipro)   & 
            & ) 
 
         CALL obs_int_comm_3d( 2, 2, ipro, kpk, igrdi, igrdj, & 
            &                  prodatqc%vdmean(:,:,:,1), zinmt ) 
         CALL obs_int_comm_3d( 2, 2, ipro, kpk, igrdi, igrdj, & 
            &                  prodatqc%vdmean(:,:,:,2), zinms ) 
 
      ENDIF 
       
      ! Return if no observations to process 
      ! Has to be done after comm commands to ensure processors 
      ! stay in sync 
      IF ( ipro == 0 ) RETURN 
 
      DO jobs = prodatqc%nprofup + 1, prodatqc%nprofup + ipro 
    
         iobs = jobs - prodatqc%nprofup 
    
         IF ( kt /= prodatqc%mstp(jobs) ) THEN 
             
            IF(lwp) THEN 
               WRITE(numout,*) 
               WRITE(numout,*) ' E R R O R : Observation',              & 
                  &            ' time step is not consistent with the', & 
                  &            ' model time step' 
               WRITE(numout,*) ' =========' 
               WRITE(numout,*) 
               WRITE(numout,*) ' Record  = ', jobs,                    & 
                  &            ' kt      = ', kt,                      & 
                  &            ' mstp    = ', prodatqc%mstp(jobs), & 
                  &            ' ntyp    = ', prodatqc%ntyp(jobs) 
            ENDIF 
            CALL ctl_stop( 'obs_pro_opt', 'Inconsistent time' ) 
         ENDIF 
          
         zlam = prodatqc%rlam(jobs) 
         zphi = prodatqc%rphi(jobs) 
          
         ! Horizontal weights 
         ! Only calculated once, for both T and S. 
         ! Masked values are calculated later.  
 
         IF ( ( prodatqc%npvend(jobs,1) > 0 ) .OR. & 
            & ( prodatqc%npvend(jobs,2) > 0 ) ) THEN 
 
            CALL obs_int_h2d_init( 1, 1, k2dint, zlam, zphi,     & 
               &                   zglam(:,:,iobs), zgphi(:,:,iobs), & 
               &                   zmask(:,:,1,iobs), zweig, zmsk_1 ) 
 
         ENDIF 
         
         ! IF zmsk_1 = 0; then ob is on land 
         IF (zmsk_1(1) < 0.1) THEN 
            WRITE(numout,*) 'WARNING (obs_oper) :- profile found within landmask' 
   
         ELSE  
             
            ! Temperature 
             
            IF ( prodatqc%npvend(jobs,1) > 0 ) THEN  
    
               zobsk(:) = obfillflt 
    
               IF ( ANY (idailyavtypes(:) == prodatqc%ntyp(jobs)) ) THEN 
    
                  IF ( idayend == 0 )  THEN 
                   
                     ! Daily averaged moored buoy (MRB) data 
                   
                     ! vertically interpolate all 4 corners 
                     ista = prodatqc%npvsta(jobs,1) 
                     iend = prodatqc%npvend(jobs,1) 
                     inum_obs = iend - ista + 1 
                     ALLOCATE(interp_corner(2,2,inum_obs),iv_indic(inum_obs)) 
      
                     DO iin=1,2 
                        DO ijn=1,2 
                                       
                                       
           
                           IF ( k1dint == 1 ) THEN 
                              CALL obs_int_z1d_spl( kpk, & 
                                 &     zinmt(iin,ijn,:,iobs), & 
                                 &     zobs2k, zgdept(iin,ijn,:,iobs), & 
                                 &     zmask(iin,ijn,:,iobs)) 
                           ENDIF 
       
                           CALL obs_level_search(kpk, & 
                              &    zgdept(iin,ijn,:,iobs), & 
                              &    inum_obs, prodatqc%var(1)%vdep(ista:iend), & 
                              &    iv_indic) 
                           CALL obs_int_z1d(kpk, iv_indic, k1dint, inum_obs, & 
                              &    prodatqc%var(1)%vdep(ista:iend), & 
                              &    zinmt(iin,ijn,:,iobs), & 
                              &    zobs2k, interp_corner(iin,ijn,:), & 
                              &    zgdept(iin,ijn,:,iobs), & 
                              &    zmask(iin,ijn,:,iobs)) 
       
                        ENDDO 
                     ENDDO 
                   
                   
                  ELSE 
                
                     CALL ctl_stop( ' A nonzero' //     & 
                        &           ' number of profile T BUOY data should' // & 
                        &           ' only occur at the end of a given day' ) 
    
                  ENDIF 
         
               ELSE  
                
                  ! Point data 
     
                  ! vertically interpolate all 4 corners 
                  ista = prodatqc%npvsta(jobs,1) 
                  iend = prodatqc%npvend(jobs,1) 
                  inum_obs = iend - ista + 1 
                  ALLOCATE(interp_corner(2,2,inum_obs), iv_indic(inum_obs)) 
                  DO iin=1,2  
                     DO ijn=1,2 
                                    
                                    
                        IF ( k1dint == 1 ) THEN 
                           CALL obs_int_z1d_spl( kpk, & 
                              &    zintt(iin,ijn,:,iobs),& 
                              &    zobs2k, zgdept(iin,ijn,:,iobs), & 
                              &    zmask(iin,ijn,:,iobs)) 
  
                        ENDIF 
       
                        CALL obs_level_search(kpk, & 
                            &        zgdept(iin,ijn,:,iobs),& 
                            &        inum_obs, prodatqc%var(1)%vdep(ista:iend), & 
                            &         iv_indic) 
                        CALL obs_int_z1d(kpk, iv_indic, k1dint, inum_obs,     & 
                            &          prodatqc%var(1)%vdep(ista:iend),     & 
                            &          zintt(iin,ijn,:,iobs),            & 
                            &          zobs2k,interp_corner(iin,ijn,:), & 
                            &          zgdept(iin,ijn,:,iobs),         & 
                            &          zmask(iin,ijn,:,iobs) )      
         
                     ENDDO 
                  ENDDO 
             
               ENDIF 
       
               !------------------------------------------------------------- 
               ! Compute the horizontal interpolation for every profile level 
               !------------------------------------------------------------- 
             
               DO ikn=1,inum_obs 
                  iend=ista+ikn-1
                  
                  l_zweig(:,:,1) = 0._wp 
   
                  ! This code forces the horizontal weights to be  
                  ! zero IF the observation is below the bottom of the  
                  ! corners of the interpolation nodes, Or if it is in  
                  ! the mask. This is important for observations are near  
                  ! steep bathymetry 
                  DO iin=1,2 
                     DO ijn=1,2 
     
                        depth_loop1: DO ik=kpk,2,-1 
                           IF(zmask(iin,ijn,ik-1,iobs ) > 0.9 )THEN   
                            
                              l_zweig(iin,ijn,1) = &  
                                 & zweig(iin,ijn,1) * & 
                                 & MAX( SIGN(1._wp,(zgdepw(iin,ijn,ik,iobs) ) & 
                                 &  - prodatqc%var(1)%vdep(iend)),0._wp) 
                            
                              EXIT depth_loop1 
                           ENDIF 
                        ENDDO depth_loop1 
     
                     ENDDO 
                  ENDDO 
   
                  CALL obs_int_h2d( 1, 1, l_zweig, interp_corner(:,:,ikn), & 
                  &          prodatqc%var(1)%vmod(iend:iend) ) 

                  ! Set QC flag for any observations found below the bottom
                  ! needed as the check here is more strict than that in obs_prep
                  IF (sum(l_zweig) == 0.0_wp) prodatqc%var(1)%nvqc(iend:iend)=4
 
               ENDDO 
 
 
               DEALLOCATE(interp_corner,iv_indic) 
          
            ENDIF 
       
 
            ! Salinity  
          
            IF ( prodatqc%npvend(jobs,2) > 0 ) THEN  
    
               zobsk(:) = obfillflt 
    
               IF ( ANY (idailyavtypes(:) == prodatqc%ntyp(jobs)) ) THEN 
    
                  IF ( idayend == 0 )  THEN 
                   
                     ! Daily averaged moored buoy (MRB) data 
                   
                     ! vertically interpolate all 4 corners 
                     ista = prodatqc%npvsta(jobs,2) 
                     iend = prodatqc%npvend(jobs,2) 
                     inum_obs = iend - ista + 1 
                     ALLOCATE(interp_corner(2,2,inum_obs),iv_indic(inum_obs)) 
      
                     DO iin=1,2 
                        DO ijn=1,2 
                                       
                                       
           
                           IF ( k1dint == 1 ) THEN 
                              CALL obs_int_z1d_spl( kpk, & 
                                 &     zinms(iin,ijn,:,iobs), & 
                                 &     zobs2k, zgdept(iin,ijn,:,iobs), & 
                                 &     zmask(iin,ijn,:,iobs)) 
                           ENDIF 
       
                           CALL obs_level_search(kpk, & 
                              &    zgdept(iin,ijn,:,iobs), & 
                              &    inum_obs, prodatqc%var(2)%vdep(ista:iend), & 
                              &    iv_indic) 
                           CALL obs_int_z1d(kpk, iv_indic, k1dint, inum_obs, & 
                              &    prodatqc%var(2)%vdep(ista:iend), & 
                              &    zinms(iin,ijn,:,iobs), & 
                              &    zobs2k, interp_corner(iin,ijn,:), & 
                              &    zgdept(iin,ijn,:,iobs), & 
                              &    zmask(iin,ijn,:,iobs)) 
       
                        ENDDO 
                     ENDDO 
                   
                   
                  ELSE 
                
                     CALL ctl_stop( ' A nonzero' //     & 
                        &           ' number of profile T BUOY data should' // & 
                        &           ' only occur at the end of a given day' ) 
    
                  ENDIF 
         
               ELSE  
                
                  ! Point data 
     
                  ! vertically interpolate all 4 corners 
                  ista = prodatqc%npvsta(jobs,2) 
                  iend = prodatqc%npvend(jobs,2) 
                  inum_obs = iend - ista + 1 
                  ALLOCATE(interp_corner(2,2,inum_obs), iv_indic(inum_obs)) 
                   
                  DO iin=1,2     
                     DO ijn=1,2  
                                 
                                 
                        IF ( k1dint == 1 ) THEN 
                           CALL obs_int_z1d_spl( kpk, & 
                              &    zints(iin,ijn,:,iobs),& 
                              &    zobs2k, zgdept(iin,ijn,:,iobs), & 
                              &    zmask(iin,ijn,:,iobs)) 
  
                        ENDIF 
       
                        CALL obs_level_search(kpk, & 
                           &        zgdept(iin,ijn,:,iobs),& 
                           &        inum_obs, prodatqc%var(2)%vdep(ista:iend), & 
                           &         iv_indic) 
                        CALL obs_int_z1d(kpk, iv_indic, k1dint, inum_obs,  & 
                           &          prodatqc%var(2)%vdep(ista:iend),     & 
                           &          zints(iin,ijn,:,iobs),               & 
                           &          zobs2k,interp_corner(iin,ijn,:),     & 
                           &          zgdept(iin,ijn,:,iobs),              & 
                           &          zmask(iin,ijn,:,iobs) )      
         
                     ENDDO 
                  ENDDO 
             
               ENDIF 
       
               !------------------------------------------------------------- 
               ! Compute the horizontal interpolation for every profile level 
               !------------------------------------------------------------- 
             
               DO ikn=1,inum_obs 
                  iend=ista+ikn-1
                  
                  l_zweig(:,:,1) = 0._wp 
   
                  ! This code forces the horizontal weights to be  
                  ! zero IF the observation is below the bottom of the  
                  ! corners of the interpolation nodes, Or if it is in  
                  ! the mask. This is important for observations are near  
                  ! steep bathymetry 
                  DO iin=1,2 
                     DO ijn=1,2 
     
                        depth_loop2: DO ik=kpk,2,-1 
                           IF(zmask(iin,ijn,ik-1,iobs ) > 0.9 )THEN   
                            
                              l_zweig(iin,ijn,1) = &  
                                 &  zweig(iin,ijn,1) * & 
                                 &  MAX( SIGN(1._wp,(zgdepw(iin,ijn,ik,iobs) ) & 
                                 &  - prodatqc%var(2)%vdep(iend)),0._wp) 
                            
                              EXIT depth_loop2 
                           ENDIF 
                        ENDDO depth_loop2 
     
                     ENDDO 
                  ENDDO 
   
                  CALL obs_int_h2d( 1, 1, l_zweig, interp_corner(:,:,ikn), & 
                  &          prodatqc%var(2)%vmod(iend:iend) ) 

                  ! Set QC flag for any observations found below the bottom
                  ! needed as the check here is more strict than that in obs_prep
                  IF (sum(l_zweig) == 0.0_wp) prodatqc%var(2)%nvqc(iend:iend)=4
 
               ENDDO 
 
 
               DEALLOCATE(interp_corner,iv_indic) 
          
            ENDIF 
          
         ENDIF 
       
      END DO 
     
      ! Deallocate the data for interpolation 
      DEALLOCATE( & 
         & igrdi, & 
         & igrdj, & 
         & zglam, & 
         & zgphi, & 
         & zmask, & 
         & zintt, & 
         & zints, & 
         & zgdept,& 
         & zgdepw & 
         & ) 
      ! At the end of the day also get interpolated means 
      IF ( idayend == 0 ) THEN 
         DEALLOCATE( & 
            & zinmt,  & 
            & zinms   & 
            & ) 
      ENDIF 
    
      prodatqc%nprofup = prodatqc%nprofup + ipro  
       
   END SUBROUTINE obs_pro_sco_opt 
 
   SUBROUTINE obs_sla_opt( sladatqc, kt, kpi, kpj, kit000, &
      &                    psshn, psshmask, k2dint )
      !!-----------------------------------------------------------------------
      !!
      !!                     ***  ROUTINE obs_sla_opt  ***
      !!
      !! ** Purpose : Compute the model counterpart of sea level anomaly
      !!              data by interpolating from the model grid to the 
      !!              observation point.
      !!
      !! ** Method  : Linearly interpolate to each observation point using 
      !!              the model values at the corners of the surrounding grid box.
      !!
      !!    The now model SSH is first computed at the obs (lon, lat) point.
      !!
      !!    Several horizontal interpolation schemes are available:
      !!        - distance-weighted (great circle) (k2dint = 0)
      !!        - distance-weighted (small angle)  (k2dint = 1)
      !!        - bilinear (geographical grid)     (k2dint = 2)
      !!        - bilinear (quadrilateral grid)    (k2dint = 3)
      !!        - polynomial (quadrilateral grid)  (k2dint = 4)
      !!  
      !!    The sea level anomaly at the observation points is then computed 
      !!    by removing a mean dynamic topography (defined at the obs. point).
      !!
      !! ** Action  :
      !!
      !! History :
      !!      ! 07-03 (A. Weaver)
      !!-----------------------------------------------------------------------
  
      !! * Modules used
      USE obs_surf_def  ! Definition of storage space for surface observations

      IMPLICIT NONE

      !! * Arguments
      TYPE(obs_surf), INTENT(INOUT) :: sladatqc     ! Subset of surface data not failing screening
      INTEGER, INTENT(IN) :: kt      ! Time step
      INTEGER, INTENT(IN) :: kpi     ! Model grid parameters
      INTEGER, INTENT(IN) :: kpj
      INTEGER, INTENT(IN) :: kit000   ! Number of the first time step 
                                      !   (kit000-1 = restart time)
      INTEGER, INTENT(IN) :: k2dint   ! Horizontal interpolation type (see header)
      REAL(KIND=wp), INTENT(IN), DIMENSION(kpi,kpj) :: &
         & psshn,  &    ! Model SSH field
         & psshmask     ! Land-sea mask
         
      !! * Local declarations
      INTEGER :: ji
      INTEGER :: jj
      INTEGER :: jobs
      INTEGER :: inrc
      INTEGER :: isla
      INTEGER :: iobs
      REAL(KIND=wp) :: zlam
      REAL(KIND=wp) :: zphi
      REAL(KIND=wp) :: zext(1), zobsmask(1)
      REAL(kind=wp), DIMENSION(2,2,1) :: &
         & zweig
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: &
         & zmask, &
         & zsshl, &
         & zglam, &
         & zgphi
      INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: &
         & igrdi, &
         & igrdj

      !------------------------------------------------------------------------
      ! Local initialization 
      !------------------------------------------------------------------------
      ! ... Record and data counters
      inrc = kt - kit000 + 2
      isla = sladatqc%nsstp(inrc)

      ! Get the data for interpolation

      ALLOCATE( &
         & igrdi(2,2,isla), &
         & igrdj(2,2,isla), &
         & zglam(2,2,isla), &
         & zgphi(2,2,isla), &
         & zmask(2,2,isla), &
         & zsshl(2,2,isla)  &
         & )
      
      DO jobs = sladatqc%nsurfup + 1, sladatqc%nsurfup + isla
         iobs = jobs - sladatqc%nsurfup
         igrdi(1,1,iobs) = sladatqc%mi(jobs)-1
         igrdj(1,1,iobs) = sladatqc%mj(jobs)-1
         igrdi(1,2,iobs) = sladatqc%mi(jobs)-1
         igrdj(1,2,iobs) = sladatqc%mj(jobs)
         igrdi(2,1,iobs) = sladatqc%mi(jobs)
         igrdj(2,1,iobs) = sladatqc%mj(jobs)-1
         igrdi(2,2,iobs) = sladatqc%mi(jobs)
         igrdj(2,2,iobs) = sladatqc%mj(jobs)
      END DO

      CALL obs_int_comm_2d( 2, 2, isla, &
         &                  igrdi, igrdj, glamt, zglam )
      CALL obs_int_comm_2d( 2, 2, isla, &
         &                  igrdi, igrdj, gphit, zgphi )
      CALL obs_int_comm_2d( 2, 2, isla, &
         &                  igrdi, igrdj, psshmask, zmask )
      CALL obs_int_comm_2d( 2, 2, isla, &
         &                  igrdi, igrdj, psshn, zsshl )

      ! Loop over observations

      DO jobs = sladatqc%nsurfup + 1, sladatqc%nsurfup + isla

         iobs = jobs - sladatqc%nsurfup

         IF ( kt /= sladatqc%mstp(jobs) ) THEN
            
            IF(lwp) THEN
               WRITE(numout,*)
               WRITE(numout,*) ' E R R O R : Observation',              &
                  &            ' time step is not consistent with the', &
                  &            ' model time step'
               WRITE(numout,*) ' ========='
               WRITE(numout,*)
               WRITE(numout,*) ' Record  = ', jobs,                &
                  &            ' kt      = ', kt,                  &
                  &            ' mstp    = ', sladatqc%mstp(jobs), &
                  &            ' ntyp    = ', sladatqc%ntyp(jobs)
            ENDIF
            CALL ctl_stop( 'obs_sla_opt', 'Inconsistent time' )
            
         ENDIF
         
         zlam = sladatqc%rlam(jobs)
         zphi = sladatqc%rphi(jobs)

         ! Get weights to interpolate the model SSH to the observation point
         CALL obs_int_h2d_init( 1, 1, k2dint, zlam, zphi,         &
            &                   zglam(:,:,iobs), zgphi(:,:,iobs), &
            &                   zmask(:,:,iobs), zweig, zobsmask )
         

         ! Interpolate the model SSH to the observation point
         CALL obs_int_h2d( 1, 1,      &
            &              zweig, zsshl(:,:,iobs),  zext )
         
         sladatqc%rext(jobs,1) = zext(1)
         ! ... Remove the MDT at the observation point
         sladatqc%rmod(jobs,1) = sladatqc%rext(jobs,1) - sladatqc%rext(jobs,2)

      END DO

      ! Deallocate the data for interpolation
      DEALLOCATE( &
         & igrdi, &
         & igrdj, &
         & zglam, &
         & zgphi, &
         & zmask, &
         & zsshl  &
         & )

      sladatqc%nsurfup = sladatqc%nsurfup + isla

   END SUBROUTINE obs_sla_opt

   SUBROUTINE obs_sst_opt( sstdatqc, kt, kpi, kpj, kit000, kdaystp, &
      &                    psstn, psstmask, k2dint, ld_nightav )
      !!-----------------------------------------------------------------------
      !!
      !!                     ***  ROUTINE obs_sst_opt  ***
      !!
      !! ** Purpose : Compute the model counterpart of surface temperature
      !!              data by interpolating from the model grid to the 
      !!              observation point.
      !!
      !! ** Method  : Linearly interpolate to each observation point using 
      !!              the model values at the corners of the surrounding grid box.
      !!
      !!    The now model SST is first computed at the obs (lon, lat) point.
      !!
      !!    Several horizontal interpolation schemes are available:
      !!        - distance-weighted (great circle) (k2dint = 0)
      !!        - distance-weighted (small angle)  (k2dint = 1)
      !!        - bilinear (geographical grid)     (k2dint = 2)
      !!        - bilinear (quadrilateral grid)    (k2dint = 3)
      !!        - polynomial (quadrilateral grid)  (k2dint = 4)
      !!
      !!
      !! ** Action  :
      !!
      !! History :
      !!        !  07-07  (S. Ricci ) : Original
      !!      
      !!-----------------------------------------------------------------------

      !! * Modules used
      USE obs_surf_def  ! Definition of storage space for surface observations
      USE sbcdcy

      IMPLICIT NONE

      !! * Arguments
      TYPE(obs_surf), INTENT(INOUT) :: &
         & sstdatqc     ! Subset of surface data not failing screening
      INTEGER, INTENT(IN) :: kt        ! Time step
      INTEGER, INTENT(IN) :: kpi       ! Model grid parameters
      INTEGER, INTENT(IN) :: kpj
      INTEGER, INTENT(IN) :: kit000    ! Number of the first time step 
                                       !   (kit000-1 = restart time)
      INTEGER, INTENT(IN) :: k2dint    ! Horizontal interpolation type (see header)
      INTEGER, INTENT(IN) :: kdaystp   ! Number of time steps per day  
      REAL(KIND=wp), INTENT(IN), DIMENSION(kpi,kpj) :: &
         & psstn,  &    ! Model SST field
         & psstmask     ! Land-sea mask

      !! * Local declarations
      INTEGER :: ji
      INTEGER :: jj
      INTEGER :: jobs
      INTEGER :: inrc
      INTEGER :: isst
      INTEGER :: iobs
      INTEGER :: idayend
      REAL(KIND=wp) :: zlam
      REAL(KIND=wp) :: zphi
      REAL(KIND=wp) :: zext(1), zobsmask(1)
      REAL(KIND=wp) :: zdaystp
      INTEGER, DIMENSION(:,:), SAVE, ALLOCATABLE :: &
         & icount_sstnight,      &
         & imask_night
      REAL(kind=wp), DIMENSION(:,:), SAVE, ALLOCATABLE :: &
         & zintmp, &
         & zouttmp, & 
         & zmeanday    ! to compute model sst in region of 24h daylight (pole)
      REAL(kind=wp), DIMENSION(2,2,1) :: &
         & zweig
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: &
         & zmask, &
         & zsstl, &
         & zsstm, &
         & zglam, &
         & zgphi
      INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: &
         & igrdi, &
         & igrdj
      LOGICAL, INTENT(IN) :: ld_nightav

      !-----------------------------------------------------------------------
      ! Local initialization 
      !-----------------------------------------------------------------------
      ! ... Record and data counters
      inrc = kt - kit000 + 2
      isst = sstdatqc%nsstp(inrc)

      IF ( ld_nightav ) THEN

      ! Initialize array for night mean

      IF ( kt .EQ. 0 ) THEN
         ALLOCATE ( icount_sstnight(kpi,kpj) )
         ALLOCATE ( imask_night(kpi,kpj) )
         ALLOCATE ( zintmp(kpi,kpj) )
         ALLOCATE ( zouttmp(kpi,kpj) )
         ALLOCATE ( zmeanday(kpi,kpj) )
         nday_qsr = -1   ! initialisation flag for nbc_dcy
      ENDIF

      ! Initialize daily mean for first timestep
      idayend = MOD( kt - kit000 + 1, kdaystp )

      ! Added kt == 0 test to catch restart case 
      IF ( idayend == 1 .OR. kt == 0) THEN
         IF (lwp) WRITE(numout,*) 'Reset sstdatqc%vdmean on time-step: ',kt
         DO jj = 1, jpj
            DO ji = 1, jpi
               sstdatqc%vdmean(ji,jj) = 0.0
               zmeanday(ji,jj) = 0.0
               icount_sstnight(ji,jj) = 0
            END DO
         END DO
      ENDIF

      zintmp(:,:) = 0.0
      zouttmp(:,:) = sbc_dcy( zintmp(:,:), .TRUE. )
      imask_night(:,:) = INT( zouttmp(:,:) )

      DO jj = 1, jpj
         DO ji = 1, jpi
            ! Increment the temperature field for computing night mean and counter
            sstdatqc%vdmean(ji,jj) = sstdatqc%vdmean(ji,jj)  &
                   &                        + psstn(ji,jj)*imask_night(ji,jj)
            zmeanday(ji,jj)        = zmeanday(ji,jj) + psstn(ji,jj)
            icount_sstnight(ji,jj) = icount_sstnight(ji,jj) + imask_night(ji,jj)
         END DO
      END DO
   
      ! Compute the daily mean at the end of day

      zdaystp = 1.0 / REAL( kdaystp )

      IF ( idayend == 0 ) THEN 
         DO jj = 1, jpj
            DO ji = 1, jpi
               ! Test if "no night" point
               IF ( icount_sstnight(ji,jj) .NE. 0 ) THEN
                  sstdatqc%vdmean(ji,jj) = sstdatqc%vdmean(ji,jj) &
                    &                        / icount_sstnight(ji,jj) 
               ELSE
                  sstdatqc%vdmean(ji,jj) = zmeanday(ji,jj) * zdaystp
               ENDIF
            END DO
         END DO
      ENDIF

      ENDIF

      ! Get the data for interpolation
      
      ALLOCATE( &
         & igrdi(2,2,isst), &
         & igrdj(2,2,isst), &
         & zglam(2,2,isst), &
         & zgphi(2,2,isst), &
         & zmask(2,2,isst), &
         & zsstl(2,2,isst)  &
         & )
      
      DO jobs = sstdatqc%nsurfup + 1, sstdatqc%nsurfup + isst
         iobs = jobs - sstdatqc%nsurfup
         igrdi(1,1,iobs) = sstdatqc%mi(jobs)-1
         igrdj(1,1,iobs) = sstdatqc%mj(jobs)-1
         igrdi(1,2,iobs) = sstdatqc%mi(jobs)-1
         igrdj(1,2,iobs) = sstdatqc%mj(jobs)
         igrdi(2,1,iobs) = sstdatqc%mi(jobs)
         igrdj(2,1,iobs) = sstdatqc%mj(jobs)-1
         igrdi(2,2,iobs) = sstdatqc%mi(jobs)
         igrdj(2,2,iobs) = sstdatqc%mj(jobs)
      END DO
      
      CALL obs_int_comm_2d( 2, 2, isst, &
         &                  igrdi, igrdj, glamt, zglam )
      CALL obs_int_comm_2d( 2, 2, isst, &
         &                  igrdi, igrdj, gphit, zgphi )
      CALL obs_int_comm_2d( 2, 2, isst, &
         &                  igrdi, igrdj, psstmask, zmask )
      CALL obs_int_comm_2d( 2, 2, isst, &
         &                  igrdi, igrdj, psstn, zsstl )

      ! At the end of the day get interpolated means
      IF ( idayend == 0 .AND. ld_nightav ) THEN

         ALLOCATE( &
            & zsstm(2,2,isst)  &
            & )

         CALL obs_int_comm_2d( 2, 2, isst, igrdi, igrdj, &
            &               sstdatqc%vdmean(:,:), zsstm )

      ENDIF

      ! Loop over observations

      DO jobs = sstdatqc%nsurfup + 1, sstdatqc%nsurfup + isst
         
         iobs = jobs - sstdatqc%nsurfup
         
         IF ( kt /= sstdatqc%mstp(jobs) ) THEN
            
            IF(lwp) THEN
               WRITE(numout,*)
               WRITE(numout,*) ' E R R O R : Observation',              &
                  &            ' time step is not consistent with the', &
                  &            ' model time step'
               WRITE(numout,*) ' ========='
               WRITE(numout,*)
               WRITE(numout,*) ' Record  = ', jobs,                &
                  &            ' kt      = ', kt,                  &
                  &            ' mstp    = ', sstdatqc%mstp(jobs), &
                  &            ' ntyp    = ', sstdatqc%ntyp(jobs)
            ENDIF
            CALL ctl_stop( 'obs_sst_opt', 'Inconsistent time' )
            
         ENDIF
         
         zlam = sstdatqc%rlam(jobs)
         zphi = sstdatqc%rphi(jobs)
         
         ! Get weights to interpolate the model SST to the observation point
         CALL obs_int_h2d_init( 1, 1, k2dint, zlam, zphi,         &
            &                   zglam(:,:,iobs), zgphi(:,:,iobs), &
            &                   zmask(:,:,iobs), zweig, zobsmask )
            
         ! Interpolate the model SST to the observation point 

         IF ( ld_nightav ) THEN

           IF ( idayend == 0 )  THEN
               ! Daily averaged/diurnal cycle of SST  data
               CALL obs_int_h2d( 1, 1,      & 
                     &              zweig, zsstm(:,:,iobs), zext )
            ELSE 
               CALL ctl_stop( ' ld_nightav is set to true: a nonzero' //     &
                     &           ' number of night SST data should' // &
                     &           ' only occur at the end of a given day' )
            ENDIF

         ELSE

            CALL obs_int_h2d( 1, 1,      &
            &              zweig, zsstl(:,:,iobs),  zext )

         ENDIF
         sstdatqc%rmod(jobs,1) = zext(1)
         
      END DO
      
      ! Deallocate the data for interpolation
      DEALLOCATE( &
         & igrdi, &
         & igrdj, &
         & zglam, &
         & zgphi, &
         & zmask, &
         & zsstl  &
         & )

      ! At the end of the day also get interpolated means
      IF ( idayend == 0 .AND. ld_nightav ) THEN
         DEALLOCATE( &
            & zsstm  &
            & )
      ENDIF
      
      sstdatqc%nsurfup = sstdatqc%nsurfup + isst

   END SUBROUTINE obs_sst_opt

   SUBROUTINE obs_sss_opt
      !!-----------------------------------------------------------------------
      !!
      !!                     ***  ROUTINE obs_sss_opt  ***
      !!
      !! ** Purpose : Compute the model counterpart of sea surface salinity
      !!              data by interpolating from the model grid to the 
      !!              observation point.
      !!
      !! ** Method  : 
      !!
      !! ** Action  :
      !!
      !! History :
      !!      ! ??-?? 
      !!-----------------------------------------------------------------------

      IMPLICIT NONE

   END SUBROUTINE obs_sss_opt

   SUBROUTINE obs_seaice_opt( seaicedatqc, kt, kpi, kpj, kit000, &
      &                    pseaicen, pseaicemask, k2dint )

      !!-----------------------------------------------------------------------
      !!
      !!                     ***  ROUTINE obs_seaice_opt  ***
      !!
      !! ** Purpose : Compute the model counterpart of surface temperature
      !!              data by interpolating from the model grid to the 
      !!              observation point.
      !!
      !! ** Method  : Linearly interpolate to each observation point using 
      !!              the model values at the corners of the surrounding grid box.
      !!
      !!    The now model sea ice is first computed at the obs (lon, lat) point.
      !!
      !!    Several horizontal interpolation schemes are available:
      !!        - distance-weighted (great circle) (k2dint = 0)
      !!        - distance-weighted (small angle)  (k2dint = 1)
      !!        - bilinear (geographical grid)     (k2dint = 2)
      !!        - bilinear (quadrilateral grid)    (k2dint = 3)
      !!        - polynomial (quadrilateral grid)  (k2dint = 4)
      !!
      !!
      !! ** Action  :
      !!
      !! History :
      !!        !  07-07  (S. Ricci ) : Original
      !!      
      !!-----------------------------------------------------------------------

      !! * Modules used
      USE obs_surf_def  ! Definition of storage space for surface observations

      IMPLICIT NONE

      !! * Arguments
      TYPE(obs_surf), INTENT(INOUT) :: seaicedatqc     ! Subset of surface data not failing screening
      INTEGER, INTENT(IN) :: kt       ! Time step
      INTEGER, INTENT(IN) :: kpi      ! Model grid parameters
      INTEGER, INTENT(IN) :: kpj
      INTEGER, INTENT(IN) :: kit000   ! Number of the first time step 
                                      !   (kit000-1 = restart time)
      INTEGER, INTENT(IN) :: k2dint   ! Horizontal interpolation type (see header)
      REAL(KIND=wp), INTENT(IN), DIMENSION(kpi,kpj) :: &
         & pseaicen,  &    ! Model sea ice field
         & pseaicemask     ! Land-sea mask
         
      !! * Local declarations
      INTEGER :: ji
      INTEGER :: jj
      INTEGER :: jobs
      INTEGER :: inrc
      INTEGER :: iseaice
      INTEGER :: iobs
       
      REAL(KIND=wp) :: zlam
      REAL(KIND=wp) :: zphi
      REAL(KIND=wp) :: zext(1), zobsmask(1)
      REAL(kind=wp), DIMENSION(2,2,1) :: &
         & zweig
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: &
         & zmask, &
         & zseaicel, &
         & zglam, &
         & zgphi
      INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: &
         & igrdi, &
         & igrdj

      !------------------------------------------------------------------------
      ! Local initialization 
      !------------------------------------------------------------------------
      ! ... Record and data counters
      inrc = kt - kit000 + 2
      iseaice = seaicedatqc%nsstp(inrc)

      ! Get the data for interpolation
      
      ALLOCATE( &
         & igrdi(2,2,iseaice), &
         & igrdj(2,2,iseaice), &
         & zglam(2,2,iseaice), &
         & zgphi(2,2,iseaice), &
         & zmask(2,2,iseaice), &
         & zseaicel(2,2,iseaice)  &
         & )
      
      DO jobs = seaicedatqc%nsurfup + 1, seaicedatqc%nsurfup + iseaice
         iobs = jobs - seaicedatqc%nsurfup
         igrdi(1,1,iobs) = seaicedatqc%mi(jobs)-1
         igrdj(1,1,iobs) = seaicedatqc%mj(jobs)-1
         igrdi(1,2,iobs) = seaicedatqc%mi(jobs)-1
         igrdj(1,2,iobs) = seaicedatqc%mj(jobs)
         igrdi(2,1,iobs) = seaicedatqc%mi(jobs)
         igrdj(2,1,iobs) = seaicedatqc%mj(jobs)-1
         igrdi(2,2,iobs) = seaicedatqc%mi(jobs)
         igrdj(2,2,iobs) = seaicedatqc%mj(jobs)
      END DO
      
      CALL obs_int_comm_2d( 2, 2, iseaice, &
         &                  igrdi, igrdj, glamt, zglam )
      CALL obs_int_comm_2d( 2, 2, iseaice, &
         &                  igrdi, igrdj, gphit, zgphi )
      CALL obs_int_comm_2d( 2, 2, iseaice, &
         &                  igrdi, igrdj, pseaicemask, zmask )
      CALL obs_int_comm_2d( 2, 2, iseaice, &
         &                  igrdi, igrdj, pseaicen, zseaicel )
      
      DO jobs = seaicedatqc%nsurfup + 1, seaicedatqc%nsurfup + iseaice
         
         iobs = jobs - seaicedatqc%nsurfup
         
         IF ( kt /= seaicedatqc%mstp(jobs) ) THEN
            
            IF(lwp) THEN
               WRITE(numout,*)
               WRITE(numout,*) ' E R R O R : Observation',              &
                  &            ' time step is not consistent with the', &
                  &            ' model time step'
               WRITE(numout,*) ' ========='
               WRITE(numout,*)
               WRITE(numout,*) ' Record  = ', jobs,                &
                  &            ' kt      = ', kt,                  &
                  &            ' mstp    = ', seaicedatqc%mstp(jobs), &
                  &            ' ntyp    = ', seaicedatqc%ntyp(jobs)
            ENDIF
            CALL ctl_stop( 'obs_seaice_opt', 'Inconsistent time' )
            
         ENDIF
         
         zlam = seaicedatqc%rlam(jobs)
         zphi = seaicedatqc%rphi(jobs)
         
         ! Get weights to interpolate the model sea ice to the observation point
         CALL obs_int_h2d_init( 1, 1, k2dint, zlam, zphi,         &
            &                   zglam(:,:,iobs), zgphi(:,:,iobs), &
            &                   zmask(:,:,iobs), zweig, zobsmask )
         
         ! ... Interpolate the model sea ice to the observation point
         CALL obs_int_h2d( 1, 1,      &
            &              zweig, zseaicel(:,:,iobs),  zext )
         
         seaicedatqc%rmod(jobs,1) = zext(1)
         
      END DO
      
      ! Deallocate the data for interpolation
      DEALLOCATE( &
         & igrdi,    &
         & igrdj,    &
         & zglam,    &
         & zgphi,    &
         & zmask,    &
         & zseaicel  &
         & )
      
      seaicedatqc%nsurfup = seaicedatqc%nsurfup + iseaice

   END SUBROUTINE obs_seaice_opt

   SUBROUTINE obs_vel_opt( prodatqc, kt, kpi, kpj, kpk, kit000, kdaystp, &
      &                    pun, pvn, pgdept, pumask, pvmask, k1dint, k2dint, &
      &                    ld_dailyav )
      !!-----------------------------------------------------------------------
      !!
      !!                     ***  ROUTINE obs_vel_opt  ***
      !!
      !! ** Purpose : Compute the model counterpart of velocity profile
      !!              data by interpolating from the model grid to the 
      !!              observation point.
      !!
      !! ** Method  : Linearly interpolate zonal and meridional components of velocity 
      !!              to each observation point using the model values at the corners of 
      !!              the surrounding grid box. The model velocity components are on a 
      !!              staggered C- grid.
      !!
      !!    For velocity data from the TAO array, the model equivalent is
      !!    a daily mean velocity field. So, we first compute
      !!    the mean, then interpolate only at the end of the day.
      !!
      !! ** Action  :
      !!
      !! History :
      !!    ! 07-03 (K. Mogensen)      : Temperature and Salinity profiles
      !!    ! 08-10 (Maria Valdivieso) : Velocity component (U,V) profiles
      !!-----------------------------------------------------------------------
    
      !! * Modules used
      USE obs_profiles_def ! Definition of storage space for profile obs.

      IMPLICIT NONE

      !! * Arguments
      TYPE(obs_prof), INTENT(INOUT) :: &
         & prodatqc        ! Subset of profile data not failing screening
      INTEGER, INTENT(IN) :: kt        ! Time step
      INTEGER, INTENT(IN) :: kpi       ! Model grid parameters
      INTEGER, INTENT(IN) :: kpj
      INTEGER, INTENT(IN) :: kpk 
      INTEGER, INTENT(IN) :: kit000    ! Number of the first time step 
                                       !   (kit000-1 = restart time)
      INTEGER, INTENT(IN) :: k1dint    ! Vertical interpolation type (see header)
      INTEGER, INTENT(IN) :: k2dint    ! Horizontal interpolation type (see header)
      INTEGER, INTENT(IN) :: kdaystp   ! Number of time steps per day                    
      REAL(KIND=wp), INTENT(IN), DIMENSION(kpi,kpj,kpk) :: &
         & pun,    &    ! Model zonal component of velocity
         & pvn,    &    ! Model meridional component of velocity
         & pumask, &    ! Land-sea mask
         & pvmask       ! Land-sea mask
      REAL(KIND=wp), INTENT(IN), DIMENSION(kpk) :: &
         & pgdept       ! Model array of depth levels
      LOGICAL, INTENT(IN) :: ld_dailyav
         
      !! * Local declarations
      INTEGER :: ji
      INTEGER :: jj
      INTEGER :: jk
      INTEGER :: jobs
      INTEGER :: inrc
      INTEGER :: ipro
      INTEGER :: idayend
      INTEGER :: ista
      INTEGER :: iend
      INTEGER :: iobs
      INTEGER, DIMENSION(imaxavtypes) :: &
         & idailyavtypes
      REAL(KIND=wp) :: zlam
      REAL(KIND=wp) :: zphi
      REAL(KIND=wp) :: zdaystp
      REAL(KIND=wp), DIMENSION(kpk) :: &
         & zobsmasku, &
         & zobsmaskv, &
         & zobsmask,  &
         & zobsk,     &
         & zobs2k
      REAL(KIND=wp), DIMENSION(2,2,kpk) :: &
         & zweigu,zweigv
      REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE :: &
         & zumask, zvmask, &
         & zintu, &
         & zintv, &
         & zinmu, &
         & zinmv
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: &
         & zglamu, zglamv, &
         & zgphiu, zgphiv
      INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: &
         & igrdiu, &
         & igrdju, &
         & igrdiv, &
         & igrdjv

      !------------------------------------------------------------------------
      ! Local initialization 
      !------------------------------------------------------------------------
      ! ... Record and data counters
      inrc = kt - kit000 + 2
      ipro = prodatqc%npstp(inrc)

      ! Initialize daily mean for first timestep
      idayend = MOD( kt - kit000 + 1, kdaystp )

      ! Added kt == 0 test to catch restart case 
      IF ( idayend == 1 .OR. kt == 0) THEN
         IF (lwp) WRITE(numout,*) 'Reset prodatqc%vdmean on time-step: ',kt
         prodatqc%vdmean(:,:,:,1) = 0.0
         prodatqc%vdmean(:,:,:,2) = 0.0
      ENDIF

      ! Increment the zonal velocity field for computing daily mean
      prodatqc%vdmean(:,:,:,1) = prodatqc%vdmean(:,:,:,1) + pun(:,:,:)
      ! Increment the meridional velocity field for computing daily mean
      prodatqc%vdmean(:,:,:,2) = prodatqc%vdmean(:,:,:,2) + pvn(:,:,:)
   
      ! Compute the daily mean at the end of day
      zdaystp = 1.0 / REAL( kdaystp )
      IF ( idayend == 0 ) THEN
         prodatqc%vdmean(:,:,:,1) = prodatqc%vdmean(:,:,:,1) * zdaystp
         prodatqc%vdmean(:,:,:,2) = prodatqc%vdmean(:,:,:,2) * zdaystp
      ENDIF

      ! Get the data for interpolation
      ALLOCATE( &
         & igrdiu(2,2,ipro),      &
         & igrdju(2,2,ipro),      &
         & igrdiv(2,2,ipro),      &
         & igrdjv(2,2,ipro),      &
         & zglamu(2,2,ipro), zglamv(2,2,ipro), &
         & zgphiu(2,2,ipro), zgphiv(2,2,ipro), &
         & zumask(2,2,kpk,ipro), zvmask(2,2,kpk,ipro), &
         & zintu(2,2,kpk,ipro),  &
         & zintv(2,2,kpk,ipro)   &
         & )

      DO jobs = prodatqc%nprofup + 1, prodatqc%nprofup + ipro
         iobs = jobs - prodatqc%nprofup
         igrdiu(1,1,iobs) = prodatqc%mi(jobs,1)-1
         igrdju(1,1,iobs) = prodatqc%mj(jobs,1)-1
         igrdiu(1,2,iobs) = prodatqc%mi(jobs,1)-1
         igrdju(1,2,iobs) = prodatqc%mj(jobs,1)
         igrdiu(2,1,iobs) = prodatqc%mi(jobs,1)
         igrdju(2,1,iobs) = prodatqc%mj(jobs,1)-1
         igrdiu(2,2,iobs) = prodatqc%mi(jobs,1)
         igrdju(2,2,iobs) = prodatqc%mj(jobs,1)
         igrdiv(1,1,iobs) = prodatqc%mi(jobs,2)-1
         igrdjv(1,1,iobs) = prodatqc%mj(jobs,2)-1
         igrdiv(1,2,iobs) = prodatqc%mi(jobs,2)-1
         igrdjv(1,2,iobs) = prodatqc%mj(jobs,2)
         igrdiv(2,1,iobs) = prodatqc%mi(jobs,2)
         igrdjv(2,1,iobs) = prodatqc%mj(jobs,2)-1
         igrdiv(2,2,iobs) = prodatqc%mi(jobs,2)
         igrdjv(2,2,iobs) = prodatqc%mj(jobs,2)
      END DO

      CALL obs_int_comm_2d( 2, 2, ipro, igrdiu, igrdju, glamu, zglamu )
      CALL obs_int_comm_2d( 2, 2, ipro, igrdiu, igrdju, gphiu, zgphiu )
      CALL obs_int_comm_3d( 2, 2, ipro, kpk, igrdiu, igrdju, pumask, zumask )
      CALL obs_int_comm_3d( 2, 2, ipro, kpk, igrdiu, igrdju, pun, zintu )

      CALL obs_int_comm_2d( 2, 2, ipro, igrdiv, igrdjv, glamv, zglamv )
      CALL obs_int_comm_2d( 2, 2, ipro, igrdiv, igrdjv, gphiv, zgphiv )
      CALL obs_int_comm_3d( 2, 2, ipro, kpk, igrdiv, igrdjv, pvmask, zvmask )
      CALL obs_int_comm_3d( 2, 2, ipro, kpk, igrdiv, igrdjv, pvn, zintv )

      ! At the end of the day also get interpolated means
      IF ( idayend == 0 ) THEN

         ALLOCATE( &
            & zinmu(2,2,kpk,ipro),  &
            & zinmv(2,2,kpk,ipro)   &
            & )

         CALL obs_int_comm_3d( 2, 2, ipro, kpk, igrdiu, igrdju, &
            &                  prodatqc%vdmean(:,:,:,1), zinmu )
         CALL obs_int_comm_3d( 2, 2, ipro, kpk, igrdiv, igrdjv, &
            &                  prodatqc%vdmean(:,:,:,2), zinmv )

      ENDIF

! loop over observations

      DO jobs = prodatqc%nprofup + 1, prodatqc%nprofup + ipro

         iobs = jobs - prodatqc%nprofup

         IF ( kt /= prodatqc%mstp(jobs) ) THEN
            
            IF(lwp) THEN
               WRITE(numout,*)
               WRITE(numout,*) ' E R R O R : Observation',              &
                  &            ' time step is not consistent with the', &
                  &            ' model time step'
               WRITE(numout,*) ' ========='
               WRITE(numout,*)
               WRITE(numout,*) ' Record  = ', jobs,                    &
                  &            ' kt      = ', kt,                      &
                  &            ' mstp    = ', prodatqc%mstp(jobs), &
                  &            ' ntyp    = ', prodatqc%ntyp(jobs)
            ENDIF
            CALL ctl_stop( 'obs_pro_opt', 'Inconsistent time' )
         ENDIF
         
         zlam = prodatqc%rlam(jobs)
         zphi = prodatqc%rphi(jobs)

         ! Initialize observation masks

         zobsmasku(:) = 0.0
         zobsmaskv(:) = 0.0
         
         ! Horizontal weights and vertical mask

         IF  ( prodatqc%npvend(jobs,1) > 0 ) THEN

            CALL obs_int_h2d_init( kpk, kpk, k2dint, zlam, zphi,     &
               &                   zglamu(:,:,iobs), zgphiu(:,:,iobs), &
               &                   zumask(:,:,:,iobs), zweigu, zobsmasku )

         ENDIF

         
         IF ( prodatqc%npvend(jobs,2) > 0 ) THEN

            CALL obs_int_h2d_init( kpk, kpk, k2dint, zlam, zphi,     &
               &                   zglamv(:,:,iobs), zgphiv(:,:,iobs), &
               &                   zvmask(:,:,:,iobs), zweigv, zobsmaskv )

         ENDIF

         ! Ensure that the vertical mask on u and v are consistent.

         zobsmask(:) = MIN( zobsmasku(:), zobsmaskv(:) )

         IF ( prodatqc%npvend(jobs,1) > 0 ) THEN

            zobsk(:) = obfillflt

	    IF ( ld_dailyav ) THEN

               IF ( idayend == 0 )  THEN
                  
                  ! Daily averaged data
                  
                  CALL obs_int_h2d( kpk, kpk,      &
                     &              zweigu, zinmu(:,:,:,iobs), zobsk )
                  
                  
               ELSE
               
                  CALL ctl_stop( ' A nonzero' //     &
                     &           ' number of U profile data should' // &
                     &           ' only occur at the end of a given day' )

               ENDIF
	       
            ELSE 
               
               ! Point data

               CALL obs_int_h2d( kpk, kpk,      &
                  &              zweigu, zintu(:,:,:,iobs), zobsk )

            ENDIF

            !-------------------------------------------------------------
            ! Compute vertical second-derivative of the interpolating 
            ! polynomial at obs points
            !-------------------------------------------------------------
            
            IF ( k1dint == 1 ) THEN
               CALL obs_int_z1d_spl( kpk, zobsk, zobs2k,   &
                  &                  pgdept, zobsmask )
            ENDIF
            
            !-----------------------------------------------------------------
            !  Vertical interpolation to the observation point
            !-----------------------------------------------------------------
            ista = prodatqc%npvsta(jobs,1)
            iend = prodatqc%npvend(jobs,1)
            CALL obs_int_z1d( kpk,                &
               & prodatqc%var(1)%mvk(ista:iend),  &
               & k1dint, iend - ista + 1,         &
               & prodatqc%var(1)%vdep(ista:iend), &
               & zobsk, zobs2k,                   &
               & prodatqc%var(1)%vmod(ista:iend), &
               & pgdept, zobsmask )

         ENDIF

         IF ( prodatqc%npvend(jobs,2) > 0 ) THEN

            zobsk(:) = obfillflt

            IF ( ld_dailyav ) THEN

               IF ( idayend == 0 )  THEN

                  ! Daily averaged data
                  
                  CALL obs_int_h2d( kpk, kpk,      &
                     &              zweigv, zinmv(:,:,:,iobs), zobsk )
                  
               ELSE

                  CALL ctl_stop( ' A nonzero' //     &
                     &           ' number of V profile data should' // &
                     &           ' only occur at the end of a given day' )

               ENDIF

            ELSE
               
               ! Point data

               CALL obs_int_h2d( kpk, kpk,      &
                  &              zweigv, zintv(:,:,:,iobs), zobsk )

            ENDIF


            !-------------------------------------------------------------
            ! Compute vertical second-derivative of the interpolating 
            ! polynomial at obs points
            !-------------------------------------------------------------
            
            IF ( k1dint == 1 ) THEN
               CALL obs_int_z1d_spl( kpk, zobsk, zobs2k, &
                  &                  pgdept, zobsmask )
            ENDIF
            
            !----------------------------------------------------------------
            !  Vertical interpolation to the observation point
            !----------------------------------------------------------------
            ista = prodatqc%npvsta(jobs,2)
            iend = prodatqc%npvend(jobs,2)
            CALL obs_int_z1d( kpk, &
               & prodatqc%var(2)%mvk(ista:iend),&
               & k1dint, iend - ista + 1, &
               & prodatqc%var(2)%vdep(ista:iend),&
               & zobsk, zobs2k, &
               & prodatqc%var(2)%vmod(ista:iend),&
               & pgdept, zobsmask )

         ENDIF

      END DO
 
      ! Deallocate the data for interpolation
      DEALLOCATE( &
         & igrdiu, &
         & igrdju, &
         & igrdiv, &
         & igrdjv, &
         & zglamu, zglamv, &
         & zgphiu, zgphiv, &
         & zumask, zvmask, &
         & zintu, &
         & zintv  &
         & )
      ! At the end of the day also get interpolated means
      IF ( idayend == 0 ) THEN
         DEALLOCATE( &
            & zinmu,  &
            & zinmv   &
            & )
      ENDIF

      prodatqc%nprofup = prodatqc%nprofup + ipro 
      
   END SUBROUTINE obs_vel_opt

   SUBROUTINE obs_logchl_opt( logchldatqc, kt, kpi, kpj, kit000, &
      &                    plogchln, plogchlmask, k2dint )

      !!-----------------------------------------------------------------------
      !!
      !!                     ***  ROUTINE obs_logchl_opt  ***
      !!
      !! ** Purpose : Compute the model counterpart of logchl
      !!              data by interpolating from the model grid to the 
      !!              observation point.
      !!
      !! ** Method  : Linearly interpolate to each observation point using 
      !!              the model values at the corners of the surrounding grid box.
      !!
      !!    The now model logchl is first computed at the obs (lon, lat) point.
      !!
      !!    Several horizontal interpolation schemes are available:
      !!        - distance-weighted (great circle) (k2dint = 0)
      !!        - distance-weighted (small angle)  (k2dint = 1)
      !!        - bilinear (geographical grid)     (k2dint = 2)
      !!        - bilinear (quadrilateral grid)    (k2dint = 3)
      !!        - polynomial (quadrilateral grid)  (k2dint = 4)
      !!
      !!
      !! ** Action  :
      !!
      !! History :
      !!      
      !!-----------------------------------------------------------------------

      !! * Modules used
      USE obs_surf_def  ! Definition of storage space for surface observations

      IMPLICIT NONE

      !! * Arguments
      TYPE(obs_surf), INTENT(INOUT) :: logchldatqc     ! Subset of surface data not failing screening
      INTEGER, INTENT(IN) :: kt       ! Time step
      INTEGER, INTENT(IN) :: kpi      ! Model grid parameters
      INTEGER, INTENT(IN) :: kpj
      INTEGER, INTENT(IN) :: kit000   ! Number of the first time step 
                                      !   (kit000-1 = restart time)
      INTEGER, INTENT(IN) :: k2dint   ! Horizontal interpolation type (see header)
      REAL(KIND=wp), INTENT(IN), DIMENSION(kpi,kpj) :: &
         & plogchln,  &    ! Model logchl field
         & plogchlmask     ! Land-sea mask
         
      !! * Local declarations
      INTEGER :: ji
      INTEGER :: jj
      INTEGER :: jobs
      INTEGER :: inrc
      INTEGER :: ilogchl
      INTEGER :: iobs
       
      REAL(KIND=wp) :: zlam
      REAL(KIND=wp) :: zphi
      REAL(KIND=wp) :: zext(1), zobsmask(1)
      REAL(kind=wp), DIMENSION(2,2,1) :: &
         & zweig
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: &
         & zmask, &
         & zlogchll, &
         & zglam, &
         & zgphi
      INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: &
         & igrdi, &
         & igrdj

      !------------------------------------------------------------------------
      ! Local initialization 
      !------------------------------------------------------------------------
      ! ... Record and data counters
      inrc = kt - kit000 + 2
      ilogchl = logchldatqc%nsstp(inrc)

      ! Get the data for interpolation
      
      ALLOCATE( &
         & igrdi(2,2,ilogchl), &
         & igrdj(2,2,ilogchl), &
         & zglam(2,2,ilogchl), &
         & zgphi(2,2,ilogchl), &
         & zmask(2,2,ilogchl), &
         & zlogchll(2,2,ilogchl)  &
         & )
      
      DO jobs = logchldatqc%nsurfup + 1, logchldatqc%nsurfup + ilogchl
         iobs = jobs - logchldatqc%nsurfup
         igrdi(1,1,iobs) = logchldatqc%mi(jobs)-1
         igrdj(1,1,iobs) = logchldatqc%mj(jobs)-1
         igrdi(1,2,iobs) = logchldatqc%mi(jobs)-1
         igrdj(1,2,iobs) = logchldatqc%mj(jobs)
         igrdi(2,1,iobs) = logchldatqc%mi(jobs)
         igrdj(2,1,iobs) = logchldatqc%mj(jobs)-1
         igrdi(2,2,iobs) = logchldatqc%mi(jobs)
         igrdj(2,2,iobs) = logchldatqc%mj(jobs)
      END DO
      
      CALL obs_int_comm_2d( 2, 2, ilogchl, &
         &                  igrdi, igrdj, glamt, zglam )
      CALL obs_int_comm_2d( 2, 2, ilogchl, &
         &                  igrdi, igrdj, gphit, zgphi )
      CALL obs_int_comm_2d( 2, 2, ilogchl, &
         &                  igrdi, igrdj, plogchlmask, zmask )
      CALL obs_int_comm_2d( 2, 2, ilogchl, &
         &                  igrdi, igrdj, plogchln, zlogchll )
      
      DO jobs = logchldatqc%nsurfup + 1, logchldatqc%nsurfup + ilogchl
         
         iobs = jobs - logchldatqc%nsurfup
         
         IF ( kt /= logchldatqc%mstp(jobs) ) THEN
            
            IF(lwp) THEN
               WRITE(numout,*)
               WRITE(numout,*) ' E R R O R : Observation',              &
                  &            ' time step is not consistent with the', &
                  &            ' model time step'
               WRITE(numout,*) ' ========='
               WRITE(numout,*)
               WRITE(numout,*) ' Record  = ', jobs,                &
                  &            ' kt      = ', kt,                  &
                  &            ' mstp    = ', logchldatqc%mstp(jobs), &
                  &            ' ntyp    = ', logchldatqc%ntyp(jobs)
            ENDIF
            CALL ctl_stop( 'obs_logchl_opt', 'Inconsistent time' )
            
         ENDIF
         
         zlam = logchldatqc%rlam(jobs)
         zphi = logchldatqc%rphi(jobs)
         
         ! Get weights to interpolate the model logchl to the observation point
         CALL obs_int_h2d_init( 1, 1, k2dint, zlam, zphi,         &
            &                   zglam(:,:,iobs), zgphi(:,:,iobs), &
            &                   zmask(:,:,iobs), zweig, zobsmask )
         
         ! ... Interpolate the model logchl to the observation point
         CALL obs_int_h2d( 1, 1,      &
            &              zweig, zlogchll(:,:,iobs),  zext )
         
         logchldatqc%rmod(jobs,1) = zext(1)
         
      END DO
      
      ! Deallocate the data for interpolation
      DEALLOCATE( &
         & igrdi,    &
         & igrdj,    &
         & zglam,    &
         & zgphi,    &
         & zmask,    &
         & zlogchll  &
         & )
      
      logchldatqc%nsurfup = logchldatqc%nsurfup + ilogchl

   END SUBROUTINE obs_logchl_opt

   SUBROUTINE obs_spm_opt( spmdatqc, kt, kpi, kpj, kit000, &
      &                    pspmn, pspmmask, k2dint )

      !!-----------------------------------------------------------------------
      !!
      !!                     ***  ROUTINE obs_spm_opt  ***
      !!
      !! ** Purpose : Compute the model counterpart of spm
      !!              data by interpolating from the model grid to the 
      !!              observation point.
      !!
      !! ** Method  : Linearly interpolate to each observation point using 
      !!              the model values at the corners of the surrounding grid box.
      !!
      !!    The now model spm is first computed at the obs (lon, lat) point.
      !!
      !!    Several horizontal interpolation schemes are available:
      !!        - distance-weighted (great circle) (k2dint = 0)
      !!        - distance-weighted (small angle)  (k2dint = 1)
      !!        - bilinear (geographical grid)     (k2dint = 2)
      !!        - bilinear (quadrilateral grid)    (k2dint = 3)
      !!        - polynomial (quadrilateral grid)  (k2dint = 4)
      !!
      !!
      !! ** Action  :
      !!
      !! History :
      !!      
      !!-----------------------------------------------------------------------

      !! * Modules used
      USE obs_surf_def  ! Definition of storage space for surface observations

      IMPLICIT NONE

      !! * Arguments
      TYPE(obs_surf), INTENT(INOUT) :: spmdatqc     ! Subset of surface data not failing screening
      INTEGER, INTENT(IN) :: kt       ! Time step
      INTEGER, INTENT(IN) :: kpi      ! Model grid parameters
      INTEGER, INTENT(IN) :: kpj
      INTEGER, INTENT(IN) :: kit000   ! Number of the first time step 
                                      !   (kit000-1 = restart time)
      INTEGER, INTENT(IN) :: k2dint   ! Horizontal interpolation type (see header)
      REAL(KIND=wp), INTENT(IN), DIMENSION(kpi,kpj) :: &
         & pspmn,  &    ! Model spm field
         & pspmmask     ! Land-sea mask
         
      !! * Local declarations
      INTEGER :: ji
      INTEGER :: jj
      INTEGER :: jobs
      INTEGER :: inrc
      INTEGER :: ispm
      INTEGER :: iobs
       
      REAL(KIND=wp) :: zlam
      REAL(KIND=wp) :: zphi
      REAL(KIND=wp) :: zext(1), zobsmask(1)
      REAL(kind=wp), DIMENSION(2,2,1) :: &
         & zweig
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: &
         & zmask, &
         & zspml, &
         & zglam, &
         & zgphi
      INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: &
         & igrdi, &
         & igrdj

      !------------------------------------------------------------------------
      ! Local initialization 
      !------------------------------------------------------------------------
      ! ... Record and data counters
      inrc = kt - kit000 + 2
      ispm = spmdatqc%nsstp(inrc)

      ! Get the data for interpolation
      
      ALLOCATE( &
         & igrdi(2,2,ispm), &
         & igrdj(2,2,ispm), &
         & zglam(2,2,ispm), &
         & zgphi(2,2,ispm), &
         & zmask(2,2,ispm), &
         & zspml(2,2,ispm)  &
         & )
      
      DO jobs = spmdatqc%nsurfup + 1, spmdatqc%nsurfup + ispm
         iobs = jobs - spmdatqc%nsurfup
         igrdi(1,1,iobs) = spmdatqc%mi(jobs)-1
         igrdj(1,1,iobs) = spmdatqc%mj(jobs)-1
         igrdi(1,2,iobs) = spmdatqc%mi(jobs)-1
         igrdj(1,2,iobs) = spmdatqc%mj(jobs)
         igrdi(2,1,iobs) = spmdatqc%mi(jobs)
         igrdj(2,1,iobs) = spmdatqc%mj(jobs)-1
         igrdi(2,2,iobs) = spmdatqc%mi(jobs)
         igrdj(2,2,iobs) = spmdatqc%mj(jobs)
      END DO
      
      CALL obs_int_comm_2d( 2, 2, ispm, &
         &                  igrdi, igrdj, glamt, zglam )
      CALL obs_int_comm_2d( 2, 2, ispm, &
         &                  igrdi, igrdj, gphit, zgphi )
      CALL obs_int_comm_2d( 2, 2, ispm, &
         &                  igrdi, igrdj, pspmmask, zmask )
      CALL obs_int_comm_2d( 2, 2, ispm, &
         &                  igrdi, igrdj, pspmn, zspml )
      
      DO jobs = spmdatqc%nsurfup + 1, spmdatqc%nsurfup + ispm
         
         iobs = jobs - spmdatqc%nsurfup
         
         IF ( kt /= spmdatqc%mstp(jobs) ) THEN
            
            IF(lwp) THEN
               WRITE(numout,*)
               WRITE(numout,*) ' E R R O R : Observation',              &
                  &            ' time step is not consistent with the', &
                  &            ' model time step'
               WRITE(numout,*) ' ========='
               WRITE(numout,*)
               WRITE(numout,*) ' Record  = ', jobs,                &
                  &            ' kt      = ', kt,                  &
                  &            ' mstp    = ', spmdatqc%mstp(jobs), &
                  &            ' ntyp    = ', spmdatqc%ntyp(jobs)
            ENDIF
            CALL ctl_stop( 'obs_spm_opt', 'Inconsistent time' )
            
         ENDIF
         
         zlam = spmdatqc%rlam(jobs)
         zphi = spmdatqc%rphi(jobs)
         
         ! Get weights to interpolate the model spm to the observation point
         CALL obs_int_h2d_init( 1, 1, k2dint, zlam, zphi,         &
            &                   zglam(:,:,iobs), zgphi(:,:,iobs), &
            &                   zmask(:,:,iobs), zweig, zobsmask )
         
         ! ... Interpolate the model spm to the observation point
         CALL obs_int_h2d( 1, 1,      &
            &              zweig, zspml(:,:,iobs),  zext )
         
         spmdatqc%rmod(jobs,1) = zext(1)
         
      END DO
      
      ! Deallocate the data for interpolation
      DEALLOCATE( &
         & igrdi,    &
         & igrdj,    &
         & zglam,    &
         & zgphi,    &
         & zmask,    &
         & zspml  &
         & )
      
      spmdatqc%nsurfup = spmdatqc%nsurfup + ispm

   END SUBROUTINE obs_spm_opt

   SUBROUTINE obs_fco2_opt( fco2datqc, kt, kpi, kpj, kit000, &
      &                    pfco2n, pfco2mask, k2dint )

      !!-----------------------------------------------------------------------
      !!
      !!                     ***  ROUTINE obs_fco2_opt  ***
      !!
      !! ** Purpose : Compute the model counterpart of fco2
      !!              data by interpolating from the model grid to the 
      !!              observation point.
      !!
      !! ** Method  : Linearly interpolate to each observation point using 
      !!              the model values at the corners of the surrounding grid box.
      !!
      !!    The now model fco2 is first computed at the obs (lon, lat) point.
      !!
      !!    Several horizontal interpolation schemes are available:
      !!        - distance-weighted (great circle) (k2dint = 0)
      !!        - distance-weighted (small angle)  (k2dint = 1)
      !!        - bilinear (geographical grid)     (k2dint = 2)
      !!        - bilinear (quadrilateral grid)    (k2dint = 3)
      !!        - polynomial (quadrilateral grid)  (k2dint = 4)
      !!
      !!
      !! ** Action  :
      !!
      !! History :
      !!      
      !!-----------------------------------------------------------------------

      !! * Modules used
      USE obs_surf_def  ! Definition of storage space for surface observations

      IMPLICIT NONE

      !! * Arguments
      TYPE(obs_surf), INTENT(INOUT) :: fco2datqc     ! Subset of surface data not failing screening
      INTEGER, INTENT(IN) :: kt       ! Time step
      INTEGER, INTENT(IN) :: kpi      ! Model grid parameters
      INTEGER, INTENT(IN) :: kpj
      INTEGER, INTENT(IN) :: kit000   ! Number of the first time step 
                                      !   (kit000-1 = restart time)
      INTEGER, INTENT(IN) :: k2dint   ! Horizontal interpolation type (see header)
      REAL(KIND=wp), INTENT(IN), DIMENSION(kpi,kpj) :: &
         & pfco2n,  &    ! Model fco2 field
         & pfco2mask     ! Land-sea mask
         
      !! * Local declarations
      INTEGER :: ji
      INTEGER :: jj
      INTEGER :: jobs
      INTEGER :: inrc
      INTEGER :: ifco2
      INTEGER :: iobs
       
      REAL(KIND=wp) :: zlam
      REAL(KIND=wp) :: zphi
      REAL(KIND=wp) :: zext(1), zobsmask(1)
      REAL(kind=wp), DIMENSION(2,2,1) :: &
         & zweig
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: &
         & zmask, &
         & zfco2l, &
         & zglam, &
         & zgphi
      INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: &
         & igrdi, &
         & igrdj

      !------------------------------------------------------------------------
      ! Local initialization 
      !------------------------------------------------------------------------
      ! ... Record and data counters
      inrc = kt - kit000 + 2
      ifco2 = fco2datqc%nsstp(inrc)

      ! Get the data for interpolation
      
      ALLOCATE( &
         & igrdi(2,2,ifco2), &
         & igrdj(2,2,ifco2), &
         & zglam(2,2,ifco2), &
         & zgphi(2,2,ifco2), &
         & zmask(2,2,ifco2), &
         & zfco2l(2,2,ifco2)  &
         & )
      
      DO jobs = fco2datqc%nsurfup + 1, fco2datqc%nsurfup + ifco2
         iobs = jobs - fco2datqc%nsurfup
         igrdi(1,1,iobs) = fco2datqc%mi(jobs)-1
         igrdj(1,1,iobs) = fco2datqc%mj(jobs)-1
         igrdi(1,2,iobs) = fco2datqc%mi(jobs)-1
         igrdj(1,2,iobs) = fco2datqc%mj(jobs)
         igrdi(2,1,iobs) = fco2datqc%mi(jobs)
         igrdj(2,1,iobs) = fco2datqc%mj(jobs)-1
         igrdi(2,2,iobs) = fco2datqc%mi(jobs)
         igrdj(2,2,iobs) = fco2datqc%mj(jobs)
      END DO
      
      CALL obs_int_comm_2d( 2, 2, ifco2, &
         &                  igrdi, igrdj, glamt, zglam )
      CALL obs_int_comm_2d( 2, 2, ifco2, &
         &                  igrdi, igrdj, gphit, zgphi )
      CALL obs_int_comm_2d( 2, 2, ifco2, &
         &                  igrdi, igrdj, pfco2mask, zmask )
      CALL obs_int_comm_2d( 2, 2, ifco2, &
         &                  igrdi, igrdj, pfco2n, zfco2l )
      
      DO jobs = fco2datqc%nsurfup + 1, fco2datqc%nsurfup + ifco2
         
         iobs = jobs - fco2datqc%nsurfup
         
         IF ( kt /= fco2datqc%mstp(jobs) ) THEN
            
            IF(lwp) THEN
               WRITE(numout,*)
               WRITE(numout,*) ' E R R O R : Observation',              &
                  &            ' time step is not consistent with the', &
                  &            ' model time step'
               WRITE(numout,*) ' ========='
               WRITE(numout,*)
               WRITE(numout,*) ' Record  = ', jobs,                &
                  &            ' kt      = ', kt,                  &
                  &            ' mstp    = ', fco2datqc%mstp(jobs), &
                  &            ' ntyp    = ', fco2datqc%ntyp(jobs)
            ENDIF
            CALL ctl_stop( 'obs_fco2_opt', 'Inconsistent time' )
            
         ENDIF
         
         zlam = fco2datqc%rlam(jobs)
         zphi = fco2datqc%rphi(jobs)
         
         ! Get weights to interpolate the model fco2 to the observation point
         CALL obs_int_h2d_init( 1, 1, k2dint, zlam, zphi,         &
            &                   zglam(:,:,iobs), zgphi(:,:,iobs), &
            &                   zmask(:,:,iobs), zweig, zobsmask )
         
         ! ... Interpolate the model fco2 to the observation point
         CALL obs_int_h2d( 1, 1,      &
            &              zweig, zfco2l(:,:,iobs),  zext )
         
         fco2datqc%rmod(jobs,1) = zext(1)
         
      END DO
      
      ! Deallocate the data for interpolation
      DEALLOCATE( &
         & igrdi,    &
         & igrdj,    &
         & zglam,    &
         & zgphi,    &
         & zmask,    &
         & zfco2l  &
         & )
      
      fco2datqc%nsurfup = fco2datqc%nsurfup + ifco2

   END SUBROUTINE obs_fco2_opt

   SUBROUTINE obs_pco2_opt( pco2datqc, kt, kpi, kpj, kit000, &
      &                    ppco2n, ppco2mask, k2dint )

      !!-----------------------------------------------------------------------
      !!
      !!                     ***  ROUTINE obs_pco2_opt  ***
      !!
      !! ** Purpose : Compute the model counterpart of pco2
      !!              data by interpolating from the model grid to the 
      !!              observation point.
      !!
      !! ** Method  : Linearly interpolate to each observation point using 
      !!              the model values at the corners of the surrounding grid box.
      !!
      !!    The now model pco2 is first computed at the obs (lon, lat) point.
      !!
      !!    Several horizontal interpolation schemes are available:
      !!        - distance-weighted (great circle) (k2dint = 0)
      !!        - distance-weighted (small angle)  (k2dint = 1)
      !!        - bilinear (geographical grid)     (k2dint = 2)
      !!        - bilinear (quadrilateral grid)    (k2dint = 3)
      !!        - polynomial (quadrilateral grid)  (k2dint = 4)
      !!
      !!
      !! ** Action  :
      !!
      !! History :
      !!      
      !!-----------------------------------------------------------------------

      !! * Modules used
      USE obs_surf_def  ! Definition of storage space for surface observations

      IMPLICIT NONE

      !! * Arguments
      TYPE(obs_surf), INTENT(INOUT) :: pco2datqc     ! Subset of surface data not failing screening
      INTEGER, INTENT(IN) :: kt       ! Time step
      INTEGER, INTENT(IN) :: kpi      ! Model grid parameters
      INTEGER, INTENT(IN) :: kpj
      INTEGER, INTENT(IN) :: kit000   ! Number of the first time step 
                                      !   (kit000-1 = restart time)
      INTEGER, INTENT(IN) :: k2dint   ! Horizontal interpolation type (see header)
      REAL(KIND=wp), INTENT(IN), DIMENSION(kpi,kpj) :: &
         & ppco2n,  &    ! Model pco2 field
         & ppco2mask     ! Land-sea mask
         
      !! * Local declarations
      INTEGER :: ji
      INTEGER :: jj
      INTEGER :: jobs
      INTEGER :: inrc
      INTEGER :: ipco2
      INTEGER :: iobs
       
      REAL(KIND=wp) :: zlam
      REAL(KIND=wp) :: zphi
      REAL(KIND=wp) :: zext(1), zobsmask(1)
      REAL(kind=wp), DIMENSION(2,2,1) :: &
         & zweig
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: &
         & zmask, &
         & zpco2l, &
         & zglam, &
         & zgphi
      INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: &
         & igrdi, &
         & igrdj

      !------------------------------------------------------------------------
      ! Local initialization 
      !------------------------------------------------------------------------
      ! ... Record and data counters
      inrc = kt - kit000 + 2
      ipco2 = pco2datqc%nsstp(inrc)

      ! Get the data for interpolation
      
      ALLOCATE( &
         & igrdi(2,2,ipco2), &
         & igrdj(2,2,ipco2), &
         & zglam(2,2,ipco2), &
         & zgphi(2,2,ipco2), &
         & zmask(2,2,ipco2), &
         & zpco2l(2,2,ipco2)  &
         & )
      
      DO jobs = pco2datqc%nsurfup + 1, pco2datqc%nsurfup + ipco2
         iobs = jobs - pco2datqc%nsurfup
         igrdi(1,1,iobs) = pco2datqc%mi(jobs)-1
         igrdj(1,1,iobs) = pco2datqc%mj(jobs)-1
         igrdi(1,2,iobs) = pco2datqc%mi(jobs)-1
         igrdj(1,2,iobs) = pco2datqc%mj(jobs)
         igrdi(2,1,iobs) = pco2datqc%mi(jobs)
         igrdj(2,1,iobs) = pco2datqc%mj(jobs)-1
         igrdi(2,2,iobs) = pco2datqc%mi(jobs)
         igrdj(2,2,iobs) = pco2datqc%mj(jobs)
      END DO
      
      CALL obs_int_comm_2d( 2, 2, ipco2, &
         &                  igrdi, igrdj, glamt, zglam )
      CALL obs_int_comm_2d( 2, 2, ipco2, &
         &                  igrdi, igrdj, gphit, zgphi )
      CALL obs_int_comm_2d( 2, 2, ipco2, &
         &                  igrdi, igrdj, ppco2mask, zmask )
      CALL obs_int_comm_2d( 2, 2, ipco2, &
         &                  igrdi, igrdj, ppco2n, zpco2l )
      
      DO jobs = pco2datqc%nsurfup + 1, pco2datqc%nsurfup + ipco2
         
         iobs = jobs - pco2datqc%nsurfup
         
         IF ( kt /= pco2datqc%mstp(jobs) ) THEN
            
            IF(lwp) THEN
               WRITE(numout,*)
               WRITE(numout,*) ' E R R O R : Observation',              &
                  &            ' time step is not consistent with the', &
                  &            ' model time step'
               WRITE(numout,*) ' ========='
               WRITE(numout,*)
               WRITE(numout,*) ' Record  = ', jobs,                &
                  &            ' kt      = ', kt,                  &
                  &            ' mstp    = ', pco2datqc%mstp(jobs), &
                  &            ' ntyp    = ', pco2datqc%ntyp(jobs)
            ENDIF
            CALL ctl_stop( 'obs_pco2_opt', 'Inconsistent time' )
            
         ENDIF
         
         zlam = pco2datqc%rlam(jobs)
         zphi = pco2datqc%rphi(jobs)
         
         ! Get weights to interpolate the model pco2 to the observation point
         CALL obs_int_h2d_init( 1, 1, k2dint, zlam, zphi,         &
            &                   zglam(:,:,iobs), zgphi(:,:,iobs), &
            &                   zmask(:,:,iobs), zweig, zobsmask )
         
         ! ... Interpolate the model pco2 to the observation point
         CALL obs_int_h2d( 1, 1,      &
            &              zweig, zpco2l(:,:,iobs),  zext )
         
         pco2datqc%rmod(jobs,1) = zext(1)
         
      END DO
      
      ! Deallocate the data for interpolation
      DEALLOCATE( &
         & igrdi,    &
         & igrdj,    &
         & zglam,    &
         & zgphi,    &
         & zmask,    &
         & zpco2l  &
         & )
      
      pco2datqc%nsurfup = pco2datqc%nsurfup + ipco2

   END SUBROUTINE obs_pco2_opt

END MODULE obs_oper

