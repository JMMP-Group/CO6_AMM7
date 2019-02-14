MODULE diaobs
   !!======================================================================
   !!                       ***  MODULE diaobs  ***
   !! Observation diagnostics: Computation of the misfit between data and
   !!                          their model equivalent 
   !!======================================================================

   !!----------------------------------------------------------------------
   !!   'key_diaobs' : Switch on the observation diagnostic computation
   !!----------------------------------------------------------------------
   !!   dia_obs_init : Reading and prepare observations
   !!   dia_obs      : Compute model equivalent to observations
   !!   dia_obs_wri  : Write observational diagnostics
   !!   ini_date     : Compute the initial date YYYYMMDD.HHMMSS
   !!   fin_date     : Compute the final date YYYYMMDD.HHMMSS
   !!----------------------------------------------------------------------
   !! * Modules used   
   USE wrk_nemo                 ! Memory Allocation
   USE par_kind                 ! Precision variables
   USE in_out_manager           ! I/O manager
   USE par_oce
   USE dom_oce                  ! Ocean space and time domain variables
   USE obs_const, ONLY: obfillflt ! Fill value
   USE obs_fbm, ONLY: ln_cl4    ! Class 4 diagnostic switch
   USE obs_read_prof            ! Reading and allocation of observations (Coriolis)
   USE obs_read_sla             ! Reading and allocation of SLA observations  
   USE obs_read_sst             ! Reading and allocation of SST observations  
   USE obs_sstbias              ! Bias correction routine for SST
   USE obs_readmdt              ! Reading and allocation of MDT for SLA.
   USE obs_read_seaice          ! Reading and allocation of Sea Ice observations  
   USE obs_read_vel             ! Reading and allocation of velocity component observations
   USE obs_read_logchl          ! Reading and allocation of logchl observations
   USE obs_read_spm             ! Reading and allocation of spm observations
   USE obs_read_fco2            ! Reading and allocation of fco2 observations
   USE obs_read_pco2            ! Reading and allocation of pco2 observations
   USE obs_prep                 ! Preparation of obs. (grid search etc).
   USE obs_oper                 ! Observation operators
   USE obs_write                ! Writing of observation related diagnostics
   USE obs_grid                 ! Grid searching
   USE obs_read_altbias         ! Bias treatment for altimeter
   USE obs_profiles_def         ! Profile data definitions
   USE obs_profiles             ! Profile data storage
   USE obs_surf_def             ! Surface data definitions
   USE obs_sla                  ! SLA data storage
   USE obs_sst                  ! SST data storage
   USE obs_seaice               ! Sea Ice data storage
   USE obs_logchl               ! logchl data storage
   USE obs_spm                  ! spm data storage
   USE obs_fco2                 ! fco2 data storage
   USE obs_pco2                 ! pco2 data storage
   USE obs_types                ! Definitions for observation types
   USE mpp_map                  ! MPP mapping
   USE lib_mpp                  ! For ctl_warn/stop

   IMPLICIT NONE

   !! * Routine accessibility
   PRIVATE
   PUBLIC dia_obs_init, &  ! Initialize and read observations
      &   dia_obs,      &  ! Compute model equivalent to observations
      &   dia_obs_wri,  &  ! Write model equivalent to observations
      &   dia_obs_dealloc  ! Deallocate dia_obs data

   !! * Shared Module variables
   LOGICAL, PUBLIC, PARAMETER :: &
#if defined key_diaobs
      & lk_diaobs = .TRUE.   !: Logical switch for observation diangostics
#else
      & lk_diaobs = .FALSE.  !: Logical switch for observation diangostics
#endif

   !! * Module variables
   LOGICAL, PUBLIC :: ln_t3d         !: Logical switch for temperature profiles
   LOGICAL, PUBLIC :: ln_s3d         !: Logical switch for salinity profiles
   LOGICAL, PUBLIC :: ln_ena         !: Logical switch for the ENACT data set
   LOGICAL, PUBLIC :: ln_cor         !: Logical switch for the Coriolis data set
   LOGICAL, PUBLIC :: ln_profb       !: Logical switch for profile feedback datafiles
   LOGICAL, PUBLIC :: ln_sla         !: Logical switch for sea level anomalies 
   LOGICAL, PUBLIC :: ln_sladt       !: Logical switch for SLA from AVISO files
   LOGICAL, PUBLIC :: ln_slafb       !: Logical switch for SLA from feedback files
   LOGICAL, PUBLIC :: ln_sst         !: Logical switch for sea surface temperature
   LOGICAL, PUBLIC :: ln_reysst      !: Logical switch for Reynolds sea surface temperature
   LOGICAL, PUBLIC :: ln_ghrsst      !: Logical switch for GHRSST data
   LOGICAL, PUBLIC :: ln_sstfb       !: Logical switch for SST from feedback files
   LOGICAL, PUBLIC :: ln_seaice      !: Logical switch for sea ice concentration
   LOGICAL, PUBLIC :: ln_vel3d       !: Logical switch for velocity component (u,v) observations
   LOGICAL, PUBLIC :: ln_velavcur    !: Logical switch for raw daily averaged netCDF current meter vel. data 
   LOGICAL, PUBLIC :: ln_velhrcur    !: Logical switch for raw high freq netCDF current meter vel. data 
   LOGICAL, PUBLIC :: ln_velavadcp   !: Logical switch for raw daily averaged netCDF ADCP vel. data 
   LOGICAL, PUBLIC :: ln_velhradcp   !: Logical switch for raw high freq netCDF ADCP vel. data 
   LOGICAL, PUBLIC :: ln_velfb       !: Logical switch for velocities from feedback files
   LOGICAL, PUBLIC :: ln_logchl      !: Logical switch for log10(chlorophyll)
   LOGICAL, PUBLIC :: ln_logchlfb    !: Logical switch for logchl from feedback files
   LOGICAL, PUBLIC :: ln_spm         !: Logical switch for spm
   LOGICAL, PUBLIC :: ln_spmfb       !: Logical switch for spm from feedback files
   LOGICAL, PUBLIC :: ln_fco2        !: Logical switch for fco2
   LOGICAL, PUBLIC :: ln_fco2fb      !: Logical switch for fco2 from feedback files
   LOGICAL, PUBLIC :: ln_pco2        !: Logical switch for pco2
   LOGICAL, PUBLIC :: ln_pco2fb      !: Logical switch for pco2 from feedback files
   LOGICAL, PUBLIC :: ln_ssh         !: Logical switch for sea surface height
   LOGICAL, PUBLIC :: ln_sss         !: Logical switch for sea surface salinity
   LOGICAL, PUBLIC :: ln_sstnight    !: Logical switch for night mean SST observations
   LOGICAL, PUBLIC :: ln_nea         !: Remove observations near land
   LOGICAL, PUBLIC :: ln_altbias     !: Logical switch for altimeter bias  
   LOGICAL, PUBLIC :: ln_ignmis      !: Logical switch for ignoring missing files
   LOGICAL, PUBLIC :: ln_s_at_t      !: Logical switch to compute model S at T observations
   LOGICAL, PUBLIC :: ln_sstbias     !: Logical switch for bias corection of SST

   REAL(KIND=dp), PUBLIC :: dobsini   !: Observation window start date YYYYMMDD.HHMMSS
   REAL(KIND=dp), PUBLIC :: dobsend   !: Observation window end date YYYYMMDD.HHMMSS
  
   INTEGER, PUBLIC :: n1dint       !: Vertical interpolation method
   INTEGER, PUBLIC :: n2dint       !: Horizontal interpolation method 

   INTEGER, DIMENSION(imaxavtypes) :: &
      & endailyavtypes !: ENACT data types which are daily average

   INTEGER, PARAMETER :: MaxNumFiles = 1000
   LOGICAL, DIMENSION(MaxNumFiles) :: &
      & ln_profb_ena, & !: Is the feedback files from ENACT data ?
   !                    !: If so use endailyavtypes
      & ln_profb_enatim !: Change tim for 820 enact data set.
   
   INTEGER, DIMENSION(MaxNumFiles), PUBLIC :: sstbias_type !SST bias type

   LOGICAL, DIMENSION(MaxNumFiles) :: &
      & ln_velfb_av   !: Is the velocity feedback files daily average?
   LOGICAL, DIMENSION(:), ALLOCATABLE :: &
      & ld_enact     !: Profile data is ENACT so use endailyavtypes
   LOGICAL, DIMENSION(:), ALLOCATABLE :: &
      & ld_velav     !: Velocity data is daily averaged
   LOGICAL, DIMENSION(:), ALLOCATABLE :: &
      & ld_sstnight  !: SST observation corresponds to night mean

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

   !! * Substitutions 
#  include "domzgr_substitute.h90"
CONTAINS

   SUBROUTINE dia_obs_init
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE dia_obs_init  ***
      !!          
      !! ** Purpose : Initialize and read observations
      !!
      !! ** Method  : Read the namelist and call reading routines
      !!
      !! ** Action  : Read the namelist and call reading routines
      !!
      !! History :
      !!        !  06-03  (K. Mogensen) Original code
      !!        !  06-05  (A. Weaver) Reformatted
      !!        !  06-10  (A. Weaver) Cleaning and add controls
      !!        !  07-03  (K. Mogensen) General handling of profiles
      !!        !  14-08  (J.While) Incorporated SST bias correction
      !!----------------------------------------------------------------------

      IMPLICIT NONE

      !! * Local declarations
      CHARACTER(len=128) :: enactfiles(MaxNumFiles)
      CHARACTER(len=128) :: coriofiles(MaxNumFiles)
      CHARACTER(len=128) :: profbfiles(MaxNumFiles)
      CHARACTER(len=128) :: sstfiles(MaxNumFiles)      
      CHARACTER(len=128) :: sstfbfiles(MaxNumFiles)
      CHARACTER(len=128) :: sstbias_files(MaxNumFiles) 
      CHARACTER(len=128) :: slafilesact(MaxNumFiles)      
      CHARACTER(len=128) :: slafilespas(MaxNumFiles)      
      CHARACTER(len=128) :: slafbfiles(MaxNumFiles)
      CHARACTER(len=128) :: seaicefiles(MaxNumFiles)           
      CHARACTER(len=128) :: velcurfiles(MaxNumFiles)  
      CHARACTER(len=128) :: veladcpfiles(MaxNumFiles)    
      CHARACTER(len=128) :: velavcurfiles(MaxNumFiles)
      CHARACTER(len=128) :: velhrcurfiles(MaxNumFiles)
      CHARACTER(len=128) :: velavadcpfiles(MaxNumFiles)
      CHARACTER(len=128) :: velhradcpfiles(MaxNumFiles)
      CHARACTER(len=128) :: velfbfiles(MaxNumFiles)
      CHARACTER(len=128) :: logchlfiles(MaxNumFiles)
      CHARACTER(len=128) :: logchlfbfiles(MaxNumFiles)
      CHARACTER(len=128) :: spmfiles(MaxNumFiles)
      CHARACTER(len=128) :: spmfbfiles(MaxNumFiles)
      CHARACTER(len=128) :: fco2files(MaxNumFiles)
      CHARACTER(len=128) :: fco2fbfiles(MaxNumFiles)
      CHARACTER(len=128) :: pco2files(MaxNumFiles)
      CHARACTER(len=128) :: pco2fbfiles(MaxNumFiles)
      CHARACTER(LEN=128) :: reysstname
      CHARACTER(LEN=12)  :: reysstfmt
      CHARACTER(LEN=128) :: bias_file
      CHARACTER(LEN=20)  :: datestr=" ", timestr=" "
      NAMELIST/namobs/ln_ena, ln_cor, ln_profb, ln_t3d, ln_s3d,       &
         &            ln_sla, ln_sladt, ln_slafb,                     &
         &            ln_ssh, ln_sst, ln_sstfb, ln_sss, ln_nea,       &
         &            ln_bound_reject,                                &
         &            enactfiles, coriofiles, profbfiles,             &
         &            slafilesact, slafilespas, slafbfiles,           &
         &            sstfiles, sstfbfiles,                           &
         &            ln_seaice, seaicefiles,                         &
         &            dobsini, dobsend, n1dint, n2dint,               &
         &            nmsshc, mdtcorr, mdtcutoff,                     &
         &            ln_reysst, ln_ghrsst, reysstname, reysstfmt,    &
         &            ln_sstnight,                                    &
         &            ln_grid_search_lookup,                          &
         &            grid_search_file, grid_search_res,              &
         &            ln_grid_global, bias_file, ln_altbias,          &
         &            endailyavtypes, ln_s_at_t, ln_profb_ena,        &
         &            ln_vel3d, ln_velavcur, velavcurfiles,           &
         &            ln_velhrcur, velhrcurfiles,                     &
         &            ln_velavadcp, velavadcpfiles,                   &
         &            ln_velhradcp, velhradcpfiles,                   &
         &            ln_velfb, velfbfiles, ln_velfb_av,              &
         &            ln_logchl, ln_logchlfb,                         &
         &            logchlfiles, logchlfbfiles,                     &
         &            ln_spm, ln_spmfb,                               &
         &            spmfiles, spmfbfiles,                           &
         &            ln_fco2, ln_fco2fb,                             &
         &            fco2files, fco2fbfiles,                         &
         &            ln_pco2, ln_pco2fb,                             &
         &            pco2files, pco2fbfiles,                         &
         &            ln_profb_enatim, ln_ignmis, ln_cl4,             &
         &            ln_sstbias, sstbias_files

      INTEGER :: jprofset
      INTEGER :: jveloset
      INTEGER :: jvar
      INTEGER :: jnumenact
      INTEGER :: jnumcorio
      INTEGER :: jnumprofb
      INTEGER :: jnumslaact
      INTEGER :: jnumslapas
      INTEGER :: jnumslafb
      INTEGER :: jnumsst
      INTEGER :: jnumsstfb
      INTEGER :: jnumsstbias
      INTEGER :: jnumseaice
      INTEGER :: jnumvelavcur
      INTEGER :: jnumvelhrcur  
      INTEGER :: jnumvelavadcp
      INTEGER :: jnumvelhradcp   
      INTEGER :: jnumvelfb
      INTEGER :: jnumlogchl
      INTEGER :: jnumlogchlfb
      INTEGER :: jnumspm
      INTEGER :: jnumspmfb
      INTEGER :: jnumfco2
      INTEGER :: jnumfco2fb
      INTEGER :: jnumpco2
      INTEGER :: jnumpco2fb
      INTEGER :: ji
      INTEGER :: jset
      INTEGER :: ios                 ! Local integer output status for namelist read
      LOGICAL :: lmask(MaxNumFiles), ll_u3d, ll_v3d

      !-----------------------------------------------------------------------
      ! Read namelist parameters
      !-----------------------------------------------------------------------

      ln_logchl   = .FALSE.
      ln_logchlfb = .FALSE.
      ln_spm      = .FALSE.
      ln_spmfb    = .FALSE.
      ln_fco2     = .FALSE.
      ln_fco2fb   = .FALSE.
      ln_pco2     = .FALSE.
      ln_pco2fb   = .FALSE.
      
      !Initalise all values in namelist arrays
      enactfiles(:) = ''
      coriofiles(:) = ''
      profbfiles(:) = ''
      slafilesact(:) = ''
      slafilespas(:) = ''
      slafbfiles(:) = ''
      sstfiles(:)   = ''
      sstfbfiles(:) = ''
      seaicefiles(:) = ''
      velcurfiles(:) = ''
      veladcpfiles(:) = ''
      velavcurfiles(:) = ''
      velhrcurfiles(:) = ''
      velavadcpfiles(:) = ''
      velhradcpfiles(:) = ''
      velfbfiles(:) = ''
      velcurfiles(:) = ''
      veladcpfiles(:) = ''
      logchlfiles(:) = ''
      logchlfbfiles(:) = ''
      spmfiles(:) = ''
      spmfbfiles(:) = ''
      fco2files(:) = ''
      fco2fbfiles(:) = ''
      pco2files(:) = ''
      pco2fbfiles(:) = ''
      sstbias_files(:) = ''
      endailyavtypes(:) = -1
      endailyavtypes(1) = 820
      ln_profb_ena(:) = .FALSE.
      ln_profb_enatim(:) = .TRUE.
      ln_velfb_av(:) = .FALSE.
      ln_ignmis = .FALSE.
      ln_bound_reject = .TRUE.

      ! Read Namelist namobs : control observation diagnostics
      REWIND( numnam_ref )              ! Namelist namobs in reference namelist : Diagnostic: control observation
      READ  ( numnam_ref, namobs, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namobs in reference namelist', lwp )

      CALL ini_date( dobsini )
      CALL fin_date( dobsend )
 
      REWIND( numnam_cfg )              ! Namelist namobs in configuration namelist : Diagnostic: control observation
      READ  ( numnam_cfg, namobs, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namobs in configuration namelist', lwp )
      IF(lwm) WRITE ( numond, namobs )

      ! Count number of files for each type
      IF (ln_ena) THEN
         lmask(:) = .FALSE.
         WHERE (enactfiles(:) /= '') lmask(:) = .TRUE.
         jnumenact = COUNT(lmask)
      ENDIF
      IF (ln_cor) THEN
         lmask(:) = .FALSE.
         WHERE (coriofiles(:) /= '') lmask(:) = .TRUE.
         jnumcorio = COUNT(lmask)
      ENDIF
      IF (ln_profb) THEN
         lmask(:) = .FALSE.
         WHERE (profbfiles(:) /= '') lmask(:) = .TRUE.
         jnumprofb = COUNT(lmask)
      ENDIF
      IF (ln_sladt) THEN
         lmask(:) = .FALSE.
         WHERE (slafilesact(:) /= '') lmask(:) = .TRUE.
         jnumslaact = COUNT(lmask)
         lmask(:) = .FALSE.
         WHERE (slafilespas(:) /= '') lmask(:) = .TRUE.
         jnumslapas = COUNT(lmask)
      ENDIF
      IF (ln_slafb) THEN
         lmask(:) = .FALSE.
         WHERE (slafbfiles(:) /= '') lmask(:) = .TRUE.
         jnumslafb = COUNT(lmask)
         lmask(:) = .FALSE.
      ENDIF
      IF (ln_ghrsst) THEN
         lmask(:) = .FALSE.
         WHERE (sstfiles(:) /= '') lmask(:) = .TRUE.
         jnumsst = COUNT(lmask)
      ENDIF      
      IF (ln_sstfb) THEN
         lmask(:) = .FALSE.
         WHERE (sstfbfiles(:) /= '') lmask(:) = .TRUE.
         jnumsstfb = COUNT(lmask)
         lmask(:) = .FALSE.
      ENDIF
      IF (ln_sstbias) THEN  
         lmask(:) = .FALSE.  
         WHERE (sstbias_files(:) /= '') lmask(:) = .TRUE.  
         jnumsstbias = COUNT(lmask)  
         lmask(:) = .FALSE.  
      ENDIF       
      IF (ln_seaice) THEN
         lmask(:) = .FALSE.
         WHERE (seaicefiles(:) /= '') lmask(:) = .TRUE.
         jnumseaice = COUNT(lmask)
      ENDIF
      IF (ln_velavcur) THEN
         lmask(:) = .FALSE.
         WHERE (velavcurfiles(:) /= '') lmask(:) = .TRUE.
         jnumvelavcur = COUNT(lmask)
      ENDIF
      IF (ln_velhrcur) THEN
         lmask(:) = .FALSE.
         WHERE (velhrcurfiles(:) /= '') lmask(:) = .TRUE.
         jnumvelhrcur = COUNT(lmask)
      ENDIF
      IF (ln_velavadcp) THEN
         lmask(:) = .FALSE.
         WHERE (velavadcpfiles(:) /= '') lmask(:) = .TRUE.
         jnumvelavadcp = COUNT(lmask)
      ENDIF
      IF (ln_velhradcp) THEN
         lmask(:) = .FALSE.
         WHERE (velhradcpfiles(:) /= '') lmask(:) = .TRUE.
         jnumvelhradcp = COUNT(lmask)
      ENDIF
      IF (ln_velfb) THEN
         lmask(:) = .FALSE.
         WHERE (velfbfiles(:) /= '') lmask(:) = .TRUE.
         jnumvelfb = COUNT(lmask)
         lmask(:) = .FALSE.
      ENDIF
      IF (ln_logchl) THEN
         lmask(:) = .FALSE.
         WHERE (logchlfiles(:) /= '') lmask(:) = .TRUE.
         jnumlogchl = COUNT(lmask)
      ENDIF
      IF (ln_logchlfb) THEN
         lmask(:) = .FALSE.
         WHERE (logchlfbfiles(:) /= '') lmask(:) = .TRUE.
         jnumlogchlfb = COUNT(lmask)
      ENDIF
      IF (ln_spm) THEN
         lmask(:) = .FALSE.
         WHERE (spmfiles(:) /= '') lmask(:) = .TRUE.
         jnumspm = COUNT(lmask)
      ENDIF
      IF (ln_spmfb) THEN
         lmask(:) = .FALSE.
         WHERE (spmfbfiles(:) /= '') lmask(:) = .TRUE.
         jnumspmfb = COUNT(lmask)
      ENDIF
      IF (ln_fco2) THEN
         lmask(:) = .FALSE.
         WHERE (fco2files(:) /= '') lmask(:) = .TRUE.
         jnumfco2 = COUNT(lmask)
      ENDIF
      IF (ln_fco2fb) THEN
         lmask(:) = .FALSE.
         WHERE (fco2fbfiles(:) /= '') lmask(:) = .TRUE.
         jnumfco2fb = COUNT(lmask)
      ENDIF
      IF (ln_pco2) THEN
         lmask(:) = .FALSE.
         WHERE (pco2files(:) /= '') lmask(:) = .TRUE.
         jnumpco2 = COUNT(lmask)
      ENDIF
      IF (ln_pco2fb) THEN
         lmask(:) = .FALSE.
         WHERE (pco2fbfiles(:) /= '') lmask(:) = .TRUE.
         jnumpco2fb = COUNT(lmask)
      ENDIF
      
      ! Control print
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'dia_obs_init : Observation diagnostic initialization'
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '          Namelist namobs : set observation diagnostic parameters' 
         WRITE(numout,*) '             Logical switch for T profile observations          ln_t3d = ', ln_t3d
         WRITE(numout,*) '             Logical switch for S profile observations          ln_s3d = ', ln_s3d
         WRITE(numout,*) '             Logical switch for ENACT insitu data set           ln_ena = ', ln_ena
         WRITE(numout,*) '             Logical switch for Coriolis insitu data set        ln_cor = ', ln_cor
         WRITE(numout,*) '             Logical switch for feedback insitu data set      ln_profb = ', ln_profb
         WRITE(numout,*) '             Logical switch for SLA observations                ln_sla = ', ln_sla
         WRITE(numout,*) '             Logical switch for AVISO SLA data                ln_sladt = ', ln_sladt
         WRITE(numout,*) '             Logical switch for feedback SLA data             ln_slafb = ', ln_slafb
         WRITE(numout,*) '             Logical switch for SSH observations                ln_ssh = ', ln_ssh
         WRITE(numout,*) '             Logical switch for SST observations                ln_sst = ', ln_sst
         WRITE(numout,*) '             Logical switch for Reynolds observations        ln_reysst = ', ln_reysst    
         WRITE(numout,*) '             Logical switch for GHRSST observations          ln_ghrsst = ', ln_ghrsst
         WRITE(numout,*) '             Logical switch for feedback SST data             ln_sstfb = ', ln_sstfb
         WRITE(numout,*) '             Logical switch for SST bias correction         ln_sstbias = ', ln_sstbias
         WRITE(numout,*) '             Logical switch for night-time SST obs         ln_sstnight = ', ln_sstnight
         WRITE(numout,*) '             Logical switch for SSS observations                ln_sss = ', ln_sss
         WRITE(numout,*) '             Logical switch for Sea Ice observations         ln_seaice = ', ln_seaice
         WRITE(numout,*) '             Logical switch for velocity observations         ln_vel3d = ', ln_vel3d
         WRITE(numout,*) '             Logical switch for velocity daily av. cur.    ln_velavcur = ', ln_velavcur
         WRITE(numout,*) '             Logical switch for velocity high freq. cur.   ln_velhrcur = ', ln_velhrcur
         WRITE(numout,*) '             Logical switch for velocity daily av. ADCP   ln_velavadcp = ', ln_velavadcp
         WRITE(numout,*) '             Logical switch for velocity high freq. ADCP  ln_velhradcp = ', ln_velhradcp
         WRITE(numout,*) '             Logical switch for feedback velocity data        ln_velfb = ', ln_velfb
         WRITE(numout,*) '             Logical switch for logchl observations          ln_logchl = ', ln_logchl
         WRITE(numout,*) '             Logical switch for feedback logchl data       ln_logchlfb = ', ln_logchlfb
         WRITE(numout,*) '             Logical switch for spm observations                ln_spm = ', ln_spm
         WRITE(numout,*) '             Logical switch for feedback spm data             ln_spmfb = ', ln_spmfb
         WRITE(numout,*) '             Logical switch for fco2 observations              ln_fco2 = ', ln_fco2
         WRITE(numout,*) '             Logical switch for pco2 observations              ln_pco2 = ', ln_pco2
         WRITE(numout,*) '             Logical switch for feedback pco2 data           ln_pco2fb = ', ln_pco2fb
         WRITE(numout,*) '             Logical switch for feedback fco2 data           ln_fco2fb = ', ln_fco2fb
         WRITE(numout,*) '             Global distribtion of observations         ln_grid_global = ',ln_grid_global
         WRITE(numout,*) &
   '             Logical switch for obs grid search w/lookup table  ln_grid_search_lookup = ',ln_grid_search_lookup
         IF (ln_grid_search_lookup) &
            WRITE(numout,*) '             Grid search lookup file header       grid_search_file = ', grid_search_file
         IF (ln_ena) THEN
            DO ji = 1, jnumenact
               WRITE(numout,'(1X,2A)') '             ENACT input observation file name          enactfiles = ', &
                  TRIM(enactfiles(ji))
            END DO
         ENDIF
         IF (ln_cor) THEN
            DO ji = 1, jnumcorio
               WRITE(numout,'(1X,2A)') '             Coriolis input observation file name       coriofiles = ', &
                  TRIM(coriofiles(ji))
            END DO
         ENDIF
         IF (ln_profb) THEN
            DO ji = 1, jnumprofb
               IF (ln_profb_ena(ji)) THEN
                  WRITE(numout,'(1X,2A)') '       Enact feedback input observation file name       profbfiles = ', &
                     TRIM(profbfiles(ji))
               ELSE
                  WRITE(numout,'(1X,2A)') '             Feedback input observation file name       profbfiles = ', &
                     TRIM(profbfiles(ji))
               ENDIF
               WRITE(numout,'(1X,2A)') '       Enact feedback input time setting switch    ln_profb_enatim = ', ln_profb_enatim(ji)
            END DO
         ENDIF
         IF (ln_sladt) THEN
            DO ji = 1, jnumslaact
               WRITE(numout,'(1X,2A)') '             Active SLA input observation file name    slafilesact = ', &
                  TRIM(slafilesact(ji))
            END DO
            DO ji = 1, jnumslapas
               WRITE(numout,'(1X,2A)') '             Passive SLA input observation file name   slafilespas = ', &
                  TRIM(slafilespas(ji))
            END DO
         ENDIF
         IF (ln_slafb) THEN
            DO ji = 1, jnumslafb
               WRITE(numout,'(1X,2A)') '             Feedback SLA input observation file name   slafbfiles = ', &
                  TRIM(slafbfiles(ji))
            END DO
         ENDIF
         IF (ln_ghrsst) THEN
            DO ji = 1, jnumsst
               WRITE(numout,'(1X,2A)') '             GHRSST input observation file name           sstfiles = ', &
                  TRIM(sstfiles(ji))
            END DO
         ENDIF
         IF (ln_sstfb) THEN
            DO ji = 1, jnumsstfb
               WRITE(numout,'(1X,2A)') '             Feedback SST input observation file name   sstfbfiles = ', &
                  TRIM(sstfbfiles(ji))
            END DO
         ENDIF
         IF (ln_seaice) THEN
            DO ji = 1, jnumseaice
               WRITE(numout,'(1X,2A)') '             Sea Ice input observation file name       seaicefiles = ', &
                  TRIM(seaicefiles(ji))
            END DO
         ENDIF
         IF (ln_velavcur) THEN
            DO ji = 1, jnumvelavcur
               WRITE(numout,'(1X,2A)') '             Vel. cur. daily av. input file name     velavcurfiles = ', &
                  TRIM(velavcurfiles(ji))
            END DO
         ENDIF
         IF (ln_velhrcur) THEN
            DO ji = 1, jnumvelhrcur
               WRITE(numout,'(1X,2A)') '             Vel. cur. high freq. input file name    velhvcurfiles = ', &
                  TRIM(velhrcurfiles(ji))
            END DO
         ENDIF
         IF (ln_velavadcp) THEN
            DO ji = 1, jnumvelavadcp
               WRITE(numout,'(1X,2A)') '             Vel. ADCP daily av. input file name    velavadcpfiles = ', &
                  TRIM(velavadcpfiles(ji))
            END DO
         ENDIF
         IF (ln_velhradcp) THEN
            DO ji = 1, jnumvelhradcp
               WRITE(numout,'(1X,2A)') '             Vel. ADCP high freq. input file name   velhvadcpfiles = ', &
                  TRIM(velhradcpfiles(ji))
            END DO
         ENDIF
         IF (ln_velfb) THEN
            DO ji = 1, jnumvelfb
               IF (ln_velfb_av(ji)) THEN
                  WRITE(numout,'(1X,2A)') '             Vel. feedback daily av. input file name    velfbfiles = ', &
                     TRIM(velfbfiles(ji))
               ELSE
                  WRITE(numout,'(1X,2A)') '             Vel. feedback input observation file name  velfbfiles = ', &
                     TRIM(velfbfiles(ji))
               ENDIF
            END DO
         ENDIF
         IF (ln_logchl) THEN
            DO ji = 1, jnumlogchl
               WRITE(numout,'(1X,2A)') '             logchl input observation file name        logchlfiles = ', &
                  TRIM(logchlfiles(ji))
            END DO
         ENDIF
         IF (ln_logchlfb) THEN
            DO ji = 1, jnumlogchlfb
               WRITE(numout,'(1X,2A)') '        Feedback logchl input observation file name  logchlfbfiles = ', &
                  TRIM(logchlfbfiles(ji))
            END DO
         ENDIF
         IF (ln_spm) THEN
            DO ji = 1, jnumspm
               WRITE(numout,'(1X,2A)') '             spm input observation file name              spmfiles = ', &
                  TRIM(spmfiles(ji))
            END DO
         ENDIF
         IF (ln_spmfb) THEN
            DO ji = 1, jnumspmfb
               WRITE(numout,'(1X,2A)') '             Feedback spm input observation file name   spmfbfiles = ', &
                  TRIM(spmfbfiles(ji))
            END DO
         ENDIF
         IF (ln_fco2) THEN
            DO ji = 1, jnumfco2
               WRITE(numout,'(1X,2A)') '             fco2 input observation file name            fco2files = ', &
                  TRIM(fco2files(ji))
            END DO
         ENDIF
         IF (ln_fco2fb) THEN
            DO ji = 1, jnumfco2fb
               WRITE(numout,'(1X,2A)') '            Feedback fco2 input observation file name  fco2fbfiles = ', &
                  TRIM(fco2fbfiles(ji))
            END DO
         ENDIF
         IF (ln_pco2) THEN
            DO ji = 1, jnumpco2
               WRITE(numout,'(1X,2A)') '             pco2 input observation file name            pco2files = ', &
                  TRIM(pco2files(ji))
            END DO
         ENDIF
         IF (ln_pco2fb) THEN
            DO ji = 1, jnumpco2fb
               WRITE(numout,'(1X,2A)') '            Feedback pco2 input observation file name  pco2fbfiles = ', &
                  TRIM(pco2fbfiles(ji))
            END DO
         ENDIF
         WRITE(numout,*) '             Initial date in window YYYYMMDD.HHMMSS        dobsini = ', dobsini
         WRITE(numout,*) '             Final date in window YYYYMMDD.HHMMSS          dobsend = ', dobsend
         WRITE(numout,*) '             Type of vertical interpolation method          n1dint = ', n1dint
         WRITE(numout,*) '             Type of horizontal interpolation method        n2dint = ', n2dint
         WRITE(numout,*) '             Rejection of observations near land swithch    ln_nea = ', ln_nea
         WRITE(numout,*) '             Rejection of obs near open bdys       ln_bound_reject = ', ln_bound_reject
         WRITE(numout,*) '             MSSH correction scheme                         nmsshc = ', nmsshc
         WRITE(numout,*) '             MDT  correction                               mdtcorr = ', mdtcorr
         WRITE(numout,*) '             MDT cutoff for computed correction          mdtcutoff = ', mdtcutoff
         WRITE(numout,*) '             Logical switch for alt bias                ln_altbias = ', ln_altbias
         WRITE(numout,*) '             Logical switch for ignoring missing files   ln_ignmis = ', ln_ignmis
         WRITE(numout,*) '             ENACT daily average types                             = ',endailyavtypes

      ENDIF
      
      IF ( ln_vel3d .AND. ( .NOT. ln_grid_global ) ) THEN
         CALL ctl_stop( 'Velocity data only works with ln_grid_global=.true.' )
         RETURN
      ENDIF

      IF ( ln_grid_global ) THEN
         CALL ctl_warn( 'ln_grid_global=T may cause memory issues when used with a large number of processors' )
      ENDIF

      CALL obs_typ_init
      
      IF ( ln_grid_global ) THEN      
         CALL mppmap_init
      ENDIF
      
      ! Parameter control
#if defined key_diaobs
      IF ( ( .NOT. ln_t3d ).AND.( .NOT. ln_s3d ).AND.( .NOT. ln_sla ).AND. &
         & ( .NOT. ln_vel3d ).AND.                                         &
         & ( .NOT. ln_ssh ).AND.( .NOT. ln_sst ).AND.( .NOT. ln_sss ).AND. &
         & ( .NOT. ln_seaice ).AND.( .NOT. ln_vel3d ).AND.( .NOT. ln_logchl ).AND. &
         & ( .NOT. ln_spm ).AND.( .NOT. ln_fco2 ).AND.( .NOT. ln_pco2 ) ) THEN
         IF(lwp) WRITE(numout,cform_war)
         IF(lwp) WRITE(numout,*) ' key_diaobs is activated but logical flags', &
            &                    ' ln_t3d, ln_s3d, ln_sla, ln_ssh, ln_sst, ln_sss, ln_seaice, ln_vel3d,', &
            &                    ' ln_logchl, ln_spm, ln_fco2, ln_pco2 are all set to .FALSE.'
         nwarn = nwarn + 1
      ENDIF
#endif

      CALL obs_grid_setup( )
      IF ( ( n1dint < 0 ).OR.( n1dint > 1 ) ) THEN
         CALL ctl_stop(' Choice of vertical (1D) interpolation method', &
            &                    ' is not available')
      ENDIF
      IF ( ( n2dint < 0 ).OR.( n2dint > 4 ) ) THEN
         CALL ctl_stop(' Choice of horizontal (2D) interpolation method', &
            &                    ' is not available')
      ENDIF

      !-----------------------------------------------------------------------
      ! Depending on switches read the various observation types
      !-----------------------------------------------------------------------
      !  - Temperature/salinity profiles

      IF ( ln_t3d .OR. ln_s3d ) THEN

         ! Set the number of variables for profiles to 2 (T and S)
         nprofvars = 2
         ! Set the number of extra variables for profiles to 1 (insitu temp).
         nprofextr = 1

         ! Count how may insitu data sets we have and allocate data.
         jprofset = 0
         IF ( ln_ena ) jprofset = jprofset + 1
         IF ( ln_cor ) jprofset = jprofset + 1
         IF ( ln_profb ) jprofset = jprofset + jnumprofb
         nprofsets = jprofset
         IF ( nprofsets > 0 ) THEN
            ALLOCATE(ld_enact(nprofsets))
            ALLOCATE(profdata(nprofsets))
            ALLOCATE(prodatqc(nprofsets))
         ENDIF

         jprofset = 0
          
         ! ENACT insitu data

         IF ( ln_ena ) THEN

            jprofset = jprofset + 1
            
            ld_enact(jprofset) = .TRUE.

            CALL obs_rea_pro_dri( 1, profdata(jprofset),          &
               &                  jnumenact, enactfiles(1:jnumenact), &
               &                  nprofvars, nprofextr,        &
               &                  nitend-nit000+2,             &
               &                  dobsini, dobsend, ln_t3d, ln_s3d, &
               &                  ln_ignmis, ln_s_at_t, .TRUE., .FALSE., &
               &                  kdailyavtypes = endailyavtypes )

            DO jvar = 1, 2

               CALL obs_prof_staend( profdata(jprofset), jvar )

            END DO

            CALL obs_pre_pro( profdata(jprofset), prodatqc(jprofset),   &
               &              ln_t3d, ln_s3d, ln_nea, &
               &              kdailyavtypes=endailyavtypes )
            
         ENDIF

         ! Coriolis insitu data

         IF ( ln_cor ) THEN
           
            jprofset = jprofset + 1

            ld_enact(jprofset) = .FALSE.

            CALL obs_rea_pro_dri( 2, profdata(jprofset),          &
               &                  jnumcorio, coriofiles(1:jnumcorio), &
               &                  nprofvars, nprofextr,        &
               &                  nitend-nit000+2,             &
               &                  dobsini, dobsend, ln_t3d, ln_s3d, &
               &                  ln_ignmis, ln_s_at_t, .FALSE., .FALSE. )

            DO jvar = 1, 2

               CALL obs_prof_staend( profdata(jprofset), jvar )

            END DO

            CALL obs_pre_pro( profdata(jprofset), prodatqc(jprofset),   &
                 &            ln_t3d, ln_s3d, ln_nea )
            
         ENDIF
 
         ! Feedback insitu data

         IF ( ln_profb ) THEN
           
            DO jset = 1, jnumprofb
               
               jprofset = jprofset + 1
               ld_enact (jprofset) = ln_profb_ena(jset)

               CALL obs_rea_pro_dri( 0, profdata(jprofset),          &
                  &                  1, profbfiles(jset:jset), &
                  &                  nprofvars, nprofextr,        &
                  &                  nitend-nit000+2,             &
                  &                  dobsini, dobsend, ln_t3d, ln_s3d, &
                  &                  ln_ignmis, ln_s_at_t, &
                  &                  ld_enact(jprofset).AND.&
                  &                  ln_profb_enatim(jset), &
                  &                  .FALSE., kdailyavtypes = endailyavtypes )
               
               DO jvar = 1, 2
                  
                  CALL obs_prof_staend( profdata(jprofset), jvar )
                  
               END DO
               
               IF ( ld_enact(jprofset) ) THEN
                  CALL obs_pre_pro( profdata(jprofset), prodatqc(jprofset),   &
                     &              ln_t3d, ln_s3d, ln_nea, &
                     &              kdailyavtypes = endailyavtypes )
               ELSE
                  CALL obs_pre_pro( profdata(jprofset), prodatqc(jprofset),   &
                     &              ln_t3d, ln_s3d, ln_nea )
               ENDIF
               
            END DO

         ENDIF

      ENDIF

      !  - Sea level anomalies
      IF ( ln_sla ) THEN
        ! Set the number of variables for sla to 1
         nslavars = 1

         ! Set the number of extra variables for sla to 2
         nslaextr = 2
         
         ! Set the number of sla data sets to 2
         nslasets = 0
         IF ( ln_sladt ) THEN
            nslasets = nslasets + 2
         ENDIF
         IF ( ln_slafb ) THEN
            nslasets = nslasets + jnumslafb
         ENDIF
         
         ALLOCATE(sladata(nslasets))
         ALLOCATE(sladatqc(nslasets))
         sladata(:)%nsurf=0
         sladatqc(:)%nsurf=0

         nslasets = 0

         ! AVISO SLA data

         IF ( ln_sladt ) THEN

            ! Active SLA observations
            
            nslasets = nslasets + 1
            
            CALL obs_rea_sla( 1, sladata(nslasets), jnumslaact, &
               &              slafilesact(1:jnumslaact), &
               &              nslavars, nslaextr, nitend-nit000+2, &
               &              dobsini, dobsend, ln_ignmis, .FALSE. )
            CALL obs_pre_sla( sladata(nslasets), sladatqc(nslasets), &
               &              ln_sla, ln_nea )
            
            ! Passive SLA observations
            
            nslasets = nslasets + 1
            
            CALL obs_rea_sla( 1, sladata(nslasets), jnumslapas, &
               &              slafilespas(1:jnumslapas), &
               &              nslavars, nslaextr, nitend-nit000+2, &
               &              dobsini, dobsend, ln_ignmis, .FALSE. )
            
            CALL obs_pre_sla( sladata(nslasets), sladatqc(nslasets), &
               &              ln_sla, ln_nea )

         ENDIF
         
         ! Feedback SLA data

         IF ( ln_slafb ) THEN

            DO jset = 1, jnumslafb
            
               nslasets = nslasets + 1
            
               CALL obs_rea_sla( 0, sladata(nslasets), 1, &
                  &              slafbfiles(jset:jset), &
                  &              nslavars, nslaextr, nitend-nit000+2, &
                  &              dobsini, dobsend, ln_ignmis, .FALSE. )
               CALL obs_pre_sla( sladata(nslasets), sladatqc(nslasets), &
                  &              ln_sla, ln_nea )

            END DO               

         ENDIF
         
         CALL obs_rea_mdt( nslasets, sladatqc, n2dint )
            
         ! read in altimeter bias
         
         IF ( ln_altbias ) THEN     
            CALL obs_rea_altbias ( nslasets, sladatqc, n2dint, bias_file )
         ENDIF
     
      ENDIF

      !  - Sea surface height
      IF ( ln_ssh ) THEN
         IF(lwp) WRITE(numout,*) ' SSH currently not available'
      ENDIF

      !  - Sea surface temperature
      IF ( ln_sst ) THEN

         ! Set the number of variables for sst to 1
         nsstvars = 1

         ! Set the number of extra variables for sst to 0
         nsstextr = 0

         nsstsets = 0

         IF (ln_reysst) nsstsets = nsstsets + 1
         IF (ln_ghrsst) nsstsets = nsstsets + 1
         IF ( ln_sstfb ) THEN
            nsstsets = nsstsets + jnumsstfb
         ENDIF

         ALLOCATE(sstdata(nsstsets))
         ALLOCATE(sstdatqc(nsstsets))
         ALLOCATE(ld_sstnight(nsstsets))
         sstdata(:)%nsurf=0
         sstdatqc(:)%nsurf=0    
         ld_sstnight(:)=.false.

         nsstsets = 0

         IF (ln_reysst) THEN

            nsstsets = nsstsets + 1

            ld_sstnight(nsstsets) = ln_sstnight

            CALL obs_rea_sst_rey( reysstname, reysstfmt, sstdata(nsstsets), &
               &                  nsstvars, nsstextr, &
               &                  nitend-nit000+2, dobsini, dobsend )
            CALL obs_pre_sst( sstdata(nsstsets), sstdatqc(nsstsets), ln_sst, &
               &              ln_nea )

        ENDIF
        
        IF (ln_ghrsst) THEN
        
            nsstsets = nsstsets + 1

            ld_sstnight(nsstsets) = ln_sstnight
          
            CALL obs_rea_sst( 1, sstdata(nsstsets), jnumsst, &
               &              sstfiles(1:jnumsst), &
               &              nsstvars, nsstextr, nitend-nit000+2, &
               &              dobsini, dobsend, ln_ignmis, .FALSE. )
            CALL obs_pre_sst( sstdata(nsstsets), sstdatqc(nsstsets), ln_sst, &
               &              ln_nea )

        ENDIF
               
         ! Feedback SST data

         IF ( ln_sstfb ) THEN

            DO jset = 1, jnumsstfb
            
               nsstsets = nsstsets + 1

               ld_sstnight(nsstsets) = ln_sstnight
            
               CALL obs_rea_sst( 0, sstdata(nsstsets), 1, &
                  &              sstfbfiles(jset:jset), &
                  &              nsstvars, nsstextr, nitend-nit000+2, &
                  &              dobsini, dobsend, ln_ignmis, .FALSE. )
               CALL obs_pre_sst( sstdata(nsstsets), sstdatqc(nsstsets), &
                  &              ln_sst, ln_nea )

            END DO               

         ENDIF
         
         !Read in bias field and correct SST. 
         IF ( ln_sstbias ) THEN 
            IF ( jnumsstbias == 0 ) CALL ctl_stop("ln_sstbias set,"// & 
                                             "  but no bias"// & 
                                             " files to read in")    
            CALL obs_app_sstbias( nsstsets, sstdatqc, n2dint, & 
                                  jnumsstbias, &  
                                  sstbias_files(1:jnumsstbias) ) 
         ENDIF 

      ENDIF

      !  - Sea surface salinity
      IF ( ln_sss ) THEN
         IF(lwp) WRITE(numout,*) ' SSS currently not available'
      ENDIF

      !  - Sea Ice Concentration
      
      IF ( ln_seaice ) THEN

         ! Set the number of variables for seaice to 1
         nseaicevars = 1

         ! Set the number of extra variables for seaice to 0
         nseaiceextr = 0
         
         ! Set the number of data sets to 1
         nseaicesets = 1

         ALLOCATE(seaicedata(nseaicesets))
         ALLOCATE(seaicedatqc(nseaicesets))
         seaicedata(:)%nsurf=0
         seaicedatqc(:)%nsurf=0

         CALL obs_rea_seaice( 1, seaicedata(nseaicesets), jnumseaice, &
            &                 seaicefiles(1:jnumseaice), &
            &                 nseaicevars, nseaiceextr, nitend-nit000+2, &
            &                 dobsini, dobsend, ln_ignmis, .FALSE. )

         CALL obs_pre_seaice( seaicedata(nseaicesets), seaicedatqc(nseaicesets), &
            &                 ln_seaice, ln_nea )
 
      ENDIF

      IF (ln_vel3d) THEN

         ! Set the number of variables for profiles to 2 (U and V)
         nvelovars = 2

         ! Set the number of extra variables for profiles to 2 to store 
         ! rotation parameters
         nveloextr = 2

         jveloset = 0
         
         IF ( ln_velavcur ) jveloset = jveloset + 1
         IF ( ln_velhrcur ) jveloset = jveloset + 1
         IF ( ln_velavadcp ) jveloset = jveloset + 1
         IF ( ln_velhradcp ) jveloset = jveloset + 1
         IF (ln_velfb) jveloset = jveloset + jnumvelfb

         nvelosets = jveloset
         IF ( nvelosets > 0 ) THEN
            ALLOCATE( velodata(nvelosets) )
            ALLOCATE( veldatqc(nvelosets) )
            ALLOCATE( ld_velav(nvelosets) )
         ENDIF
         
         jveloset = 0
         
         ! Daily averaged data

         IF ( ln_velavcur ) THEN
            
            jveloset = jveloset + 1
            
            ld_velav(jveloset) = .TRUE.
            
            CALL obs_rea_vel_dri( 1, velodata(jveloset), jnumvelavcur, &
               &                  velavcurfiles(1:jnumvelavcur), &
               &                  nvelovars, nveloextr, &
               &                  nitend-nit000+2,              &
               &                  dobsini, dobsend, ln_ignmis, &
               &                  ld_velav(jveloset), &
               &                  .FALSE. )
            
            DO jvar = 1, 2
               CALL obs_prof_staend( velodata(jveloset), jvar )
            END DO
            
            CALL obs_pre_vel( velodata(jveloset), veldatqc(jveloset), &
               &              ln_vel3d, ln_nea, ld_velav(jveloset) )
            
         ENDIF

         ! High frequency data

         IF ( ln_velhrcur ) THEN
            
            jveloset = jveloset + 1
            
            ld_velav(jveloset) = .FALSE.
               
            CALL obs_rea_vel_dri( 1, velodata(jveloset), jnumvelhrcur, &
               &                  velhrcurfiles(1:jnumvelhrcur), &
               &                  nvelovars, nveloextr, &
               &                  nitend-nit000+2,              &
               &                  dobsini, dobsend, ln_ignmis, &
               &                  ld_velav(jveloset), &
               &                  .FALSE. )
            
            DO jvar = 1, 2
               CALL obs_prof_staend( velodata(jveloset), jvar )
            END DO
            
            CALL obs_pre_vel( velodata(jveloset), veldatqc(jveloset), &
               &              ln_vel3d, ln_nea, ld_velav(jveloset) )
            
         ENDIF

         ! Daily averaged data

         IF ( ln_velavadcp ) THEN
            
            jveloset = jveloset + 1
            
            ld_velav(jveloset) = .TRUE.
            
            CALL obs_rea_vel_dri( 1, velodata(jveloset), jnumvelavadcp, &
               &                  velavadcpfiles(1:jnumvelavadcp), &
               &                  nvelovars, nveloextr, &
               &                  nitend-nit000+2,              &
               &                  dobsini, dobsend, ln_ignmis, &
               &                  ld_velav(jveloset), &
               &                  .FALSE. )
            
            DO jvar = 1, 2
               CALL obs_prof_staend( velodata(jveloset), jvar )
            END DO
            
            CALL obs_pre_vel( velodata(jveloset), veldatqc(jveloset), &
               &              ln_vel3d, ln_nea, ld_velav(jveloset) )
            
         ENDIF

         ! High frequency data

         IF ( ln_velhradcp ) THEN
            
            jveloset = jveloset + 1
            
            ld_velav(jveloset) = .FALSE.
               
            CALL obs_rea_vel_dri( 1, velodata(jveloset), jnumvelhradcp, &
               &                  velhradcpfiles(1:jnumvelhradcp), &
               &                  nvelovars, nveloextr, &
               &                  nitend-nit000+2,              &
               &                  dobsini, dobsend, ln_ignmis, &
               &                  ld_velav(jveloset), &
               &                  .FALSE. )
            
            DO jvar = 1, 2
               CALL obs_prof_staend( velodata(jveloset), jvar )
            END DO
            
            CALL obs_pre_vel( velodata(jveloset), veldatqc(jveloset), &
               &              ln_vel3d, ln_nea, ld_velav(jveloset) )
            
         ENDIF

         IF ( ln_velfb ) THEN

            DO jset = 1, jnumvelfb
            
               jveloset = jveloset + 1

               ld_velav(jveloset) = ln_velfb_av(jset)
               
               CALL obs_rea_vel_dri( 0, velodata(jveloset), 1, &
                  &                  velfbfiles(jset:jset), &
                  &                  nvelovars, nveloextr, &
                  &                  nitend-nit000+2,              &
                  &                  dobsini, dobsend, ln_ignmis, &
                  &                  ld_velav(jveloset), &
                  &                  .FALSE. )
               
               DO jvar = 1, 2
                  CALL obs_prof_staend( velodata(jveloset), jvar )
               END DO
               
               CALL obs_pre_vel( velodata(jveloset), veldatqc(jveloset), &
                  &              ln_vel3d, ln_nea, ld_velav(jveloset) )


            END DO
            
         ENDIF

      ENDIF

      !  - log10(chlorophyll)
      
      IF ( ln_logchl ) THEN

         ! Set the number of variables for logchl to 1
         nlogchlvars = 1

         ! Set the number of extra variables for logchl to 0
         nlogchlextr = 0
         
         IF ( ln_logchlfb ) THEN
            nlogchlsets = jnumlogchlfb
         ELSE
            nlogchlsets = 1
         ENDIF

         ALLOCATE(logchldata(nlogchlsets))
         ALLOCATE(logchldatqc(nlogchlsets))
         logchldata(:)%nsurf=0
         logchldatqc(:)%nsurf=0

         nlogchlsets = 0

         IF ( ln_logchlfb ) THEN             ! Feedback file format

            DO jset = 1, jnumlogchlfb
            
               nlogchlsets = nlogchlsets + 1

               CALL obs_rea_logchl( 0, logchldata(nlogchlsets), 1, &
                  &                 logchlfbfiles(jset:jset), &
                  &                 nlogchlvars, nlogchlextr, nitend-nit000+2, &
                  &                 dobsini, dobsend, ln_ignmis, .FALSE. )

               CALL obs_pre_logchl( logchldata(nlogchlsets), logchldatqc(nlogchlsets), &
                  &                 ln_logchl, ln_nea )
            
            ENDDO

         ELSE                              ! Original file format

            nlogchlsets = nlogchlsets + 1

            CALL obs_rea_logchl( 1, logchldata(nlogchlsets), jnumlogchl, &
               &                 logchlfiles(1:jnumlogchl), &
               &                 nlogchlvars, nlogchlextr, nitend-nit000+2, &
               &                 dobsini, dobsend, ln_ignmis, .FALSE. )

            CALL obs_pre_logchl( logchldata(nlogchlsets), logchldatqc(nlogchlsets), &
               &                 ln_logchl, ln_nea )

         ENDIF
 
      ENDIF

      !  - spm
      
      IF ( ln_spm ) THEN

         ! Set the number of variables for spm to 1
         nspmvars = 1

         ! Set the number of extra variables for spm to 0
         nspmextr = 0
         
         IF ( ln_spmfb ) THEN
            nspmsets = jnumspmfb
         ELSE
            nspmsets = 1
         ENDIF

         ALLOCATE(spmdata(nspmsets))
         ALLOCATE(spmdatqc(nspmsets))
         spmdata(:)%nsurf=0
         spmdatqc(:)%nsurf=0

         nspmsets = 0

         IF ( ln_spmfb ) THEN             ! Feedback file format

            DO jset = 1, jnumspmfb
            
               nspmsets = nspmsets + 1

               CALL obs_rea_spm( 0, spmdata(nspmsets), 1, &
                  &                 spmfbfiles(jset:jset), &
                  &                 nspmvars, nspmextr, nitend-nit000+2, &
                  &                 dobsini, dobsend, ln_ignmis, .FALSE. )

               CALL obs_pre_spm( spmdata(nspmsets), spmdatqc(nspmsets), &
                  &                 ln_spm, ln_nea )
            
            ENDDO

         ELSE                              ! Original file format

            nspmsets = nspmsets + 1

            CALL obs_rea_spm( 1, spmdata(nspmsets), jnumspm, &
               &                 spmfiles(1:jnumspm), &
               &                 nspmvars, nspmextr, nitend-nit000+2, &
               &                 dobsini, dobsend, ln_ignmis, .FALSE. )

            CALL obs_pre_spm( spmdata(nspmsets), spmdatqc(nspmsets), &
               &                 ln_spm, ln_nea )

         ENDIF
 
      ENDIF

      !  - fco2
      
      IF ( ln_fco2 ) THEN

         ! Set the number of variables for fco2 to 1
         nfco2vars = 1

         ! Set the number of extra variables for fco2 to 0
         nfco2extr = 0
         
         IF ( ln_fco2fb ) THEN
            nfco2sets = jnumfco2fb
         ELSE
            nfco2sets = 1
         ENDIF

         ALLOCATE(fco2data(nfco2sets))
         ALLOCATE(fco2datqc(nfco2sets))
         fco2data(:)%nsurf=0
         fco2datqc(:)%nsurf=0

         nfco2sets = 0

         IF ( ln_fco2fb ) THEN             ! Feedback file format

            DO jset = 1, jnumfco2fb
            
               nfco2sets = nfco2sets + 1

               CALL obs_rea_fco2( 0, fco2data(nfco2sets), 1, &
                  &                 fco2fbfiles(jset:jset), &
                  &                 nfco2vars, nfco2extr, nitend-nit000+2, &
                  &                 dobsini, dobsend, ln_ignmis, .FALSE. )

               CALL obs_pre_fco2( fco2data(nfco2sets), fco2datqc(nfco2sets), &
                  &                 ln_fco2, ln_nea )
            
            ENDDO

         ELSE                              ! Original file format

            nfco2sets = nfco2sets + 1

            CALL obs_rea_fco2( 1, fco2data(nfco2sets), jnumfco2, &
               &                 fco2files(1:jnumfco2), &
               &                 nfco2vars, nfco2extr, nitend-nit000+2, &
               &                 dobsini, dobsend, ln_ignmis, .FALSE. )

            CALL obs_pre_fco2( fco2data(nfco2sets), fco2datqc(nfco2sets), &
               &                 ln_fco2, ln_nea )

         ENDIF
 
      ENDIF

      !  - pco2
      
      IF ( ln_pco2 ) THEN

         ! Set the number of variables for pco2 to 1
         npco2vars = 1

         ! Set the number of extra variables for pco2 to 0
         npco2extr = 0
         
         IF ( ln_pco2fb ) THEN
            npco2sets = jnumpco2fb
         ELSE
            npco2sets = 1
         ENDIF

         ALLOCATE(pco2data(npco2sets))
         ALLOCATE(pco2datqc(npco2sets))
         pco2data(:)%nsurf=0
         pco2datqc(:)%nsurf=0

         npco2sets = 0

         IF ( ln_pco2fb ) THEN             ! Feedback file format

            DO jset = 1, jnumpco2fb
            
               npco2sets = npco2sets + 1

               CALL obs_rea_pco2( 0, pco2data(npco2sets), 1, &
                  &                 pco2fbfiles(jset:jset), &
                  &                 npco2vars, npco2extr, nitend-nit000+2, &
                  &                 dobsini, dobsend, ln_ignmis, .FALSE. )

               CALL obs_pre_pco2( pco2data(npco2sets), pco2datqc(npco2sets), &
                  &                 ln_pco2, ln_nea )
            
            ENDDO

         ELSE                              ! Original file format

            npco2sets = npco2sets + 1

            CALL obs_rea_pco2( 1, pco2data(npco2sets), jnumpco2, &
               &                 pco2files(1:jnumpco2), &
               &                 npco2vars, npco2extr, nitend-nit000+2, &
               &                 dobsini, dobsend, ln_ignmis, .FALSE. )

            CALL obs_pre_pco2( pco2data(npco2sets), pco2datqc(npco2sets), &
               &                 ln_pco2, ln_nea )

         ENDIF
 
      ENDIF
     
   END SUBROUTINE dia_obs_init

   SUBROUTINE dia_obs( kstp )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE dia_obs  ***
      !!          
      !! ** Purpose : Call the observation operators on each time step
      !!
      !! ** Method  : Call the observation operators on each time step to
      !!              compute the model equivalent of the following date:
      !!               - T profiles
      !!               - S profiles
      !!               - Sea surface height (referenced to a mean)
      !!               - Sea surface temperature
      !!               - Sea surface salinity
      !!               - Velocity component (U,V) profiles
      !!               - Sea surface log10(chlorophyll)
      !!               - Sea surface spm
      !!               - Sea surface fco2
      !!               - Sea surface pco2
      !!
      !! ** Action  : 
      !!
      !! History :
      !!        !  06-03  (K. Mogensen) Original code
      !!        !  06-05  (K. Mogensen) Reformatted
      !!        !  06-10  (A. Weaver) Cleaning
      !!        !  07-03  (K. Mogensen) General handling of profiles
      !!        !  07-04  (G. Smith) Generalized surface operators
      !!        !  08-10  (M. Valdivieso) obs operator for velocity profiles
      !!        !  14-08  (J. While) observation operator for profiles in 
      !!                             generalised vertical coordinates
      !!----------------------------------------------------------------------
      !! * Modules used
      USE dom_oce, ONLY : &             ! Ocean space and time domain variables
         & rdt,           &                       
         & gdept_1d,       &             
#if defined key_vvl 
         & gdept_n,       &
#else 
         & gdept_1d,      &
#endif                                        
         & tmask, umask, vmask                            
      USE phycst, ONLY : &              ! Physical constants
         & rday, &
         & rt0
      USE oce, ONLY : &                 ! Ocean dynamics and tracers variables
         & tsn,  &             
         & un, vn,  &
         & sshn
#if defined  key_lim3
      USE ice, ONLY : &                     ! LIM Ice model variables
         & frld
#endif
#if defined key_lim2
      USE ice_2, ONLY : &                     ! LIM Ice model variables
         & frld
#endif
#if defined key_hadocc
      USE trc, ONLY :  &                ! HadOCC chlorophyll, fCO2 and pCO2
         & HADOCC_CHL, &
         & HADOCC_FCO2, &
         & HADOCC_PCO2, &
         & HADOCC_FILL_FLT
#elif defined key_medusa && defined key_foam_medusa
      USE trc, ONLY :  &                ! MEDUSA chlorophyll, fCO2 and pCO2
         & MEDUSA_CHL, &
         & MEDUSA_FCO2, &
         & MEDUSA_PCO2, &
         & MEDUSA_FILL_FLT
#elif defined key_fabm
      USE fabm
      USE par_fabm
#endif
#if defined key_spm
      USE par_spm, ONLY: &              ! ERSEM/SPM sediments
         & jp_spm
      USE trc, ONLY :  &
         & trn
#endif
      IMPLICIT NONE

      !! * Arguments
      INTEGER, INTENT(IN) :: kstp                         ! Current timestep
      !! * Local declarations
      INTEGER :: idaystp                ! Number of timesteps per day
      INTEGER :: jprofset               ! Profile data set loop variable
      INTEGER :: jslaset                ! SLA data set loop variable
      INTEGER :: jsstset                ! SST data set loop variable
      INTEGER :: jseaiceset             ! sea ice data set loop variable
      INTEGER :: jveloset               ! velocity profile data loop variable
      INTEGER :: jlogchlset             ! logchl data set loop variable
      INTEGER :: jspmset                ! spm data set loop variable
      INTEGER :: jfco2set               ! fco2 data set loop variable
      INTEGER :: jpco2set               ! pco2 data set loop variable
      INTEGER :: jvar                   ! Variable number    
#if ! defined key_lim2 && ! defined key_lim3
      REAL(wp), POINTER, DIMENSION(:,:) :: frld   
#endif
      REAL(wp) :: tiny                  ! small number
      REAL(wp), DIMENSION(jpi,jpj) :: &
         logchl                         ! array for log chlorophyll
      REAL(wp), DIMENSION(jpi,jpj) :: &
         maskchl                        ! array for special chlorophyll mask
      REAL(wp), DIMENSION(jpi,jpj) :: &
         spm                            ! array for spm
      REAL(wp), DIMENSION(jpi,jpj) :: &
         fco2                           ! array for fco2
      REAL(wp), DIMENSION(jpi,jpj) :: &
         maskfco2                       ! array for special fco2 mask
      REAL(wp), DIMENSION(jpi,jpj) :: &
         pco2                           ! array for pco2
      REAL(wp), DIMENSION(jpi,jpj) :: &
         maskpco2                       ! array for special pco2 mask
      INTEGER :: jn                     ! loop index
#if defined key_fabm
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: logchl_3d
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: pco2_3d
#endif
      CHARACTER(LEN=20) :: datestr=" ",timestr=" "
 
#if ! defined key_lim2 && ! defined key_lim3
      CALL wrk_alloc(jpi,jpj,frld) 
#endif

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'dia_obs : Call the observation operators', kstp
         WRITE(numout,*) '~~~~~~~'
      ENDIF

      idaystp = NINT( rday / rdt )

      !-----------------------------------------------------------------------
      ! No LIM => frld == 0.0_wp
      !-----------------------------------------------------------------------
#if ! defined key_lim2 && ! defined key_lim3
      frld(:,:) = 0.0_wp
#endif
      !-----------------------------------------------------------------------
      ! Depending on switches call various observation operators
      !-----------------------------------------------------------------------

      !  - Temperature/salinity profiles
      IF ( ln_t3d .OR. ln_s3d ) THEN
         DO jprofset = 1, nprofsets
            IF( (.NOT. lk_vvl) .AND. (ln_zco .OR. ln_zps) ) THEN 
               IF(lwp) THEN
                  WRITE(numout,*) 'dia_obs : calling obs_pro_opt'
               ENDIF
               IF ( ld_enact(jprofset) ) THEN 
                  CALL obs_pro_opt( prodatqc(jprofset),                     & 
                     &              kstp, jpi, jpj, jpk, nit000, idaystp,   & 
                     &              tsn(:,:,:,jp_tem), tsn(:,:,:,jp_sal),   & 
                     &              gdept_1d, tmask, n1dint, n2dint,        & 
                     &              kdailyavtypes = endailyavtypes ) 
               ELSE 
                  CALL obs_pro_opt( prodatqc(jprofset),                     & 
                     &              kstp, jpi, jpj, jpk, nit000, idaystp,   & 
                     &              tsn(:,:,:,jp_tem), tsn(:,:,:,jp_sal),   & 
                     &              gdept_1d, tmask, n1dint, n2dint               ) 
               ENDIF 
            ELSE
               IF(lwp) THEN
                  WRITE(numout,*) 'dia_obs : calling obs_pro_sco_opt'
               ENDIF
               IF ( ld_enact(jprofset) ) THEN 
                  CALL obs_pro_sco_opt( prodatqc(jprofset),                 & 
                     &              kstp, jpi, jpj, jpk, nit000, idaystp,   & 
                     &              tsn(:,:,:,jp_tem), tsn(:,:,:,jp_sal),   & 
                     &              fsdept(:,:,:), fsdepw(:,:,:),           & 
                     &              tmask, n1dint, n2dint,                  & 
                     &              kdailyavtypes = endailyavtypes ) 
               ELSE 
                  CALL obs_pro_sco_opt( prodatqc(jprofset),                 & 
                     &              kstp, jpi, jpj, jpk, nit000, idaystp,   & 
                     &              tsn(:,:,:,jp_tem), tsn(:,:,:,jp_sal),   & 
                     &              fsdept(:,:,:), fsdepw(:,:,:),           &
                     &              tmask, n1dint, n2dint ) 
               ENDIF 
            ENDIF
         END DO
      ENDIF

      !  - Sea surface anomaly
      IF ( ln_sla ) THEN
         DO jslaset = 1, nslasets
            CALL obs_sla_opt( sladatqc(jslaset),            &
               &              kstp, jpi, jpj, nit000, sshn, &
               &              tmask(:,:,1), n2dint )
         END DO         
      ENDIF

      !  - Sea surface temperature
      IF ( ln_sst ) THEN
         DO jsstset = 1, nsstsets
            CALL obs_sst_opt( sstdatqc(jsstset),                &
               &              kstp, jpi, jpj, nit000, idaystp,  &
               &              tsn(:,:,1,jp_tem), tmask(:,:,1),  &
               &              n2dint, ld_sstnight(jsstset) )
         END DO
      ENDIF

      !  - Sea surface salinity
      IF ( ln_sss ) THEN
         IF(lwp) WRITE(numout,*) ' SSS currently not available'
      ENDIF

#if defined key_lim2 || defined key_lim3
      IF ( ln_seaice ) THEN
         DO jseaiceset = 1, nseaicesets
            CALL obs_seaice_opt( seaicedatqc(jseaiceset),      &
               &              kstp, jpi, jpj, nit000, 1.-frld, &
               &              tmask(:,:,1), n2dint )
         END DO
      ENDIF      
#endif

      !  - Velocity profiles
      IF ( ln_vel3d ) THEN
         DO jveloset = 1, nvelosets
           ! zonal component of velocity
           CALL obs_vel_opt( veldatqc(jveloset), kstp, jpi, jpj, jpk, &
              &              nit000, idaystp, un, vn, gdept_1d, umask, vmask, &
                             n1dint, n2dint, ld_velav(jveloset) )
         END DO
      ENDIF

      IF ( ln_logchl ) THEN

#if defined key_hadocc
         logchl(:,:) = HADOCC_CHL(:,:,1)    ! (not log) chlorophyll from HadOCC
#elif defined key_medusa && defined key_foam_medusa
         logchl(:,:) = MEDUSA_CHL(:,:,1)    ! (not log) chlorophyll from HadOCC
#elif defined key_fabm
         logchl_3d(:,:,:) = fabm_get_bulk_diagnostic_data(model, jp_fabmdia_chltot)
         logchl(:,:) = logchl_3d(:,:,1)
#else
         CALL ctl_stop( ' Trying to run logchl observation operator', &
            &           ' but no biogeochemical model appears to have been defined' )
#endif

         maskchl(:,:) = tmask(:,:,1)         ! create a special mask to exclude certain things

         ! Take the log10 where we can, otherwise exclude
         tiny = 1.0e-20
         WHERE(logchl(:,:) > tiny .AND. logchl(:,:) /= obfillflt )
            logchl(:,:)  = LOG10(logchl(:,:))
         ELSEWHERE
            logchl(:,:)  = obfillflt
            maskchl(:,:) = 0
         END WHERE

         DO jlogchlset = 1, nlogchlsets
             CALL obs_logchl_opt( logchldatqc(jlogchlset),             &
               &                  kstp, jpi, jpj, nit000, logchl(:,:), &
               &                  maskchl(:,:), n2dint )
         END DO         
      ENDIF 

      IF ( ln_spm ) THEN
#if defined key_spm
         spm(:,:) = 0.0
         DO jn = 1, jp_spm
            spm(:,:) = spm(:,:) + trn(:,:,1,jn)   ! sum SPM sizes
         END DO
#else
         CALL ctl_stop( ' Trying to run spm observation operator', &
            &           ' but no spm model appears to have been defined' )
#endif

         DO jspmset = 1, nspmsets
             CALL obs_spm_opt( spmdatqc(jspmset),                &
               &               kstp, jpi, jpj, nit000, spm(:,:), &
               &               tmask(:,:,1), n2dint )
         END DO         
      ENDIF

      IF ( ln_fco2 ) THEN
         maskfco2(:,:) = tmask(:,:,1)         ! create a special mask to exclude certain things
#if defined key_hadocc
         fco2(:,:) = HADOCC_FCO2(:,:)    ! fCO2 from HadOCC
         IF ( ( MINVAL( HADOCC_FCO2 ) == HADOCC_FILL_FLT ).AND.( MAXVAL( HADOCC_FCO2 ) == HADOCC_FILL_FLT ) ) THEN
            fco2(:,:) = obfillflt
            maskfco2(:,:) = 0
            CALL ctl_warn( ' HadOCC fCO2 values masked out for observation operator', &
               &           ' on timestep ' // TRIM(STR(kstp)),                              &
               &           ' as HADOCC_FCO2(:,:) == HADOCC_FILL_FLT' )
         ENDIF
#elif defined key_medusa && defined key_foam_medusa
         fco2(:,:) = MEDUSA_FCO2(:,:)    ! fCO2 from MEDUSA
         IF ( ( MINVAL( MEDUSA_FCO2 ) == MEDUSA_FILL_FLT ).AND.( MAXVAL( MEDUSA_FCO2 ) == MEDUSA_FILL_FLT ) ) THEN
            fco2(:,:) = obfillflt
            maskfco2(:,:) = 0
            CALL ctl_warn( ' MEDUSA fCO2 values masked out for observation operator', &
               &           ' on timestep ' // TRIM(STR(kstp)),                              &
               &           ' as MEDUSA_FCO2(:,:) == MEDUSA_FILL_FLT' )
         ENDIF
#elif defined key_fabm
         ! First, get pCO2 from FABM
         pco2_3d(:,:,:) = fabm_get_bulk_diagnostic_data(model, jp_fabm_o3pc)
         pco2(:,:) = pco2_3d(:,:,1)
         ! Now, convert pCO2 to fCO2, based on SST in K. This follows the standard methodology of:
         ! Pierrot et al. (2009), Recommendations for autonomous underway pCO2 measuring systems
         ! and data reduction routines, Deep-Sea Research II, 56: 512-522.
         ! and
         ! Weiss (1974), Carbon dioxide in water and seawater: the solubility of a non-ideal gas,
         ! Marine Chemistry, 2: 203-215.
         ! In the implementation below, atmospheric pressure has been assumed to be 1 atm and so
         ! not explicitly included - atmospheric pressure is not necessarily available so this is
         ! the best assumption.
         ! Further, the (1-xCO2)^2 term has been neglected. This is common practice
         ! (see e.g. Zeebe and Wolf-Gladrow (2001), CO2 in Seawater: Equilibrium, Kinetics, Isotopes)
         ! because xCO2 in atm is ~0, and so this term will only affect the result to the 3rd decimal
         ! place for typical values, and xCO2 would need to be approximated from pCO2 anyway.
         fco2(:,:) = pco2(:,:) * EXP((-1636.75                                                                               + &
            &                         12.0408      * (tsn(:,:,1,jp_tem)+rt0)                                                 - &
            &                         0.0327957    * (tsn(:,:,1,jp_tem)+rt0)*(tsn(:,:,1,jp_tem)+rt0)                         + &
            &                         0.0000316528 * (tsn(:,:,1,jp_tem)+rt0)*(tsn(:,:,1,jp_tem)+rt0)*(tsn(:,:,1,jp_tem)+rt0) + &
            &                         2.0 * (57.7 - 0.118 * (tsn(:,:,1,jp_tem)+rt0)))                                        / &
            &                        (82.0578 * (tsn(:,:,1,jp_tem)+rt0)))
#else
         CALL ctl_stop( ' Trying to run fco2 observation operator', &
            &           ' but no biogeochemical model appears to have been defined' )
#endif

         DO jfco2set = 1, nfco2sets
             CALL obs_fco2_opt( fco2datqc(jfco2set),                      &
               &                kstp, jpi, jpj, nit000, fco2(:,:), &
               &                maskfco2(:,:), n2dint )
         END DO
      ENDIF

      IF ( ln_pco2 ) THEN
         maskpco2(:,:) = tmask(:,:,1)         ! create a special mask to exclude certain things
#if defined key_hadocc
         pco2(:,:) = HADOCC_PCO2(:,:)    ! pCO2 from HadOCC
         IF ( ( MINVAL( HADOCC_PCO2 ) == HADOCC_FILL_FLT ).AND.( MAXVAL( HADOCC_PCO2 ) == HADOCC_FILL_FLT ) ) THEN
            pco2(:,:) = obfillflt
            maskpco2(:,:) = 0
            CALL ctl_warn( ' HadOCC pCO2 values masked out for observation operator', &
               &           ' on timestep ' // TRIM(STR(kstp)),                              &
               &           ' as HADOCC_PCO2(:,:) == HADOCC_FILL_FLT' )
         ENDIF
#elif defined key_medusa && defined key_foam_medusa
         pco2(:,:) = MEDUSA_PCO2(:,:)    ! pCO2 from MEDUSA
         IF ( ( MINVAL( MEDUSA_PCO2 ) == MEDUSA_FILL_FLT ).AND.( MAXVAL( MEDUSA_PCO2 ) == MEDUSA_FILL_FLT ) ) THEN
            pco2(:,:) = obfillflt
            maskpco2(:,:) = 0
            CALL ctl_warn( ' MEDUSA pCO2 values masked out for observation operator', &
               &           ' on timestep ' // TRIM(STR(kstp)),                              &
               &           ' as MEDUSA_PCO2(:,:) == MEDUSA_FILL_FLT' )
         ENDIF
#elif defined key_fabm
         pco2_3d(:,:,:) = fabm_get_bulk_diagnostic_data(model, jp_fabm_o3pc)
         pco2(:,:) = pco2_3d(:,:,1)
#else
         CALL ctl_stop( ' Trying to run pCO2 observation operator', &
            &           ' but no biogeochemical model appears to have been defined' )
#endif

         DO jpco2set = 1, npco2sets
             CALL obs_pco2_opt( pco2datqc(jpco2set),                      &
               &                kstp, jpi, jpj, nit000, pco2(:,:), &
               &                maskpco2(:,:), n2dint )
         END DO
      ENDIF

#if ! defined key_lim2 && ! defined key_lim3
      CALL wrk_dealloc(jpi,jpj,frld) 
#endif

   END SUBROUTINE dia_obs
  
   SUBROUTINE dia_obs_wri 
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE dia_obs_wri  ***
      !!          
      !! ** Purpose : Call observation diagnostic output routines
      !!
      !! ** Method  : Call observation diagnostic output routines
      !!
      !! ** Action  : 
      !!
      !! History :
      !!        !  06-03  (K. Mogensen) Original code
      !!        !  06-05  (K. Mogensen) Reformatted
      !!        !  06-10  (A. Weaver) Cleaning
      !!        !  07-03  (K. Mogensen) General handling of profiles
      !!        !  08-09  (M. Valdivieso) Velocity component (U,V) profiles
      !!----------------------------------------------------------------------
      IMPLICIT NONE

      !! * Local declarations

      INTEGER :: jprofset                 ! Profile data set loop variable
      INTEGER :: jveloset                 ! Velocity data set loop variable
      INTEGER :: jslaset                  ! SLA data set loop variable
      INTEGER :: jsstset                  ! SST data set loop variable
      INTEGER :: jseaiceset               ! Sea Ice data set loop variable
      INTEGER :: jlogchlset               ! logchl data set loop variable
      INTEGER :: jspmset                  ! spm data set loop variable
      INTEGER :: jfco2set                 ! fco2 data set loop variable
      INTEGER :: jpco2set                 ! pco2 data set loop variable
      INTEGER :: jset
      INTEGER :: jfbini
      CHARACTER(LEN=20) :: datestr=" ",timestr=" "
      CHARACTER(LEN=20) :: cdtmp
      !-----------------------------------------------------------------------
      ! Depending on switches call various observation output routines
      !-----------------------------------------------------------------------

      !  - Temperature/salinity profiles

      IF( ln_t3d .OR. ln_s3d ) THEN

         ! Copy data from prodatqc to profdata structures
         DO jprofset = 1, nprofsets

            CALL obs_prof_decompress( prodatqc(jprofset), &
                 &                    profdata(jprofset), .TRUE., numout )

         END DO

         ! Write the profiles.

         jprofset = 0

         ! ENACT insitu data

         IF ( ln_ena ) THEN
           
            jprofset = jprofset + 1

            CALL obs_wri_p3d( 'enact', profdata(jprofset) )

         ENDIF

         ! Coriolis insitu data

         IF ( ln_cor ) THEN
            
            jprofset = jprofset + 1

            CALL obs_wri_p3d( 'corio', profdata(jprofset) )
            
         ENDIF
         
         ! Feedback insitu data

         IF ( ln_profb ) THEN

            jfbini = jprofset + 1

            DO jprofset = jfbini, nprofsets
               
               jset = jprofset - jfbini + 1
               WRITE(cdtmp,'(A,I2.2)')'profb_',jset
               CALL obs_wri_p3d( cdtmp, profdata(jprofset) )

            END DO

         ENDIF

      ENDIF

      !  - Sea surface anomaly
      IF ( ln_sla ) THEN

         ! Copy data from sladatqc to sladata structures
         DO jslaset = 1, nslasets

              CALL obs_surf_decompress( sladatqc(jslaset), &
                 &                    sladata(jslaset), .TRUE., numout )

         END DO

         jslaset = 0 

         ! Write the AVISO SLA data

         IF ( ln_sladt ) THEN
            
            jslaset = 1
            CALL obs_wri_sla( 'aviso_act', sladata(jslaset) )
            jslaset = 2
            CALL obs_wri_sla( 'aviso_pas', sladata(jslaset) )

         ENDIF

         IF ( ln_slafb ) THEN
            
            jfbini = jslaset + 1

            DO jslaset = jfbini, nslasets
               
               jset = jslaset - jfbini + 1
               WRITE(cdtmp,'(A,I2.2)')'slafb_',jset
               CALL obs_wri_sla( cdtmp, sladata(jslaset) )

            END DO

         ENDIF

      ENDIF

      !  - Sea surface temperature
      IF ( ln_sst ) THEN

         ! Copy data from sstdatqc to sstdata structures
         DO jsstset = 1, nsstsets
     
              CALL obs_surf_decompress( sstdatqc(jsstset), &
                 &                    sstdata(jsstset), .TRUE., numout )

         END DO

         jsstset = 0 

         ! Write the AVISO SST data

         IF ( ln_reysst ) THEN
            
            jsstset = jsstset + 1
            CALL obs_wri_sst( 'reynolds', sstdata(jsstset) )

         ENDIF

         IF ( ln_ghrsst ) THEN
            
            jsstset = jsstset + 1
            CALL obs_wri_sst( 'ghr', sstdata(jsstset) )

         ENDIF

         IF ( ln_sstfb ) THEN
            
            jfbini = jsstset + 1

            DO jsstset = jfbini, nsstsets
               
               jset = jsstset - jfbini + 1
               WRITE(cdtmp,'(A,I2.2)')'sstfb_',jset
               CALL obs_wri_sst( cdtmp, sstdata(jsstset) )

            END DO

         ENDIF

      ENDIF

      !  - Sea surface salinity
      IF ( ln_sss ) THEN
         IF(lwp) WRITE(numout,*) ' SSS currently not available'
      ENDIF

      !  - Sea Ice Concentration
      IF ( ln_seaice ) THEN

         ! Copy data from seaicedatqc to seaicedata structures
         DO jseaiceset = 1, nseaicesets

              CALL obs_surf_decompress( seaicedatqc(jseaiceset), &
                 &                    seaicedata(jseaiceset), .TRUE., numout )

         END DO

         ! Write the Sea Ice data
         DO jseaiceset = 1, nseaicesets
      
            WRITE(cdtmp,'(A,I2.2)')'seaicefb_',jseaiceset
            CALL obs_wri_seaice( cdtmp, seaicedata(jseaiceset) )

         END DO

      ENDIF
      
      ! Velocity data
      IF( ln_vel3d ) THEN

         ! Copy data from veldatqc to velodata structures
         DO jveloset = 1, nvelosets

            CALL obs_prof_decompress( veldatqc(jveloset), &
                 &                    velodata(jveloset), .TRUE., numout )

         END DO

         ! Write the profiles.

         jveloset = 0

         ! Daily averaged data

         IF ( ln_velavcur ) THEN
            
            jveloset = jveloset + 1

            CALL obs_wri_vel( 'velavcurr', velodata(jveloset), n2dint )

         ENDIF

         ! High frequency data

         IF ( ln_velhrcur ) THEN
            
            jveloset = jveloset + 1

            CALL obs_wri_vel( 'velhrcurr', velodata(jveloset), n2dint )

         ENDIF

         ! Daily averaged data

         IF ( ln_velavadcp ) THEN
            
            jveloset = jveloset + 1

            CALL obs_wri_vel( 'velavadcp', velodata(jveloset), n2dint )

         ENDIF

         ! High frequency data

         IF ( ln_velhradcp ) THEN
            
            jveloset = jveloset + 1
            
            CALL obs_wri_vel( 'velhradcp', velodata(jveloset), n2dint )
               
         ENDIF

         ! Feedback velocity data

         IF ( ln_velfb ) THEN

            jfbini = jveloset + 1

            DO jveloset = jfbini, nvelosets
               
               jset = jveloset - jfbini + 1
               WRITE(cdtmp,'(A,I2.2)')'velfb_',jset
               CALL obs_wri_vel( cdtmp, velodata(jveloset), n2dint )

            END DO

         ENDIF
         
      ENDIF

      !  - log10(chlorophyll)
      IF ( ln_logchl ) THEN

         ! Copy data from logchldatqc to logchldata structures
         DO jlogchlset = 1, nlogchlsets

            CALL obs_surf_decompress( logchldatqc(jlogchlset), &
                 &                    logchldata(jlogchlset), .TRUE., numout )

         END DO
         
         ! Mark as bad observations with no valid model counterpart due to activities in dia_obs
         ! Seem to need to set to fill value rather than marking as bad to be effective, so do both
         DO jlogchlset = 1, nlogchlsets
            WHERE ( logchldata(jlogchlset)%rmod(:,1) == obfillflt )
               logchldata(jlogchlset)%nqc(:)    = 1
               logchldata(jlogchlset)%robs(:,1) = obfillflt
            END WHERE
         END DO

         ! Write the logchl data
         DO jlogchlset = 1, nlogchlsets
      
            WRITE(cdtmp,'(A,I2.2)')'logchlfb_',jlogchlset
            CALL obs_wri_logchl( cdtmp, logchldata(jlogchlset) )

         END DO

      ENDIF

      !  - spm
      IF ( ln_spm ) THEN

         ! Copy data from spmdatqc to spmdata structures
         DO jspmset = 1, nspmsets

            CALL obs_surf_decompress( spmdatqc(jspmset), &
                 &                    spmdata(jspmset), .TRUE., numout )

         END DO

         ! Write the spm data
         DO jspmset = 1, nspmsets
      
            WRITE(cdtmp,'(A,I2.2)')'spmfb_',jspmset
            CALL obs_wri_spm( cdtmp, spmdata(jspmset) )

         END DO

      ENDIF

      !  - fco2
      IF ( ln_fco2 ) THEN

         ! Copy data from fco2datqc to fco2data structures
         DO jfco2set = 1, nfco2sets

            CALL obs_surf_decompress( fco2datqc(jfco2set), &
                 &                    fco2data(jfco2set), .TRUE., numout )

         END DO
         
         ! Mark as bad observations with no valid model counterpart due to fCO2 not being in the restart
         ! Seem to need to set to fill value rather than marking as bad to be effective, so do both
         DO jfco2set = 1, nfco2sets
            WHERE ( fco2data(jfco2set)%rmod(:,1) == obfillflt )
               fco2data(jfco2set)%nqc(:)    = 1
               fco2data(jfco2set)%robs(:,1) = obfillflt
            END WHERE
         END DO

         ! Write the fco2 data
         DO jfco2set = 1, nfco2sets
      
            WRITE(cdtmp,'(A,I2.2)')'fco2fb_',jfco2set
            CALL obs_wri_fco2( cdtmp, fco2data(jfco2set) )

         END DO

      ENDIF

      !  - pco2
      IF ( ln_pco2 ) THEN

         ! Copy data from pco2datqc to pco2data structures
         DO jpco2set = 1, npco2sets

            CALL obs_surf_decompress( pco2datqc(jpco2set), &
                 &                    pco2data(jpco2set), .TRUE., numout )

         END DO
         
         ! Mark as bad observations with no valid model counterpart due to pco2 not being in the restart
         ! Seem to need to set to fill value rather than marking as bad to be effective, so do both
         DO jpco2set = 1, npco2sets
            WHERE ( pco2data(jpco2set)%rmod(:,1) == obfillflt )
               pco2data(jpco2set)%nqc(:)    = 1
               pco2data(jpco2set)%robs(:,1) = obfillflt
            END WHERE
         END DO

         ! Write the pco2 data
         DO jpco2set = 1, npco2sets
      
            WRITE(cdtmp,'(A,I2.2)')'pco2fb_',jpco2set
            CALL obs_wri_pco2( cdtmp, pco2data(jpco2set) )

         END DO

      ENDIF

   END SUBROUTINE dia_obs_wri

   SUBROUTINE dia_obs_dealloc
      IMPLICIT NONE
      !!----------------------------------------------------------------------
      !!                    *** ROUTINE dia_obs_dealloc ***
      !!
      !!  ** Purpose : To deallocate data to enable the obs_oper online loop.
      !!               Specifically: dia_obs_init --> dia_obs --> dia_obs_wri
      !!
      !!  ** Method : Clean up various arrays left behind by the obs_oper.
      !!
      !!  ** Action :
      !!
      !!----------------------------------------------------------------------
      !! obs_grid deallocation
      CALL obs_grid_deallocate

      !! diaobs deallocation
      IF ( nprofsets > 0 ) THEN
          DEALLOCATE(ld_enact, &
                  &  profdata, &
                  &  prodatqc)
      END IF
      IF ( ln_sla ) THEN
          DEALLOCATE(sladata, &
                  &  sladatqc)
      END IF
      IF ( ln_seaice ) THEN
          DEALLOCATE(sladata, &
                  &  sladatqc)
      END IF
      IF ( ln_sst ) THEN
          DEALLOCATE(sstdata, &
                  &  sstdatqc)
      END IF
      IF ( ln_vel3d ) THEN
          DEALLOCATE(ld_velav, &
                  &  velodata, &
                  &  veldatqc)
      END IF
   END SUBROUTINE dia_obs_dealloc

   SUBROUTINE ini_date( ddobsini )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE ini_date  ***
      !!          
      !! ** Purpose : Get initial data in double precision YYYYMMDD.HHMMSS format
      !!
      !! ** Method  : Get initial data in double precision YYYYMMDD.HHMMSS format
      !!
      !! ** Action  : Get initial data in double precision YYYYMMDD.HHMMSS format
      !!
      !! History :
      !!        !  06-03  (K. Mogensen)  Original code
      !!        !  06-05  (K. Mogensen)  Reformatted
      !!        !  06-10  (A. Weaver) Cleaning
      !!        !  06-10  (G. Smith) Calculates initial date the same as method for final date
      !!        !  10-05  (D. Lea) Update to month length calculation for NEMO vn3.2
      !!----------------------------------------------------------------------
      USE phycst, ONLY : &            ! Physical constants
         & rday
!      USE daymod, ONLY : &            ! Time variables
!         & nmonth_len           
      USE dom_oce, ONLY : &           ! Ocean space and time domain variables
         & rdt

      IMPLICIT NONE

      !! * Arguments
      REAL(KIND=dp), INTENT(OUT) :: ddobsini                         ! Initial date in YYYYMMDD.HHMMSS

      !! * Local declarations
      INTEGER :: iyea        ! date - (year, month, day, hour, minute)
      INTEGER :: imon
      INTEGER :: iday
      INTEGER :: ihou
      INTEGER :: imin
      INTEGER :: imday         ! Number of days in month.
      REAL(KIND=wp) :: zdayfrc ! Fraction of day

      INTEGER, DIMENSION(12) ::   imonth_len    !: length in days of the months of the current year

      !!----------------------------------------------------------------------
      !! Initial date initialization (year, month, day, hour, minute)
      !! (This assumes that the initial date is for 00z))
      !!----------------------------------------------------------------------
      iyea =   ndate0 / 10000
      imon = ( ndate0 - iyea * 10000 ) / 100
      iday =   ndate0 - iyea * 10000 - imon * 100
      ihou = 0
      imin = 0

      !!----------------------------------------------------------------------
      !! Compute number of days + number of hours + min since initial time
      !!----------------------------------------------------------------------
      iday = iday + ( nit000 -1 ) * rdt / rday
      zdayfrc = ( nit000 -1 ) * rdt / rday
      zdayfrc = zdayfrc - aint(zdayfrc)
      ihou = int( zdayfrc * 24 )
      imin = int( (zdayfrc * 24 - ihou) * 60 )

      !!-----------------------------------------------------------------------
      !! Convert number of days (iday) into a real date
      !!----------------------------------------------------------------------

      CALL calc_month_len( iyea, imonth_len )
      
      DO WHILE ( iday > imonth_len(imon) )
         iday = iday - imonth_len(imon)
         imon = imon + 1 
         IF ( imon > 12 ) THEN
            imon = 1
            iyea = iyea + 1
            CALL calc_month_len( iyea, imonth_len )  ! update month lengths
         ENDIF
      END DO

      !!----------------------------------------------------------------------
      !! Convert it into YYYYMMDD.HHMMSS format.
      !!----------------------------------------------------------------------
      ddobsini = iyea * 10000_dp + imon * 100_dp + &
         &       iday + ihou * 0.01_dp + imin * 0.0001_dp


   END SUBROUTINE ini_date

   SUBROUTINE fin_date( ddobsfin )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE fin_date  ***
      !!          
      !! ** Purpose : Get final data in double precision YYYYMMDD.HHMMSS format
      !!
      !! ** Method  : Get final data in double precision YYYYMMDD.HHMMSS format
      !!
      !! ** Action  : Get final data in double precision YYYYMMDD.HHMMSS format
      !!
      !! History :
      !!        !  06-03  (K. Mogensen)  Original code
      !!        !  06-05  (K. Mogensen)  Reformatted
      !!        !  06-10  (A. Weaver) Cleaning
      !!        !  10-05  (D. Lea) Update to month length calculation for NEMO vn3.2
      !!----------------------------------------------------------------------
      USE phycst, ONLY : &            ! Physical constants
         & rday
!      USE daymod, ONLY : &            ! Time variables
!         & nmonth_len                
      USE dom_oce, ONLY : &           ! Ocean space and time domain variables
         & rdt

      IMPLICIT NONE

      !! * Arguments
      REAL(KIND=dp), INTENT(OUT) :: ddobsfin                   ! Final date in YYYYMMDD.HHMMSS

      !! * Local declarations
      INTEGER :: iyea        ! date - (year, month, day, hour, minute)
      INTEGER :: imon
      INTEGER :: iday
      INTEGER :: ihou
      INTEGER :: imin
      INTEGER :: imday         ! Number of days in month.
      REAL(KIND=wp) :: zdayfrc       ! Fraction of day
         
      INTEGER, DIMENSION(12) ::   imonth_len    !: length in days of the months of the current year
            
      !-----------------------------------------------------------------------
      ! Initial date initialization (year, month, day, hour, minute)
      ! (This assumes that the initial date is for 00z)
      !-----------------------------------------------------------------------
      iyea =   ndate0 / 10000
      imon = ( ndate0 - iyea * 10000 ) / 100
      iday =   ndate0 - iyea * 10000 - imon * 100
      ihou = 0
      imin = 0
      
      !-----------------------------------------------------------------------
      ! Compute number of days + number of hours + min since initial time
      !-----------------------------------------------------------------------
      iday    = iday +  nitend  * rdt / rday
      zdayfrc =  nitend  * rdt / rday
      zdayfrc = zdayfrc - AINT( zdayfrc )
      ihou    = INT( zdayfrc * 24 )
      imin    = INT( ( zdayfrc * 24 - ihou ) * 60 )

      !-----------------------------------------------------------------------
      ! Convert number of days (iday) into a real date
      !----------------------------------------------------------------------

      CALL calc_month_len( iyea, imonth_len )
      
      DO WHILE ( iday > imonth_len(imon) )
         iday = iday - imonth_len(imon)
         imon = imon + 1 
         IF ( imon > 12 ) THEN
            imon = 1
            iyea = iyea + 1
            CALL calc_month_len( iyea, imonth_len )  ! update month lengths
         ENDIF
      END DO

      !-----------------------------------------------------------------------
      ! Convert it into YYYYMMDD.HHMMSS format
      !-----------------------------------------------------------------------
      ddobsfin = iyea * 10000_dp + imon * 100_dp    + iday &
         &     + ihou * 0.01_dp  + imin * 0.0001_dp

    END SUBROUTINE fin_date
    
END MODULE diaobs
