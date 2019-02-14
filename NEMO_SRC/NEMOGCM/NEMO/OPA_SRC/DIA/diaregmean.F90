MODULE diaregmean 
   !!======================================================================
   !!                       ***  MODULE  diaharm  ***
   !! Timeseries of Regional Means 
   !!======================================================================
   !! History :  3.6  !  11/2016  (J Tinker)  Original code
   !!----------------------------------------------------------------------
#if defined key_diareg
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain
   USE in_out_manager  ! I/O units
   USE iom             ! I/0 library
   USE wrk_nemo        ! working arrays
   USE diatmb          ! Top,middle,bottom output
   USE diapea          ! Top,middle,bottom output
   USE zdfmxl          ! MLD
   USE sbc_oce
   USE sbcapr
   USE sbcflx
#if defined key_diaar5
   USE diaar5
#endif

   IMPLICIT NONE
   PRIVATE

   LOGICAL , PUBLIC ::   ln_diaregmean  ! region mean calculation
   LOGICAL, PUBLIC, PARAMETER ::   lk_diareg = .TRUE.   !: model-data diagnostics flag
   PUBLIC   dia_regmean_init            ! routine called by nemogcm.F90
   PUBLIC   dia_regmean                 ! routine called by diawri.F90
   
   LOGICAL :: ln_diaregmean_ascii  ! region mean calculation ascii output
   LOGICAL :: ln_diaregmean_bin  ! region mean calculation binary output
   LOGICAL :: ln_diaregmean_nc  ! region mean calculation netcdf output
   LOGICAL :: ln_diaregmean_diaar5   ! region mean calculation including AR5 SLR terms
   LOGICAL :: ln_diaregmean_diasbc   ! region mean calculation including Surface BC
   LOGICAL :: ln_diaregmean_karamld  ! region mean calculation including kara mld terms
   LOGICAL :: ln_diaregmean_pea  ! region mean calculation including pea terms
   
   REAL(wp), SAVE, ALLOCATABLE,   DIMENSION(:,:,:)  ::   tmp_region_mask_real   ! tempory region_mask of reals
   INTEGER, SAVE, ALLOCATABLE,   DIMENSION(:,:,:)   ::   region_mask            ! region_mask matrix
   INTEGER                                          ::   nmasks                 ! Number of mask files in region_mask.nc file - 
   INTEGER, SAVE, ALLOCATABLE,   DIMENSION(:)       ::   nreg_mat               ! Number of regions in each mask
   
   REAL(wp),  ALLOCATABLE,   DIMENSION(:,:,:) ::   tmp_field_mat !: temporary region_mask
   REAL(wp),  ALLOCATABLE,   DIMENSION(:,:,:) ::   tmp_field_AR5_mat !: temporary region_mask
   REAL(wp),  ALLOCATABLE,   DIMENSION(:,:,:) ::   tmp_field_SBC_mat !: temporary region_mask
   INTEGER  ::   tmp_field_cnt                                   ! tmp_field_cnt integer
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.6 , NEMO Consortium (2014)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dia_regmean_init 
      !!---------------------------------------------------------------------------
      !!                  ***  ROUTINE dia_regmean_init  ***
      !!     
      !! ** Purpose: Initialization of region mask namelist 
      !!        
      !! ** Method : Read namelist
      !!   History
      !!   3.6  !  11-16  (J Tinker) Routine to initialize dia_regmean
      !!---------------------------------------------------------------------------
      !!
      INTEGER  ::   ios                  ! Local integer output status for namelist read
      INTEGER  ::   inum                 ! temporary logical unit ! copied from DOM/domzgr.F90
      INTEGER  ::   ierr                ! error integer for IOM_get
      INTEGER  ::   idmaskvar           ! output of iom_varid
      INTEGER  ::   maskno              ! counter for number of masks
      INTEGER  ::   jj,ji               ! i and j index
      INTEGER  ::   tmpint              ! temporary integer
      REAL(wp),  ALLOCATABLE,   DIMENSION(:,:) ::   tmpregion !: temporary region_mask
      INTEGER, DIMENSION(3) ::   zdimsz   ! number of elements in each of the 3 dimensions (i.e., lon, lat, no of masks, 297,  375,  4) for an array
      INTEGER               ::   zndims   ! number of dimensions in an array (i.e. 3, )
      !
      NAMELIST/nam_diaregmean/ ln_diaregmean,ln_diaregmean_ascii,ln_diaregmean_bin,ln_diaregmean_nc,&
        & ln_diaregmean_karamld, ln_diaregmean_pea,ln_diaregmean_diaar5,ln_diaregmean_diasbc
      
      
      ! read in Namelist. 
      !!----------------------------------------------------------------------
      !
      REWIND ( numnam_ref )              ! Read Namelist nam_diatmb in referdiatmbence namelist : TMB diagnostics
      READ   ( numnam_ref, nam_diaregmean, IOSTAT=ios, ERR= 901 )
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nam_diaregmean in reference namelist', lwp )

      REWIND( numnam_cfg )              ! Namelist nam_diatmb in configuration namelist  TMB diagnostics
      READ  ( numnam_cfg, nam_diaregmean, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nam_diaregmean in configuration namelist', lwp )
      IF(lwm) WRITE ( numond, nam_diaregmean )

      IF(lwp) THEN                   ! Control print
          WRITE(numout,*)
          WRITE(numout,*) 'dia_regmean_init : Output regional mean Diagnostics'
          WRITE(numout,*) '~~~~~~~~~~~~'
          WRITE(numout,*) 'Namelist nam_regmean : set regmeanoutputs '
          WRITE(numout,*) 'Switch for regmean diagnostics (T) or not (F)  ln_diaregmean  = ', ln_diaregmean
          WRITE(numout,*) 'Switch for regmean ascii output (T) or not (F)  ln_diaregmean_ascii  = ', ln_diaregmean_ascii
          WRITE(numout,*) 'Switch for regmean binary output (T) or not (F)  ln_diaregmean_bin  = ', ln_diaregmean_bin
          WRITE(numout,*) 'Switch for regmean netcdf output (T) or not (F)  ln_diaregmean_nc  = ', ln_diaregmean_nc
          WRITE(numout,*) 'Switch for regmean kara mld terms (T) or not (F)  ln_diaregmean_karamld  = ', ln_diaregmean_karamld
          WRITE(numout,*) 'Switch for regmean PEA terms (T) or not (F)  ln_diaregmean_pea  = ', ln_diaregmean_pea
          WRITE(numout,*) 'Switch for regmean AR5 SLR terms (T) or not (F)  ln_diaregmean_diaar5  = ', ln_diaregmean_diaar5
          WRITE(numout,*) 'Switch for regmean Surface forcing terms (T) or not (F)  ln_diaregmean_diasbc  = ', ln_diaregmean_diasbc
      ENDIF
      
      
      !ALLOCATE( tmp_field_mat(jpi,jpj,7),  STAT= ierr ) !SS/NB/DT T/S, SSH, MLD, PEA, PEAT, PEAS
      ALLOCATE( tmp_field_mat(jpi,jpj,11),  STAT= ierr ) !SS/NB/DT T/S, SSH, MLD, PEA, PEAT, PEAS
          IF( ierr /= 0 )   CALL ctl_stop( 'tmp_field_mat: failed to allocate tmp_region_mask_real array' )
      tmp_field_mat(:,:,:) = 0.
      tmp_field_cnt = 0
      
      IF(ln_diaregmean_diaar5) THEN   
        ALLOCATE( tmp_field_AR5_mat(jpi,jpj,4),  STAT= ierr ) !SS/NB/DT T/S, SSH, MLD, PEA, PEAT, PEAS
            IF( ierr /= 0 )   CALL ctl_stop( 'tmp_field_AR5_mat: failed to allocate tmp_region_mask_real array' )
        tmp_field_AR5_mat(:,:,:) = 0.
      ENDIF
      
      IF(ln_diaregmean_diasbc) THEN   
        ALLOCATE( tmp_field_SBC_mat(jpi,jpj,7),  STAT= ierr ) !SS/NB/DT T/S, SSH, MLD, PEA, PEAT, PEAS
            IF( ierr /= 0 )   CALL ctl_stop( 'tmp_field_SBC_mat: failed to allocate tmp_region_mask_real array' )
        tmp_field_SBC_mat(:,:,:) = 0.
      ENDIF
      
      IF (ln_diaregmean) THEN
      
          ! Open region mask for region means, and retrieve the size of the mask (number of levels)          
          CALL iom_open ( 'region_mask.nc', inum )
          idmaskvar = iom_varid( inum, 'mask', kdimsz=zdimsz, kndims=zndims, ldstop = .FALSE.)          
          nmasks = zdimsz(3)
          
          ! read in the region mask (which contains floating point numbers) into a temporary array of reals.
          ALLOCATE( tmp_region_mask_real(jpi,jpj,nmasks),  STAT= ierr )
          IF( ierr /= 0 )   CALL ctl_stop( 'dia_regmean_init: failed to allocate tmp_region_mask_real array' )
          
          ! Use jpdom_unknown to read in a n layer mask.
          tmp_region_mask_real(:,:,:) = 0
          CALL iom_get( inum, jpdom_unknown, 'mask', tmp_region_mask_real(1:nlci,1:nlcj,1:nmasks),   &
              &          kstart = (/ mig(1),mjg(1),1 /), kcount = (/ nlci,nlcj,nmasks /) )
          
          CALL iom_close( inum )
          
          !Convert the region mask of reals into one of integers. 
          
          ALLOCATE( region_mask(jpi,jpj,nmasks),  STAT= ierr )
          IF( ierr /= 0 )   CALL ctl_stop( 'dia_regmean_init: failed to allocate region_mask array' )
          region_mask(:,:,:) = 0
          region_mask = int(tmp_region_mask_real(:,:,:))
          DEALLOCATE( tmp_region_mask_real)
          
          
          ALLOCATE( nreg_mat(nmasks),  STAT= ierr )
          IF( ierr /= 0 )   CALL ctl_stop( 'dia_regmean_init: failed to allocate nreg_mat array' )

          ! work out the number of regions in each mask, asssuming land is 0, and the regions are consectively numbered, 
          ! without missing any number, so the number of regions is the maximum number + 1 (for land). mpp_max across the 
          ! processors to get the global maxima
          DO maskno = 1,nmasks
              tmpint = maxval(region_mask(:,:,maskno))
              CALL mpp_max( tmpint )
              nreg_mat(maskno) = tmpint + 1
          END DO
        
          IF(lwp) THEN 
              ! if writing out as binary and text, open the files. 
              IF ( ln_diaregmean_bin ) THEN
                  ! Open binary for region means
                  !OPEN( UNIT=73, FILE='region_mean_timeseries.dat', FORM='UNFORMATTED', STATUS='REPLACE' )
                  
                  
                  CALL ctl_opn( numdct_reg_bin  ,'region_mean_timeseries.dat'  , 'NEW', 'UNFORMATTED', 'SEQUENTIAL', -1, numout,  .TRUE. )
                  
                  
              ENDIF
              
              IF ( ln_diaregmean_ascii ) THEN
                  ! Open text files for region means
                  !OPEN( UNIT=37, FILE='region_mean_timeseries.txt', FORM='FORMATTED', STATUS='REPLACE' )
                  CALL ctl_opn( numdct_reg_txt  ,'region_mean_timeseries.txt'  , 'NEW', 'FORMATTED', 'SEQUENTIAL', -1, numout,  .TRUE. )
              ENDIF
          ENDIF
     ENDIF

   END SUBROUTINE dia_regmean_init

   SUBROUTINE dia_calctmb_region_mean( pinfield,pouttmb )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE dia_calctmb_region_mean  ***
      !!                   
      !! ** Purpose :    Find the Top, Mid, Bottom and Top minus Bottom fields of water Column
      !!
      !! ** Method  :   
      !!      use mbathy to find surface, mid and bottom of model levels
      !!
      !! History :
      !!   3.6  !  08-14  (E. O'Dea) Routine based on dia_wri_foam
      !!----------------------------------------------------------------------
      !! * Modules used

      ! Routine to map 3d field to top, middle, bottom
      IMPLICIT NONE


      ! Routine arguments
      REAL(wp), DIMENSION(jpi, jpj, jpk), INTENT(IN   ) :: pinfield    ! Input 3d field and mask
      !REAL(wp), DIMENSION(jpi, jpj, 4  ), INTENT(  OUT) :: pouttmb     ! Output top, middle, bottom and surface minus bed
      REAL(wp), DIMENSION(jpi, jpj, 3  ), INTENT(  OUT) :: pouttmb     ! Output top, middle, bottom and surface minus bed

      ! Local variables
      INTEGER :: ji,jj,jk  ! Dummy loop indices

      ! Local Real
      REAL(wp)                         ::   zmdi  !  set masked values

      zmdi=1.e+20 !missing data indicator for masking

      ! Calculate top
      pouttmb(:,:,1) = pinfield(:,:,1)*tmask(:,:,1)  + zmdi*(1.0-tmask(:,:,1))

     ! Calculate middle
      !DO jj = 1,jpj
      !    DO ji = 1,jpi
      !        jk              = max(1,mbathy(ji,jj)/2)
      !        pouttmb(ji,jj,2) = pinfield(ji,jj,jk)*tmask(ji,jj,jk)  + zmdi*(1.0-tmask(ji,jj,jk))
      !    END DO
      !END DO

      ! Calculate bottom, and top minus bottom
      DO jj = 1,jpj
          DO ji = 1,jpi
              jk              = max(1,mbathy(ji,jj))
              !pouttmb(ji,jj,3) = pinfield(ji,jj,jk)*tmask(ji,jj,jk)  + zmdi*(1.0-tmask(ji,jj,jk))
              !pouttmb(ji,jj,4) = (pouttmb(ji,jj,1) - pouttmb(ji,jj,3))*tmask(ji,jj,1)  + zmdi*(1.0-tmask(ji,jj,1))
              pouttmb(ji,jj,2) = pinfield(ji,jj,jk)*tmask(ji,jj,jk)  + zmdi*(1.0-tmask(ji,jj,jk))
              pouttmb(ji,jj,3) = (pouttmb(ji,jj,1) - pouttmb(ji,jj,2))*tmask(ji,jj,1)  + zmdi*(1.0-tmask(ji,jj,1))
          END DO
      END DO

   END SUBROUTINE dia_calctmb_region_mean


   SUBROUTINE dia_regmean( kt ) 
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE dia_regmean  ***
      !! ** Purpose :   Produce regional mean diagnostics
      !!
      !! ** Method  :   calls dia_wri_region_mean to calculate and write the regional means for a number of variables, 
      !!                (calling dia_calctmb_region_mean where necessary).
      !!                
      !!                Closes all text and binary files on last time step
      !!                
      !!      
      !!      
      !!
      !! History :
      !!   3.6  !  11-16  (J. Tinker) 
      !!         
      !!--------------------------------------------------------------------
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zwtmbT    ! temporary T workspace 
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zwtmbS    ! temporary S workspace 
      REAL(wp)                         ::   zmdi      ! set masked values
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
      
      REAL(wp)                         ::   zdt  ! temporary reals
      INTEGER                          ::   i_steps, ierr         ! no of timesteps per hour, allocation error index
      INTEGER                          ::   maskno,jj,ji,jm,nreg ! indices of mask, i and j, and number of regions
      
      zmdi=1.e+20 !missing data indicator for maskin

      IF (ln_diaregmean) THEN
        ! If regional mean calculations required by namelist
        ! -----------------
        ! identify hourly time steps (not used)
        zdt = rdt
        IF( nacc == 1 ) zdt = rdtmin

        IF( MOD( 3600,INT(zdt) ) == 0 ) THEN
            i_steps = 3600/INT(zdt)
        ELSE
            CALL ctl_stop('STOP', 'dia_regmean: timestep must give MOD(3600,rdt) = 0 otherwise no hourly values are possible')
        ENDIF
        
        !i_steps = 1
        
        !Extract 2d fields from 3d T and S with dia_calctmb_region_mean
        CALL wrk_alloc( jpi , jpj, 3 , zwtmbT )
        CALL wrk_alloc( jpi , jpj, 3 , zwtmbS )
            
        CALL dia_calctmb_region_mean(  tsn(:,:,:,jp_tem),zwtmbT)
        CALL dia_calctmb_region_mean(  tsn(:,:,:,jp_sal),zwtmbS)
            
          tmp_field_mat(:,:,1) = tmp_field_mat(:,:,1) + (zwtmbT(:,:,1)*tmask(:,:,1))
          tmp_field_mat(:,:,2) = tmp_field_mat(:,:,2) + (zwtmbT(:,:,2)*tmask(:,:,1))
          tmp_field_mat(:,:,3) = tmp_field_mat(:,:,3) + (zwtmbT(:,:,3)*tmask(:,:,1))
          tmp_field_mat(:,:,4) = tmp_field_mat(:,:,4) + (zwtmbS(:,:,1)*tmask(:,:,1))
          tmp_field_mat(:,:,5) = tmp_field_mat(:,:,5) + (zwtmbS(:,:,2)*tmask(:,:,1))
          tmp_field_mat(:,:,6) = tmp_field_mat(:,:,6) + (zwtmbS(:,:,3)*tmask(:,:,1))
          tmp_field_mat(:,:,7) = tmp_field_mat(:,:,7) + (sshn(:,:)*tmask(:,:,1))
        
        
        IF( ln_diaregmean_karamld  ) THEN
          tmp_field_mat(:,:,8) = tmp_field_mat(:,:,8) + (hmld_kara(:,:)*tmask(:,:,1)) !hmlp(:,:)
        ENDIF
        IF( ln_diaregmean_pea  ) THEN
          tmp_field_mat(:,:,9) = tmp_field_mat(:,:,9) + (pea(:,:)*tmask(:,:,1))
          tmp_field_mat(:,:,10) = tmp_field_mat(:,:,10) + (peat(:,:)*tmask(:,:,1))
          tmp_field_mat(:,:,11) = tmp_field_mat(:,:,11) + (peas(:,:)*tmask(:,:,1))
        ENDIF

#if defined key_diaar5
        IF( ln_diaregmean_diaar5  ) THEN
            tmp_field_AR5_mat(:,:,1) = tmp_field_AR5_mat(:,:,1) + (sshsteric_mat(:,:)*tmask(:,:,1))
            tmp_field_AR5_mat(:,:,2) = tmp_field_AR5_mat(:,:,2) + (sshthster_mat(:,:)*tmask(:,:,1))
            tmp_field_AR5_mat(:,:,3) = tmp_field_AR5_mat(:,:,3) + (sshhlster_mat(:,:)*tmask(:,:,1))
            tmp_field_AR5_mat(:,:,4) = tmp_field_AR5_mat(:,:,4) + (zbotpres_mat(:,:)*tmask(:,:,1))

        ENDIF
#endif
        
        IF( ln_diaregmean_diasbc  ) THEN
        
            tmp_field_SBC_mat(:,:,1) = tmp_field_SBC_mat(:,:,1) + ((qsr  + qns)*tmask(:,:,1))
            tmp_field_SBC_mat(:,:,2) = tmp_field_SBC_mat(:,:,2) + (qsr*tmask(:,:,1))
            tmp_field_SBC_mat(:,:,3) = tmp_field_SBC_mat(:,:,3) + (qns*tmask(:,:,1))
            tmp_field_SBC_mat(:,:,4) = tmp_field_SBC_mat(:,:,4) + (emp*tmask(:,:,1))
            tmp_field_SBC_mat(:,:,5) = tmp_field_SBC_mat(:,:,5) + (wndm*tmask(:,:,1))
            IF ( ln_apr_dyn  .AND. .NOT. ln_shelf_flx ) THEN
               WRITE(numout,*) 'PETE - using apr for pressure'
               tmp_field_SBC_mat(:,:,6) = tmp_field_SBC_mat(:,:,6) + (apr(:,:)*tmask(:,:,1))
            ELSE
               tmp_field_SBC_mat(:,:,6) = tmp_field_SBC_mat(:,:,6) + (pressnow*tmask(:,:,1))
            ENDIF
            tmp_field_SBC_mat(:,:,7) = tmp_field_SBC_mat(:,:,7) + (rnf*tmask(:,:,1))


        ENDIF
        
        tmp_field_cnt = tmp_field_cnt + 1
        
        IF( MOD( kt, i_steps ) == 0  .and. kt .ne. nn_it000 ) THEN

          
          CALL dia_wri_region_mean(kt, "sst" , tmp_field_mat(:,:,1)/real(tmp_field_cnt,wp))
          CALL dia_wri_region_mean(kt, "nbt" , tmp_field_mat(:,:,2)/real(tmp_field_cnt,wp))
          CALL dia_wri_region_mean(kt, "dft" , tmp_field_mat(:,:,3)/real(tmp_field_cnt,wp))

          CALL dia_wri_region_mean(kt, "sss" , tmp_field_mat(:,:,4)/real(tmp_field_cnt,wp))
          CALL dia_wri_region_mean(kt, "nbs" , tmp_field_mat(:,:,5)/real(tmp_field_cnt,wp))
          CALL dia_wri_region_mean(kt, "dfs" , tmp_field_mat(:,:,6)/real(tmp_field_cnt,wp))

          CALL dia_wri_region_mean(kt, "ssh" , tmp_field_mat(:,:,7)/real(tmp_field_cnt,wp))
          
          
        IF( ln_diaregmean_karamld  ) THEN
          
          CALL dia_wri_region_mean(kt, "mldkara" , tmp_field_mat(:,:,8)/real(tmp_field_cnt,wp)) ! tm
        ENDIF
        IF( ln_diaregmean_pea  ) THEN
          
          CALL dia_wri_region_mean(kt, "pea"  , tmp_field_mat(:,:,9)/real(tmp_field_cnt,wp))
          CALL dia_wri_region_mean(kt, "peat" , tmp_field_mat(:,:,10)/real(tmp_field_cnt,wp))
          CALL dia_wri_region_mean(kt, "peas" , tmp_field_mat(:,:,11)/real(tmp_field_cnt,wp)) ! tmb
        ENDIF
          
          tmp_field_mat(:,:,:) = 0.

#if defined key_diaar5        
          IF( ln_diaregmean_diaar5  ) THEN
          
            CALL dia_wri_region_mean(kt, "ssh_steric" ,      tmp_field_AR5_mat(:,:,1)/real(tmp_field_cnt,wp))
            CALL dia_wri_region_mean(kt, "ssh_thermosteric", tmp_field_AR5_mat(:,:,2)/real(tmp_field_cnt,wp))
            CALL dia_wri_region_mean(kt, "ssh_halosteric" ,  tmp_field_AR5_mat(:,:,3)/real(tmp_field_cnt,wp))
            CALL dia_wri_region_mean(kt, "bot_pres" ,        tmp_field_AR5_mat(:,:,4)/real(tmp_field_cnt,wp))
            tmp_field_AR5_mat(:,:,:) = 0.
          ENDIF
#endif
          
          IF( ln_diaregmean_diasbc  ) THEN
          
            CALL dia_wri_region_mean(kt, "qt"   , tmp_field_SBC_mat(:,:,1)/real(tmp_field_cnt,wp))
            CALL dia_wri_region_mean(kt, "qsr"  , tmp_field_SBC_mat(:,:,2)/real(tmp_field_cnt,wp))
            CALL dia_wri_region_mean(kt, "qns"  , tmp_field_SBC_mat(:,:,3)/real(tmp_field_cnt,wp))
            CALL dia_wri_region_mean(kt, "emp"  , tmp_field_SBC_mat(:,:,4)/real(tmp_field_cnt,wp))
            CALL dia_wri_region_mean(kt, "wspd" , tmp_field_SBC_mat(:,:,5)/real(tmp_field_cnt,wp))
            CALL dia_wri_region_mean(kt, "mslp" , tmp_field_SBC_mat(:,:,6)/real(tmp_field_cnt,wp))
            CALL dia_wri_region_mean(kt, "rnf"  , tmp_field_SBC_mat(:,:,7)/real(tmp_field_cnt,wp))
            tmp_field_SBC_mat(:,:,:) = 0.
          ENDIF
          
          tmp_field_cnt = 0
  
        ENDIF
        
        
        ! If on the last time step, close binary and ascii files. 
        IF( kt == nitend ) THEN
          IF(lwp) THEN
            IF ( ln_diaregmean_bin ) THEN
              !Closing binary files for regional mean time series.
              CLOSE(numdct_reg_bin)
            ENDIF
            IF ( ln_diaregmean_ascii ) THEN
              !Closing text files for regional mean time series.
              CLOSE(numdct_reg_txt)
            ENDIF
            
            DEALLOCATE( region_mask, nreg_mat, tmp_field_mat)
            IF( ln_diaregmean_diaar5  ) DEALLOCATE( tmp_field_AR5_mat)
            IF( ln_diaregmean_diasbc  ) DEALLOCATE( tmp_field_SBC_mat)
          ENDIF
        ENDIF
          
          
      ELSE
        CALL ctl_warn('dia_regmean: regmean diagnostic is set to false you should not have seen this')
      ENDIF
      
   END SUBROUTINE dia_regmean
   
   
   SUBROUTINE dia_wri_region_mean(kt, name,         infield  )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE dia_tmb  ***
      !!                   
      !! ** Purpose :   Calculate and write region mean time series for 2d arrays
      !!
      !! ** Method  :   
      !!      use 
      !!
      !! History :
      !!   ??  !  15/10/2015  (JTinker) Routine taken from old dia_wri_foam
      !!----------------------------------------------------------------------
      !! * Modules used
      !use lib_mpp
      !use lib_fortr
      IMPLICIT NONE
      
      INTEGER, INTENT(in) ::   kt
      CHARACTER (len=60) , INTENT(IN   ) ::    name
      REAL(wp), DIMENSION(jpi, jpj), INTENT(IN   ) :: infield    ! Input 3d field and mask
      
      ! Local variables
      INTEGER, DIMENSION(jpi, jpj) :: internal_region_mask    ! Input 3d field and mask
      REAL(wp), ALLOCATABLE, DIMENSION(:) ::   zrmet_ave,zrmet_tot,zrmet_var,zrmet_cnt,zrmet_mask_id,zrmet_reg_id
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::   zrmet_out
      REAL(wp), ALLOCATABLE,   DIMENSION(:) ::   ave_mat,tot_mat,num_mat,var_mat,ssq_mat,cnt_mat,reg_id_mat,mask_id_mat   !: region_mask
      
      REAL(wp)                         ::   zmdi      ! set masked values
      INTEGER :: maskno,nreg  ! ocean time-step indexocean time step            
      INTEGER :: ji,jj,jk,ind,jm ! Dummy loop indices
      INTEGER :: reg_ind_cnt ! Dummy loop indices
      
      INTEGER  ::   ierr      
      REAL(wp)  :: tmpreal
      CHARACTER(LEN=180) :: FormatString,nreg_string
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::   dummy_zrmet
      zmdi=1.e+20 !missing data indicator for maskin
      
      !Allocate output arrays for iomput, set to zmdi, and set a region counter = 1
      ALLOCATE( zrmet_ave(n_regions_output),  STAT= ierr )
        IF( ierr /= 0 )   CALL ctl_stop( 'dia_wri_region_mean: failed to allocate zrmet_ave array' )
      ALLOCATE( zrmet_tot(n_regions_output),  STAT= ierr )
        IF( ierr /= 0 )   CALL ctl_stop( 'dia_wri_region_mean: failed to allocate zrmet_tot array' )
      ALLOCATE( zrmet_var(n_regions_output),  STAT= ierr )
        IF( ierr /= 0 )   CALL ctl_stop( 'dia_wri_region_mean: failed to allocate zrmet_var array' )
      ALLOCATE( zrmet_cnt(n_regions_output),  STAT= ierr )
        IF( ierr /= 0 )   CALL ctl_stop( 'dia_wri_region_mean: failed to allocate zrmet_cnt array' )
      ALLOCATE( zrmet_mask_id(n_regions_output),  STAT= ierr )
        IF( ierr /= 0 )   CALL ctl_stop( 'dia_wri_region_mean: failed to allocate zrmet_mask_id array' )
      ALLOCATE( zrmet_reg_id(n_regions_output),  STAT= ierr )
        IF( ierr /= 0 )   CALL ctl_stop( 'dia_wri_region_mean: failed to allocate zrmet_reg_id array' )
      
      ALLOCATE( zrmet_out(jpi,jpj,n_regions_output),  STAT= ierr )
        IF( ierr /= 0 )   CALL ctl_stop( 'dia_wri_region_mean: failed to allocate zrmet_reg_id array' )
      
      
      zrmet_ave(:) = zmdi
      zrmet_tot(:) = zmdi
      zrmet_var(:) = zmdi
      zrmet_cnt(:) = zmdi
      zrmet_mask_id(:) = zmdi
      zrmet_reg_id(:) = zmdi
      
      reg_ind_cnt = 1
      
      
      ! loop though the masks
      DO maskno = 1,nmasks
          
          
          ! For each mask, get the number of regions (nreg), and a local copy of the region. 
          nreg = nreg_mat(maskno)
          internal_region_mask = region_mask(:,:,maskno)
          
          ! allocate temporary stat arrays, and set to zero
          ALLOCATE( ave_mat(nreg),  STAT= ierr )
          IF( ierr /= 0 )   CALL ctl_stop( 'dia_wri_region_mean: failed to allocate ave_mat array' )
          ALLOCATE( tot_mat(nreg),      STAT= ierr )
          IF( ierr /= 0 )   CALL ctl_stop( 'dia_wri_region_mean: failed to allocate tot_mat array' )
          ALLOCATE( num_mat(nreg),  STAT= ierr )
          IF( ierr /= 0 )   CALL ctl_stop( 'dia_wri_region_mean: failed to allocate num_mat array' )
          ALLOCATE( var_mat(nreg),  STAT= ierr )
          IF( ierr /= 0 )   CALL ctl_stop( 'dia_wri_region_mean: failed to allocate var_mat array' )
          ALLOCATE( ssq_mat(nreg),  STAT= ierr )
          IF( ierr /= 0 )   CALL ctl_stop( 'dia_wri_region_mean: failed to allocate ssq_mat array' )
          ALLOCATE( cnt_mat(nreg),  STAT= ierr )
          IF( ierr /= 0 )   CALL ctl_stop( 'dia_wri_region_mean: failed to allocate cnt_mat array' )
          ALLOCATE( reg_id_mat(nreg),  STAT= ierr )
          IF( ierr /= 0 )   CALL ctl_stop( 'dia_wri_region_mean: failed to allocate reg_id_mat array' )
          ALLOCATE( mask_id_mat(nreg),  STAT= ierr )
          IF( ierr /= 0 )   CALL ctl_stop( 'dia_wri_region_mean: failed to allocate mask_id_mat array' )
          
          ave_mat(:) = 0.
          tot_mat(:) = 0.
          num_mat(:) = 0.
          var_mat(:) = 0.
          cnt_mat(:) = 0.
          ssq_mat(:) = 0.
          reg_id_mat(:) = 0.
          mask_id_mat(:) = 0.
          
          ! loop though the array. for each sea grid box where tmask == 1), 
          ! read which region the grid box is in, add the value of the gridbox (and its square) 
          ! to the total for that region, and then increment the counter for that region.
          !CALL cpu_time(start_reg_mean_loop)
          !WRITE(numout,*) kt,start_reg_mean_loop
          DO ji = 1,jpi
              DO jj = 1,jpj
                    IF ( tmask(ji,jj,1) == 1.0_wp ) THEN
                        ind = internal_region_mask(ji,jj)+1
                        tot_mat(ind) = tot_mat(ind) + (infield(ji,jj))
                        ssq_mat(ind) = ssq_mat(ind) + ( infield(ji,jj) *  infield(ji,jj))
                        cnt_mat(ind) = cnt_mat(ind) + 1.
                    ENDIF
              END DO
          END DO
          ! sum the totals, the counts, and the squares across the processors          
          CALL mpp_sum( tot_mat,nreg )
          CALL mpp_sum( ssq_mat,nreg )
          CALL mpp_sum( cnt_mat,nreg )
          
          
          !calculate the mean and variance from the total, sum of squares and the count. 
          
          ave_mat = tot_mat(:)/cnt_mat(:)
          var_mat = ssq_mat(:)/cnt_mat(:) - (ave_mat(:)*ave_mat(:))
          
          
          !mask array of mask and region number. 
          DO jj = 1,nreg
              reg_id_mat(jj) = real(jj-1)
              mask_id_mat(jj) = real(maskno)
          END DO
          
          
          !write text and binary, and note region statistics for current mask for later iom_put
          IF( lwp ) THEN 
          
              !Write out ascii and binary if requred
              IF ( ln_diaregmean_bin ) THEN
                  !Writing out regional mean time series to binary files
                  WRITE(numdct_reg_bin) name,kt,maskno,n_regions_output
                  WRITE(numdct_reg_bin) ave_mat
                  WRITE(numdct_reg_bin) tot_mat
                  WRITE(numdct_reg_bin) var_mat
                  WRITE(numdct_reg_bin) ssq_mat
                  WRITE(numdct_reg_bin) cnt_mat
              ENDIF
              
              IF ( ln_diaregmean_ascii  ) THEN
                  !Writing out regional mean time series to text files

                  WRITE(nreg_string, "(I5)") nreg
                  FormatString = "(A17,"//trim(nreg_string)//"F15.3)"
                  WRITE(numdct_reg_txt, FMT="(A17,I6,I6)") name,kt,maskno            
                  WRITE(numdct_reg_txt, FMT=trim(FormatString)) trim(name)//" "//"ave_mat:", ave_mat
                  WRITE(numdct_reg_txt, FMT=trim(FormatString)) trim(name)//" "//"tot_mat:", tot_mat
                  WRITE(numdct_reg_txt, FMT=trim(FormatString)) trim(name)//" "//"var_mat:", var_mat
                  WRITE(numdct_reg_txt, FMT=trim(FormatString)) trim(name)//" "//"ssq_mat:", ssq_mat
                  WRITE(numdct_reg_txt, FMT=trim(FormatString)) trim(name)//" "//"cnt_mat:", cnt_mat
                  WRITE(numdct_reg_txt, FMT=trim(FormatString)) trim(name)//" "//"reg_mat:", reg_id_mat
                  WRITE(numdct_reg_txt, FMT=trim(FormatString)) trim(name)//" "//"msk_mat:", mask_id_mat

              ENDIF
              
              DO jm = 1,nreg
                  zrmet_ave(    reg_ind_cnt) =     ave_mat(jm)
                  zrmet_tot(    reg_ind_cnt) =     tot_mat(jm)
                  zrmet_var(    reg_ind_cnt) =     var_mat(jm)
                  zrmet_cnt(    reg_ind_cnt) =     cnt_mat(jm)
                  zrmet_reg_id( reg_ind_cnt) =  reg_id_mat(jm)
                  zrmet_mask_id(reg_ind_cnt) = mask_id_mat(jm)
                
                  reg_ind_cnt = reg_ind_cnt + 1 
              END DO
          
          ENDIF
        
          DEALLOCATE(ave_mat,tot_mat,num_mat,var_mat,ssq_mat,cnt_mat,reg_id_mat,mask_id_mat)
          !DEALLOCATE(tot_mat_2)
      END DO
      
      
      !With current field_def.xml and iodef.xml, these fields must be output, so set to dummy values if not required.
      
      IF ( ln_diaregmean_nc ) THEN
      
      
      
          DO jm = 1,n_regions_output
            zrmet_out(:,:,jm) = zrmet_ave(jm)
          END DO      
          CALL iom_put( "reg_" // trim(name) // '_ave', zrmet_out )
          zrmet_out(:,:,:) = 0
      
          DO jm = 1,n_regions_output
            zrmet_out(:,:,jm) = zrmet_tot(jm)
          END DO          
          CALL iom_put( "reg_" // trim(name) // '_tot', zrmet_out )
          zrmet_out(:,:,:) = 0
      
          DO jm = 1,n_regions_output
            zrmet_out(:,:,jm) = zrmet_var(jm)
          END DO          
          CALL iom_put( "reg_" // trim(name) // '_var', zrmet_out )
          zrmet_out(:,:,:) = 0
      
          DO jm = 1,n_regions_output
            zrmet_out(:,:,jm) = zrmet_cnt(jm)
          END DO          
          CALL iom_put( "reg_" // trim(name) // '_cnt', zrmet_out )
          zrmet_out(:,:,:) = 0
      
          DO jm = 1,n_regions_output
            zrmet_out(:,:,jm) = zrmet_reg_id(jm)
          END DO          
          CALL iom_put( "reg_" // trim(name) // '_reg_id', zrmet_out )
          zrmet_out(:,:,:) = 0
      
          DO jm = 1,n_regions_output
            zrmet_out(:,:,jm) = zrmet_mask_id(jm)
          END DO          
          CALL iom_put( "reg_" // trim(name) // '_mask_id', zrmet_out )
          zrmet_out(:,:,:) = 0
      ELSE
        
          ALLOCATE( dummy_zrmet(jpi,jpj,n_regions_output),  STAT= ierr )
            IF( ierr /= 0 )   CALL ctl_stop( 'dia_wri_region_mean: failed to allocate dummy_zrmet array' )

          DO jm = 1,n_regions_output
              dummy_zrmet(:,:,jm) =     real(jm,wp)
          END DO

          DO jm = 1,9
              CALL iom_put( "reg_" // trim(name) // '_ave', dummy_zrmet )
              CALL iom_put( "reg_" // trim(name) // '_tot', dummy_zrmet )
              CALL iom_put( "reg_" // trim(name) // '_var', dummy_zrmet )
              CALL iom_put( "reg_" // trim(name) // '_cnt', dummy_zrmet )
              CALL iom_put( "reg_" // trim(name) // '_reg_id', dummy_zrmet )
              CALL iom_put( "reg_" // trim(name) // '_mask_id', dummy_zrmet )
          END DO
    
          DEALLOCATE( dummy_zrmet)
      ENDIF
      
      DEALLOCATE(zrmet_ave,zrmet_tot,zrmet_var,zrmet_cnt,zrmet_mask_id,zrmet_reg_id,zrmet_out)
   END SUBROUTINE dia_wri_region_mean

#else
   !!----------------------------------------------------------------------
   !!   Default option :                                       Dummy module
   !!----------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_diareg = .FALSE.    !: diareg flag
   PUBLIC 
   !! $Id$
CONTAINS

   SUBROUTINE dia_regmean_init          ! Dummy routine
      WRITE(*,*) 'dia_regmean_init: You should not have seen this print! error?', kt
   END SUBROUTINE dia_regmean_init

   SUBROUTINE dia_regmean( kt )         ! Dummy routine
      INTEGER, INTENT( in ) :: kt   ! ocean time-step index
      WRITE(*,*) 'dia_regmean: You should not have seen this print! error?', kt
   END SUBROUTINE dia_regmean
#endif

   !!======================================================================
END MODULE diaregmean
