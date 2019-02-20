MODULE diatmb 
   !!======================================================================
   !!                       ***  MODULE  diaharm  ***
   !! Harmonic analysis of tidal constituents 
   !!======================================================================
   !! History :  3.6  !  2014  (E O'Dea)  Original code
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain
   USE in_out_manager  ! I/O units
   USE iom             ! I/0 library
   USE wrk_nemo        ! working arrays
#if defined key_fabm
   USE trc, ONLY: trn
   USE par_fabm
   USE fabm, ONLY: fabm_get_interior_diagnostic_data
#endif


   IMPLICIT NONE
   PRIVATE

   LOGICAL , PUBLIC ::   ln_diatmb     !: Top Middle and Bottom output 
   PUBLIC   dia_tmb_init            ! routine called by nemogcm.F90
   PUBLIC   dia_tmb                 ! routine called by diawri.F90
   PUBLIC   dia_calctmb             ! routine called by dia25h.F90

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.6 , NEMO Consortium (2014)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dia_tmb_init 
      !!---------------------------------------------------------------------------
      !!                  ***  ROUTINE dia_tmb_init  ***
      !!     
      !! ** Purpose: Initialization of tmb namelist 
      !!        
      !! ** Method : Read namelist
      !!   History
      !!   3.6  !  08-14  (E. O'Dea) Routine to initialize dia_tmb
      !!---------------------------------------------------------------------------
      !!
      INTEGER ::   ios                 ! Local integer output status for namelist read
      !
      NAMELIST/nam_diatmb/ ln_diatmb
      !!----------------------------------------------------------------------
      !
      REWIND ( numnam_ref )              ! Read Namelist nam_diatmb in reference namelist : TMB diagnostics
      READ   ( numnam_ref, nam_diatmb, IOSTAT=ios, ERR= 901 )
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nam_diatmb in reference namelist', lwp )
 
      REWIND( numnam_cfg )              ! Namelist nam_diatmb in configuration namelist  TMB diagnostics
      READ  ( numnam_cfg, nam_diatmb, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nam_diatmb in configuration namelist', lwp )
      IF(lwm) WRITE ( numond, nam_diatmb )

      IF(lwp) THEN                   ! Control print
         WRITE(numout,*)
         WRITE(numout,*) 'dia_tmb_init : Output Top, Middle, Bottom Diagnostics'
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) 'Namelist nam_diatmb : set tmb outputs '
         WRITE(numout,*) 'Switch for TMB diagnostics (T) or not (F)  ln_diatmb  = ', ln_diatmb
      ENDIF

   END SUBROUTINE dia_tmb_init

   SUBROUTINE dia_calctmb( pinfield,pouttmb )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE dia_tmb  ***
      !!                   
      !! ** Purpose :    Find the Top, Mid and Bottom fields of water Column
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
      REAL(wp), DIMENSION(jpi, jpj, 3  ), INTENT(  OUT) :: pouttmb     ! Output top, middle, bottom



      ! Local variables
      INTEGER :: ji,jj,jk  ! Dummy loop indices

      ! Local Real
      REAL(wp)                         ::   zmdi  !  set masked values

      zmdi=1.e+20 !missing data indicator for masking

      ! Calculate top
      pouttmb(:,:,1) = pinfield(:,:,1)*tmask(:,:,1)  + zmdi*(1.0-tmask(:,:,1))

      ! Calculate middle
      DO jj = 1,jpj
         DO ji = 1,jpi
            jk              = max(1,mbathy(ji,jj)/2)
            pouttmb(ji,jj,2) = pinfield(ji,jj,jk)*tmask(ji,jj,jk)  + zmdi*(1.0-tmask(ji,jj,jk))
         END DO
      END DO

      ! Calculate bottom
      DO jj = 1,jpj
         DO ji = 1,jpi
            jk              = max(1,mbathy(ji,jj) )
            pouttmb(ji,jj,3) = pinfield(ji,jj,jk)*tmask(ji,jj,jk)  + zmdi*(1.0-tmask(ji,jj,jk))
         END DO
      END DO

   END SUBROUTINE dia_calctmb



   SUBROUTINE dia_tmb
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE dia_tmb  ***
      !! ** Purpose :   Write diagnostics for Top, Mid and Bottom of water Column
      !!
      !! ** Method  :   
      !!      use mbathy to find surface, mid and bottom of model levels
      !!      calls calctmb to retrieve TMB values before sending to iom_put
      !!
      !! History :
      !!   3.6  !  08-14  (E. O'Dea) 
      !!         
      !!--------------------------------------------------------------------
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zwtmb    ! temporary workspace 
      REAL(wp)                         ::   zmdi      ! set masked values
      INTEGER                          ::   jn        ! loop counter

      zmdi=1.e+20 !missing data indicator for maskin

      IF (ln_diatmb) THEN
         CALL wrk_alloc( jpi , jpj, 3 , zwtmb )
         CALL dia_calctmb(  tsn(:,:,:,jp_tem),zwtmb )
         !ssh already output but here we output it masked
         CALL iom_put( "sshnmasked" , sshn(:,:)*tmask(:,:,1) + zmdi*(1.0 - tmask(:,:,1)) )   ! tmb Temperature
         CALL iom_put( "top_temp" , zwtmb(:,:,1) )    ! tmb Temperature
         CALL iom_put( "mid_temp" , zwtmb(:,:,2) )    ! tmb Temperature
         CALL iom_put( "bot_temp" , zwtmb(:,:,3) )    ! tmb Temperature
!         CALL iom_put( "sotrefml" , hmld_tref(:,:) )    ! "T criterion Mixed Layer Depth

         CALL dia_calctmb(  tsn(:,:,:,jp_sal),zwtmb )
         CALL iom_put( "top_sal" , zwtmb(:,:,1) )    ! tmb Salinity 
         CALL iom_put( "mid_sal" , zwtmb(:,:,2) )    ! tmb Salinity
         CALL iom_put( "bot_sal" , zwtmb(:,:,3) )    ! tmb Salinity

         CALL dia_calctmb(  un(:,:,:),zwtmb )
         CALL iom_put( "top_u" , zwtmb(:,:,1) )    ! tmb  U Velocity
         CALL iom_put( "mid_u" , zwtmb(:,:,2) )    ! tmb  U Velocity
         CALL iom_put( "bot_u" , zwtmb(:,:,3) )    ! tmb  U Velocity
!Called in  dynspg_ts.F90        CALL iom_put( "baro_u" , un_b )    ! Barotropic  U Velocity

         CALL dia_calctmb(  vn(:,:,:),zwtmb )
         CALL iom_put( "top_v" , zwtmb(:,:,1) )    ! tmb  V Velocity
         CALL iom_put( "mid_v" , zwtmb(:,:,2) )    ! tmb  V Velocity
         CALL iom_put( "bot_v" , zwtmb(:,:,3) )    ! tmb  V Velocity
!Called in  dynspg_ts.F90       CALL iom_put( "baro_v" , vn_b )    ! Barotropic  V Velocity

#if defined key_fabm
         DO jn = 1, jp_fabm
            CALL dia_calctmb( trn(:,:,:,jp_fabm_m1+jn), zwtmb )
            CALL iom_put( "top_"//TRIM(model%state_variables(jn)%name) , zwtmb(:,:,1) )
            CALL iom_put( "mid_"//TRIM(model%state_variables(jn)%name) , zwtmb(:,:,2) )
            CALL iom_put( "bot_"//TRIM(model%state_variables(jn)%name) , zwtmb(:,:,3) )
         END DO
         DO jn = 1, jp_fabm_3d
            CALL dia_calctmb( fabm_get_interior_diagnostic_data(model, jn), zwtmb )
            CALL iom_put( "top_"//TRIM(model%diagnostic_variables(jn)%name) , zwtmb(:,:,1) )
            CALL iom_put( "mid_"//TRIM(model%diagnostic_variables(jn)%name) , zwtmb(:,:,2) )
            CALL iom_put( "bot_"//TRIM(model%diagnostic_variables(jn)%name) , zwtmb(:,:,3) )
         END DO
#endif
      ELSE
         CALL ctl_warn('dia_tmb: tmb diagnostic is set to false you should not have seen this')
      ENDIF

   END SUBROUTINE dia_tmb
   !!======================================================================
END MODULE diatmb
