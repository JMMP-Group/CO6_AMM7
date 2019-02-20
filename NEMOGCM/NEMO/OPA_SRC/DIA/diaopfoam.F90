MODULE diaopfoam
   !!======================================================================
   !!                       ***  MODULE  diaopfoam  ***
   !! Output stream for operational use
   !!======================================================================
   !! History :  3.6  !  2016  (P Sykes)  Original code
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain
   USE diainsitutem, ONLY: rinsitu_t, theta2t
   USE in_out_manager  ! I/O units
   USE iom             ! I/0 library
   USE wrk_nemo        ! working arrays
   USE diatmb
   USE diurnal_bulk
   USE cool_skin


   IMPLICIT NONE
   PRIVATE

   LOGICAL , PUBLIC ::   ln_diaopfoam     !: Diaopfoam output 
   PUBLIC   dia_diaopfoam_init            ! routine called by nemogcm.F90
   PUBLIC   dia_diaopfoam                 ! routine called by diawri.F90
   PUBLIC   calc_max_cur                  ! routine called by diaopfoam.F90

   !! * Substitutions
#  include "domzgr_substitute.h90"

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.6 , NEMO Consortium (2014)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dia_diaopfoam_init
      !!---------------------------------------------------------------------------
      !!                  ***  ROUTINE dia_wri_diaop_init  ***
      !!     
      !! ** Purpose: Initialization of diaopfoam namelist 
      !!        
      !! ** Method : Read namelist
      !!   History
      !!   3.4  !  03-14  (P. Sykes) Routine to initialize dia_wri_diaop
      !!---------------------------------------------------------------------------
      !!
      INTEGER            ::   ios      ! Local integer output status for namelist read
      INTEGER            ::   ierror   ! local integer
      !!
      NAMELIST/nam_diadiaop/ ln_diaopfoam
      !!
      !!----------------------------------------------------------------------
      !
      ln_diaopfoam = .false.         ! default value for diaopfoam stream
      REWIND ( numnam_ref )              ! Read Namelist nam_diadiaop in reference namelist : 3D hourly diagnostics
      READ   ( numnam_ref, nam_diadiaop, IOSTAT=ios, ERR= 901 )
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nam_diadiaop in reference namelist', lwp )

      REWIND( numnam_cfg )              ! Namelist nam_diadiaop in configuration namelist  3D hourly diagnostics
      READ  ( numnam_cfg, nam_diadiaop, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nam_diadiaop in configuration namelist', lwp )
      IF(lwm) WRITE ( numond, nam_diadiaop )
      !
      IF(lwp) THEN                   ! Control print
         WRITE(numout,*)
         WRITE(numout,*) 'dia_diaopfoam_init : Output Diaopfoam Diagnostics'
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist nam_diadiaop : set diaopfoam outputs '
         WRITE(numout,*) '      Switch for diaopfoam diagnostics (T) or not (F)  ln_diaopfoam  = ', ln_diaopfoam
      ENDIF
    END SUBROUTINE dia_diaopfoam_init

    SUBROUTINE dia_diaopfoam
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE dia_diaopfoam  ***
      !! ** Purpose :   Write 3D hourly diagnostics for operational use
      !!
      !!
      !! History :
      !!   3.6  !  11-16  (P. Sykes) 
      !!         
      !!--------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj)   :: zw2d       ! 2D workspace
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zw3d       ! 3D workspace
      REAL(wp) :: zmdi
      REAL(wp), POINTER, DIMENSION(:,:)   :: zwu
      REAL(wp), POINTER, DIMENSION(:,:)   :: zwv
      REAL(wp), POINTER, DIMENSION(:,:)   :: zwz

      CALL wrk_alloc( jpi , jpj      , zwu )
      CALL wrk_alloc( jpi , jpj      , zwv )
      CALL wrk_alloc( jpi , jpj      , zwz )

      zmdi=1.e+20                               !  missing data indicator for masking
      ! Diaopfoam stream if needed
      IF (ln_diaopfoam) THEN
         CALL theta2t
         CALL iom_put( "insitut_op" , rinsitu_t(:,:,:)                      )    ! insitu temperature
         CALL iom_put( "toce_op"   , tsn(:,:,:,jp_tem)                     )    ! temperature
         CALL iom_put( "soce_op"   , tsn(:,:,:,jp_sal)                     )    ! salinity
         IF (ln_diurnal) THEN
            CALL iom_put( "sst_wl_op"   , x_dsst                    )    ! warm layer 
            CALL iom_put( "sst_cs_op"   , x_csdsst                  )    ! cool skin 
         ENDIF
         zw2d(:,:)=sshn(:,:)*tmask(:,:,1) + zmdi*(1.0-tmask(:,:,1))
         CALL iom_put( "ssh_op"  , zw2d(:,:)                          ) ! sea surface height
         CALL iom_put( "uoce_op"   , un                                    )    ! i-current      
         CALL iom_put( "voce_op"   , vn                                    )    ! j-current
         !CALL iom_put( "woce_op"   , wn                                    )    ! k-current
         CALL calc_max_cur(zwu,zwv,zwz,zmdi)
         CALL iom_put( "maxu" , zwu                                     ) ! max u current
         CALL iom_put( "maxv" , zwv                                     ) ! max v current
         CALL iom_put( "maxz" , zwz                                     ) ! max current depth
      ENDIF
   END SUBROUTINE dia_diaopfoam

   SUBROUTINE calc_max_cur(zmax_u, zmax_v, zmax_z, inmdi)
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE calc_max_cur  ***
      !!
      !! ** Purpose :   To locate within the water column the magnitude and
      !!                vertical location of the strongest horizontal
      !!                current. The vertical component is ignored since it
      !!                is an order of magnitude smaller than the horizontal
      !!                flow, in general. 
      !!
      !! ** Method :    A. Map U,V to T-grid.
      !!                B. Calculate the magnitude of the current for every
      !!                   grid cell.
      !!                C. Locate the vertical index using FORTRAN's builtin
      !!                   MAXLOC function.
      !!                D. Copy the U,V,Z components of the relevant
      !!                   indices.
      !!
      !! ** Returns :   Value of u, v component and depth of maximum
      !!                horizontal current on T-grid.
      !!
      !!---------------------------------------------------------------------
      IMPLICIT NONE
      REAL(wp), DIMENSION(jpi, jpj), INTENT(  OUT) :: zmax_u, zmax_v, zmax_z
      REAL(wp), DIMENSION(jpi, jpj, jpk)           :: zmax_u_t, zmax_v_t
      REAL(wp), DIMENSION(jpi, jpj, jpk)           :: zcmag
      REAL(wp), DIMENSION(jpi, jpj)                :: zmaxk
      REAL(wp),                      INTENT(IN   ) :: inmdi
      INTEGER :: ji, jj, jk

      ! Initialise output arrays
      zmax_u(:, :) = inmdi
      zmax_v(:, :) = inmdi
      zmax_z(:, :) = inmdi
      ! Map to T-grid
      zmax_u_t(:, :, :) = 0._wp
      zmax_v_t(:, :, :) = 0._wp
      DO jk = 1,jpk
         DO jj = 2,jpj
            DO ji = 2,jpi
               zmax_u_t(ji, jj, jk) = 0.5 * (un(ji, jj, jk) + un(ji-1, jj, jk))
               zmax_v_t(ji, jj, jk) = 0.5 * (vn(ji, jj, jk) + vn(ji, jj-1, jk))
            END DO
         END DO
      END DO
      ! Calculate absolute velocity
      zcmag = sqrt((zmax_u_t)**2 + (zmax_v_t)**2)
      ! Find max. current
      zmaxk = maxloc(zcmag, dim=3)
      ! Output values
      DO jj = 1,jpj
         DO ji = 1,jpi
             zmax_u(ji, jj) = zmax_u_t(ji, jj, INT(zmaxk(ji, jj)))
             zmax_v(ji, jj) = zmax_v_t(ji, jj, INT(zmaxk(ji, jj)))
             zmax_z(ji, jj) = fsdept(ji, jj, INT(zmaxk(ji, jj)))
         END DO
      END DO      
   END SUBROUTINE calc_max_cur

END MODULE diaopfoam
