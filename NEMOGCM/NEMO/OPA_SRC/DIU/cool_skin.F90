MODULE cool_skin
   !!======================================================================
   !!                    ***  MODULE  cool_skin  ***
   !!     Cool skin thickness and delta T correction using Artele et al. (2002)
   !!     [see also Tu and Tsuang (2005)]
   !!
   !!=====================================================================
   !! History :        !  2012-01  (P. Sykes)  Original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   diurnal_sst_coolskin_step  : timestep the cool skin corrections
   !!----------------------------------------------------------------------
   USE par_kind
   USE phycst
   USE dom_oce
   USE in_out_manager
   USE sbc_oce
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   
   IMPLICIT NONE
   
   !Namelist parameters

   !Parameters
   REAL(wp), PRIVATE, PARAMETER :: pp_k = 0.596_wp          !Thermal conductivity of seawater
   REAL(wp), PRIVATE, PARAMETER :: pp_v = 1.05e-6_wp        !Kinematic viscosity of seawater
   REAL(wp), PRIVATE, PARAMETER :: pp_C = 86400             !seconds [see Tu and Tsuang (2005)]
   REAL(wp), PRIVATE, PARAMETER :: pp_cw = 3993._wp         !specific heat capacity of seawater
   REAL(wp), PRIVATE, PARAMETER :: pp_h = 10._wp            !reference depth [using 10m from Artale et al. (2002)]
   REAL(wp), PRIVATE, PARAMETER :: pp_rhoa = 1.20421_wp     !density of air (at 20C)
   REAL(wp), PRIVATE, PARAMETER :: pp_cda = 1.45e-3_wp      !assumed air-sea drag coefficient for calculating wind speed
   
   !Key variables
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: x_csdsst     !Cool skin delta SST
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: x_csthick  !Cool skin thickness

   PRIVATE
   PUBLIC diurnal_sst_coolskin_step, diurnal_sst_coolskin_init

      !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   
   CONTAINS 
   
   SUBROUTINE diurnal_sst_coolskin_init
      !!----------------------------------------------------------------------
      !! *** ROUTINE diurnal_sst_coolskin_init ***
      !!
      !! ** Purpose :   initalise the coolskin model
      !!
      !! ** Method : 
      !!
      !! ** Reference :
      !! 
      !!----------------------------------------------------------------------
      
      IMPLICIT NONE
      
      ALLOCATE( x_csdsst(jpi,jpj), x_csthick(jpi,jpj) )
      x_csdsst = 0.
      x_csthick = 0.
      
   END SUBROUTINE diurnal_sst_coolskin_init
 
   SUBROUTINE diurnal_sst_coolskin_step(psqflux, pstauflux, psrho, rdt)
      !!----------------------------------------------------------------------
      !! *** ROUTINE diurnal_sst_takaya_step ***
      !!
      !! ** Purpose :   Timestep the Artale cool skin model
      !!
      !! ** Method : 
      !!
      !! ** Reference : 
      !!----------------------------------------------------------------------
     
      IMPLICIT NONE
     
      !Dummy variables
      REAL(wp), INTENT(IN), DIMENSION(jpi,jpj) :: psqflux     !Heat (non-solar)(Watts)
      REAL(wp), INTENT(IN), DIMENSION(jpi,jpj) :: pstauflux   !Wind stress (kg/ m s^2)
      REAL(wp), INTENT(IN), DIMENSION(jpi,jpj) :: psrho       !Water density (kg/m^3)
      REAL(wp), INTENT(IN) :: rdt                             !Timestep
     
      !Local variables
      REAL(wp), DIMENSION(jpi,jpj) :: z_fv                    !Friction velocity     
      REAL(wp), DIMENSION(jpi,jpj) :: z_gamma                 !Dimensionless function of windspeed
      REAL(wp), DIMENSION(jpi,jpj) :: z_lamda                 !Sauders (dimensionless) proportionality constant
      REAL(wp), DIMENSION(jpi,jpj) :: z_wspd                  !Windspeed (m/s)
      REAL(wp) :: z_ztx                                       !Temporary u wind stress
      REAL(wp) :: z_zty                                       !Temporary v wind stress
      REAL(wp) :: z_zmod                                      !Temporary total wind stress
     
      INTEGER :: ji,jj
      
      DO jj = 1,jpj
         DO ji = 1,jpi
            
            !Calcualte wind speed from wind stress and friction velocity
            IF( tmask(ji,jj,1) == 1. .AND. pstauflux(ji,jj) /= 0 .AND. psrho(ji,jj) /=0 ) THEN
               z_fv(ji,jj) = SQRT( pstauflux(ji,jj) / psrho(ji,jj) )
               z_wspd(ji,jj) = SQRT( pstauflux(ji,jj) / ( pp_cda * pp_rhoa ) )
            ELSE
               z_fv(ji,jj) = 0.
               z_wspd(ji,jj) = 0.     
            ENDIF

 
            !Calculate gamma function which is dependent upon wind speed
            IF( tmask(ji,jj,1) == 1. ) THEN
               IF( ( z_wspd(ji,jj) <= 7.5 ) ) z_gamma(ji,jj) = ( 0.2 * z_wspd(ji,jj) ) + 0.5
               IF( ( z_wspd(ji,jj) > 7.5 ) .AND. ( z_wspd(ji,jj) < 10. ) ) z_gamma(ji,jj) = ( 1.6 * z_wspd(ji,jj) ) - 10.
               IF( ( z_wspd(ji,jj) >= 10. ) ) z_gamma(ji,jj) = 6.
            ENDIF


            !Calulate lamda function
            IF( tmask(ji,jj,1) == 1. .AND. z_fv(ji,jj) /= 0 ) THEN
               z_lamda(ji,jj) = ( z_fv(ji,jj) * pp_k * pp_C ) / ( z_gamma(ji,jj) * psrho(ji,jj) * pp_cw * pp_h * pp_v )
            ELSE
               z_lamda(ji,jj) = 0.
            ENDIF



            !Calculate the cool skin thickness - only when heat flux is out of the ocean
            IF( tmask(ji,jj,1) == 1. .AND. z_fv(ji,jj) /= 0 .AND. psqflux(ji,jj) < 0 ) THEN
                x_csthick(ji,jj) = ( z_lamda(ji,jj) * pp_v ) / z_fv(ji,jj)
            ELSE
                x_csthick(ji,jj) = 0.
            ENDIF



            !Calculate the cool skin correction - only when the heat flux is out of the ocean
            IF( tmask(ji,jj,1) == 1. .AND. x_csthick(ji,jj) /= 0. .AND. psqflux(ji,jj) < 0. ) THEN
               x_csdsst(ji,jj) = ( psqflux(ji,jj) * x_csthick(ji,jj) ) / pp_k
             ELSE
               x_csdsst(ji,jj) = 0.
            ENDIF

         ENDDO
      ENDDO

   END SUBROUTINE diurnal_sst_coolskin_step


END MODULE cool_skin
