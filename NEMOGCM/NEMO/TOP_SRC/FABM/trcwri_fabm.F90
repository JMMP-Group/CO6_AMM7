MODULE trcwri_fabm
   !!======================================================================
   !!                       *** MODULE trcwri_fabm ***
   !!    fabm :   Output of FABM tracers
   !!======================================================================
   !! History :   1.0  !  2009-05 (C. Ethe)  Original code
   !!----------------------------------------------------------------------
#if defined key_top && key_fabm && defined key_iomput
   !!----------------------------------------------------------------------
   !!   'key_fabm'                                           FABM model
   !!----------------------------------------------------------------------
   !! trc_wri_fabm   :  outputs of concentration fields
   !!----------------------------------------------------------------------
   USE trc         ! passive tracers common variables 
   USE iom         ! I/O manager
   USE trdtrc_oce
   USE trcsms_fabm, only: trc_sms_fabm_check_mass
   USE par_fabm
   USE st2d_fabm
   USE fabm, only: fabm_get_bulk_diagnostic_data, fabm_get_horizontal_diagnostic_data

   IMPLICIT NONE
   PRIVATE

#if defined key_tracer_budget
   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:,:), SAVE :: tr_temp
   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:), SAVE :: fabm_st2d_temp
#endif

   INTERFACE trc_wri_fabm
       MODULE PROCEDURE wri_fabm,wri_fabm_fl
   END INTERFACE trc_wri_fabm


   PUBLIC trc_wri_fabm 

#  include "top_substitute.h90"
CONTAINS

   SUBROUTINE wri_fabm_fl(kt,fl)
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE trc_wri_trc  ***
      !!
      !! ** Purpose :   output passive tracers fields 
      !!---------------------------------------------------------------------
      INTEGER, INTENT( in )               :: fl
      INTEGER, INTENT( in )               :: kt

#if defined key_tracer_budget
      INTEGER              :: jn
      CHARACTER (len=20)   :: cltra
      REAL(wp), DIMENSION(jpi,jpj,jpk)    :: trpool !temporary storage pool 3D
      REAL(wp), DIMENSION(jpi,jpj)    :: st2dpool !temporary storage pool 2D
      !!---------------------------------------------------------------------
 
      ! write the tracer concentrations in the file
      ! ---------------------------------------
! depth integrated
! for strict budgetting write this out at end of timestep as an average between 'now' and 'after' at kt
      DO jn = 1, jp_fabm1
        IF(ln_trdtrc (jn))THEN
         trpool(:,:,:) = 0.5 * ( trn(:,:,:,jp_fabm0+jn-1)*fse3t_a(:,:,:) + &
                             tr_temp(:,:,:,jn)*fse3t(:,:,:) )
         cltra = TRIM( model%state_variables(jn)%name )//"_e3t"     ! depth integrated output
         IF( kt == nittrc000 ) write(6,*)'output pool ',cltra
         CALL iom_put( cltra, trpool)
        ENDIF
      END DO
#else
      CONTINUE
#endif

   END SUBROUTINE wri_fabm_fl

   SUBROUTINE wri_fabm(kt)
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE trc_wri_trc  ***
      !!
      !! ** Purpose :   output passive tracers fields 
      !!---------------------------------------------------------------------
      INTEGER, INTENT( in )               :: kt
      INTEGER              :: jn

#if defined key_tracer_budget
      IF( kt == nittrc000 ) THEN
         ALLOCATE(tr_temp(jpi,jpj,jpk,jp_fabm),fabm_st2d_temp(jpi,jpj,jp_fabm_surface+jp_fabm_bottom))
      ENDIF
      tr_temp(:,:,:,:)=trn(:,:,:,jp_fabm0:jp_fabm1) ! slwa save for tracer budget (unfiltered trn)
      fabm_st2d_temp(:,:,:)=fabm_st2dn(:,:,:)
#endif
      DO jn = 1, jp_fabm
         CALL iom_put( model%state_variables(jn)%name, trn(:,:,:,jp_fabm0+jn-1) )
      END DO
      DO jn = 1, jp_fabm_surface
         CALL iom_put( model%surface_state_variables(jn)%name, fabm_st2dn(:,:,jn) )
      END DO
      DO jn = 1, jp_fabm_bottom
         CALL iom_put( model%bottom_state_variables(jn)%name, fabm_st2dn(:,:,jp_fabm_surface+jn) )
      END DO

      ! write 3D diagnostics in the file
      ! ---------------------------------------
      DO jn = 1, size(model%diagnostic_variables)
         IF (model%diagnostic_variables(jn)%save) &
             CALL iom_put( model%diagnostic_variables(jn)%name, fabm_get_bulk_diagnostic_data(model,jn))
      END DO

      ! write 2D diagnostics in the file
      ! ---------------------------------------
      DO jn = 1, size(model%horizontal_diagnostic_variables)
         IF (model%horizontal_diagnostic_variables(jn)%save) &
             CALL iom_put( model%horizontal_diagnostic_variables(jn)%name, fabm_get_horizontal_diagnostic_data(model,jn))
      END DO
      !
      CALL trc_sms_fabm_check_mass

   END SUBROUTINE wri_fabm

#else
   !!----------------------------------------------------------------------
   !!  Dummy module :                                     No passive tracer
   !!----------------------------------------------------------------------
   PUBLIC trc_wri_fabm
CONTAINS
   SUBROUTINE trc_wri_fabm                     ! Empty routine  
   END SUBROUTINE trc_wri_fabm
#endif

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!======================================================================
END MODULE trcwri_fabm
