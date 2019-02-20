MODULE par_fabm

   USE fabm

   IMPLICIT NONE

   TYPE (type_model) :: model !FABM model instance

   INTEGER, PUBLIC :: jp_fabm0, jp_fabm1, jp_fabm, &
                      jp_fabm_surface, jp_fabm_bottom, &
                      jp_fabm_m1, jp_fabm_2d, jp_fabm_3d

   ! Variables needed for OBS/ASM
   INTEGER, PUBLIC :: jp_fabm_chl1, jp_fabm_chl2, &
                      jp_fabm_chl3, jp_fabm_chl4, &
                      jp_fabm_p1c,  jp_fabm_p1n,  &
                      jp_fabm_p1p,  jp_fabm_p1s,  &
                      jp_fabm_p2c,  jp_fabm_p2n,  &
                      jp_fabm_p2p,  jp_fabm_p3c,  &
                      jp_fabm_p3n,  jp_fabm_p3p,  &
                      jp_fabm_p4c,  jp_fabm_p4n,  &
                      jp_fabm_p4p,  jp_fabm_z4c,  &
                      jp_fabm_z5c,  jp_fabm_z5n,  &
                      jp_fabm_z5p,  jp_fabm_z6c,  &
                      jp_fabm_z6n,  jp_fabm_z6p,  &
                      jp_fabm_n1p,  jp_fabm_n3n,  &
                      jp_fabm_n4n,  jp_fabm_n5s,  &
                      jp_fabm_o2o,  jp_fabm_o3c,  &
                      jp_fabm_o3a,  jp_fabm_o3ph, &
                      jp_fabm_o3pc

#if defined key_fabm
   !!---------------------------------------------------------------------
   !!   'key_fabm'                     FABM tracers
   !!---------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_fabm     = .TRUE.   !: FABM flag 
   LOGICAL, PUBLIC, ALLOCATABLE, DIMENSION(:) ::   lk_rad_fabm !: FABM negativity correction flag array 
#else
   !!---------------------------------------------------------------------
   !!   Default                           No user defined tracers (FABM)
   !!---------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_fabm     = .FALSE.  !: FABM flag 
#endif

   !!======================================================================
END MODULE par_fabm
