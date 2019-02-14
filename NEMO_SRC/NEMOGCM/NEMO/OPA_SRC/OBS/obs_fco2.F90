MODULE obs_fco2
   !!=====================================================================
   !!                       ***  MODULE  obs_fco2  ***
   !! Observation diagnostics: Storage space for fco2 observations
   !!                          arrays and additional flags etc.
   !!=====================================================================
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

   
   !! * Modules used 
   USE obs_surf_def ! Definition of surface data types and tools

   IMPLICIT NONE
   
   SAVE

   !! * Routine accessibility
   PRIVATE

   PUBLIC nfco2vars, nfco2extr, nfco2sets, fco2data, fco2datqc

   !! * Shared Module variables
   INTEGER :: nfco2vars                               ! Number of fco2data variables
   INTEGER :: nfco2extr                               ! Number of fco2data extra 
                                                      ! variables
   INTEGER :: nfco2sets                               ! Number of fco2data sets
   TYPE(obs_surf), POINTER, DIMENSION(:) :: fco2data  ! Initial fco2 data
   TYPE(obs_surf), POINTER, DIMENSION(:) :: fco2datqc ! Sea ice data after quality control

END MODULE obs_fco2

