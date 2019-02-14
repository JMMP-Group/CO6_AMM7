MODULE obs_pco2
   !!=====================================================================
   !!                       ***  MODULE  obs_pco2  ***
   !! Observation diagnostics: Storage space for pco2 observations
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

   PUBLIC npco2vars, npco2extr, npco2sets, pco2data, pco2datqc

   !! * Shared Module variables
   INTEGER :: npco2vars                               ! Number of pco2data variables
   INTEGER :: npco2extr                               ! Number of pco2data extra 
                                                      ! variables
   INTEGER :: npco2sets                               ! Number of pco2data sets
   TYPE(obs_surf), POINTER, DIMENSION(:) :: pco2data  ! Initial pco2 data
   TYPE(obs_surf), POINTER, DIMENSION(:) :: pco2datqc ! Sea ice data after quality control

END MODULE obs_pco2

