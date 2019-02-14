MODULE obs_logchl
   !!=====================================================================
   !!                       ***  MODULE  obs_logchl  ***
   !! Observation diagnostics: Storage space for logchl observations
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

   PUBLIC nlogchlvars, nlogchlextr, nlogchlsets, logchldata, logchldatqc

   !! * Shared Module variables
   INTEGER :: nlogchlvars                               ! Number of logchldata variables
   INTEGER :: nlogchlextr                               ! Number of logchldata extra 
                                                        ! variables
   INTEGER :: nlogchlsets                               ! Number of logchldata sets
   TYPE(obs_surf), POINTER, DIMENSION(:) :: logchldata  ! Initial logchl data
   TYPE(obs_surf), POINTER, DIMENSION(:) :: logchldatqc ! Sea ice data after quality control

END MODULE obs_logchl

