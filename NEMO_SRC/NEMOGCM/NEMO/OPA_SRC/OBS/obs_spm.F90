MODULE obs_spm
   !!=====================================================================
   !!                       ***  MODULE  obs_spm  ***
   !! Observation diagnostics: Storage space for spm observations
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

   PUBLIC nspmvars, nspmextr, nspmsets, spmdata, spmdatqc

   !! * Shared Module variables
   INTEGER :: nspmvars                               ! Number of spmdata variables
   INTEGER :: nspmextr                               ! Number of spmdata extra 
                                                     ! variables
   INTEGER :: nspmsets                               ! Number of spmdata sets
   TYPE(obs_surf), POINTER, DIMENSION(:) :: spmdata  ! Initial spm data
   TYPE(obs_surf), POINTER, DIMENSION(:) :: spmdatqc ! Sea ice data after quality control

END MODULE obs_spm

