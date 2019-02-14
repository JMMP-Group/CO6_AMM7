MODULE obs_spm_io
   !!======================================================================
   !!                       ***  MODULE obs_spm_io  ***
   !! Observation operators : I/O for spm files
   !!======================================================================
   !! History : 
   !!             !  09-01  (K. Mogensen) Initial version
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   read_spmfile    :  Read a obfbdata structure from a spm file
   !!----------------------------------------------------------------------
   USE par_kind
   USE obs_utils
   USE obs_fbm
   USE julian
   USE netcdf

   IMPLICIT NONE

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

#include "obsspm_io.h90"

END MODULE obs_spm_io
