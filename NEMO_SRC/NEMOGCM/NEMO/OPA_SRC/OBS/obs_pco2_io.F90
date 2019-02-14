MODULE obs_pco2_io
   !!======================================================================
   !!                       ***  MODULE obs_pco2_io  ***
   !! Observation operators : I/O for pco2 files
   !!======================================================================
   !! History : 
   !!             !  09-01  (K. Mogensen) Initial version
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   read_pco2file    :  Read a obfbdata structure from a pco2 file
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

#include "obspco2_io.h90"

END MODULE obs_pco2_io
