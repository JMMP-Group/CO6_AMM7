MODULE obs_fco2_io
   !!======================================================================
   !!                       ***  MODULE obs_fco2_io  ***
   !! Observation operators : I/O for fco2 files
   !!======================================================================
   !! History : 
   !!             !  09-01  (K. Mogensen) Initial version
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   read_fco2file    :  Read a obfbdata structure from a fco2 file
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

#include "obsfco2_io.h90"

END MODULE obs_fco2_io
