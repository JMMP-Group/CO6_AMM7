MODULE obs_logchl_io
   !!======================================================================
   !!                       ***  MODULE obs_logchl_io  ***
   !! Observation operators : I/O for logchl files
   !!======================================================================
   !! History : 
   !!             !  09-01  (K. Mogensen) Initial version
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   read_logchlfile    :  Read a obfbdata structure from a logchl file
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

#include "obslogchl_io.h90"

END MODULE obs_logchl_io
