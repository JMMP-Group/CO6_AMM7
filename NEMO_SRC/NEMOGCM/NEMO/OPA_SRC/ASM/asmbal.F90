MODULE asmbal
   !!======================================================================
   !!                       ***  MODULE asmbal  ***
   !! Assimilation balancing interface: Write to file the balancing increments
   !!                                   calculated for biogeochemistry
   !!======================================================================
   !!----------------------------------------------------------------------
   !!   'key_asminc' : Switch on the assimilation increment interface
   !!----------------------------------------------------------------------
   !!   asm_bal_wri  : Write out the background state
   !!----------------------------------------------------------------------
   !! * Modules used
   USE par_kind, ONLY: &   ! Precision variables
      & wp
   USE iom                 ! I/O module
   USE asminc              ! Main assimilation increments module
   USE asmpar              ! Parameters for the assimilation interface
#if defined key_fabm
   USE par_fabm
#elif defined key_medusa && defined key_foam_medusa
   USE par_medusa
#elif defined key_hadocc
   USE par_hadocc
#endif

   IMPLICIT NONE

   !! * Routine accessibility
   PRIVATE
   PUBLIC asm_bal_wri   !: Write out the background state

CONTAINS

   SUBROUTINE asm_bal_wri( kt )
      !!-----------------------------------------------------------------------
      !!
      !!                  ***  ROUTINE asm_bal_wri ***
      !!
      !! ** Purpose : Write to file the balancing increments
      !!              calculated for biogeochemistry
      !!
      !! ** Method  : Write to file the balancing increments
      !!              calculated for biogeochemistry
      !!
      !! ** Action  :
      !!                   
      !! References :
      !!
      !! History    :
      !!        ! 2014-08 (D. Ford)  Initial version, based on asm_bkg_wri
      !!-----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( IN ) :: kt        ! Current time-step
      !! * Local declarations
      CHARACTER(LEN=50) :: cl_asmbal     ! Filename (with extension)
      LOGICAL           :: llok          ! Check if file exists
      INTEGER           :: inum          ! File unit number
      REAL(wp)          :: zdate         ! Date
      !!-----------------------------------------------------------------------
     
      ! Set things up
      zdate = REAL( ndastp )
      WRITE(cl_asmbal, FMT='(A,".nc")' ) TRIM( c_asmbal )

      ! Check if file exists
      INQUIRE( FILE = TRIM( cl_asmbal ), EXIST = llok )
      IF( .NOT. llok ) THEN
         IF(lwp) WRITE(numout,*) ' Setting up assimilation balancing increments file '// &
            &                    TRIM( c_asmbal ) // ' at timestep = ', kt

         ! Define the output file       
         CALL iom_open( c_asmbal, inum, ldwrt = .TRUE., kiolib = jprstlib)

         ! Write the information
         CALL iom_rstput( kt, kt, inum, 'rdastp' , zdate   )

         IF ( ln_logchltotinc .OR. ln_logchlpftinc ) THEN
#if defined key_fabm
            CALL iom_rstput( kt, kt, inum, 'logchl_balinc_chl1', logchl_balinc(:,:,:,jp_fabm_m1+jp_fabm_chl1) )
            CALL iom_rstput( kt, kt, inum, 'logchl_balinc_chl2', logchl_balinc(:,:,:,jp_fabm_m1+jp_fabm_chl2) )
            CALL iom_rstput( kt, kt, inum, 'logchl_balinc_chl3', logchl_balinc(:,:,:,jp_fabm_m1+jp_fabm_chl3) )
            CALL iom_rstput( kt, kt, inum, 'logchl_balinc_chl4', logchl_balinc(:,:,:,jp_fabm_m1+jp_fabm_chl4) )
#elif defined key_medusa && defined key_foam_medusa
            CALL iom_rstput( kt, kt, inum, 'logchl_balinc_chn', logchl_balinc(:,:,:,jpchn) )
            CALL iom_rstput( kt, kt, inum, 'logchl_balinc_chd', logchl_balinc(:,:,:,jpchd) )
            CALL iom_rstput( kt, kt, inum, 'logchl_balinc_phn', logchl_balinc(:,:,:,jpphn) )
            CALL iom_rstput( kt, kt, inum, 'logchl_balinc_phd', logchl_balinc(:,:,:,jpphd) )
            CALL iom_rstput( kt, kt, inum, 'logchl_balinc_pds', logchl_balinc(:,:,:,jppds) )
            CALL iom_rstput( kt, kt, inum, 'logchl_balinc_zmi', logchl_balinc(:,:,:,jpzmi) )
            CALL iom_rstput( kt, kt, inum, 'logchl_balinc_zme', logchl_balinc(:,:,:,jpzme) )
            CALL iom_rstput( kt, kt, inum, 'logchl_balinc_din', logchl_balinc(:,:,:,jpdin) )
            CALL iom_rstput( kt, kt, inum, 'logchl_balinc_sil', logchl_balinc(:,:,:,jpsil) )
            CALL iom_rstput( kt, kt, inum, 'logchl_balinc_fer', logchl_balinc(:,:,:,jpfer) )
            CALL iom_rstput( kt, kt, inum, 'logchl_balinc_det', logchl_balinc(:,:,:,jpdet) )
            CALL iom_rstput( kt, kt, inum, 'logchl_balinc_dtc', logchl_balinc(:,:,:,jpdtc) )
            CALL iom_rstput( kt, kt, inum, 'logchl_balinc_dic', logchl_balinc(:,:,:,jpdic) )
            CALL iom_rstput( kt, kt, inum, 'logchl_balinc_alk', logchl_balinc(:,:,:,jpalk) )
            CALL iom_rstput( kt, kt, inum, 'logchl_balinc_oxy', logchl_balinc(:,:,:,jpoxy) )
#elif defined key_hadocc
            CALL iom_rstput( kt, kt, inum, 'logchl_balinc_nut', logchl_balinc(:,:,:,jp_had_nut) )
            CALL iom_rstput( kt, kt, inum, 'logchl_balinc_phy', logchl_balinc(:,:,:,jp_had_phy) )
            CALL iom_rstput( kt, kt, inum, 'logchl_balinc_zoo', logchl_balinc(:,:,:,jp_had_zoo) )
            CALL iom_rstput( kt, kt, inum, 'logchl_balinc_det', logchl_balinc(:,:,:,jp_had_det) )
            CALL iom_rstput( kt, kt, inum, 'logchl_balinc_dic', logchl_balinc(:,:,:,jp_had_dic) )
            CALL iom_rstput( kt, kt, inum, 'logchl_balinc_alk', logchl_balinc(:,:,:,jp_had_alk) )
#endif
         ENDIF

         CALL iom_close( inum )
      ELSE
         CALL ctl_warn( TRIM( cl_asmbal ) // ' already exists ', &
            &           ' Therefore not writing out balancing increments at this timestep', &
            &           ' Something has probably gone wrong somewhere' )
         IF(lwp) WRITE(numout,*) ' Failed to set up assimilation balancing increments file '// &
            &                    TRIM( c_asmbal ) // ' at timestep = ', kt
      ENDIF
                                 
   END SUBROUTINE asm_bal_wri
END MODULE asmbal
