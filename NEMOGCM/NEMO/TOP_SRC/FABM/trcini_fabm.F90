MODULE trcini_fabm
   !!======================================================================
   !!                         ***  MODULE trcini_fabm  ***
   !! TOP :   initialisation of the FABM tracers
   !!======================================================================
   !! History :   2.0  !  2007-12  (C. Ethe, G. Madec) Original code
   !!----------------------------------------------------------------------
#if defined key_fabm
   !!----------------------------------------------------------------------
   !!   'key_fabm'                                               FABM tracers
   !!----------------------------------------------------------------------
   !! trc_ini_fabm   : FABM model initialisation
   !!----------------------------------------------------------------------
   USE par_trc         ! TOP parameters
   USE oce_trc
   USE trc
   USE par_fabm
   USE trcsms_fabm
   USE fabm_config,ONLY: fabm_create_model_from_yaml_file
   USE fabm,ONLY: type_external_variable, fabm_initialize_library
   USE inputs_fabm,ONLY: initialize_inputs,link_inputs, &
     type_input_variable,type_input_data,type_river_data, &
     first_input_data,first_river_data
#if defined key_git_version
   USE fabm_version,ONLY: fabm_commit_id=>git_commit_id, &
                          fabm_branch_name=>git_branch_name
   USE fabm_types,ONLY: type_version,first_module_version
#endif


   IMPLICIT NONE
   PRIVATE

#if defined key_git_version
!#include "gitversion.h90"
   CHARACTER(len=*),parameter :: git_commit_id = _NEMO_COMMIT_ID_
   CHARACTER(len=*),parameter :: git_branch_name = _NEMO_BRANCH_
#endif

   PUBLIC   trc_ini_fabm   ! called by trcini.F90 module
   PUBLIC   nemo_fabm_init

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE nemo_fabm_init()
      INTEGER :: jn
      INTEGER, PARAMETER :: xml_unit = 1979
      TYPE (type_input_data),POINTER :: input_data
      TYPE (type_river_data),POINTER :: river_data
      CLASS (type_input_variable),POINTER :: input_pointer

      ! Allow FABM to parse fabm.yaml. This ensures numbers of variables are known.
      call fabm_create_model_from_yaml_file(model)

      jp_fabm = size(model%state_variables)
      jp_fabm_bottom = size(model%bottom_state_variables)
      jp_fabm_surface = size(model%surface_state_variables)
      jp_fabm0 = jptra + 1
      jp_fabm1 = jptra + jp_fabm
      jp_fabm_m1=jptra
      jptra = jptra + jp_fabm
      jp_fabm_2d = size(model%horizontal_diagnostic_variables)
      jp_fabm_3d = size(model%diagnostic_variables)
      jpdia2d = jpdia2d + jp_fabm_2d
      jpdia3d = jpdia3d + jp_fabm_3d
      jpdiabio = jpdiabio + jp_fabm

      !Initialize input data structures.
      call initialize_inputs

      ! Get indexes for select state variables
      jp_fabm_chl1 = fabm_state_index( 'P1_Chl' )
      jp_fabm_chl2 = fabm_state_index( 'P2_Chl' )
      jp_fabm_chl3 = fabm_state_index( 'P3_Chl' )
      jp_fabm_chl4 = fabm_state_index( 'P4_Chl' )
      jp_fabm_p1c  = fabm_state_index( 'P1_c' )
      jp_fabm_p1n  = fabm_state_index( 'P1_n' )
      jp_fabm_p1p  = fabm_state_index( 'P1_p' )
      jp_fabm_p1s  = fabm_state_index( 'P1_s' )
      jp_fabm_p2c  = fabm_state_index( 'P2_c' )
      jp_fabm_p2n  = fabm_state_index( 'P2_n' )
      jp_fabm_p2p  = fabm_state_index( 'P2_p' )
      jp_fabm_p3c  = fabm_state_index( 'P3_c' )
      jp_fabm_p3n  = fabm_state_index( 'P3_n' )
      jp_fabm_p3p  = fabm_state_index( 'P3_p' )
      jp_fabm_p4c  = fabm_state_index( 'P4_c' )
      jp_fabm_p4n  = fabm_state_index( 'P4_n' )
      jp_fabm_p4p  = fabm_state_index( 'P4_p' )
      jp_fabm_z4c  = fabm_state_index( 'Z4_c' )
      jp_fabm_z5c  = fabm_state_index( 'Z5_c' )
      jp_fabm_z5n  = fabm_state_index( 'Z5_n' )
      jp_fabm_z5p  = fabm_state_index( 'Z5_p' )
      jp_fabm_z6c  = fabm_state_index( 'Z6_c' )
      jp_fabm_z6n  = fabm_state_index( 'Z6_n' )
      jp_fabm_z6p  = fabm_state_index( 'Z6_p' )
      jp_fabm_n1p  = fabm_state_index( 'N1_p' )
      jp_fabm_n3n  = fabm_state_index( 'N3_n' )
      jp_fabm_n4n  = fabm_state_index( 'N4_n' )
      jp_fabm_n5s  = fabm_state_index( 'N5_s' )
      jp_fabm_o2o  = fabm_state_index( 'O2_o' )
      jp_fabm_o3c  = fabm_state_index( 'O3_c' )
      jp_fabm_o3a  = fabm_state_index( 'O3_bioalk' )

      ! Get indexes for select diagnostic variables
      jp_fabm_o3ph = fabm_diag_index( 'O3_pH' )
      jp_fabm_o3pc = fabm_diag_index( 'O3_pCO2' )

      IF (lwp) THEN
         ! write field_def_fabm.xml on lead process
         OPEN(UNIT=xml_unit,FILE='field_def_fabm.xml',ACTION='WRITE',STATUS='REPLACE')

         WRITE (xml_unit,1000) '<field_definition level="1" prec="4" operation="average" enabled=".TRUE." default_value="1.e20" >'

         WRITE (xml_unit,1000) ' <field_group id="ptrc_T" grid_ref="grid_T_3D">'
         DO jn=1,jp_fabm
            CALL write_variable_xml(xml_unit,model%state_variables(jn))
#if defined key_trdtrc
            CALL write_trends_xml(xml_unit,model%state_variables(jn))
#endif
            CALL write_25hourm_xml(xml_unit,model%state_variables(jn))
            CALL write_tmb_xml(xml_unit,model%state_variables(jn))
         END DO
         WRITE (xml_unit,1000) ' </field_group>'

         WRITE (xml_unit,1000) ' <field_group id="sf_T" grid_ref="grid_T_2D">'
         DO jn=1,jp_fabm_surface
            CALL write_variable_xml(xml_unit,model%surface_state_variables(jn))
            CALL write_25hourm_xml(xml_unit,model%surface_state_variables(jn))
         END DO
         DO jn=1,jp_fabm_bottom
            CALL write_variable_xml(xml_unit,model%bottom_state_variables(jn))
            CALL write_25hourm_xml(xml_unit,model%bottom_state_variables(jn))
         END DO
         WRITE (xml_unit,1000) ' </field_group>'

         WRITE (xml_unit,1000) ' <field_group id="diad_T" grid_ref="grid_T_2D">'
         DO jn=1,size(model%diagnostic_variables)
            CALL write_variable_xml(xml_unit,model%diagnostic_variables(jn),3)
            CALL write_25hourm_xml(xml_unit,model%diagnostic_variables(jn),3)
            CALL write_tmb_xml(xml_unit,model%diagnostic_variables(jn))
         END DO
         DO jn=1,size(model%horizontal_diagnostic_variables)
            CALL write_variable_xml(xml_unit,model%horizontal_diagnostic_variables(jn))
            CALL write_25hourm_xml(xml_unit,model%horizontal_diagnostic_variables(jn))
         END DO
         WRITE (xml_unit,1000) ' </field_group>'

         WRITE (xml_unit,1000) ' <field_group id="fabm_scalar" grid_ref="grid_0">'
         DO jn=1,size(model%conserved_quantities)
            CALL write_variable_xml(xml_unit,model%conserved_quantities(jn))
         END DO
         WRITE (xml_unit,1000) ' </field_group>'

         WRITE (xml_unit,1000) ' <field_group id="fabm_input" grid_ref="grid_T_2D">'
         input_data => first_input_data
         DO WHILE (ASSOCIATED(input_data))
           input_pointer => input_data
           CALL write_input_xml(xml_unit,input_pointer)
            input_data => input_data%next
         END DO
         river_data => first_river_data
         DO WHILE (ASSOCIATED(river_data))
           input_pointer => river_data
           CALL write_input_xml(xml_unit,input_pointer,3)
            river_data => river_data%next
         END DO
         WRITE (xml_unit,1000) ' </field_group>'

         WRITE (xml_unit,1000) '</field_definition>'

         CLOSE(xml_unit)
      END IF
      IF( lk_mpp )   CALL mppsync !Ensure field_def_fabm is ready.

1000 FORMAT (A)

   END SUBROUTINE nemo_fabm_init

   SUBROUTINE write_variable_xml(xml_unit,variable,flag_grid_ref)
      INTEGER,INTENT(IN) :: xml_unit
      INTEGER,INTENT(IN),OPTIONAL :: flag_grid_ref
      CLASS (type_external_variable),INTENT(IN) :: variable

      CHARACTER(LEN=20) :: missing_value,string_dimensions
      INTEGER :: number_dimensions

      ! Check variable dimension for grid_ref specificaiton.
      ! Default is to not specify the grid_ref in the field definition.
      IF (present(flag_grid_ref)) THEN
          number_dimensions=flag_grid_ref
      ELSE
          number_dimensions=-1 !default, don't specify grid_ref
      ENDIF

      WRITE (missing_value,'(E9.3)') variable%missing_value
      WRITE (string_dimensions,'(I1)') number_dimensions
      SELECT CASE (number_dimensions)
      CASE (3)
         WRITE (xml_unit,'(A)') '  <field id="'//TRIM(variable%name)//'" long_name="'//TRIM(variable%long_name)//'" unit="'//TRIM(variable%units)//'" default_value="'//TRIM(ADJUSTL(missing_value))//'" grid_ref="grid_T_3D" />'
      CASE (2)
         WRITE (xml_unit,'(A)') '  <field id="'//TRIM(variable%name)//'" long_name="'//TRIM(variable%long_name)//'" unit="'//TRIM(variable%units)//'" default_value="'//TRIM(ADJUSTL(missing_value))//'" grid_ref="grid_T_2D"/>'
      CASE (0)
         WRITE (xml_unit,'(A)') '  <field id="'//TRIM(variable%name)//'" long_name="'//TRIM(variable%long_name)//'" unit="'//TRIM(variable%units)//'" default_value="'//TRIM(ADJUSTL(missing_value))//'" grid_ref="1point"/>'
      CASE (-1)
         WRITE (xml_unit,'(A)') '  <field id="'//TRIM(variable%name)//'" long_name="'//TRIM(variable%long_name)//'" unit="'//TRIM(variable%units)//'" default_value="'//TRIM(ADJUSTL(missing_value))//'" />'
      CASE default
         IF(lwp) WRITE(numout,*) ' trc_ini_fabm: Failing to initialise output of variable '//TRIM(variable%name)//': Output of '//TRIM(ADJUSTL(string_dimensions))//'-dimensional variables not supported!!!'
      END SELECT

   END SUBROUTINE write_variable_xml

   SUBROUTINE write_25hourm_xml(xml_unit,variable,flag_grid_ref)
      INTEGER,INTENT(IN) :: xml_unit
      INTEGER,INTENT(IN),OPTIONAL :: flag_grid_ref
      CLASS (type_external_variable),INTENT(IN) :: variable

      CHARACTER(LEN=20) :: missing_value,string_dimensions
      INTEGER :: number_dimensions

      ! Check variable dimension for grid_ref specificaiton.
      ! Default is to not specify the grid_ref in the field definition.
      IF (present(flag_grid_ref)) THEN
          number_dimensions=flag_grid_ref
      ELSE
          number_dimensions=-1 !default, don't specify grid_ref
      ENDIF

      WRITE (missing_value,'(E9.3)') 1.e+20
      WRITE (string_dimensions,'(I1)') number_dimensions
      SELECT CASE (number_dimensions)
      CASE (3)
         WRITE (xml_unit,'(A)') '  <field id="'//TRIM(variable%name)//'25h'//'" long_name="'//TRIM(variable%long_name)//' 25-hour mean'//'" unit="'//TRIM(variable%units)//'" default_value="'//TRIM(ADJUSTL(missing_value))//'" grid_ref="grid_T_3D" />'
      CASE (2)
         WRITE (xml_unit,'(A)') '  <field id="'//TRIM(variable%name)//'25h'//'" long_name="'//TRIM(variable%long_name)//' 25-hour mean'//'" unit="'//TRIM(variable%units)//'" default_value="'//TRIM(ADJUSTL(missing_value))//'" grid_ref="grid_T_2D"/>'
      CASE (0)
         WRITE (xml_unit,'(A)') '  <field id="'//TRIM(variable%name)//'25h'//'" long_name="'//TRIM(variable%long_name)//' 25-hour mean'//'" unit="'//TRIM(variable%units)//'" default_value="'//TRIM(ADJUSTL(missing_value))//'" grid_ref="1point"/>'
      CASE (-1)
         WRITE (xml_unit,'(A)') '  <field id="'//TRIM(variable%name)//'25h'//'" long_name="'//TRIM(variable%long_name)//' 25-hour mean'//'" unit="'//TRIM(variable%units)//'" default_value="'//TRIM(ADJUSTL(missing_value))//'" />'
      CASE default
         IF(lwp) WRITE(numout,*) ' trc_ini_fabm: Failing to initialise output of variable '//TRIM(variable%name)//'25h'//': Output of '//TRIM(ADJUSTL(string_dimensions))//'-dimensional variables not supported!!!'
      END SELECT

   END SUBROUTINE write_25hourm_xml

   SUBROUTINE write_tmb_xml(xml_unit,variable)
      INTEGER,INTENT(IN) :: xml_unit
      CLASS (type_external_variable),INTENT(IN) :: variable

      CHARACTER(LEN=20) :: missing_value

      WRITE (missing_value,'(E9.3)') 1.e+20
      WRITE (xml_unit,'(A)') '  <field id="top_'//TRIM(variable%name)//'" long_name="Top-level '//TRIM(variable%long_name)//'" unit="'//TRIM(variable%units)//'" default_value="'//TRIM(ADJUSTL(missing_value))//'" grid_ref="grid_T_2D"/>'
      WRITE (xml_unit,'(A)') '  <field id="mid_'//TRIM(variable%name)//'" long_name="Middle-level '//TRIM(variable%long_name)//'" unit="'//TRIM(variable%units)//'" default_value="'//TRIM(ADJUSTL(missing_value))//'" grid_ref="grid_T_2D"/>'
      WRITE (xml_unit,'(A)') '  <field id="bot_'//TRIM(variable%name)//'" long_name="Bottom-level '//TRIM(variable%long_name)//'" unit="'//TRIM(variable%units)//'" default_value="'//TRIM(ADJUSTL(missing_value))//'" grid_ref="grid_T_2D"/>'

   END SUBROUTINE write_tmb_xml

   SUBROUTINE write_trends_xml(xml_unit,variable,flag_grid_ref)
      INTEGER,INTENT(IN) :: xml_unit
      INTEGER,INTENT(IN),OPTIONAL :: flag_grid_ref
      CLASS (type_external_variable),INTENT(IN) :: variable

      INTEGER :: number_dimensions,i
      CHARACTER(LEN=20) :: missing_value,string_dimensions
#if defined key_tracer_budget
      CHARACTER(LEN=3),DIMENSION(10),PARAMETER :: trd_tags = (/ &
        'LDF','BBL','FOR','ZDF','DMP','SMS','ATF', &
        'RDB','RDN','VMV' /)
      CHARACTER(LEN=3),DIMENSION(3),PARAMETER :: trd_e3t_tags = (/ &
        'XAD','YAD','ZAD' /)
#else
      CHARACTER(LEN=3),DIMENSION(13),PARAMETER :: trd_tags = (/ &
        'XAD','YAD','ZAD','LDF','BBL','FOR','ZDF','DMP','SMS','ATF', &
        'RDB','RDN','VMV' /)
#endif

      ! Check variable dimension for grid_ref specificaiton.
      ! Default is to not specify the grid_ref in the field definition.
      IF (present(flag_grid_ref)) THEN
          number_dimensions=flag_grid_ref
      ELSE
          number_dimensions=-1 !default, don't specify grid_ref
      ENDIF

      WRITE (missing_value,'(E9.3)') -2.E20
      WRITE (string_dimensions,'(I1)') number_dimensions
      SELECT CASE (number_dimensions)
      CASE (3)
        DO i=1,size(trd_tags)
         WRITE (xml_unit,'(A)') '  <field id="'//TRIM(trd_tags(i))//'_'//TRIM(variable%name)//'" long_name="'// &
            &                   TRIM(variable%long_name)//' '//TRIM(trd_tags(i))//' trend" unit="'//            &
            &                   TRIM(variable%units)//'/s" default_value="'//                                   &
            &                   TRIM(ADJUSTL(missing_value))//'" grid_ref="grid_T_3D" />'
        END DO
#if defined key_tracer_budget
        DO i=1,size(trd_e3t_tags)
         WRITE (xml_unit,'(A)') '  <field id="'//TRIM(trd_e3t_tags(i))//'_'//TRIM(variable%name)//                     &
            &                   '_e3t" long_name="'//TRIM(variable%long_name)//' cell depth integrated '//             &
            &                   TRIM(trd_e3t_tags(i))//' trend" unit="'//TRIM(variable%units)//'/s" default_value="'// &
            &                   TRIM(ADJUSTL(missing_value))//'" grid_ref="grid_T_3D" />'
        END DO
        WRITE (xml_unit,'(A)') '  <field id="'//TRIM(variable%name)//'_e3t" long_name="'//TRIM(variable%long_name)// &
           &                   ' cell depth integrated" unit="'//TRIM(variable%units)//'/s" default_value="'//       &
           &                   TRIM(ADJUSTL(missing_value))//'" grid_ref="grid_T_3D" />'
#endif
      CASE (-1)
        DO i=1,size(trd_tags)
         WRITE (xml_unit,'(A)') '  <field id="'//TRIM(trd_tags(i))//'_'//TRIM(variable%name)//'" long_name="'//   &
            &                   TRIM(variable%long_name)//' '//TRIM(trd_tags(i))//' trend" unit="'//              &
            &                   TRIM(variable%units)//'/s" default_value="'//TRIM(ADJUSTL(missing_value))//'" />'
        END DO
#if defined key_tracer_budget
        DO i=1,size(trd_e3t_tags)
         WRITE (xml_unit,'(A)') '  <field id="'//TRIM(trd_e3t_tags(i))//'_'//TRIM(variable%name)//                     &
            &                   '_e3t" long_name="'//TRIM(variable%long_name)//' cell depth integrated '//             &
            &                   TRIM(trd_e3t_tags(i))//' trend" unit="'//TRIM(variable%units)//'/s" default_value="'// &
            &                   TRIM(ADJUSTL(missing_value))//'" />'
        END DO
        WRITE (xml_unit,'(A)') '  <field id="'//TRIM(variable%name)//'_e3t" long_name="'//TRIM(variable%long_name)// &
           &                   ' cell depth integrated" unit="'//TRIM(variable%units)//'/s" default_value="'//       &
           &                   TRIM(ADJUSTL(missing_value))//'" />'
#endif
      CASE default
         IF(lwp) WRITE(numout,*) ' trc_ini_fabm: Failing to initialise trends of variable '//TRIM(variable%name)// &
            &                    ': Output of '//TRIM(ADJUSTL(string_dimensions))//'-dimensional trends not supported!!!'
      END SELECT

   END SUBROUTINE write_trends_xml

   SUBROUTINE write_input_xml(xml_unit,variable,flag_grid_ref)
      INTEGER,INTENT(IN) :: xml_unit
      INTEGER,INTENT(IN),OPTIONAL :: flag_grid_ref
      CLASS(type_input_variable),POINTER,INTENT(IN) :: variable

      INTEGER :: number_dimensions,i
      CHARACTER(LEN=20) :: missing_value,string_dimensions

      ! Check variable dimension for grid_ref specificaiton.
      ! Default is to not specify the grid_ref in the field definition.
      IF (present(flag_grid_ref)) THEN
          number_dimensions=flag_grid_ref
      ELSE
          number_dimensions=-1 !default, don't specify grid_ref
      ENDIF

      WRITE (missing_value,'(E9.3)') -2.E20
      WRITE (string_dimensions,'(I1)') number_dimensions
      SELECT CASE (number_dimensions)
      CASE (3)
        WRITE (xml_unit,'(A)') '  <field id="'//'INP_'//TRIM(variable%sf(1)%clvar)//'" long_name="'//TRIM(variable%sf(1)%clvar)//' input" unit="" default_value="'//TRIM(ADJUSTL(missing_value))//'" grid_ref="grid_T_3D" />'
      CASE (-1)
        WRITE (xml_unit,'(A)') '  <field id="'//'INP_'//TRIM(variable%sf(1)%clvar)//'" long_name="'//TRIM(variable%sf(1)%clvar)//' input" unit="" default_value="'//TRIM(ADJUSTL(missing_value))//'" />'
      CASE default
         IF(lwp) WRITE(numout,*) ' trc_ini_fabm: Failing to initialise input diagnostic of variable '//TRIM(variable%sf(1)%clvar)//': Output of '//TRIM(ADJUSTL(string_dimensions))//'-dimensional diagnostic not supported!!!'
      END SELECT

   END SUBROUTINE write_input_xml

   SUBROUTINE trc_ini_fabm
      !!----------------------------------------------------------------------
      !!                     ***  trc_ini_fabm  ***
      !!
      !! ** Purpose :   initialization for FABM model
      !!
      !! ** Method  : - Read the namcfc namelist and check the parameter values
      !!----------------------------------------------------------------------
#if defined key_git_version
      TYPE (type_version),POINTER :: version
#endif
      INTEGER :: jn

      !                       ! Allocate FABM arrays
      IF( trc_sms_fabm_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'trc_ini_fabm: unable to allocate FABM arrays' )

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_ini_fabm: initialisation of FABM model'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'
#if defined key_git_version
      IF(lwp) WRITE(numout,*) ' NEMO version:   ',git_commit_id,' (',git_branch_name,' branch)'
      IF(lwp) WRITE(numout,*) ' FABM version:   ',fabm_commit_id,' (',fabm_branch_name,' branch)'
#endif

      call fabm_initialize_library()
#if defined key_git_version
      version => first_module_version

      do while (associated(version))
         IF(lwp) WRITE(numout,*)  ' '//trim(version%module_name)//' version:   ',trim(version%version_string)
         version => version%next
      end do
#endif

      ! Log mapping of FABM states:
      IF (lwp) THEN
         IF (jp_fabm.gt.0) WRITE(numout,*) " FABM tracers:"
         DO jn=1,jp_fabm
            WRITE(numout,*) "   State",jn,":",trim(model%state_variables(jn)%name), &
               " (",trim(model%state_variables(jn)%long_name), &
               ") [",trim(model%state_variables(jn)%units),"]"
         ENDDO
         IF (jp_fabm_surface.gt.0) WRITE(numout,*) "FABM seasurface states:"
         DO jn=1,jp_fabm_surface
            WRITE(numout,*) "   State",jn,":",trim(model%surface_state_variables(jn)%name), &
               " (",trim(model%surface_state_variables(jn)%long_name), &
               ") [",trim(model%surface_state_variables(jn)%units),"]"
         ENDDO
         IF (jp_fabm_bottom.gt.0) WRITE(numout,*) "FABM seafloor states:"
         DO jn=1,jp_fabm_bottom
            WRITE(numout,*) "   State",jn,":",trim(model%bottom_state_variables(jn)%name), &
               " (",trim(model%bottom_state_variables(jn)%long_name), &
               ") [",trim(model%bottom_state_variables(jn)%units),"]"
         ENDDO
      ENDIF

   END SUBROUTINE trc_ini_fabm

   INTEGER FUNCTION fabm_state_index( state_name )
      !!----------------------------------------------------------------------
      !!                     ***  fabm_state_index  ***  
      !!
      !! ** Purpose :   return index of a given FABM state variable
      !!
      !! ** Method  : - loop through state variables until found
      !!----------------------------------------------------------------------
      
      IMPLICIT NONE
      
      CHARACTER(LEN=256), INTENT(IN) :: state_name
      
      INTEGER                        :: jn

      !!----------------------------------------------------------------------
      
      fabm_state_index = -1
      DO jn=1,jp_fabm
         IF (TRIM(model%state_variables(jn)%name) == TRIM(state_name)) THEN
            fabm_state_index = jn
            EXIT
         ENDIF
      END DO
      IF (fabm_state_index == -1) THEN
         CALL ctl_stop( 'Could not find '//TRIM(state_name)//' state variable' )
      ELSE
         IF (lwp) WRITE(numout,*) 'Index for '//TRIM(state_name)//' is: ', fabm_state_index
      ENDIF
   
   END FUNCTION fabm_state_index

   INTEGER FUNCTION fabm_diag_index( diag_name )
      !!----------------------------------------------------------------------
      !!                     ***  fabm_state_index  ***  
      !!
      !! ** Purpose :   return index of a given FABM diagnostic variable
      !!
      !! ** Method  : - loop through diagnostic variables until found
      !!----------------------------------------------------------------------
      
      IMPLICIT NONE
      
      CHARACTER(LEN=256), INTENT(IN) :: diag_name
      
      INTEGER                        :: jn

      !!----------------------------------------------------------------------
      
      fabm_diag_index = -1
      DO jn = 1, SIZE(model%diagnostic_variables)
         IF (TRIM(model%diagnostic_variables(jn)%name) == TRIM(diag_name)) THEN
            fabm_diag_index = jn
            EXIT
         ENDIF
      END DO
      IF (fabm_diag_index == -1) THEN
         CALL ctl_stop( 'Could not find '//TRIM(diag_name)//' diagnostic' )
      ELSE
         IF (lwp) WRITE(numout,*) 'Index for '//TRIM(diag_name)//' is: ', fabm_diag_index
      ENDIF
   
   END FUNCTION fabm_diag_index

#else
   !!----------------------------------------------------------------------
   !!   Dummy module                                        No FABM model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE nemo_fabm_init
   END SUBROUTINE nemo_fabm_init

   SUBROUTINE trc_ini_fabm            ! Empty routine
   END SUBROUTINE trc_ini_fabm
#endif

   !!======================================================================
END MODULE trcini_fabm
