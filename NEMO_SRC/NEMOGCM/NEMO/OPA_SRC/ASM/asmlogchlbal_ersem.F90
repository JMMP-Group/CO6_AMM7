MODULE asmlogchlbal_ersem
   !!======================================================================
   !!                       ***  MODULE asmlogchlbal_ersem  ***
   !! Calculate increments to ERSEM based on surface logchl increments
   !!======================================================================
   !! History :  3.6  ! 2016-09  (D. Ford)     Original code
   !!----------------------------------------------------------------------
#if defined key_asminc && defined key_fabm
   !!----------------------------------------------------------------------
   !!   'key_asminc'         : assimilation increment interface
   !!   'key_fabm'           : FABM-ERSEM coupling
   !!----------------------------------------------------------------------
   !!   asm_logchl_bal_ersem : routine to calculate increments to ERSEM
   !!----------------------------------------------------------------------
   USE par_kind,    ONLY: wp             ! kind parameters
   USE par_oce,     ONLY: jpi, jpj, jpk  ! domain array sizes
   USE dom_oce,     ONLY: gdepw_n        ! domain information
   USE zdfmxl                            ! mixed layer depth
   USE iom                               ! i/o
   USE trc,         ONLY: trn            ! ERSEM state variables
   USE par_fabm                          ! FABM parameters
   USE par_trc,     ONLY: jptra          ! Tracer parameters

   IMPLICIT NONE
   PRIVATE                   

   PUBLIC asm_logchl_bal_ersem

CONTAINS

   SUBROUTINE asm_logchl_bal_ersem( ld_logchlpftinc, npfts, mld_choice_bgc, &
      &                             k_maxchlinc, logchl_bkginc, logchl_balinc )
      !!---------------------------------------------------------------------------
      !!                    ***  ROUTINE asm_logchl_bal_ersem  ***
      !!
      !! ** Purpose :   calculate increments to ERSEM from logchl increments
      !!
      !! ** Method  :   convert logchl increments to chl increments
      !!                split between the ERSEM PFTs
      !!                spread through the mixed layer
      !!                [forthcoming: calculate increments to nutrients and zooplankton]
      !!
      !! ** Action  :   populate logchl_balinc
      !!
      !! References :   forthcoming...
      !!---------------------------------------------------------------------------
      !!
      LOGICAL,  INTENT(in   )                               :: ld_logchlpftinc
      INTEGER,  INTENT(in   )                               :: npfts
      INTEGER,  INTENT(in   )                               :: mld_choice_bgc
      REAL(wp), INTENT(in   )                               :: k_maxchlinc
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj,npfts)     :: logchl_bkginc
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj,jpk,jptra) :: logchl_balinc
      !!
      INTEGER                      :: ji, jj, jk
      INTEGER                      :: jkmax
      REAL(wp)                     :: chl_tot, chl_inc
      REAL(wp), DIMENSION(jpi,jpj) :: zmld
      !!---------------------------------------------------------------------------

      ! Split surface logchl incs into surface Chl1-4 incs
      !
      ! In order to transform logchl incs to chl incs, need to account for the background,
      ! cannot simply do 10^logchl_bkginc. Need to:
      ! 1) Add logchl inc to log10(background) to get log10(analysis)
      ! 2) Take 10^log10(analysis) to get analysis
      ! 3) Subtract background from analysis to get chl incs
      ! If k_maxchlinc > 0 then cap total absolute chlorophyll increment at that value
      !
      ! Only apply increments if all of Chl1-4 background values are > 0
      ! In theory, they always will be, and if any are not that's a sign
      ! that something's going wrong which the assimilation might make worse
      !
      IF ( ld_logchlpftinc ) THEN
         !
         ! Assimilating separate PFTs, so separately transform each from LogChl to Chl
         !
         IF ( npfts /= 4 ) THEN
            CALL ctl_stop( 'If assimilating PFTs into ERSEM, nn_asmpfts must be 4' )
         ENDIF
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF ( ( trn(ji,jj,1,jp_fabm_m1+jp_fabm_chl1) > 0.0 ) .AND. &
                  & ( trn(ji,jj,1,jp_fabm_m1+jp_fabm_chl2) > 0.0 ) .AND. &
                  & ( trn(ji,jj,1,jp_fabm_m1+jp_fabm_chl3) > 0.0 ) .AND. &
                  & ( trn(ji,jj,1,jp_fabm_m1+jp_fabm_chl4) > 0.0 ) ) THEN
                  IF ( logchl_bkginc(ji,jj,1) /= 0.0 ) THEN
                     logchl_balinc(ji,jj,1,jp_fabm_m1+jp_fabm_chl1) = 10**( LOG10( trn(ji,jj,1,jp_fabm_m1+jp_fabm_chl1) ) + &
                        &                                                   logchl_bkginc(ji,jj,1) )                      - &
                        &                                             trn(ji,jj,1,jp_fabm_m1+jp_fabm_chl1)
                  ENDIF
                  IF ( logchl_bkginc(ji,jj,2) /= 0.0 ) THEN
                     logchl_balinc(ji,jj,1,jp_fabm_m1+jp_fabm_chl2) = 10**( LOG10( trn(ji,jj,1,jp_fabm_m1+jp_fabm_chl2) ) + &
                        &                                                   logchl_bkginc(ji,jj,2) )                      - &
                        &                                             trn(ji,jj,1,jp_fabm_m1+jp_fabm_chl2)
                  ENDIF
                  IF ( logchl_bkginc(ji,jj,3) /= 0.0 ) THEN
                     logchl_balinc(ji,jj,1,jp_fabm_m1+jp_fabm_chl3) = 10**( LOG10( trn(ji,jj,1,jp_fabm_m1+jp_fabm_chl3) ) + &
                        &                                                   logchl_bkginc(ji,jj,3) )                      - &
                        &                                             trn(ji,jj,1,jp_fabm_m1+jp_fabm_chl3)
                  ENDIF
                  IF ( logchl_bkginc(ji,jj,4) /= 0.0 ) THEN
                     logchl_balinc(ji,jj,1,jp_fabm_m1+jp_fabm_chl4) = 10**( LOG10( trn(ji,jj,1,jp_fabm_m1+jp_fabm_chl4) ) + &
                        &                                                   logchl_bkginc(ji,jj,4) )                      - &
                        &                                             trn(ji,jj,1,jp_fabm_m1+jp_fabm_chl4)
                  ENDIF
                  IF (k_maxchlinc > 0.0) THEN
                     chl_inc = logchl_balinc(ji,jj,1,jp_fabm_m1+jp_fabm_chl1) + &
                               logchl_balinc(ji,jj,1,jp_fabm_m1+jp_fabm_chl2) + &
                               logchl_balinc(ji,jj,1,jp_fabm_m1+jp_fabm_chl3) + &
                               logchl_balinc(ji,jj,1,jp_fabm_m1+jp_fabm_chl4)
                     IF ( ABS(chl_inc) > k_maxchlinc ) THEN
                        chl_tot = ABS(chl_inc) / k_maxchlinc
                        logchl_balinc(ji,jj,1,jp_fabm_m1+jp_fabm_chl1) = logchl_balinc(ji,jj,1,jp_fabm_m1+jp_fabm_chl1) / chl_tot
                        logchl_balinc(ji,jj,1,jp_fabm_m1+jp_fabm_chl2) = logchl_balinc(ji,jj,1,jp_fabm_m1+jp_fabm_chl2) / chl_tot
                        logchl_balinc(ji,jj,1,jp_fabm_m1+jp_fabm_chl3) = logchl_balinc(ji,jj,1,jp_fabm_m1+jp_fabm_chl3) / chl_tot
                        logchl_balinc(ji,jj,1,jp_fabm_m1+jp_fabm_chl4) = logchl_balinc(ji,jj,1,jp_fabm_m1+jp_fabm_chl4) / chl_tot
                     ENDIF
                  ENDIF
               ENDIF
            END DO
         END DO
      ELSE
         !
         ! Assimilating total Chl, so transform total from LogChl to Chl
         ! and split between PFTs according to the existing background ratios
         !
         IF ( npfts /= 1 ) THEN
            CALL ctl_stop( 'If assimilating total chlorophyll, nn_asmpfts must be 1' )
         ENDIF
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF ( ( logchl_bkginc(ji,jj,1)               /= 0.0 ) .AND. &
                  & ( trn(ji,jj,1,jp_fabm_m1+jp_fabm_chl1) >  0.0 ) .AND. &
                  & ( trn(ji,jj,1,jp_fabm_m1+jp_fabm_chl2) >  0.0 ) .AND. &
                  & ( trn(ji,jj,1,jp_fabm_m1+jp_fabm_chl3) >  0.0 ) .AND. &
                  & ( trn(ji,jj,1,jp_fabm_m1+jp_fabm_chl4) >  0.0 ) ) THEN
                  chl_tot = trn(ji,jj,1,jp_fabm_m1+jp_fabm_chl1) + trn(ji,jj,1,jp_fabm_m1+jp_fabm_chl2) + &
                     &      trn(ji,jj,1,jp_fabm_m1+jp_fabm_chl3) + trn(ji,jj,1,jp_fabm_m1+jp_fabm_chl4)
                  chl_inc = 10**( LOG10( chl_tot ) + logchl_bkginc(ji,jj,1) ) - chl_tot
                  IF (k_maxchlinc > 0.0) THEN
                     chl_inc = MAX( -1.0 * k_maxchlinc, MIN( chl_inc, k_maxchlinc ) )
                  ENDIF
                  logchl_balinc(ji,jj,1,jp_fabm_m1+jp_fabm_chl1) = chl_inc * trn(ji,jj,1,jp_fabm_m1+jp_fabm_chl1) / chl_tot
                  logchl_balinc(ji,jj,1,jp_fabm_m1+jp_fabm_chl2) = chl_inc * trn(ji,jj,1,jp_fabm_m1+jp_fabm_chl2) / chl_tot
                  logchl_balinc(ji,jj,1,jp_fabm_m1+jp_fabm_chl3) = chl_inc * trn(ji,jj,1,jp_fabm_m1+jp_fabm_chl3) / chl_tot
                  logchl_balinc(ji,jj,1,jp_fabm_m1+jp_fabm_chl4) = chl_inc * trn(ji,jj,1,jp_fabm_m1+jp_fabm_chl4) / chl_tot
               ENDIF
            END DO
         END DO
      ENDIF

      ! Propagate surface Chl1-4 incs through mixed layer
      ! First, choose mixed layer definition
      !
      SELECT CASE( mld_choice_bgc )
      CASE ( 1 )                   ! Turbocline/mixing depth [W points]
         zmld(:,:) = hmld(:,:)
      CASE ( 2 )                   ! Density criterion (0.01 kg/m^3 change from 10m) [W points]
         zmld(:,:) = hmlp(:,:)
      CASE ( 3 )                   ! Kara MLD [Interpolated]
#if defined key_karaml
         IF ( ln_kara ) THEN
            zmld(:,:) = hmld_kara(:,:)
         ELSE
            CALL ctl_stop( ' Kara mixed layer requested for LogChl assimilation,', &
               &           ' but ln_kara=.false.' )
         ENDIF
#else
         CALL ctl_stop( ' Kara mixed layer requested for LogChl assimilation,', &
            &           ' but is not defined' )
#endif
      CASE ( 4 )                   ! Temperature criterion (0.2 K change from surface) [T points]
         zmld(:,:) = hmld_tref(:,:)
      CASE ( 5 )                   ! Density criterion (0.01 kg/m^3 change from 10m) [T points]
         zmld(:,:) = hmlpt(:,:)
      END SELECT
      !
      ! Now set MLD to bottom of a level and propagate incs equally through mixed layer
      !
      DO jj = 1, jpj
         DO ji = 1, jpi
            !
            jkmax = jpk-1
            DO jk = jpk-1, 1, -1
               IF ( ( zmld(ji,jj) >  gdepw_n(ji,jj,jk)   ) .AND. &
                  & ( zmld(ji,jj) <= gdepw_n(ji,jj,jk+1) ) ) THEN
                  zmld(ji,jj) = gdepw_n(ji,jj,jk+1)
                  jkmax = jk
               ENDIF
            END DO
            !
            DO jk = 2, jkmax
               logchl_balinc(ji,jj,jk,jp_fabm_m1+jp_fabm_chl1) = logchl_balinc(ji,jj,1,jp_fabm_m1+jp_fabm_chl1)
               logchl_balinc(ji,jj,jk,jp_fabm_m1+jp_fabm_chl2) = logchl_balinc(ji,jj,1,jp_fabm_m1+jp_fabm_chl2)
               logchl_balinc(ji,jj,jk,jp_fabm_m1+jp_fabm_chl3) = logchl_balinc(ji,jj,1,jp_fabm_m1+jp_fabm_chl3)
               logchl_balinc(ji,jj,jk,jp_fabm_m1+jp_fabm_chl4) = logchl_balinc(ji,jj,1,jp_fabm_m1+jp_fabm_chl4)
            END DO
            !
         END DO
      END DO

      ! Multivariate balancing forthcoming...

   END SUBROUTINE asm_logchl_bal_ersem

#else
   !!----------------------------------------------------------------------
   !!   Default option : Empty routine
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE asm_logchl_bal_ersem( ld_logchlpftinc, npfts, mld_choice_bgc, &
      &                             k_maxchlinc, logchl_bkginc, logchl_balinc )
      LOGICAL :: ld_logchlpftinc
      INTEGER :: npfts
      INTEGER :: mld_choice_bgc
      REAL    :: k_maxchlinc
      REAL    :: logchl_bkginc(:,:,:)
      REAL    :: logchl_balinc(:,:,:,:)
      WRITE(*,*) 'asm_logchl_bal_ersem: You should not have seen this print! error?'
   END SUBROUTINE asm_logchl_bal_ersem
#endif

   !!======================================================================
END MODULE asmlogchlbal_ersem
