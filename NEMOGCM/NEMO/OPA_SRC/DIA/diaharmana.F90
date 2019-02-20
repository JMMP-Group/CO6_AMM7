MODULE harmonic_analysis
   !!======================================================================
   !!                       ***  MODULE  example  ***
   !! Ocean physics:  On line harmonic analyser
   !!                 
   !!=====================================================================

#if defined key_harm_ana
   !!----------------------------------------------------------------------
   !!   'key_harm_ana'  :                Calculate harmonic analysis
   !!----------------------------------------------------------------------
   !!   harm_ana        :
   !!   harm_ana_init   :
   !!----------------------------------------------------------------------

   USE oce             ! ocean dynamics and tracers 
   USE dom_oce         ! ocean space and time domain
   USE iom
   USE in_out_manager  ! I/O units
   USE phycst          ! physical constants
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE bdy_par         ! Unstructured boundary parameters
   USE bdy_oce         ! ocean open boundary conditions
   USE bdytides        ! tidal bdy forcing
   USE daymod          ! calendar
   USE tideini
   USE restart
   USE ioipsl, ONLY : ju2ymds    ! for calendar

   IMPLICIT NONE
   PRIVATE

   !! *  Routine accessibility
   PUBLIC harm_ana    ! routine called in step.F90 module

   !! * Module variables
   INTEGER, PARAMETER ::  nharm_max  =   jptides_max ! max number of harmonics to be analysed 
   INTEGER, PARAMETER ::  nhm_max    =   2*jptides_max+1 
   INTEGER, PARAMETER ::  nvab       = 2 ! number of 3D variables
   INTEGER            ::  nharm
   INTEGER            ::  nhm 
   INTEGER ::                 & !!! ** toto namelist (namtoto) **
      nflag  =  1                ! default value of nflag 
   REAL(wp), DIMENSION(nharm_max) ::                & 
      om_tide                     ! tidal frequencies ( rads/sec)
   REAL(wp), ALLOCATABLE,SAVE,DIMENSION(:)   ::                & 
      bzz,c,x    ! work arrays
   REAL(wp) :: cca,ssa,zm,bt
   REAL(wp) :: ccau,ccav,ssau,ssav
   REAL(wp), PUBLIC, ALLOCATABLE,DIMENSION(:) :: anau, anav, anaf   ! nodel/phase corrections used by diaharmana
!
   REAL(wp), ALLOCATABLE,SAVE,DIMENSION(:,:,:) ::   &
      bssh,bubar,bvbar        ! work array for ssh anaylsis
   REAL(wp), ALLOCATABLE,SAVE, DIMENSION(:,:,:) ::   &
      cosampz, cosampu,cosampv,      &
      sinampz, sinampu,sinampv
   REAL(WP), ALLOCATABLE,SAVE,DIMENSION(:,:) :: cc,a
!
   REAL(wp), ALLOCATABLE,SAVE,DIMENSION(:,:) ::   &
      gout,hout        ! arrays for output
   REAL(wp), ALLOCATABLE,SAVE,DIMENSION(:,:) ::   &
      huout,hvout,guout,gvout        ! arrays for output
   REAL(wp), PUBLIC ::   fjulday_startharm       !: Julian Day since start of harmonic analysis

   !! * Substitutions

   !!----------------------------------------------------------------------
   !!  OPA 9.0 , LOCEAN-IPSL (2005) 
   !! or LIM 2.0 , UCL-LOCEAN-IPSL (2005)
   !! or  TOP 1.0 , LOCEAN-IPSL (2005)
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/module_example,v 1.3 2005/03/27 18:34:47 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE harm_ana( kt )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE harm_ana  ***
      !!
      !! ** Purpose :   Harmonic analyser
      !!
      !! ** Method  :   
      !!
      !! ** Action  : - first action (share memory array/varible modified
      !!                in this routine
      !!              - second action .....
      !!              - .....
      !!
      !! References :
      !!   Give references if exist otherwise suppress these lines
      !!
      !! History :
      !!   9.0  !  03-08  (Autor Names)  Original code
      !!        !  02-08  (Author names)  brief description of modifications
      !!----------------------------------------------------------------------
      !! * Modules used
      
      !! * arguments
      INTEGER, INTENT( in  ) ::   &  
         kt                          ! describe it!!!

      !! * local declarations
      INTEGER ::   ji, jj, jk,  &        ! dummy loop arguments
                   ih,i1,i2,iv,jgrid


      !!--------------------------------------------------------------------



      IF( kt == nit000  )   CALL harm_ana_init    ! Initialization (first time-step only)

! this bit done every time step
        !nharm=ncpt_anal
        nharm=nb_harmo
        nhm=2*nb_harmo+1
       
        c(1)=1.0

        do ih=1,nb_harmo
             c(2*ih)=    cos((fjulday-fjulday_startharm)*86400._wp*om_tide(ih)  )
             c(2*ih+1)=  sin((fjulday-fjulday_startharm)*86400._wp*om_tide(ih)  )
        enddo


         do ji=1,jpi
          do jj=1,jpj
            do ih=1,nhm
            bssh (ih,ji,jj)=bssh (ih,ji,jj)+c(ih)*sshn(ji,jj)
            bubar(ih,ji,jj)=bubar(ih,ji,jj)+c(ih)*un_b(ji,jj)
            bvbar(ih,ji,jj)=bvbar(ih,ji,jj)+c(ih)*vn_b(ji,jj)
            enddo
          enddo
        enddo
!  
!        IF(lwp) WRITE(666,'(4E15.5)') adatrj*86400._wp, sshn(10,10),bssh(2,10,10),bssh(3,10,10)
!
         do i1=1,nhm
           do i2=1,nhm
             cc(i1,i2)=cc(i1,i2)+c(i1)*c(i2)
            enddo
         enddo
!
! At End of run

      IF (kt ==  nitend  ) then

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'harm_ana : harmonic analysis of tides at end of run'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~'
      CALL harm_rst_write(kt)     ! Dump out data for a restarted run


      if(ln_harm_ana_compute) then
               WRITE(numout,*) "Computing harmonics at last step"
!
        cosampu=0.0_wp
        sinampu=0.0_wp
        cosampv=0.0_wp
        sinampv=0.0_wp
        cosampz=0.0_wp
        sinampz=0.0_wp
!
     do jgrid=1,3 ! elevation, Ubar,Vbar
        do ji=1,jpi
            do jj=1,jpj
               bt=1.0
               do ih=1,nhm
                   if(jgrid .eq. 1) then
                     bzz(ih)=bssh(ih,ji,jj)
                   endif
                   if(jgrid .eq. 2) then
                     bzz(ih)=bubar(ih,ji,jj)
                   endif
                   if(jgrid .eq. 3) then
                     bzz(ih)=bvbar(ih,ji,jj)
                   endif
                   bt=bt*bzz(ih)
               enddo

               do i1=1,nhm
                   do i2=1,nhm
                       a(i1,i2)=cc(i1,i2)
                   enddo
               enddo
!	now do gaussian elimination of the system
!	a * x = b
! 	the matrix x is (a0,a1,b1,a2,b2 ...)
!	the matrix a and rhs b solved here for x
              x=0.0d0
              if(bt.ne.0.) then
                  if(lwp .and. ji == 10 .and. ji == 10) then
              endif
              call gelim(a,bzz,x,nhm)

              do ih=1,nb_harmo
                if(jgrid .eq. 1) then
                   cosampz(ih,ji,jj)=x(ih*2)
                   sinampz(ih,ji,jj)=x(ih*2+1)
                endif
                if(jgrid .eq. 2) then
                   cosampu(ih,ji,jj)=x(ih*2)
                   sinampu(ih,ji,jj)=x(ih*2+1)
                endif
                if(jgrid .eq. 3) then
                   cosampv(ih,ji,jj)=x(ih*2)
                   sinampv(ih,ji,jj)=x(ih*2+1)
                endif
              enddo
                if(jgrid .eq. 1) then
                  cosampz(0,ji,jj)=x(1)
                  sinampz(0,ji,jj)=0.0_wp
                endif
                if(jgrid .eq. 2) then
                  cosampu(0,ji,jj)=x(1)
                  sinampu(0,ji,jj)=0.0_wp
                endif
                if(jgrid .eq. 3) then
                  cosampv(0,ji,jj)=x(1)
                  sinampv(0,ji,jj)=0.0_wp
                endif
              endif
           enddo
        enddo
      enddo ! jgrid
!
!
        CALL harm_ana_out     ! output analysis (last time step)
     ELSE ! ln_harmana_compute 
         WRITE(numout,*) "Skipping Computing harmonics at last step"
     ENDIF
       
    ENDIF
   END SUBROUTINE harm_ana

   SUBROUTINE harm_ana_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE harm_ana_init  ***
      !!                   
      !! ** Purpose :   initialization of ....
      !!
      !! ** Method  :   blah blah blah ...
      !!
      !! ** input   :   Namlist namexa
      !!
      !! ** Action  :   ...  
      !!
      !! history :
      !!   9.0  !  03-08  (Autor Names)  Original code
      !!----------------------------------------------------------------------
      !! * local declarations
      INTEGER ::   ji, jj, jk, jit,ih   ! dummy loop indices

      !!----------------------------------------------------------------------
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'harm_init : initialization of harmonic analysis of tides'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~'
     ! DO ALLOCATIONS

      ALLOCATE( bubar(nb_harmo*2+1,jpi,jpj) )
      ALLOCATE( bvbar(nb_harmo*2+1,jpi,jpj) )
      ALLOCATE( bssh (nb_harmo*2+1,jpi,jpj) )

      ALLOCATE( cosampz(0:nb_harmo*2+1,jpi,jpj))
      ALLOCATE( cosampu(0:nb_harmo*2+1,jpi,jpj))
      ALLOCATE( cosampv(0:nb_harmo*2+1,jpi,jpj))

      ALLOCATE( sinampz(0:nb_harmo*2+1,jpi,jpj))
      ALLOCATE( sinampu(0:nb_harmo*2+1,jpi,jpj))
      ALLOCATE( sinampv(0:nb_harmo*2+1,jpi,jpj))

      ALLOCATE( cc(nb_harmo*2+1,nb_harmo*2+1) )

      ALLOCATE( a(nb_harmo*2+1,nb_harmo*2+1) )

      ALLOCATE( anau(nb_harmo) )
      ALLOCATE( anav(nb_harmo) )
      ALLOCATE( anaf(nb_harmo) )

      ALLOCATE( bzz(nb_harmo*2+1) )
      ALLOCATE( x(nb_harmo*2+1) )
      ALLOCATE( c(nb_harmo*2+1) )

      ALLOCATE( gout(jpi,jpj))
      ALLOCATE( hout(jpi,jpj))

      ALLOCATE( guout(jpi,jpj))
      ALLOCATE( huout(jpi,jpj))

      ALLOCATE( gvout(jpi,jpj))
      ALLOCATE( hvout(jpi,jpj))

        ! END ALLOCATE 

                 do ih=1,nb_harmo
                  om_tide(ih)= omega_tide(ih) 
                 enddo
       if(ln_harmana_read) then
               WRITE(numout,*) "Reading previous harmonic data from previous run"
               ! Need to read in  bssh bz, cc anau anav and anaf 
               call harm_rst_read  ! This reads in from the previous day
                                   ! Currrently the data in in assci format
       else
               WRITE(numout,*) "Starting harmonic analysis from Fresh "
         bssh(:,:,:)  = 0.0_wp
         bubar(:,:,:) = 0.0_wp
         bvbar(:,:,:) = 0.0_wp

         cc=0.0_wp

         anau(:) = 0.0_wp
         anav(:) = 0.0_wp
         anaf(:) = 0.0_wp
         bzz(:) = 0.0_wp
         x(:) = 0.0_wp
         c(:) = 0.0_wp

         DO ih = 1, nb_harmo
            anau(ih) = utide(ih)
            anav(ih) = v0tide(ih)
            anaf(ih) = ftide(ih)

         END DO
         fjulday_startharm=fjulday !Set this at very start and store
      endif

   END SUBROUTINE harm_ana_init
!
   SUBROUTINE gelim (a,b,x,n)
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE harm_ana  ***
      !!
      !! ** Purpose :   Guassian elimination
      !!
      !! ** Method  :   
      !!
      !! ** Action  : - first action (share memory array/varible modified
      !!                in this routine
      !!              - second action .....
      !!              - .....
      !!
      !! References :
      !!   Give references if exist otherwise suppress these lines
      !!
      !! History :
        implicit none
!
        integer  :: n
   REAL(WP) :: b(nb_harmo*2+1),a(nb_harmo*2+1,nb_harmo*2+1)
   REAL(WP) :: x(nb_harmo*2+1)
        INTEGER  :: row,col,prow,pivrow,rrow,ntemp
   REAL(WP) ::  atemp
   REAL(WP) ::  pivot
   REAL(WP) ::  m
   REAL(WP) ::  tempsol(nb_harmo*2+1)


        do row=1,n-1
           pivrow=row
           pivot=a(row,n-row+1)
           do prow=row+1,n
              if (abs(a(prow,n-row+1)).gt.abs(pivot)  ) then
                 pivot=a(prow,n-row+1)
                 pivrow=prow
              endif
           enddo
!	swap row and prow
           if ( pivrow .ne. row ) then
              atemp=b(pivrow)
              b(pivrow)=b(row)
              b(row)=atemp
              do col=1,n
                 atemp=a(pivrow,col)
                 a(pivrow,col)=a(row,col)
                 a(row,col)=atemp
              enddo
           endif

           do rrow=row+1,n
              if (a(row,row).ne.0) then
   
                 m=-a(rrow,n-row+1)/a(row,n-row+1)
                 do col=1,n
                    a(rrow,col)=m*a(row,col)+a(rrow,col)
                 enddo
                 b(rrow)=m*b(row)+b(rrow)
              endif
           enddo
        enddo
!	back substitution now

        x(1)=b(n)/a(n,1)
        do row=n-1,1,-1
           x(n-row+1)=b(row)
           do col=1,(n-row)
              x(n-row+1)=(x(n-row+1)-a(row,col)*x(col)) 
           enddo

           x(n-row+1)=(x(n-row+1)/a(row,(n-row)+1))
        enddo

        return
        END SUBROUTINE gelim

   SUBROUTINE harm_ana_out
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE harm_ana_init  ***
      !!                   
      !! ** Purpose :   initialization of ....
      !!
      !! ** Method  :   blah blah blah ...
      !!
      !! ** input   :   Namlist namexa
      !!
      !! ** Action  :   ...  
      !!
      !! history :
      !!   9.0  !  03-08  (Autor Names)  Original code
      !!----------------------------------------------------------------------
        USE dianam          ! build name of file (routine)
 
      !! * local declarations
      INTEGER ::   ji, jj, jk, jgrid,jit, ih, jjk   ! dummy loop indices
      INTEGER :: nh_T
      INTEGER :: nid_harm
      CHARACTER (len=40) ::           &
         clhstnamt, clop1, clop2          ! temporary names 
      CHARACTER (len=40) ::           &
         clhstnamu, clhstnamv   ! temporary names 
      REAL(wp) ::   &
         zsto1, zsto2, zout, zmax,            &  ! temporary scalars
         zjulian, zdt, zmdi  
         
        do jgrid=1,3
          do ih=1,nb_harmo
            hout = 0.0
            gout = 0.0
            do jj=1,nlcj
                do ji=1,nlci
                    if(jgrid .eq. 1) then ! SSH 
                       cca=cosampz(ih,ji,jj)
                       ssa=sinampz(ih,ji,jj)
                    endif
                    if(jgrid .eq. 2) then ! UBAR
                       cca=cosampu(ih,ji,jj)
                       ssa=sinampu(ih,ji,jj)
                    endif
                    if(jgrid .eq. 3) then ! VBAR
                       cca=cosampv(ih,ji,jj)
                       ssa=sinampv(ih,ji,jj)
                    endif

                    hout(ji,jj)=sqrt(cca**2+ssa**2)
        
                    if (cca.ne.0.0) then
                         gout(ji,jj)=(180.0/rpi)*atan(ssa/cca)
                    else
                        if (ssa.ne.0.0) then
                            if (ssa.gt. 0.0) then
                                 gout(ji,jj)=90.0
                            else
                                 gout(ji,jj)=270.0
                            endif
                        else
                            gout(ji,jj)=0.0
                        endif
                    endif

                    if (gout(ji,jj).lt.0) then
                        gout(ji,jj)=gout(ji,jj)+180.0
                    endif
                    if (ssa.lt.0) then
                        gout(ji,jj)=gout(ji,jj)+180.0
                    endif
                    
                    if (hout(ji,jj).ne.0) then
                        hout(ji,jj)=hout(ji,jj)/anaf(ih)
                    endif
                    if (gout(ji,jj).ne.0) then
                                                !Convert Rad to degree and take
                                                !modulus
                         gout(ji,jj)=gout(ji,jj)+MOD( (anau(ih)+anav(ih))/rad , 360.0)
                        if (gout(ji,jj).gt.360.0) then
                            gout(ji,jj)=gout(ji,jj)-360.0
                        else if (gout(ji,jj).lt.0) then
                            gout(ji,jj)=gout(ji,jj)+360.0
                        endif
                    endif
                enddo
            enddo
!NETCDF OUTPUT
            if(jgrid==1) then
             WRITE(numout,*) TRIM(Wave(ntide(ih))%cname_tide)//'amp'
             WRITE(numout,*) TRIM(Wave(ntide(ih))%cname_tide)//'phase'
             CALL iom_put( TRIM(Wave(ntide(ih))%cname_tide)//'amp', hout(:,:) )
             CALL iom_put( TRIM(Wave(ntide(ih))%cname_tide)//'phase', gout(:,:) )
            endif
            if(jgrid==2) then
             CALL iom_put( TRIM(Wave(ntide(ih))%cname_tide)//'amp_u', hout(:,:) )
             CALL iom_put( TRIM(Wave(ntide(ih))%cname_tide)//'phase_u', gout(:,:) )
            endif
            if(jgrid==3) then
             CALL iom_put( TRIM(Wave(ntide(ih))%cname_tide)//'amp_v', hout(:,:) )
             CALL iom_put( TRIM(Wave(ntide(ih))%cname_tide)//'phase_v', gout(:,:) )
            endif
          enddo
        enddo
 
!
   END SUBROUTINE harm_ana_out

   SUBROUTINE harm_rst_write(kt)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE harm_ana_init  ***
      !!                   
      !! ** Purpose :  To write out cummulated Tidal Harmomnic data to file for
      !!               restarting
      !!
      !! ** Method  :   restart files will be dated by default
      !!
      !! ** input   :   
      !!
      !! ** Action  :   ...  
      !!
      !! history :
      !!   0.0  !  01-16  (Enda O'Dea)  Original code
      !! ASSUMES  dated file for rose  , can change later to be more generic
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt     ! ocean time-step
      INTEGER             ::   iyear, imonth, iday,ih
      REAL (wp)           ::   zsec
      !!
      CHARACTER(LEN=20)   ::   clkt     ! ocean time-step define as a character
      CHARACTER(LEN=50)   ::   clname   ! ocean output restart file name
      CHARACTER(LEN=150)  ::   clpath   ! full path to ocean output restart file
      CHARACTER(LEN=250)  ::   clfinal   ! full name

      CALL ju2ymds( fjulday , iyear, imonth, iday, zsec)

      WRITE(clkt, '(i4.4,2i2.2)') iyear, imonth, iday
      ! create the file
      clname = TRIM(cexper)//"_"//TRIM(ADJUSTL(clkt))//"_restart_harm_ana"
      clpath = TRIM(cn_ocerst_outdir)
      IF( clpath(LEN_TRIM(clpath):) /= '/' ) clpath = TRIM(clpath) // '/'
      WRITE(numout,*) 'Open tidal harmonics restart file for writing: ',TRIM(clpath)//clname




!!      clhstnam=TRIM(cexper)//'.restart_harm_ana'
      ! Open a file to write the data into
         write (clfinal,'(a,''_'',i4.4)') trim(clpath)//trim(clname),nproc
         open(66,file=TRIM(clfinal))
         write(66,'(1e20.9)') cc
         !These Three which contain the most data can be moved to a regular
         !restart file

            CALL iom_rstput( kt, nitrst, numrow, 'Mean_bssh'     , bssh(1,:,:)       )
            CALL iom_rstput( kt, nitrst, numrow, 'Mean_bubar'    , bubar(1,:,:)       )
            CALL iom_rstput( kt, nitrst, numrow, 'Mean_bvbar'    , bvbar(1,:,:)       )
         do ih=1,nb_harmo
            CALL iom_rstput( kt, nitrst, numrow, TRIM(Wave(ntide(ih))%cname_tide)//'bssh_cos'    , bssh (ih*2  , : , : )     )
            CALL iom_rstput( kt, nitrst, numrow, TRIM(Wave(ntide(ih))%cname_tide)//'bssh_sin'    , bssh (ih*2+1, : , : )     )

            CALL iom_rstput( kt, nitrst, numrow, TRIM(Wave(ntide(ih))%cname_tide)//'bubar_cos'   , bubar(ih*2  , : , : )    )
            CALL iom_rstput( kt, nitrst, numrow, TRIM(Wave(ntide(ih))%cname_tide)//'bubar_sin'   , bubar(ih*2+1, : , : )    )

            CALL iom_rstput( kt, nitrst, numrow, TRIM(Wave(ntide(ih))%cname_tide)//'bvbar_cos'   , bvbar(ih*2  , : , : )    )
            CALL iom_rstput( kt, nitrst, numrow, TRIM(Wave(ntide(ih))%cname_tide)//'bvbar_sin'   , bvbar(ih*2+1, : , : )    )
         enddo

         write(66,'(1e20.9)') anau
         write(66,'(1e20.9)') anav
         write(66,'(1e20.9)') anaf
         write(66,'(1e25.20)') fjulday_startharm
         close(66)
 
   END SUBROUTINE harm_rst_write

   SUBROUTINE harm_rst_read
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE harm_ana_init  ***
      !!                   
      !! ** Purpose :  To read in  cummulated Tidal Harmomnic data to file for
      !!               restarting
      !!
      !! ** Method  :   
      !!
      !! ** input   :   
      !!
      !! ** Action  :   ...  
      !!
      !! history :
      !!   0.0  !  01-16  (Enda O'Dea)  Original code
      !! ASSUMES  dated file for rose  , can change later to be more generic
      !!----------------------------------------------------------------------
      CHARACTER(LEN=50)   ::   clname   ! ocean output restart file name
      CHARACTER(LEN=150)  ::   clpath   ! full path to ocean output restart file
      CHARACTER(LEN=250)  ::   clfinal   ! full name
      INTEGER             ::   ih

      ! create the file
      clname = "restart_harm_ana"
      clpath = TRIM(cn_ocerst_indir)
      IF( clpath(LEN_TRIM(clpath):) /= '/' ) clpath = TRIM(clpath) // '/'
      WRITE(numout,*) 'Open tidal harmonics restart file for reading: ',TRIM(clpath)//clname




!!      clhstnam=TRIM(cexper)//'.restart_harm_ana'
      ! Open a file to read the data ifrom
         write (clfinal,'(a,''_'',i4.4)') trim(clpath)//trim(clname),nproc
         open(66,file=TRIM(clfinal),status='old')
         read(66,'(1e20.9)') cc
!2D regular fields can be read from normal restart this saves space and handy to
!view in netcdf format also.
         CALL iom_get( numror,jpdom_autoglo, 'Mean_bssh'     , bssh(1,:,:)       )
         CALL iom_get( numror,jpdom_autoglo, 'Mean_bubar'    , bubar(1,:,:)       )
         CALL iom_get( numror,jpdom_autoglo, 'Mean_bvbar'    , bvbar(1,:,:)       )
          
         do ih=1,nb_harmo

            CALL iom_get(   numror,jpdom_autoglo, TRIM(Wave(ntide(ih))%cname_tide)//'bssh_cos'    , bssh (ih*2  , : , : )     )
            CALL iom_get(   numror,jpdom_autoglo, TRIM(Wave(ntide(ih))%cname_tide)//'bssh_sin'    , bssh (ih*2+1, : , : )     )

            CALL iom_get(   numror,jpdom_autoglo, TRIM(Wave(ntide(ih))%cname_tide)//'bubar_cos'   , bubar(ih*2  , : , : )    )
            CALL iom_get(   numror,jpdom_autoglo, TRIM(Wave(ntide(ih))%cname_tide)//'bubar_sin'   , bubar(ih*2+1, : , : )    )

            CALL iom_get(   numror,jpdom_autoglo, TRIM(Wave(ntide(ih))%cname_tide)//'bvbar_cos'   , bvbar(ih*2  , : , : )    )
            CALL iom_get(   numror,jpdom_autoglo, TRIM(Wave(ntide(ih))%cname_tide)//'bvbar_sin'   , bvbar(ih*2+1, : , : )    )
         enddo
         read(66,'(1e20.9)') anau
         read(66,'(1e20.9)') anav
         read(66,'(1e20.9)') anaf
         read(66,'(1e25.20)') fjulday_startharm
      WRITE(numout,*) 'Checking anaf is correct'
      WRITE(numout,*) anaf
      WRITE(numout,*) '---end of anaf checking----'

         close(66)

   END SUBROUTINE harm_rst_read

   !!======================================================================
#else
!!---------------------------------------------------------------------------------
!!   Dummy module                                   NO harmonic Analysis
!!---------------------------------------------------------------------------------
        CONTAINS
           SUBROUTINE harm_rst_write(kt)     ! Dummy routine
           END SUBROUTINE harm_rst_write
           SUBROUTINE harm_rst_read    ! Dummy routine
           END SUBROUTINE harm_rst_read
           SUBROUTINE harm_ana_out      ! Dummy routine
           END SUBROUTINE harm_ana_out
           SUBROUTINE harm_ana_init
           END SUBROUTINE harm_ana_init
           SUBROUTINE harm_ana( kt )
           END SUBROUTINE harm_ana
           SUBROUTINE gelim (a,b,x,n)
           END SUBROUTINE gelim
           
#endif

END MODULE harmonic_analysis
