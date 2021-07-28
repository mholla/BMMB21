************************************************************************
!  
! A umat for manuscript: 'Numerical investigation of biomechanically-coulped
! growth in cortical folding'
! 
! 
! Shuolun Wang, August 2019, Implemented in Abaqus 6.19
! 
************************************************************************
!     State Variables
!     --------------------------------------------------------------
!      statev(1) = growth parameter 
!      statev(2) = coordx     
!      statev(3) = coordy 
!      statev(4) = coordz 
!      statev(5) = curvature
!     Material Properties Vector
!     --------------------------------------------------------------
!      !!!!white matters!!! 
!      mu        = props(1)
!      lambda    = props(2)
!      Gaxn      = props(3)
!      lam0      = props(4)
!      !!!!!gray matter!!!! 
!      mu_g      = props(5)
!      lambda_g  = props(6)
!      Gctx      = props(7)
!      Aconst    = props(8)
!      charlength = props(9)
!**********************************************************************
      module GlobalStorage
        

          ! set total number of gray matter elements 
          ! max - min + 1
          parameter(NElements=150) 
          parameter(NElements2=150)
          ! set the element number where gray matter starts from 
          ! min - 1
          parameter(BaseElement=0)     
          parameter(BaseElement2=150)        

          real*8 StorageOld(NElements,2)
          real*8 StorageOld2(NElements2,2)



          real*8 curvOld(NElements,1)
          real*8 curvOld2(NElements2,1)

          real*8 inicoord(2327,3)
 
      end module    
***********************************************************************
C************************************************************************
C subroutine for calculating the curvature 
C
      SUBROUTINE GetCurvature(coordAx,coordBx,coordCx,
     +                        coordAy,coordBy,coordCy, 
     +                        curvature)

          real*8 coordAx,coordBx,coordCx
          real*8 coordAy,coordBy,coordCy 
          real*8 curvature,x0,y0,x1,y1,x2,y2
          real*8 nx,ny,Qmat(2,2)
          real*8 globalarray1(2,1),globalarray2(2,1)
          real*8 localarray1(2,1),localarray2(2,1)
          real*8 x_p1,x_p2,y_p1,y_p2,a

          ! read the global coordinate of each point in space (2D)
          x0 = coordAx
          y0 = coordAy

          x1 = coordBx
          y1 = coordBy

          x2 = coordCx
          y2 = coordCy


          ! calculate the local outward surface normal 
          ! in the global coordinates
          !
          nx = -(y2 - y1)/sqrt((y2 - y1)**2.0 + (x2 - x1)**2.0)
          ny = (x2 - x1)/sqrt((y2 - y1)**2.0 + (x2 - x1)**2.0)

          ! calculate the rotation matrix
          !
          Qmat(1,1) = ny
          Qmat(1,2) = -nx
          Qmat(2,1) = nx
          Qmat(2,2) = ny


          ! Calculate the coordinates of the nearest neighbors
          ! in the primed (local) coordinate system
          !
          globalarray1(1,1) = x1 - x0
          globalarray1(2,1) = y1 - y0
          
          globalarray2(1,1) = x2 - x0
          globalarray2(2,1) = y2 - y0
          
          localarray1 = matmul(Qmat,globalarray1)
          localarray2 = matmul(Qmat,globalarray2)


          x_p1 = localarray1(1,1)
          y_p1 = localarray1(2,1)

          x_p2 = localarray2(1,1)
          y_p2 = localarray2(2,1)

          ! obtain the curvature by fitting parabola to three points 
          ! 
          a = ((x_p1**2.0) * y_p1 + (x_p2**2.0) * y_p2) / (x_p1**4.0 + x_p2**4.0)
          curvature = a

 
      END SUBROUTINE 
***********************************************************************
      subroutine vumat (
c Read only -
     +     jblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     +     stepTime, totalTime, dt, cmname, coordMp, charLength,
     +     props, density, strainInc, relSpinInc,
     +     tempOld, stretchOld, defgradOld, fieldOld,
     +     stressOld, stateOld, enerInternOld, enerInelasOld,
     +     tempNew, stretchNew, defgradNew, fieldNew,
c Write only -
     +     stressNew, stateNew, enerInternNew, enerInelasNew )
c
      include 'vaba_param.inc'
c
      dimension jblock(*), props(nprops),density(*), coordMp(*),
     +     charLength(*), strainInc(*),
     +     relSpinInc(*), tempOld(*),
     +     stretchOld(*),
     +     defgradOld(*),
     +     fieldOld(*), stressOld(*),
     +     stateOld(*), enerInternOld(*),
     +     enerInelasOld(*), tempNew(*),
     +     stretchNew(*),
     +     defgradNew(*),
     +     fieldNew(*),
     +     stressNew(*), stateNew(*),
     +     enerInternNew(*), enerInelasNew(*)
c
      character*80 cmname
      character*256 WHIT,GRAY

      parameter (     
     +     i_umt_nblock = 1,
     +     i_umt_npt    = 2,
     +     i_umt_layer  = 3,
     +     i_umt_kspt   = 4,
     +     i_umt_noel   = 5 )


      !--------------------------------------------------------
      ! 
      ! call particular user material to perform the analysis 
      ! 
      IF (CMNAME(1:4) .EQ. 'WHIT') THEN

      !
      ! this is white matter 
      !
      call  vumatXtrArg_white (jblock(i_umt_nblock),
     +     ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     +     stepTime, totalTime, dt, cmname, coordMp, charLength,
     +     props, density, strainInc, relSpinInc,
     +     tempOld, stretchOld, defgradOld, fieldOld,
     +     stressOld, stateOld, enerInternOld, enerInelasOld,
     +     tempNew, stretchNew, defgradNew, fieldNew,
     +     stressNew, stateNew, enerInternNew, enerInelasNew,
     +     jblock(i_umt_noel), jblock(i_umt_npt),
     +     jblock(i_umt_layer), jblock(i_umt_kspt))
      !
      !
      ELSE IF(CMNAME(1:7) .EQ. 'GRAYONE') THEN
      !
      ! this is gray matter 
      !
      call  vumatXtrArg_grayone (jblock(i_umt_nblock),
     +     ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     +     stepTime, totalTime, dt, cmname, coordMp, charLength,
     +     props, density, strainInc, relSpinInc,
     +     tempOld, stretchOld, defgradOld, fieldOld,
     +     stressOld, stateOld, enerInternOld, enerInelasOld,
     +     tempNew, stretchNew, defgradNew, fieldNew,
     +     stressNew, stateNew, enerInternNew, enerInelasNew,
     +     jblock(i_umt_noel), jblock(i_umt_npt),
     +     jblock(i_umt_layer), jblock(i_umt_kspt))
      !
      !
      !
      !
      ELSE IF(CMNAME(1:7) .EQ. 'GRAYTWO') THEN
      !
      ! this is gray matter 
      !
      call  vumatXtrArg_graytwo (jblock(i_umt_nblock),
     +     ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     +     stepTime, totalTime, dt, cmname, coordMp, charLength,
     +     props, density, strainInc, relSpinInc,
     +     tempOld, stretchOld, defgradOld, fieldOld,
     +     stressOld, stateOld, enerInternOld, enerInelasOld,
     +     tempNew, stretchNew, defgradNew, fieldNew,
     +     stressNew, stateNew, enerInternNew, enerInelasNew,
     +     jblock(i_umt_noel), jblock(i_umt_npt),
     +     jblock(i_umt_layer), jblock(i_umt_kspt))
      !
      !
      Endif

     
      end subroutine vumat
***********************************************************************
      subroutine vumatXtrArg_white (
c Read only -
     +     nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     +     stepTime, totalTime, dt, cmname, coordMp, charLength,
     +     props, density, strainInc, relSpinInc,
     +     tempOld, stretchOld, defgradOld, fieldOld,
     +     stressOld, stateOld, enerInternOld, enerInelasOld,
     +     tempNew, stretchNew, defgradNew, fieldNew,
c Write only -
     +     stressNew, stateNew, enerInternNew, enerInelasNew,
c Read only extra arguments -
     +     nElement, nMatPoint, nLayer, nSecPoint )

      use GlobalStorage


c$$$        implicit none ! This is used during compilation testing to make
      include 'vaba_param.inc'

c$$$!        When implicit none is used during compilation testing, all the following
c$$$!        variables need to be defined.
c$$$        integer nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal
c$$$        real*8  stepTime, totalTime, dt, coordMp, charLength, props,
c$$$     +          density, strainInc, relSpinInc, tempOld, stretchOld,
c$$$     +          defgradOld, fieldOld, stressOld, stateOld, 
c$$$     +          enerInternOld, enerInelasOld, tempNew, stretchNew,
c$$$     +          defgradNew, fieldNew, stressNew, stateNew, 
c$$$     +          enerInternNew, enerInelasNew

      dimension props(nprops), density(nblock), coordMp(nblock,*),
     +     charLength(nblock), strainInc(nblock,ndir+nshr),
     +     relSpinInc(nblock,nshr), tempOld(nblock),
     +     stretchOld(nblock,ndir+nshr),
     +     defgradOld(nblock,ndir+nshr+nshr),
     +     fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr),
     +     stateOld(nblock,nstatev), enerInternOld(nblock),
     +     enerInelasOld(nblock), tempNew(nblock),
     +     stretchNew(nblock,ndir+nshr),
     +     defgradNew(nblock,ndir+nshr+nshr),
     +     fieldNew(nblock,nfieldv),
     +     stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     +     enerInternNew(nblock), enerInelasNew(nblock)

c
c Documentation of extra arguments:
c  nElement: Array of internal element numbers
      dimension nElement(nblock)
c  nMatPoint: Integration point number
c  nLayer   : Layer number for composite shells and layered solids
c  nSecPoint: Section point number within the current layer
c
      character*80 cmname

      integer i,j,l,i1,j1,ii,jj,kk,ll,km,ifail

      real*8 Iden(3,3),F_t(3,3),F_tau(3,3),U_t(3,3),U_tau(3,3),Fp_t(3,3)
      real*8 Fp_tau(3,3),Me_t(3,3),Me_tau(3,3),nuP_t,nuP_tau,Y_t,Y_tau
      real*8 gBarP_t,gBarP_tau,T_tau(3,3),R_tau(3,3),U_inv(3,3),detF
      real*8 Fp_inv(3,3),Ee_tau(3,3),Re_tau(3,3),Ue_tau(3,3),Fe_tau(3,3)
      real*8 pnu0,damage_t,damage_tau,mag_Dp_tau,pwrinct,stress_power
      real*8 nu1,nu3,nu5
      real*8 lamg_t,lamg_tau,matProps(nprops)
      real*8 E_tau,E_t
      real*8 coordx,coordy,coordz



      ! Parameters
      !
      real*8 zero,one,two,three,half,third,four,Pi,two_third
      parameter(zero=0.d0,one=1.d0,two=2.d0,three=3.d0,half=0.5d0,
     +     third=1.d0/3.d0,two_third=2.d0/3.d0,four=4.d0,Pi=3.1415926d0)



      ! Pour initial coordinates into the global variable matrix 
      !
      if (totaltime.lt. 0.001) then
          do km=1,nblock
             inicoord(nElement(km),1) = coordMp(km,1)
             inicoord(nElement(km),2) = coordMp(km,2)
             inicoord(nElement(km),3) = coordMp(km,3)
          enddo
      end if 

c      open(unit=80, file='C:\Abaqus_Working_Dir\FeFp\aaMSGS')

      ! Identity matrix for later use.
      ! 
      call onem(Iden)


      ! 
      ! START LOOP OVER MATERIAL POINTS:
      ! 
      do km=1,nblock
      
           
         ! Copy old and new deformation gradients
         !
         F_t(1,1) = defgradOld(km,1)
         F_t(2,2) = defgradOld(km,2)
         F_t(3,3) = defgradOld(km,3)
         F_t(1,2) = defgradOld(km,4)
         F_tau(1,1) = defgradNew(km,1)
         F_tau(2,2) = defgradNew(km,2)
         F_tau(3,3) = defgradNew(km,3)
         F_tau(1,2) = defgradNew(km,4)
         U_tau(1,1) = stretchNew(km,1)
         U_tau(2,2) = stretchNew(km,2)
         U_tau(3,3) = stretchNew(km,3)
         U_tau(1,2) = stretchNew(km,4)
         if(nshr .lt. 2) then
            ! 2D case
            F_t(2,1) = defgradOld(km,5)
            F_t(1,3) = zero
            F_t(2,3) = zero
            F_t(3,1) = zero
            F_t(3,2) = zero
            F_tau(2,1) = defgradNew(km,5)
            F_tau(1,3) = zero
            F_tau(2,3) = zero
            F_tau(3,1) = zero
            F_tau(3,2) = zero
            U_tau(2,1) = U_tau(1,2)
            U_tau(1,3) = zero
            U_tau(2,3) = zero
            U_tau(3,1) = zero
            U_tau(3,2) = zero
         else
            ! 3D case
            F_t(2,3) = defgradOld(km,5)
            F_t(3,1) = defgradOld(km,6)
            F_t(2,1) = defgradOld(km,7)
            F_t(3,2) = defgradOld(km,8)
            F_t(1,3) = defgradOld(km,9)
            F_tau(2,3) = defgradNew(km,5)
            F_tau(3,1) = defgradNew(km,6)
            F_tau(2,1) = defgradNew(km,7)
            F_tau(3,2) = defgradNew(km,8)
            F_tau(1,3) = defgradNew(km,9)
            U_tau(2,3) = stretchNew(km,5)
            U_tau(3,1) = stretchNew(km,6)
            U_tau(2,1) = U_tau(1,2)
            U_tau(3,2) = U_tau(2,3)
            U_tau(1,3) = U_tau(3,1)
         end if


         if((totalTime.eq.zero).and.(stepTime.eq.zero)) then
            ! Dummy step, initalize state variables
            
            stateOld(km,1)   = one ! growth parameter at t=0
         endif

         ! Read old state variables
         
         lamg_t = stateOld(km,1) ! growth parameter at time t




          coordx = inicoord(nElement(km),1)
          coordy = inicoord(nElement(km),2)
          coordz = inicoord(nElement(km),3)
c           coordx = 1.0
c           coordy = 2.0
c           coordz = 3.0

         !---------------------------------------------------------------
         ! Perform the time integration and compute the 
         !  constitutive response based on the material model.
         
         matProps = props

         if((totalTime.eq.zero).and.(stepTime.eq.zero)) then
            !
            ! dummy step, call elastic response, note dt=-1.0 is sent
            !  into the integ subroutine
            !
            call integ_white(matProps,nprops,F_tau,-1.0,T_tau,lamg_t,lamg_tau,
     +                       coord,coordy,coordz)

         else
            !
            ! Perform explicit time integration procedure
            !
            call integ_white(matProps,nprops,F_tau,dt,T_tau,lamg_t,lamg_tau,
     +                       coord,coordy,coordz)

         endif
         !---------------------------------------------------------------


         ! Update state variables
         !
         stateNew(km,1) = lamg_tau ! growth parameter at time tau

         stateNew(km,2) = coordx ! 
         stateNew(km,3) = coordy ! 
         stateNew(km,4) = coordz ! 


         ! ABAQUS/Explicit uses stress measure (transpose(R) T R)
         !
         call m3inv(U_tau,U_inv)
         R_tau = matmul(F_tau,U_inv)
         T_tau = matmul(transpose(R_tau),matmul(T_tau,R_tau))

         do i=1,ndir
            stressNew(km,i) = T_tau(i,i)
         end do
         if(nshr.ne.0) then
            stressNew(km,ndir+1) = T_tau(1,2)
            if(nshr.ne.1) then
               stressNew(km, ndir+2) = T_tau(2,3)
               if(nshr.ne.2) then
                  stressNew(km,ndir+3) = T_tau(1,3)
               endif
            endif
         endif


         ! Update the specific internal energy
         !
         stress_power = 0.d0
         do i = 1,ndir
            stress_power = stress_power +
     +           0.5*((StressOld(km,i)+StressNew(km,i))*
     +           StrainInc(km,i))
         enddo
         
         select case (nshr)
         case(1)
            stress_power = stress_power + 
     +           0.5*((StressOld(km,ndir+1)+StressNew(km,ndir+1))*
     +           StrainInc(km,ndir+1))
         case(3)
            stress_power = stress_power + 
     +           0.5*(((StressOld(km,ndir+1) + StressNew(km,ndir+1))*
     +           StrainInc(km,ndir+1)) +
     +           ((StressOld(km,ndir+2)+ StressNew(km,ndir+2)) *
     +           StrainInc(km,ndir+2))+
     +           ((StressOld(km,ndir+3) + StressNew(km,ndir+3))*
     +           StrainInc(km,ndir+3)))
         end select
           
         enerInternNew(km) = enerInternOld(km) + 
     +        stress_power/density(km)
           
         enerInelasNew(km) = enerInelasOld(km) + 
     +        pwrinct/density(km)
           
           
      enddo ! end loop over material points

      end subroutine vumatXtrArg_white
***********************************************************************
      subroutine vumatXtrArg_grayone (
c Read only -
     +     nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     +     stepTime, totalTime, dt, cmname, coordMp, charLength,
     +     props, density, strainInc, relSpinInc,
     +     tempOld, stretchOld, defgradOld, fieldOld,
     +     stressOld, stateOld, enerInternOld, enerInelasOld,
     +     tempNew, stretchNew, defgradNew, fieldNew,
c Write only -
     +     stressNew, stateNew, enerInternNew, enerInelasNew,
c Read only extra arguments -
     +     nElement, nMatPoint, nLayer, nSecPoint )

c$$$        implicit none ! This is used during compilation testing to make
      use GlobalStorage
      include 'vaba_param.inc'

c$$$!        When implicit none is used during compilation testing, all the following
c$$$!        variables need to be defined.
c$$$        integer nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal
c$$$        real*8  stepTime, totalTime, dt, coordMp, charLength, props,
c$$$     +          density, strainInc, relSpinInc, tempOld, stretchOld,
c$$$     +          defgradOld, fieldOld, stressOld, stateOld, 
c$$$     +          enerInternOld, enerInelasOld, tempNew, stretchNew,
c$$$     +          defgradNew, fieldNew, stressNew, stateNew, 
c$$$     +          enerInternNew, enerInelasNew

      dimension props(nprops), density(nblock), coordMp(nblock,*),
     +     charLength(nblock), strainInc(nblock,ndir+nshr),
     +     relSpinInc(nblock,nshr), tempOld(nblock),
     +     stretchOld(nblock,ndir+nshr),
     +     defgradOld(nblock,ndir+nshr+nshr),
     +     fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr),
     +     stateOld(nblock,nstatev), enerInternOld(nblock),
     +     enerInelasOld(nblock), tempNew(nblock),
     +     stretchNew(nblock,ndir+nshr),
     +     defgradNew(nblock,ndir+nshr+nshr),
     +     fieldNew(nblock,nfieldv),
     +     stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     +     enerInternNew(nblock), enerInelasNew(nblock)
c
c Documentation of extra arguments:
c  nElement: Array of internal element numbers
      dimension nElement(nblock)
c  nMatPoint: Integration point number
c  nLayer   : Layer number for composite shells and layered solids
c  nSecPoint: Section point number within the current layer
c




      character*80 cmname

      integer i,j,l,i1,j1,ii,jj,kk,ll,km,ifail

      real*8 Iden(3,3),F_t(3,3),F_tau(3,3),U_t(3,3),U_tau(3,3),Fp_t(3,3)
      real*8 Fp_tau(3,3),Me_t(3,3),Me_tau(3,3),nuP_t,nuP_tau,Y_t,Y_tau
      real*8 gBarP_t,gBarP_tau,T_tau(3,3),R_tau(3,3),U_inv(3,3),detF
      real*8 Fp_inv(3,3),Ee_tau(3,3),Re_tau(3,3),Ue_tau(3,3),Fe_tau(3,3)
      real*8 pnu0,damage_t,damage_tau,mag_Dp_tau,pwrinct,stress_power
      real*8 nu1,nu3,nu5
      real*8 lamg_t,lamg_tau,matProps(nprops)
      real*8 E_tau,E_t,thetag_t,thetag_tau
      real*8 coordx,coordy,coordz
      real*8 curvature

      ! Parameters
      !
      real*8 zero,one,two,three,half,third,four,Pi,two_third
      parameter(zero=0.d0,one=1.d0,two=2.d0,three=3.d0,half=0.5d0,
     +     third=1.d0/3.d0,two_third=2.d0/3.d0,four=4.d0,Pi=3.1415926d0)

c      open(unit=80, file='C:\Abaqus_Working_Dir\FeFp\aaMSGS')



      ! pour initial coordinates into the global variable
      if (totaltime.lt. 0.001) then
          do km=1,nblock
             inicoord(nElement(km),1) = coordMp(km,1)
             inicoord(nElement(km),2) = coordMp(km,2)
             inicoord(nElement(km),3) = coordMp(km,3)
          enddo
      end if

      ! Identity matrix for later use.
      !
      call onem(Iden)


      !
      ! START LOOP OVER MATERIAL POINTS:
      !
      do km=1,nblock

           
         ! Copy old and new deformation gradients
         !
         F_t(1,1) = defgradOld(km,1)
         F_t(2,2) = defgradOld(km,2)
         F_t(3,3) = defgradOld(km,3)
         F_t(1,2) = defgradOld(km,4)
         F_tau(1,1) = defgradNew(km,1)
         F_tau(2,2) = defgradNew(km,2)
         F_tau(3,3) = defgradNew(km,3)
         F_tau(1,2) = defgradNew(km,4)
         U_tau(1,1) = stretchNew(km,1)
         U_tau(2,2) = stretchNew(km,2)
         U_tau(3,3) = stretchNew(km,3)
         U_tau(1,2) = stretchNew(km,4)
         if(nshr .lt. 2) then
            ! 2D case
            F_t(2,1) = defgradOld(km,5)
            F_t(1,3) = zero
            F_t(2,3) = zero
            F_t(3,1) = zero
            F_t(3,2) = zero
            F_tau(2,1) = defgradNew(km,5)
            F_tau(1,3) = zero
            F_tau(2,3) = zero
            F_tau(3,1) = zero
            F_tau(3,2) = zero
            U_tau(2,1) = U_tau(1,2)
            U_tau(1,3) = zero
            U_tau(2,3) = zero
            U_tau(3,1) = zero
            U_tau(3,2) = zero
         else
            ! 3D case
            F_t(2,3) = defgradOld(km,5)
            F_t(3,1) = defgradOld(km,6)
            F_t(2,1) = defgradOld(km,7)
            F_t(3,2) = defgradOld(km,8)
            F_t(1,3) = defgradOld(km,9)
            F_tau(2,3) = defgradNew(km,5)
            F_tau(3,1) = defgradNew(km,6)
            F_tau(2,1) = defgradNew(km,7)
            F_tau(3,2) = defgradNew(km,8)
            F_tau(1,3) = defgradNew(km,9)
            U_tau(2,3) = stretchNew(km,5)
            U_tau(3,1) = stretchNew(km,6)
            U_tau(2,1) = U_tau(1,2)
            U_tau(3,2) = U_tau(2,3)
            U_tau(1,3) = U_tau(3,1)
         end if


         if((totalTime.eq.zero).and.(stepTime.eq.zero)) then
            ! Dummy step, initalize state variables
            
            stateOld(km,1)   = one ! growth parameter at t=0
         endif

         ! Read old state variables
         
         thetag_t = stateOld(km,1) ! growth parameter at time t




         ! reads in the original coordinates 
          coordx = inicoord(nElement(km),1)
          coordy = inicoord(nElement(km),2)
          coordz = inicoord(nElement(km),3)




         ! store all the coordinates 
         StorageOld(nElement(km)-BaseElement,1) = coordMp(km,1)
         StorageOld(nElement(km)-BaseElement,2) = coordMp(km,2)

         ! here is trick I play here, so for full size simulation,
         ! the curvature will be flipped in (-x,y) and (x-y) region,
         ! here I just mirror it...
c         if (StorageOld(nElement(km)-BaseElement,1) * StorageOld(nElement(km)-BaseElement,2).lt.0.0) then
c             StorageOld(nElement(km)-BaseElement,2) = -StorageOld(nElement(km)-BaseElement,2)
c         endif



      do i=1,nElement(km)-BaseElement

         if (i.eq.1) then 
            
               call GetCurvature(StorageOld(2,1),StorageOld(1,1),StorageOld(3,1),
     +                        StorageOld(2,2),StorageOld(1,2),StorageOld(3,2), 
     +                        curvOld(1,1))

         elseif(i.lt.NElements) then


               call GetCurvature(StorageOld(i,1),StorageOld(i-1,1),StorageOld(i+1,1),
     +                        StorageOld(i,2),StorageOld(i-1,2),StorageOld(i+1,2), 
     +                        curvOld(i,1))

         elseif(i.eq.NElements) then

               call GetCurvature(StorageOld(NElements - BaseElement -1,1),StorageOld(NElements - BaseElement -2,1),
     +                           StorageOld(NElements - BaseElement,1),StorageOld(NElements - BaseElement -1,2),
     +                           StorageOld(NElements - BaseElement -2,2),StorageOld(NElements - BaseElement,2), 
     +                           curvOld(NElements - BaseElement,1))
         endif


      enddo


       curvature = curvOld(nElement(km)-BaseElement,1)


         !---------------------------------------------------------------
         ! Perform the time integration and compute the 
         !  constitutive response based on the material model.
         
         matProps = props

         if((totalTime.eq.zero).and.(stepTime.eq.zero)) then
            !
            ! dummy step, call elastic response, note dt=-1.0 is sent
            !  into the integ subroutine
            !
            call integ_gray(matProps,nprops,F_tau,-1.0,T_tau,thetag_t,thetag_tau,
     +                     coordx,coordy,coordz,curvature)

         else
            !
            ! Perform explicit time integration procedure
            !
            call integ_gray(matProps,nprops,F_tau,dt,T_tau,thetag_t,thetag_tau,
     +                     coordx,coordy,coordz,curvature)

         endif
         !---------------------------------------------------------------
         

         ! Update state variables
         !
         stateNew(km,1) = thetag_tau ! growth parameter at time tau

         stateNew(km,2) = coordx
         stateNew(km,3) = coordy
         stateNew(km,4) = coordz

         stateNew(km,5) = curvature

         ! ABAQUS/Explicit uses stress measure (transpose(R) T R)
         !
         call m3inv(U_tau,U_inv)
         R_tau = matmul(F_tau,U_inv)
         T_tau = matmul(transpose(R_tau),matmul(T_tau,R_tau))

         do i=1,ndir
            stressNew(km,i) = T_tau(i,i)
         end do
         if(nshr.ne.0) then
            stressNew(km,ndir+1) = T_tau(1,2)
            if(nshr.ne.1) then
               stressNew(km, ndir+2) = T_tau(2,3)
               if(nshr.ne.2) then
                  stressNew(km,ndir+3) = T_tau(1,3)
               endif
            endif
         endif


         ! Update the specific internal energy
         !
         stress_power = 0.d0
         do i = 1,ndir
            stress_power = stress_power +
     +           0.5*((StressOld(km,i)+StressNew(km,i))*
     +           StrainInc(km,i))
         enddo
         
         select case (nshr)
         case(1)
            stress_power = stress_power + 
     +           0.5*((StressOld(km,ndir+1)+StressNew(km,ndir+1))*
     +           StrainInc(km,ndir+1))
         case(3)
            stress_power = stress_power + 
     +           0.5*(((StressOld(km,ndir+1) + StressNew(km,ndir+1))*
     +           StrainInc(km,ndir+1)) +
     +           ((StressOld(km,ndir+2)+ StressNew(km,ndir+2)) *
     +           StrainInc(km,ndir+2))+
     +           ((StressOld(km,ndir+3) + StressNew(km,ndir+3))*
     +           StrainInc(km,ndir+3)))
         end select
           
         enerInternNew(km) = enerInternOld(km) + 
     +        stress_power/density(km)
           
         enerInelasNew(km) = enerInelasOld(km) + 
     +        pwrinct/density(km)
           
           
      enddo ! end loop over material points

      end subroutine vumatXtrArg_grayone
***********************************************************************
      subroutine vumatXtrArg_graytwo (
c Read only -
     +     nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     +     stepTime, totalTime, dt, cmname, coordMp, charLength,
     +     props, density, strainInc, relSpinInc,
     +     tempOld, stretchOld, defgradOld, fieldOld,
     +     stressOld, stateOld, enerInternOld, enerInelasOld,
     +     tempNew, stretchNew, defgradNew, fieldNew,
c Write only -
     +     stressNew, stateNew, enerInternNew, enerInelasNew,
c Read only extra arguments -
     +     nElement, nMatPoint, nLayer, nSecPoint )

c$$$        implicit none ! This is used during compilation testing to make
      use GlobalStorage
      include 'vaba_param.inc'

c$$$!        When implicit none is used during compilation testing, all the following
c$$$!        variables need to be defined.
c$$$        integer nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal
c$$$        real*8  stepTime, totalTime, dt, coordMp, charLength, props,
c$$$     +          density, strainInc, relSpinInc, tempOld, stretchOld,
c$$$     +          defgradOld, fieldOld, stressOld, stateOld, 
c$$$     +          enerInternOld, enerInelasOld, tempNew, stretchNew,
c$$$     +          defgradNew, fieldNew, stressNew, stateNew, 
c$$$     +          enerInternNew, enerInelasNew

      dimension props(nprops), density(nblock), coordMp(nblock,*),
     +     charLength(nblock), strainInc(nblock,ndir+nshr),
     +     relSpinInc(nblock,nshr), tempOld(nblock),
     +     stretchOld(nblock,ndir+nshr),
     +     defgradOld(nblock,ndir+nshr+nshr),
     +     fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr),
     +     stateOld(nblock,nstatev), enerInternOld(nblock),
     +     enerInelasOld(nblock), tempNew(nblock),
     +     stretchNew(nblock,ndir+nshr),
     +     defgradNew(nblock,ndir+nshr+nshr),
     +     fieldNew(nblock,nfieldv),
     +     stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     +     enerInternNew(nblock), enerInelasNew(nblock)
c
c Documentation of extra arguments:
c  nElement: Array of internal element numbers
      dimension nElement(nblock)
c  nMatPoint: Integration point number
c  nLayer   : Layer number for composite shells and layered solids
c  nSecPoint: Section point number within the current layer
c




      character*80 cmname

      integer i,j,l,i1,j1,ii,jj,kk,ll,km,ifail

      real*8 Iden(3,3),F_t(3,3),F_tau(3,3),U_t(3,3),U_tau(3,3),Fp_t(3,3)
      real*8 Fp_tau(3,3),Me_t(3,3),Me_tau(3,3),nuP_t,nuP_tau,Y_t,Y_tau
      real*8 gBarP_t,gBarP_tau,T_tau(3,3),R_tau(3,3),U_inv(3,3),detF
      real*8 Fp_inv(3,3),Ee_tau(3,3),Re_tau(3,3),Ue_tau(3,3),Fe_tau(3,3)
      real*8 pnu0,damage_t,damage_tau,mag_Dp_tau,pwrinct,stress_power
      real*8 nu1,nu3,nu5
      real*8 lamg_t,lamg_tau,matProps(nprops)
      real*8 E_tau,E_t,thetag_t,thetag_tau
      real*8 coordx,coordy,coordz
      real*8 curvature

      ! Parameters
      !
      real*8 zero,one,two,three,half,third,four,Pi,two_third
      parameter(zero=0.d0,one=1.d0,two=2.d0,three=3.d0,half=0.5d0,
     +     third=1.d0/3.d0,two_third=2.d0/3.d0,four=4.d0,Pi=3.1415926d0)

c      open(unit=80, file='C:\Abaqus_Working_Dir\FeFp\aaMSGS')



      ! pour initial coordinates into the global variable
      if (totaltime.lt. 0.001) then
          do km=1,nblock
             inicoord(nElement(km),1) = coordMp(km,1)
             inicoord(nElement(km),2) = coordMp(km,2)
             inicoord(nElement(km),3) = coordMp(km,3)
          enddo
      end if

      ! Identity matrix for later use.
      !
      call onem(Iden)


      !
      ! START LOOP OVER MATERIAL POINTS:
      !
      do km=1,nblock

           
         ! Copy old and new deformation gradients
         !
         F_t(1,1) = defgradOld(km,1)
         F_t(2,2) = defgradOld(km,2)
         F_t(3,3) = defgradOld(km,3)
         F_t(1,2) = defgradOld(km,4)
         F_tau(1,1) = defgradNew(km,1)
         F_tau(2,2) = defgradNew(km,2)
         F_tau(3,3) = defgradNew(km,3)
         F_tau(1,2) = defgradNew(km,4)
         U_tau(1,1) = stretchNew(km,1)
         U_tau(2,2) = stretchNew(km,2)
         U_tau(3,3) = stretchNew(km,3)
         U_tau(1,2) = stretchNew(km,4)
         if(nshr .lt. 2) then
            ! 2D case
            F_t(2,1) = defgradOld(km,5)
            F_t(1,3) = zero
            F_t(2,3) = zero
            F_t(3,1) = zero
            F_t(3,2) = zero
            F_tau(2,1) = defgradNew(km,5)
            F_tau(1,3) = zero
            F_tau(2,3) = zero
            F_tau(3,1) = zero
            F_tau(3,2) = zero
            U_tau(2,1) = U_tau(1,2)
            U_tau(1,3) = zero
            U_tau(2,3) = zero
            U_tau(3,1) = zero
            U_tau(3,2) = zero
         else
            ! 3D case
            F_t(2,3) = defgradOld(km,5)
            F_t(3,1) = defgradOld(km,6)
            F_t(2,1) = defgradOld(km,7)
            F_t(3,2) = defgradOld(km,8)
            F_t(1,3) = defgradOld(km,9)
            F_tau(2,3) = defgradNew(km,5)
            F_tau(3,1) = defgradNew(km,6)
            F_tau(2,1) = defgradNew(km,7)
            F_tau(3,2) = defgradNew(km,8)
            F_tau(1,3) = defgradNew(km,9)
            U_tau(2,3) = stretchNew(km,5)
            U_tau(3,1) = stretchNew(km,6)
            U_tau(2,1) = U_tau(1,2)
            U_tau(3,2) = U_tau(2,3)
            U_tau(1,3) = U_tau(3,1)
         end if


         if((totalTime.eq.zero).and.(stepTime.eq.zero)) then
            ! Dummy step, initalize state variables
            
            stateOld(km,1)   = one ! growth parameter at t=0
         endif

         ! Read old state variables
         
         thetag_t = stateOld(km,1) ! growth parameter at time t




         ! reads in the original coordinates 
          coordx = inicoord(nElement(km),1)
          coordy = inicoord(nElement(km),2)
          coordz = inicoord(nElement(km),3)




         ! store all the coordinates 
         StorageOld2(nElement(km)-BaseElement2,1) = coordMp(km,1)
         StorageOld2(nElement(km)-BaseElement2,2) = coordMp(km,2)

         ! here is trick I play here, so for full size simulation,
         ! the curvature will be flipped in (-x,y) and (x-y) region,
         ! here I just mirror it...
c         if (StorageOld(nElement(km)-BaseElement,1) * StorageOld(nElement(km)-BaseElement,2).lt.0.0) then
c             StorageOld(nElement(km)-BaseElement,2) = -StorageOld(nElement(km)-BaseElement,2)
c         endif



      do i=1,nElement(km)-BaseElement2

         if (i.eq.1) then 
            
               call GetCurvature(StorageOld2(2,1),StorageOld2(1,1),StorageOld2(3,1),
     +                        StorageOld2(2,2),StorageOld2(1,2),StorageOld2(3,2), 
     +                        curvOld2(1,1))

         elseif(i.lt.NElements2) then


               call GetCurvature(StorageOld2(i,1),StorageOld2(i-1,1),StorageOld2(i+1,1),
     +                        StorageOld2(i,2),StorageOld2(i-1,2),StorageOld2(i+1,2), 
     +                        curvOld2(i,1))

         elseif(i.eq.NElements2) then

               call GetCurvature(StorageOld2(NElements2 -1,1),StorageOld2(NElements2 -2,1),
     +                           StorageOld2(NElements2,1),StorageOld2(NElements2 -1,2),
     +                           StorageOld2(NElements2 -2,2),StorageOld2(NElements2,2), 
     +                           curvOld2(NElements2,1))
         endif


      enddo


       curvature = -curvOld2(nElement(km)-BaseElement2,1)


         !---------------------------------------------------------------
         ! Perform the time integration and compute the 
         !  constitutive response based on the material model.
         
         matProps = props

         if((totalTime.eq.zero).and.(stepTime.eq.zero)) then
            !
            ! dummy step, call elastic response, note dt=-1.0 is sent
            !  into the integ subroutine
            !
            call integ_graytwo(matProps,nprops,F_tau,-1.0,T_tau,thetag_t,thetag_tau,
     +                     coordx,coordy,coordz,curvature)

         else
            !
            ! Perform explicit time integration procedure
            !
            call integ_graytwo(matProps,nprops,F_tau,dt,T_tau,thetag_t,thetag_tau,
     +                     coordx,coordy,coordz,curvature)

         endif
         !---------------------------------------------------------------
         

         ! Update state variables
         !
         stateNew(km,1) = thetag_tau ! growth parameter at time tau

         stateNew(km,2) = coordx
         stateNew(km,3) = coordy
         stateNew(km,4) = coordz

         stateNew(km,5) = curvature

         ! ABAQUS/Explicit uses stress measure (transpose(R) T R)
         !
         call m3inv(U_tau,U_inv)
         R_tau = matmul(F_tau,U_inv)
         T_tau = matmul(transpose(R_tau),matmul(T_tau,R_tau))

         do i=1,ndir
            stressNew(km,i) = T_tau(i,i)
         end do
         if(nshr.ne.0) then
            stressNew(km,ndir+1) = T_tau(1,2)
            if(nshr.ne.1) then
               stressNew(km, ndir+2) = T_tau(2,3)
               if(nshr.ne.2) then
                  stressNew(km,ndir+3) = T_tau(1,3)
               endif
            endif
         endif


         ! Update the specific internal energy
         !
         stress_power = 0.d0
         do i = 1,ndir
            stress_power = stress_power +
     +           0.5*((StressOld(km,i)+StressNew(km,i))*
     +           StrainInc(km,i))
         enddo
         
         select case (nshr)
         case(1)
            stress_power = stress_power + 
     +           0.5*((StressOld(km,ndir+1)+StressNew(km,ndir+1))*
     +           StrainInc(km,ndir+1))
         case(3)
            stress_power = stress_power + 
     +           0.5*(((StressOld(km,ndir+1) + StressNew(km,ndir+1))*
     +           StrainInc(km,ndir+1)) +
     +           ((StressOld(km,ndir+2)+ StressNew(km,ndir+2)) *
     +           StrainInc(km,ndir+2))+
     +           ((StressOld(km,ndir+3) + StressNew(km,ndir+3))*
     +           StrainInc(km,ndir+3)))
         end select
           
         enerInternNew(km) = enerInternOld(km) + 
     +        stress_power/density(km)
           
         enerInelasNew(km) = enerInelasOld(km) + 
     +        pwrinct/density(km)
           
           
      enddo ! end loop over material points

      end subroutine vumatXtrArg_graytwo
***********************************************************************
      subroutine integ_white(Props,nprops,F_tau,dtime,T_tau,
     +                       lamg_t,lamg_tau,coordx,coordy,coordz)

      implicit none


      character*80 cmname,file1
      character*256 jobName,outDir,fileName

      integer i,j,k,l,a,b,c,d,iterError,lenJobName,lenOutDir,nargs,nprops
      parameter(nargs=5)

      real*8 Iden(3,3),F_t(3,3),F_tau(3,3),theta_t,theta_tau,T_tau(3,3)
      real*8 Gshear,effStr,detF,dTRdF(3,3,3,3),SpTanMod(3,3,3,3),trB
      real*8 Finv(3,3),Bdis(3,3),trBdis,Bdis0(3,3),B_tau(3,3)
      real*8 a0(3,1),ac(3,1),growthvec(3,1),thetatotal
      real*8 lambda,mu,thetamax,tau,gammavar,thetacrit
      real*8 Be_tau(3,3),Fg_tau(3,3),Fe_tau(3,3),Je
      real*8 thetag_tau,args(nargs),fac,thetag_t
      real*8 props(nprops),dtime,Jg
      real*8 lamg_t,lamg_tau
      real*8 H0(3,3),C_tau(3,3),HC(3,3),lam_tau,Gaxn,dE,lam0,trHC
      real*8 Ce(3,3,3,3),Cg(3,3,3,3),Cs(3,3,3,3),Klam,Kirk(3,3)
      real*8 jac(3,3,3,3),tmp
      real*8 Fginv(3,3),Cinv(3,3),Piola_tau(3,3)
      real*8 Cg_mat1(3,3,3,3),Cg_mat2(3,3,3,3),Cg_mat3(3,3,3,3),Cg_mat(3,3,3,3)
      real*8 coordx,coordy,coordz


      ! Parameters
      !
      real*8 zero,one,two,half,three,third,nine
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,three=3.d0,
     +     third=1.d0/3.d0,nine=9.d0)


      ! Identity matrix
      !
      call onem(Iden)


      ! Compute the relative volume change
      !
      call mdet(F_tau,detF)


      ! Compute the inverse of the deformation gradient
      !
      call m3inv(F_tau,Finv)


      ! compute the right Cauchy Green tensor 
      ! 
      C_tau = matmul(transpose(F_tau),F_tau)



      ! obtain the referential unit vector 
      ! 
      a0(1,1) = 0.0
      a0(2,1) = 1.0 
      a0(3,1) = 0.0 
c      a0(1,1) = 2.0*coordx/1.0**2.0
c      a0(2,1) = 2.0*coordy/1.2**2.0
c      a0(3,1) = 0.0 


      tmp = sqrt(a0(1,1)*a0(1,1) + a0(2,1)*a0(2,1) + a0(3,1)*a0(3,1))
      a0 = a0/tmp

      ! map the referential unit vector to the current one 
      ! 
      ac = matmul(F_tau,a0)


 
      ! structure tensor 
      H0 = zero
      do i=1,3
         do j=1,3
             H0(i,j) = H0(i,j) + a0(i,1)*a0(j,1)
         enddo
      enddo
      !HC tensor 
      !          
      HC = matmul(H0,C_tau)

      ! its traace
      ! 
      trHC = HC(1,1) + HC(2,2) + HC(3,3)

      ! total stretch
      ! 
      lam_tau = dsqrt(trHC)



      ! Obtain material properties 
      !
       mu        = props(1)
       lambda    = props(2)
       Gaxn      = props(3)
       lam0      = props(4)

      !!!!!!!!!!!!!!!!!!!!!!!!!!! dummy step !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(dtime.lt.zero) then

      lamg_tau = lamg_t
      Fg_tau  = Iden + (lamg_tau - 1.0)*matmul(a0,transpose(a0))
      B_tau = matmul(F_tau, transpose(F_tau))
      Be_tau = B_tau + ((1.0 - lamg_tau**2.0)/lamg_tau**2.0)*matmul(ac,transpose(ac))

 
      call mdet(Fg_tau,Jg)

      Je = detF/Jg
      

      call m3inv(Fg_tau,Fginv)     
      call m3inv(C_tau,Cinv)

      ! compute Cauchy stress 
      ! 
!      T_tau = ((lambda*dlog(Je) - mu)*Iden  + mu*Be_tau)/detF
      T_tau = ((lambda*dlog(Je) - mu)*Iden  + mu*Be_tau)/Je
         return
      endif  
      !!!!!!!!!!!!!!!!!!!!!!!!!!! dummy step !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      ! here we should do a locally iteration to get current area growth
      ! 
      args(1) = dtime 
      args(2) = lamg_t 
      args(3) = Gaxn
      args(4) = lam_tau
      args(5) = lam0



c      call solvegrowth(lamg_tau,args,nargs) 

      lamg_tau = 1.0 ! debug


      ! update  kinematics 
      ! 

      Fg_tau  = Iden + (lamg_tau - 1.0)*matmul(a0,transpose(a0))
      B_tau = matmul(F_tau, transpose(F_tau))
      Be_tau = B_tau + ((1.0 - lamg_tau**2.0)/lamg_tau**2.0)*matmul(ac,transpose(ac))

 
      call mdet(Fg_tau,Jg)

      Je = detF/Jg


      call m3inv(Fg_tau,Fginv)     
      call m3inv(C_tau,Cinv)

      

      ! compute Cauchy stress 
      ! 
!      T_tau = ((lambda*dlog(Je) - mu)*Iden  + mu*Be_tau)/detF
      T_tau = ((lambda*dlog(Je) - mu)*Iden  + mu*Be_tau)/Je


      end subroutine integ_white
****************************************************************************
      subroutine integ_gray(Props,nprops,F_tau,dtime,T_tau,
     +                       thetag_t,thetag_tau,coordx,coordy,coordz,
     +                      curvature)

      implicit none


      character*80 cmname,file1
      character*256 jobName,outDir,fileName

      integer i,j,k,l,iterError,lenJobName,lenOutDir,nargs,nprops
      parameter(nargs=5)

      real*8 Iden(3,3),F_t(3,3),F_tau(3,3),theta_t,theta_tau,T_tau(3,3)
      real*8 Gshear,effStr,detF,dTRdF(3,3,3,3),SpTanMod(3,3,3,3),trB
      real*8 Finv(3,3),Bdis(3,3),trBdis,Bdis0(3,3),B_tau(3,3)
      real*8 a0(3,1),ac(3,1),growthvec(3,1),thetatotal
      real*8 lambda,mu,thetamax,tau,gammavar,thetacrit
      real*8 Be_tau(3,3),Fg_tau(3,3),Fe_tau(3,3),Je
      real*8 thetag_tau,args(nargs),fac,thetag_t
      real*8 props(nprops),dtime,Jg
      real*8 lamg_t,lamg_tau
      real*8 H0(3,3),C_tau(3,3),HC(3,3),lam_tau,Gaxn,dE,lam0,trHC
      real*8 Ce(3,3,3,3),Cg(3,3,3,3),Cs(3,3,3,3),Klam,Kirk(3,3)
      real*8 jac(3,3,3,3),n0(3,1),nc(3,1),mu_g,lambda_g,Gctx
      real*8 coordx,coordy,coordz,tmp,curvature
      real*8 Aconst,charlength,Gcur
      real*8 threshold,NormCurve,G_total
      real*8 Fginv(3,3)



      ! Parameters
      !
      real*8 zero,one,two,half,three,third,nine
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,three=3.d0,
     +     third=1.d0/3.d0,nine=9.d0)


      ! Identity matrix
      !
      call onem(Iden)


      ! Compute the relative volume change
      !
      call mdet(F_tau,detF)


      ! Compute the inverse of the deformation gradient
      !
      call m3inv(F_tau,Finv)
    

      ! Obtain material properties 
      !
       mu_g      = props(5)
       lambda_g  = props(6)
       Gctx      = props(7)
       Aconst    = props(8)
       charlength = props(9)


      ! obtain referential surface outnormal
      a0(1,1) = 2.0*coordx/1.0**2.0
      a0(2,1) = 2.0*coordy/1.2**2.0
      a0(3,1) = 0.0

      tmp = sqrt(a0(1,1)**2.0 + a0(2,1)**2.0 + a0(3,1)**2.0)

      a0 = a0/tmp   

      

      !!!!!!!!!!!!!!!!!!!!!!!!!!! dummy step !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(dtime.lt.zero) then

      
       thetag_tau = thetag_t


      ! update  kinematics 
      ! 
      ! iso growth  
c      Fg_tau  = thetag_tau*Iden
      ! area growth 
      Fg_tau  = dsqrt(thetag_tau)*Iden 
     +          +(1.0 - dsqrt(thetag_tau))*matmul(a0,transpose(a0))

      ! inverse of the growth Fg
      ! 
      call m3inv(Fg_tau,Fginv)


      ! elastic Fe
      ! 
      Fe_tau = matmul(F_tau,Fginv)


      ! Left Cauchy Green tensor  
      ! 
      Be_tau = matmul(Fe_tau,transpose(Fe_tau)) 
 

      ! Jacobian of the Fg
      ! 
      call mdet(Fg_tau,Jg)

      Je = detF/Jg
      
      ! compute Cauchy stress 
      ! 
      T_tau = ((lambda_g*dlog(Je) - mu_g)*Iden  + mu_g*Be_tau)/Je




         return
      endif    
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
      !!!!!!!!!!!!!!!!!!!!!!!!!!! dummy step !!!!!!!!!!!!!!!!!!!!!!!!!!!!






      ! do the time integration here 
      ! 
       NormCurve = curvature*charlength

       Gcur = - Aconst * NormCurve

c       if (NormCurve.ge.0.0) then
c           G_total = Gctx
c       elseif (NormCurve.lt.0.0) then 
c           G_total = Gctx + Gcur
c           G_total = Gctx*(1.0 - Aconst*tanh(NormCurve))
c       endif

       G_total = Gctx*(1.0 - Aconst*tanh(NormCurve))


       thetag_tau = thetag_t + (G_total)*dtime 


c      ! do the time integration here 
c      ! reverse 
c       NormCurve = curvature*charlength
c
c       Gcur =  Aconst * NormCurve
c
c       if (NormCurve.ge.0.0) then
c           G_total = Gctx + Gcur
c       elseif (NormCurve.lt.0.0) then 
c           G_total = Gctx
c       endif
c
c       thetag_tau = thetag_t + (G_total)*dtime 



      ! update  kinematics 
      ! 
      ! iso growth  
c      Fg_tau  = thetag_tau*Iden
      ! area growth 
      Fg_tau  = dsqrt(thetag_tau)*Iden 
     +          +(1.0 - dsqrt(thetag_tau))*matmul(a0,transpose(a0))

      ! inverse of the growth Fg
      ! 
      call m3inv(Fg_tau,Fginv)


      ! elastic Fe
      ! 
      Fe_tau = matmul(F_tau,Fginv)


      ! Left Cauchy Green tensor  
      ! 
      Be_tau = matmul(Fe_tau,transpose(Fe_tau)) 
 

      ! Jacobian of the Fg
      ! 
      call mdet(Fg_tau,Jg)

      Je = detF/Jg
      
      ! compute Cauchy stress 
      ! 
      T_tau = ((lambda_g*dlog(Je) - mu_g)*Iden  + mu_g*Be_tau)/Je


      end subroutine integ_gray
****************************************************************************
      subroutine integ_graytwo(Props,nprops,F_tau,dtime,T_tau,
     +                       thetag_t,thetag_tau,coordx,coordy,coordz,
     +                      curvature)

      implicit none


      character*80 cmname,file1
      character*256 jobName,outDir,fileName

      integer i,j,k,l,iterError,lenJobName,lenOutDir,nargs,nprops
      parameter(nargs=5)

      real*8 Iden(3,3),F_t(3,3),F_tau(3,3),theta_t,theta_tau,T_tau(3,3)
      real*8 Gshear,effStr,detF,dTRdF(3,3,3,3),SpTanMod(3,3,3,3),trB
      real*8 Finv(3,3),Bdis(3,3),trBdis,Bdis0(3,3),B_tau(3,3)
      real*8 a0(3,1),ac(3,1),growthvec(3,1),thetatotal
      real*8 lambda,mu,thetamax,tau,gammavar,thetacrit
      real*8 Be_tau(3,3),Fg_tau(3,3),Fe_tau(3,3),Je
      real*8 thetag_tau,args(nargs),fac,thetag_t
      real*8 props(nprops),dtime,Jg
      real*8 lamg_t,lamg_tau
      real*8 H0(3,3),C_tau(3,3),HC(3,3),lam_tau,Gaxn,dE,lam0,trHC
      real*8 Ce(3,3,3,3),Cg(3,3,3,3),Cs(3,3,3,3),Klam,Kirk(3,3)
      real*8 jac(3,3,3,3),n0(3,1),nc(3,1),mu_g,lambda_g,Gctx
      real*8 coordx,coordy,coordz,tmp,curvature
      real*8 Aconst,charlength,Gcur
      real*8 threshold,NormCurve,G_total
      real*8 Fginv(3,3)



      ! Parameters
      !
      real*8 zero,one,two,half,three,third,nine
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,three=3.d0,
     +     third=1.d0/3.d0,nine=9.d0)


      ! Identity matrix
      !
      call onem(Iden)


      ! Compute the relative volume change
      !
      call mdet(F_tau,detF)


      ! Compute the inverse of the deformation gradient
      !
      call m3inv(F_tau,Finv)
    

      ! Obtain material properties 
      !
       mu_g      = props(5)
       lambda_g  = props(6)
       Gctx      = props(7)
       Aconst    = props(8)
       charlength = props(9)

      ! obtain referential surface outnormal
      a0(1,1) = 2.0*coordx/1.0**2.0
      a0(2,1) = 2.0*coordy/1.2**2.0
      a0(3,1) = 0.0

      tmp = sqrt(a0(1,1)**2.0 + a0(2,1)**2.0 + a0(3,1)**2.0)

      a0 = a0/tmp       

      !!!!!!!!!!!!!!!!!!!!!!!!!!! dummy step !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(dtime.lt.zero) then

      
       thetag_tau = thetag_t


      ! update  kinematics 
      ! 
      ! iso growth  
c      Fg_tau  = thetag_tau*Iden
      ! area growth 
      Fg_tau  = dsqrt(thetag_tau)*Iden 
     +          +(1.0 - dsqrt(thetag_tau))*matmul(a0,transpose(a0))

      ! inverse of the growth Fg
      ! 
      call m3inv(Fg_tau,Fginv)


      ! elastic Fe
      ! 
      Fe_tau = matmul(F_tau,Fginv)


      ! Left Cauchy Green tensor  
      ! 
      Be_tau = matmul(Fe_tau,transpose(Fe_tau)) 
 

      ! Jacobian of the Fg
      ! 
      call mdet(Fg_tau,Jg)

      Je = detF/Jg
      
      ! compute Cauchy stress 
      ! 
      T_tau = ((lambda_g*dlog(Je) - mu_g)*Iden  + mu_g*Be_tau)/Je




         return
      endif    
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
      !!!!!!!!!!!!!!!!!!!!!!!!!!! dummy step !!!!!!!!!!!!!!!!!!!!!!!!!!!!






      ! do the time integration here 
      ! 
       NormCurve = curvature*charlength

       Gcur = - Aconst * NormCurve

c       if (NormCurve.ge.0.0) then
c           G_total = Gctx
c       elseif (NormCurve.lt.0.0) then 
c           G_total = Gctx + Gcur
c           G_total = Gctx*(1.0 - Aconst*tanh(NormCurve))
c       endif

       G_total = Gctx*(1.0 - Aconst*tanh(NormCurve))

       thetag_tau = thetag_t + (G_total)*dtime



c      ! do the time integration here 
c      ! reverse 
c       NormCurve = curvature*charlength
c
c       Gcur =  Aconst * NormCurve
c
c       if (NormCurve.ge.0.0) then
c           G_total = Gctx + Gcur
c       elseif (NormCurve.lt.0.0) then 
c           G_total = Gctx
c       endif
c
c       thetag_tau = thetag_t + (G_total)*dtime 



      ! update  kinematics 
      ! 
      ! iso growth  
c      Fg_tau  = thetag_tau*Iden
      ! area growth 
      Fg_tau  = dsqrt(thetag_tau)*Iden 
     +          +(1.0 - dsqrt(thetag_tau))*matmul(a0,transpose(a0))

      ! inverse of the growth Fg
      ! 
      call m3inv(Fg_tau,Fginv)


      ! elastic Fe
      ! 
      Fe_tau = matmul(F_tau,Fginv)


      ! Left Cauchy Green tensor  
      ! 
      Be_tau = matmul(Fe_tau,transpose(Fe_tau)) 
 

      ! Jacobian of the Fg
      ! 
      call mdet(Fg_tau,Jg)

      Je = detF/Jg
      
      ! compute Cauchy stress 
      ! 
      T_tau = ((lambda_g*dlog(Je) - mu_g)*Iden  + mu_g*Be_tau)/Je


      end subroutine integ_graytwo
****************************************************************************
      subroutine solvegrowth(root,args,nargs)

      ! This subroutine will numerically solve for current axonal growth

      implicit none

      integer maxit,j,nargs

      real*8 xacc,f,df,fl,fh,xl,xh,x1,x2,swap,root,dxold,one
      real*8 dx,args(nargs),zero,rootOld,temp,rootMax,rootMin
      real*8 LambdaBar_tau,LambdaBar_t

      parameter(maxit=100)
      parameter(xacc=1.d-3,zero=0.d0,one=1.d0)
c      parameter(zero=0.d0,one=1.d0)


      
!      write(*,*) 'dlam=',xacc

!     if(.GE.) then

!      else


!      endif
  



      rootMax =  10.0
      rootMin =  0.0

      x1 = rootMin
      x2 = rootMax
      call gfunc(x1,FL,DF,args,nargs)
      call gfunc(x2,FH,DF,args,nargs)

      if(fl*fh.ge.zero) then
         write(*,*) 'FYI, root not bracketed on lamg'
         write(*,*) 'fl=',fl
         write(*,*) 'x1=',x1
         write(*,*) 'fh=',fh
         write(*,*) 'x2=',x2
         call xit
         return
      endif

C
C		ORIENT THE SEARCH SO THAT F(XL) < 0.
C
      IF( FL .LT. 0.D0 ) THEN
         XL = X1
         XH = X2
      ELSE
         XH = X1
         XL = X2
         SWAP = FL
         FL = FH
         FH = SWAP
      END IF
C
C		INITIALIZE THE GUESS FOR THE ROOT, THE ''STEP SIZE
C		BEFORE LAST'', AND THE LAST STEP
C
      if(rootOld.lt.rootMin) rootOld = rootMin
      if(rootOld.gt.rootMax) rootOld = rootMax
      ROOT = 0.5D0 *( X1 + X2)
      DXOLD = DABS(X2 - X1)
      DX = DXOLD
      
      call gfunc(root,F,DF,args,nargs)

C
C			LOOP OVER ALLOWED ITERATIONS
C
      DO 10 J = 1,MAXIT
C
C			BISECT IF NEWTON OUT OF RANGE, OR NOT DECREASING
C			FAST ENOUGH.
C
         IF( ((ROOT-XH)*DF - F)*((ROOT - XL)*DF -F) .GE. 0.D0
     +        .OR. DABS(2.D0*F) .GT. DABS(DXOLD*DF) ) THEN

            DXOLD = DX
            DX = 0.5D0*(XH-XL)
            ROOT = XL + DX
            IF( XL .EQ. ROOT ) THEN
C
C			CHANGE IN ROOT IS NEGLIGIBLE
C
               RETURN
            END IF

         ELSE
C
C			NEWTON STEP IS ACCEPTABLE. TAKE IT.
C
            DXOLD = DX
            DX = F/DF
            TEMP = ROOT
            ROOT = ROOT - DX
            IF( TEMP .EQ. ROOT) THEN
C
C			 CHANGE IN ROOT IS NEGLIGIBLE
C
               RETURN
            END IF

         END IF
C
C		CONVERVEGENCE CRITERION
C
         IF( DABS(DX) .LT. XACC) RETURN

C
C			THE ONE NEW FUNCTION EVALUATION PER ITERATION
C
         call gfunc(root,F,DF,args,nargs)

C
C		MAINTAIN THE BRACKET ON THE ROOT
C
         IF( F .LT. 0.D0) THEN
            XL = ROOT
            FL = F
         ELSE
            XH = ROOT
            FH = F
         END IF

 10   CONTINUE

      WRITE(*,'(/1X,A)') 'solveRate EXCEEDING MAXIMUM ITERATIONS'
      WRITE(80,'(/1X,A)') 'solveRate EXCEEDING MAXIMUM ITERATIONS'

      return
      end subroutine solvegrowth

****************************************************************************
      subroutine gfunc(lamg_tau,f,df,args,nargs)

      implicit none

      integer nargs

      real*8 args(nargs),f,df
      real*8 dtime,thetag_t,thetamax,thetatotal,gammavar,tau,thetacrit
      real*8 k_tau,Phig_tau,dkdtheta,dphidtheta,thetag_tau
      real*8 lamg_t,Gaxn,lam_tau,lam0,lamg_tau


      real*8 zero,one,two,three,third
      parameter(zero=0.d0,one=1.d0,two=2.d0,three=3.d0,third=1.d0/3.d0)



      
      ! Obtain relevant quantities
      !
      dtime   = args(1)
      lamg_t  = args(2)
      Gaxn    = args(3)
      lam_tau = args(4)
      lam0    = args(5)


      ! Compute the residual
      !


      f = lamg_tau - lamg_t - Gaxn*(lam_tau/lamg_tau - lam0)*dtime


      ! Compute the tangent
      !


      df = 1.0 + Gaxn*(lam_tau/(lamg_tau**2.0))*dtime



      return
      end subroutine gfunc

***********************************************************************

      subroutine jac2D(SpTanMod,ddsdde)

      real*8 SpTanMod(3,3,3,3),ddsdde(4,4)

      ddsdde(1,1) = SpTanMod(1,1,1,1)
      ddsdde(1,2) = SpTanMod(1,1,2,2)
      ddsdde(1,3) = SpTanMod(1,1,3,3)
      ddsdde(1,4) = SpTanMod(1,1,1,2)

      ddsdde(2,1) = SpTanMod(2,2,1,1)
      ddsdde(2,2) = SpTanMod(2,2,2,2)
      ddsdde(2,3) = SpTanMod(2,2,3,3)
      ddsdde(2,4) = SpTanMod(2,2,1,2)

      ddsdde(3,1) = SpTanMod(3,3,1,1)
      ddsdde(3,2) = SpTanMod(3,3,2,2)
      ddsdde(3,3) = SpTanMod(3,3,3,3)
      ddsdde(3,4) = SpTanMod(3,3,1,2)

      ddsdde(4,1) = SpTanMod(1,2,1,1)
      ddsdde(4,2) = SpTanMod(1,2,2,2)
      ddsdde(4,3) = SpTanMod(1,2,3,3)
      ddsdde(4,4) = SpTanMod(1,2,1,2)

      end subroutine jac2D

***********************************************************************

      subroutine jac3D(SpTanMod,ddsdde)

      real*8 SpTanMod(3,3,3,3),ddsdde(6,6)

      ddsdde(1,1) = SpTanMod(1,1,1,1)
      ddsdde(1,2) = SpTanMod(1,1,2,2)
      ddsdde(1,3) = SpTanMod(1,1,3,3)
      ddsdde(1,4) = SpTanMod(1,1,1,2)
      ddsdde(1,5) = SpTanMod(1,1,1,3)
      ddsdde(1,6) = SpTanmod(1,1,2,3)

      ddsdde(2,1) = SpTanMod(2,2,1,1)
      ddsdde(2,2) = SpTanMod(2,2,2,2)
      ddsdde(2,3) = SpTanMod(2,2,3,3)
      ddsdde(2,4) = SpTanMod(2,2,1,2)
      ddsdde(2,5) = SpTanMod(2,2,1,3)
      ddsdde(2,6) = SpTanmod(2,2,2,3)

      ddsdde(3,1) = SpTanMod(3,3,1,1)
      ddsdde(3,2) = SpTanMod(3,3,2,2)
      ddsdde(3,3) = SpTanMod(3,3,3,3)
      ddsdde(3,4) = SpTanMod(3,3,1,2)
      ddsdde(3,5) = SpTanMod(3,3,1,3)
      ddsdde(3,6) = SpTanmod(3,3,2,3)

      ddsdde(4,1) = SpTanMod(1,2,1,1)
      ddsdde(4,2) = SpTanMod(1,2,2,2)
      ddsdde(4,3) = SpTanMod(1,2,3,3)
      ddsdde(4,4) = SpTanMod(1,2,1,2)
      ddsdde(4,5) = SpTanMod(1,2,1,3)
      ddsdde(4,6) = SpTanmod(1,2,2,3)

      ddsdde(5,1) = SpTanMod(1,3,1,1)
      ddsdde(5,2) = SpTanMod(1,3,2,2)
      ddsdde(5,3) = SpTanMod(1,3,3,3)
      ddsdde(5,4) = SpTanMod(1,3,1,2)
      ddsdde(5,5) = SpTanMod(1,3,1,3)
      ddsdde(5,6) = SpTanmod(1,3,2,3)

      ddsdde(6,1) = SpTanMod(2,3,1,1)
      ddsdde(6,2) = SpTanMod(2,3,2,2)
      ddsdde(6,3) = SpTanMod(2,3,3,3)
      ddsdde(6,4) = SpTanMod(2,3,1,2)
      ddsdde(6,5) = SpTanMod(2,3,1,3)
      ddsdde(6,6) = SpTanmod(2,3,2,3)

      end subroutine jac3D

C**********************************************************************
C	THE FOLLOWING SUBROUTINES ARE UTILITY ROUTINES
C**********************************************************************
	SUBROUTINE ZEROV(V,SIZE)

C	THIS SUBROUTINE STORES THE ZERO VECTOR IN A VECTOR V
C	OF SIZE SIZE.
C**********************************************************************

	INTEGER SIZE
	REAL*8 V(0:SIZE-1)

	DO 1 I=0,SIZE
	  V(I) = 0.D0
1	CONTINUE
	
	RETURN
	END

C**********************************************************************
      	SUBROUTINE ZEROM(A)
C
C	THIS SUBROUTINE SETS ALL ENTRIES OF A 3 BY 3 MATRIX TO 0.D0.
C**********************************************************************

        REAL*8 A(3,3)

	DO 1 I=1,3
	  DO 1 J=1,3
	    A(I,J) = 0.D0
1	CONTINUE
C	
	RETURN
	END

C**********************************************************************
	SUBROUTINE ONEM(A)

C	THIS SUBROUTINE STORES THE IDENTITY MATRIX IN THE 
C	3 BY 3 MATRIX [A]
C**********************************************************************

        REAL*8 A(3,3)
        DATA ZERO/0.D0/
        DATA ONE/1.D0/

	DO 1 I=1,3
	  DO 1 J=1,3
	    IF (I .EQ. J) THEN
              A(I,J) = 1.0
            ELSE
              A(I,J) = 0.0
            ENDIF
1       CONTINUE

	RETURN
	END

C**********************************************************************
	SUBROUTINE MTRANS(A,ATRANS)
 
C	THIS SUBROUTINE CALCULATES THE TRANSPOSE OF AN 3 BY 3 
C	MATRIX [A], AND PLACES THE RESULT IN ATRANS. 
C**********************************************************************

	REAL*8 A(3,3),ATRANS(3,3)

	DO 1 I=1,3
 	  DO 1 J=1,3
	    ATRANS(J,I) = A(I,J)
1	CONTINUE

	RETURN
	END

C**********************************************************************
	SUBROUTINE MPROD(A,B,C)
 
C 	THIS SUBROUTINE MULTIPLIES TWO 3 BY 3 MATRICES [A] AND [B],
C 	AND PLACE THEIR PRODUCT IN MATRIX [C]. 
C**********************************************************************

	REAL*8 A(3,3),B(3,3),C(3,3)

	DO 2 I = 1, 3
	  DO 2 J = 1, 3
	    C(I,J) = 0.D0
	    DO 1 K = 1, 3
	      C(I,J) = C(I,J) + A(I,K) * B(K,J)                       
1	    CONTINUE
2	CONTINUE
C
	RETURN
	END

C**********************************************************************
	SUBROUTINE MPROD4(A,B,C)
 
C	THIS SUBROUTINE MULTIPLIES TWO 3 BY 3 MATRICES [A] AND [B],
C 	AND PLACE THEIR PRODUCT IN MATRIX [C]. 
C**********************************************************************

	REAL*8 A(4,4),B(4,4),C(4,4)

	DO 2 I = 1, 4
   	  DO 2 J = 1, 4
	    C(I,J) = 0.D0
	    DO 1 K = 1, 4
	      C(I,J) = C(I,J) + A(I,K) * B(K,J)                       
1	    CONTINUE
2	CONTINUE

	RETURN
	END

C**********************************************************************
	SUBROUTINE DOTPM(A,B,C)

C	THIS SUBROUTINE CALCULATES THE SCALAR PRODUCT OF TWO
C	3 BY 3 MATRICES [A] AND [B] AND STORES THE RESULT IN THE
C	SCALAR C.
C**********************************************************************

	REAL*8 A(3,3),B(3,3),C

	C = 0.D0
	DO 1 I = 1,3
	  DO 1 J = 1,3
            C = C + A(I,J)*B(I,J)
1	CONTINUE
C
	RETURN
	END

C**********************************************************************
	SUBROUTINE MDET(A,DET)
 
C 	THIS SUBROUTINE CALCULATES THE DETERMINANT
C 	OF A 3 BY 3 MATRIX [A].
C**********************************************************************

	REAL*8  A(3,3), DET

	DET =	  A(1,1)*A(2,2)*A(3,3) 
     +	        + A(1,2)*A(2,3)*A(3,1)
     +	        + A(1,3)*A(2,1)*A(3,2)
     +		- A(3,1)*A(2,2)*A(1,3)
     +		- A(3,2)*A(2,3)*A(1,1)
     +		- A(3,3)*A(2,1)*A(1,2)

	RETURN
	END

C**********************************************************************
	SUBROUTINE M3INV(A,AINV)

C 	THIS SUBROUTINE CALCULATES THE THE INVERSE OF A 3 BY 3 MATRIX
C	[A] AND PLACES THE RESULT IN [AINV]. 
C 	IF DET(A) IS ZERO, THE CALCULATION
C 	IS TERMINATED AND A DIAGNOSTIC STATEMENT IS PRINTED.
C**********************************************************************

	REAL*8  A(3,3), AINV(3,3), DET, ACOFAC(3,3), AADJ(3,3)

C	A(3,3)	        -- THE MATRIX WHOSE INVERSE IS DESIRED.
C	DET		-- THE COMPUTED DETERMINANT OF [A].
C	ACOFAC(3,3)	-- THE MATRIX OF COFACTORS OF A(I,J).
C			   THE SIGNED MINOR (-1)**(I+J)*M_IJ
C			   IS CALLED THE COFACTOR OF A(I,J).
C	AADJ(3,3)	-- THE ADJOINT OF [A]. IT IS THE MATRIX
C			   OBTAINED BY REPLACING EACH ELEMENT OF
C			   [A] BY ITS COFACTOR, AND THEN TAKING
C			   TRANSPOSE OF THE RESULTING MATRIX.
C	AINV(3,3)	-- RETURNED AS INVERSE OF [A].
C			   [AINV] = [AADJ]/DET.
C----------------------------------------------------------------------

	CALL MDET(A,DET)
	IF ( DET .EQ. 0.D0 ) THEN
	  write(*,10)
	  STOP
	ENDIF
	CALL MCOFAC(A,ACOFAC)
	CALL MTRANS(ACOFAC,AADJ)
	DO 1 I = 1,3
	DO 1 J = 1,3
	     AINV(I,J) = AADJ(I,J)/DET
1	CONTINUE
10	FORMAT(5X,'--ERROR IN M3INV--- THE MATRIX IS SINGULAR',/,
     +         10X,'PROGRAM TERMINATED')

	RETURN
	END

C**********************************************************************
	SUBROUTINE MCOFAC(A,ACOFAC)
 
C 	THIS SUBROUTINE CALCULATES THE COFACTOR OF A 3 BY 3 MATRIX [A],
C 	AND PLACES THE RESULT IN [ACOFAC]. 
C**********************************************************************

	REAL*8  A(3,3), ACOFAC(3,3)

	ACOFAC(1,1) = A(2,2)*A(3,3) - A(3,2)*A(2,3)
	ACOFAC(1,2) = -(A(2,1)*A(3,3) - A(3,1)*A(2,3))
	ACOFAC(1,3) = A(2,1)*A(3,2) - A(3,1)*A(2,2)
	ACOFAC(2,1) = -(A(1,2)*A(3,3) - A(3,2)*A(1,3))
	ACOFAC(2,2) = A(1,1)*A(3,3) - A(3,1)*A(1,3)
	ACOFAC(2,3) = -(A(1,1)*A(3,2) - A(3,1)*A(1,2))
	ACOFAC(3,1) = A(1,2)*A(2,3)  - A(2,2)*A(1,3)
	ACOFAC(3,2) = -(A(1,1)*A(2,3) - A(2,1)*A(1,3))
	ACOFAC(3,3) = A(1,1)*A(2,2) - A(2,1)*A(1,2)

	RETURN
	END

C**********************************************************************
        SUBROUTINE INVAR(A,IA,IIA,IIIA)

C	THIS SUBROUTINE CALCULATES THE PRINCIPAL INVARIANTS 
C	IA, IIA, IIIA OF A TENSOR [A].
C**********************************************************************

        REAL*8 A(3,3), AD(3,3),AD2(3,3), DETA, IA,IIA,IIIA

        DO 1 I=1,3
          DO 1 J=1,3
            AD(I,J) = A(I,J)
1       CONTINUE
        IA = AD(1,1) + AD(2,2) + AD(3,3)

C	CALCULATE THE SQUARE OF [AD]

        CALL MPROD(AD,AD,AD2)
        IIA =0.5D0 * ( IA*IA - ( AD2(1,1) + AD2(2,2) + AD2(3,3) ) )

        CALL  MDET(AD,DETA)
        IIIA = DETA

        RETURN
        END

C**********************************************************************
	SUBROUTINE TRACEM(A,TRA)

C	THIS SUBROUTINE CALCULATES THE TRACE OF A 3 BY 3 MATRIX [A]
C	AND STORES THE RESULT IN THE SCALAR TRA
C**********************************************************************

	REAL*8 A(3,3),TRA

	TRA = A(1,1) + A(2,2) + A(3,3)

	RETURN 
	END

C**********************************************************************
	SUBROUTINE DEVM(A,ADEV)

C	THIS SUBROUTINE CALCULATES THE DEVIATORIC PART OF A
C	3 BY 3 MATRIX [A]
C**********************************************************************

	REAL*8 A(3,3),TRA,ADEV(3,3),IDEN(3,3)

	CALL TRACEM(A,TRA)
	CALL ONEM(IDEN)
	CALL ZEROM(ADEV)

	DO 1 I = 1,3
	  DO 1 J = 1,3
	    ADEV(I,J) = A(I,J) - (1.D0/3.D0)*TRA*IDEN(I,J)
1	CONTINUE

	RETURN
	END

C**********************************************************************
	SUBROUTINE EQUIVS(S,SB)

C	THIS SUBROUTINE CALCULATES THE EQUIVALENT SHEAR STRESS SB
C	CORRESPONDING TO A 3 BY 3 STRESS MATRIX [S]
C**********************************************************************

	REAL*8 S(3,3),SDEV(3,3),SDOTS,SB

	SB = 0.D0
	SDOTS = 0.D0

	CALL DEVM(S,SDEV)
	CALL DOTPM(SDEV,SDEV,SDOTS)
	SB = DSQRT(1.5D0*SDOTS)

	RETURN
	END
C **********************************************************************
        SUBROUTINE PRESSURE(A,PRA)
C
C       THIS SUBROUTINE CALCULATES THE MEAN NORMAL PRESSURE
C       OF A 3 BY 3 MATRIX [A]
C       AND STORES THE RESULT IN THE SCALAR PRA
C ----------------------------------------------------------------------
C       VARIABLES
C
        REAL*8 A(3,3),PRA

        PRA = -(1.D0 / 3.D0)*( A(1,1) + A(2,2) + A(3,3) )

        RETURN 
        END
C**********************************************************************
	SUBROUTINE PRTMAT(A,M,N)
C**********************************************************************

	INTEGER M,N
	REAL*8 A(M,N)	  

	DO 10 K=1,M
	  WRITE(80,'(2X,6E12.4,2X)') (A(K,L), L=1,N)
10      CONTINUE

        RETURN
        END

C**********************************************************************
	SUBROUTINE PRTVEC(A,M)
C**********************************************************************

	INTEGER M
	REAL*8 A(M)	  

	WRITE(80,'(2X,6E12.4,2X)') (A(K), K=1,M)

        RETURN
	END
C*************************************************************************	  
