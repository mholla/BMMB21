************************************************************************
!  
! A umat for manuscript: 'Numerical investigation of biomechanically-coulped
! growth in cortical folding'
! 
! 
! Shuolun Wang, August 2019, Implemented in Abaqus 6.19
! 
! NOTE: the suboutine should be used along with table.for 
************************************************************************
!     State Variables
!     --------------------------------------------------------------
!      statev(1) = growth parameter 
!      statev(2) = coordx     
!      statev(3) = coordy 
!      statev(4) = coordz 
!      statev(5) = mean curvature 
!      statev(6) = SI index         
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
          parameter(NElements=1600) 
          ! set the element number where gray matter starts from 
          ! min - 1
          parameter(BaseElement=0)    
          parameter(BaseElement2=1600)                  

          real*8 StorageOld(NElements,3)
          real*8 StorageOld2(NElements,3)

          real*8 Storagesnap(NElements,3)
          real*8 Storagesnap2(NElements,3)

          real*8 inicoord(11767,3)
 
      end module    
***********************************************************************
      subroutine calc_normal(point,neighbor1,neighbor2,neighbor3,
     +                       neighbor4,neighbor5,neighbor6,
     +                       neighbor7,neighbor8,normal)
      
      implicit none
      real*8 point(3,1),normal(3),neighbor1(3,1),neighbor2(3,1),
     +  neighbor3(3,1),neighbor4(3,1),neighbor5(3,1),neighbor6(3,1),
     +  neighbor7(3,1),neighbor8(3,1)
     
      real*8 A(9,3),b(9,1),ATA(3,3),RHS(3,1),ATA_INV(3,3),sol(3,1),mag
     
      A(1,1) = neighbor1(1,1)
      A(1,2) = neighbor1(2,1)
      A(1,3) = 1.d0
      A(2,1) = neighbor2(1,1)
      A(2,2) = neighbor2(2,1)
      A(2,3) = 1.d0
      A(3,1) = neighbor3(1,1)
      A(3,2) = neighbor3(2,1)
      A(3,3) = 1.d0
      A(4,1) = neighbor4(1,1)
      A(4,2) = neighbor4(2,1)
      A(4,3) = 1.d0
      A(5,1) = neighbor5(1,1)
      A(5,2) = neighbor5(2,1)
      A(5,3) = 1.d0
      A(6,1) = neighbor6(1,1)
      A(6,2) = neighbor6(2,1)
      A(6,3) = 1.d0
      A(7,1) = neighbor7(1,1)
      A(7,2) = neighbor7(2,1)
      A(7,3) = 1.d0
      A(8,1) = neighbor8(1,1)
      A(8,2) = neighbor8(2,1)
      A(8,3) = 1.d0
      A(9,1) = point(1,1)
      A(9,2) = point(2,1)
      A(9,3) = 1.d0
      
      b(1,1) = neighbor1(3,1)
      b(2,1) = neighbor2(3,1)
      b(3,1) = neighbor3(3,1)
      b(4,1) = neighbor4(3,1)
      b(5,1) = neighbor5(3,1)
      b(6,1) = neighbor6(3,1)
      b(7,1) = neighbor7(3,1)
      b(8,1) = neighbor8(3,1)
      b(9,1) = point(3,1)
      
C
C       The optimal solution of this overconstrained problem is
C       gradz = (A^T A)^(-1) A^T b
C      
      ATA = matmul(transpose(A),A)
      RHS = matmul(transpose(A),b)
      call m3inv(ATA,ATA_INV)
      sol = matmul(ATA_INV,RHS)
      
      mag = dsqrt(sol(1,1)**2.d0 + sol(2,1)**2.d0 + 1.d0)
      normal(1) = sol(1,1)/mag
      normal(2) = sol(2,1)/mag
      normal(3) = -1.d0/mag
      
      return
      end
C************************************************************************
      subroutine calc_curvature(point,normal,neighbor1,neighbor2,
     +                          neighbor3,neighbor4,neighbor5,
     +                          neighbor6,neighbor7,neighbor8,
     +                          curvature,curve1,curve2)
     
      implicit none
      real*8 point(3,1),normal(3),neighbor1(3,1),neighbor2(3,1),
     +       neighbor3(3,1),neighbor4(3,1),neighbor5(3,1),
     +       neighbor6(3,1),neighbor7(3,1),neighbor8(3,1),
     +       curvature,rot(3,3),neighbor1_prime(3,1),
     +       neighbor2_prime(3,1),neighbor3_prime(3,1),
     +       neighbor4_prime(3,1),neighbor5_prime(3,1),
     +       neighbor6_prime(3,1),neighbor7_prime(3,1),
     +       neighbor8_prime(3,1),gradz(3,1),
     +       curve1,curve2
      
      if ((normal(1).eq.0.d0).and.(normal(2).eq.0.d0)) then
        call onem(rot)
      else
        rot(1,1) =  (normal(3)*normal(1)**2.d0 + normal(2)**2.d0)/
     +              (normal(1)**2.d0 + normal(2)**2.d0)
        rot(1,2) = -normal(1)*normal(2)*(1.d0 - normal(3))/
     +              (normal(1)**2.d0 + normal(2)**2.d0)
        rot(1,3) = -normal(1)
        rot(2,1) = -normal(1)*normal(2)*(1.d0 - normal(3))/
     +              (normal(1)**2.d0 + normal(2)**2.d0)
        rot(2,2) =  (normal(1)**2.d0 + normal(3)*normal(2)**2.d0)/
     +              (normal(1)**2.d0 + normal(2)**2.d0)
        rot(2,3) = -normal(2)
        rot(3,1) =  normal(1)
        rot(3,2) =  normal(2)
        rot(3,3) =  normal(3)
      end if
      
      neighbor1_prime = matmul(rot,neighbor1-point)
      neighbor2_prime = matmul(rot,neighbor2-point)
      neighbor3_prime = matmul(rot,neighbor3-point)
      neighbor4_prime = matmul(rot,neighbor4-point)
      neighbor5_prime = matmul(rot,neighbor5-point)
      neighbor6_prime = matmul(rot,neighbor6-point)
      neighbor7_prime = matmul(rot,neighbor7-point)
      neighbor8_prime = matmul(rot,neighbor8-point)
      
      call calc_grad(neighbor1_prime,neighbor2_prime,neighbor3_prime,
     +                neighbor4_prime,neighbor5_prime,neighbor6_prime,
     +                neighbor7_prime,neighbor8_prime,gradz)
     
      curvature = -(gradz(1,1) + gradz(3,1))

      curve1 = -gradz(1,1)
      curve2 = -gradz(3,1)
      
      return
      end

C************************************************************************
C
C       Calculate first and second partial derivatives
C      
      subroutine calc_grad(neighbor1_prime,neighbor2_prime,
     +                      neighbor3_prime,neighbor4_prime,
     +                      neighbor5_prime,neighbor6_prime,
     +                      neighbor7_prime,neighbor8_prime,
     +                      gradz)
      implicit none
      real*8 neighbor1_prime(3,1), neighbor2_prime(3,1), 
     +       neighbor3_prime(3,1), neighbor4_prime(3,1),
     +       neighbor5_prime(3,1), neighbor6_prime(3,1),
     +       neighbor7_prime(3,1), neighbor8_prime(3,1),
     +       gradz(3,1)
      
      real*8 A(8,3), b(8,1), ATA(3,3), RHS(3,1), ATA_inv(3,3)
C       
C       Quadratic interpolation for function z = f(x,y)
C      
C       df = gradz(1) * 0.5*dx^2 + 
C            gradz(2) * dx*dy + gradz(3) * 0.5*dy^2
C      
C       The eight nearest points give eight equations for five
C       unknowns.  Find the optimal solution.
C       
C              A(j,1) = 0.5*dx^2
C              A(j,2) = dx*dy
C              A(j,3) = 0.5*dy^2
C      
C              b(j,1)   = df
C      
C       Solve for system of linear equations:  A gradz = b
C      
      A(1,1) = 0.5d0*neighbor1_prime(1,1)**2.d0
      A(1,2) = neighbor1_prime(1,1)*neighbor1_prime(2,1)
      A(1,3) = 0.5d0*neighbor1_prime(2,1)**2.d0

      A(2,1) = 0.5d0*neighbor2_prime(1,1)**2.d0
      A(2,2) = neighbor2_prime(1,1)*neighbor2_prime(2,1)
      A(2,3) = 0.5d0*neighbor2_prime(2,1)**2.d0

      A(3,1) = 0.5d0*neighbor3_prime(1,1)**2.d0
      A(3,2) = neighbor3_prime(1,1)*neighbor3_prime(2,1)
      A(3,3) = 0.5d0*neighbor3_prime(2,1)**2.d0

      A(4,1) = 0.5d0*neighbor4_prime(1,1)**2.d0
      A(4,2) = neighbor4_prime(1,1)*neighbor4_prime(2,1)
      A(4,3) = 0.5d0*neighbor4_prime(2,1)**2.d0

      A(5,1) = 0.5d0*neighbor5_prime(1,1)**2.d0
      A(5,2) = neighbor5_prime(1,1)*neighbor5_prime(2,1)
      A(5,3) = 0.5d0*neighbor5_prime(2,1)**2.d0

      A(6,1) = 0.5d0*neighbor6_prime(1,1)**2.d0
      A(6,2) = neighbor6_prime(1,1)*neighbor6_prime(2,1)
      A(6,3) = 0.5d0*neighbor6_prime(2,1)**2.d0

      A(7,1) = 0.5d0*neighbor7_prime(1,1)**2.d0
      A(7,2) = neighbor7_prime(1,1)*neighbor7_prime(2,1)
      A(7,3) = 0.5d0*neighbor7_prime(2,1)**2.d0

      A(8,1) = 0.5d0*neighbor8_prime(1,1)**2.d0
      A(8,2) = neighbor8_prime(1,1)*neighbor8_prime(2,1)
      A(8,3) = 0.5d0*neighbor8_prime(2,1)**2.d0
      
      b(1,1) = neighbor1_prime(3,1)
      b(2,1) = neighbor2_prime(3,1)
      b(3,1) = neighbor3_prime(3,1)
      b(4,1) = neighbor4_prime(3,1)
      b(5,1) = neighbor5_prime(3,1)
      b(6,1) = neighbor6_prime(3,1)
      b(7,1) = neighbor7_prime(3,1)
      b(8,1) = neighbor8_prime(3,1)
C
C       The optimal solution of this overconstrained problem is
C       gradz = (A^T A)^(-1) A^T b
C      
      ATA = matmul(transpose(A),A)
      RHS = matmul(transpose(A),b)
      call m3inv(ATA,ATA_inv)
      gradz = matmul(ATA_inv,RHS)
C
      end subroutine calc_grad
************************************************************************
C subroutine for calculatin the shape index 
      subroutine SIFunc(curve1,curve2,SI)

      real*8 curve1, curve2, SI
      real*8 H,K
      ! Parameters
      !
      real*8 zero,one,two,three,half,third,four,Pi,two_third
      parameter(zero=0.d0,one=1.d0,two=2.d0,three=3.d0,half=0.5d0,
     +     third=1.d0/3.d0,two_third=2.d0/3.d0,four=4.d0,Pi=3.1415926d0)


      ! mean curvature 
      H = (curve1 + curve2)/2.0
      
      ! Gaussian curvature
      K = curve1 * curve2


      ! shape index 
      SI = - (2.0/pi) * atan(H/sqrt(H**2.0 - K))



c      if (curve1.gt.curve2) then 

c        SI = (2.0/pi)*atan((curve1 + curve2)/(curve2 - curve1))

c      elseif (curve2.gt.curve1) then 

c        SI = (2.0/pi)*atan((curve1 + curve2)/(curve1 - curve2))

c      endif

c      SI = (2.0/pi)*atan((curve1 + curve2)/(curve2 - curve1))


      End subroutine SIFunc
************************************************************************
C************************************************************************
C subroutine for searching the centroid points of eight nearest neighbors
C around the current point 
C 
      SUBROUTINE GoSearch(ElementCentroidArry,CurrentElementNum,
     +                    n1,n2,n3,n4,n5,n6,n7,n8)


      use GlobalStorage


      integer CurrentElementNum
      integer i,j,n,indexx
      
      real*8 n1,n2,n3,n4,n5,n6,n7,n8
      real*8 ElementCentroidArry(NElements,3)
      real*8 coordx,coordy,coordz
      real*8 distance(2,NElements),temp,temp2
  

c      write(*,*),'current element=',CurrentElementNum

   
      ! obtain the centroid coordinate of the current element 
      !
      coordx = ElementCentroidArry(CurrentElementNum,1)
      coordy = ElementCentroidArry(CurrentElementNum,2)
      coordz = ElementCentroidArry(CurrentElementNum,3)


 
      ! compare the distance against all elements
      ! 
      do i = 1,NElements 
             distance(1,i) = sqrt((ElementCentroidArry(i,1) - coordx)**2.0
     +                           +(ElementCentroidArry(i,2) - coordy)**2.0
     +                           +(ElementCentroidArry(i,3) - coordz)**2.0)

             ! note here I label the distance, so I could 
             ! keep track its index later, sweet... 
             distance(2,i) = i 


      enddo

 
      ! here I  perform the inserton sort to put the distance array in a acending
      ! order 

      n = NElements

      do i = 2,n
         
         j = i - 1
         indexx = i

            do while(j.ge.1) 
            

                if(distance(1,indexx) .lt. distance(1,j)) then 
                  
                  temp = distance(1,j)      
                  temp2 = distance(2,j)

                  distance(1,j) = distance(1,indexx)
                  distance(2,j) = distance(2,indexx)

                 distance(1,indexx) = temp
                 distance(2,indexx) = temp2
                 
                 indexx = j
                endif

             j=j-1  

            enddo      
             
      enddo 

      n1 = distance(2,2)
      n2 = distance(2,3) 
      n3 = distance(2,4) 
      n4 = distance(2,5) 
      n5 = distance(2,6) 
      n6 = distance(2,7) 
      n7 = distance(2,8) 
      n8 = distance(2,9) 



      END SUBROUTINE
***********************************************************************
      subroutine SearchTable1 (table1)

      use GlobalStorage
      real*8 table1(NElements,9)
      
      include '3Dblock_table1.for'
     
   
      end subroutine SearchTable1
***********************************************************************
      subroutine SearchTable2 (table2)

      use GlobalStorage
      real*8 table2(NElements,9)
      
      include '3Dblock_table2.for'
     
   
      end subroutine SearchTable2
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
      ELSE IF(CMNAME(1:4) .EQ. 'GRAY') THEN
      !
      ! this is gray matter 
      !
      call  vumatXtrArg_gray (jblock(i_umt_nblock),
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
      ELSE IF(CMNAME(1:4) .EQ. 'GRAZ') THEN
      !
      ! this is gray matter 
      !
      call  vumatXtrArg_graz (jblock(i_umt_nblock),
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
      if (totaltime.lt. 10.0) then
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
      subroutine vumatXtrArg_gray (
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

      integer i,j,l,i1,j1,ii,jj,kk,ll,km,ifail,CurrentElement

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
      real*8 n1,n2,n3,n4,n5,n6,n7,n8
      real*8 point(3,1),neighbor1(3,1),neighbor2(3,1),neighbor3(3,1)
      real*8 neighbor4(3,1),neighbor5(3,1),neighbor6(3,1),neighbor7(3,1)
      real*8 neighbor8(3,1),normal
      real*8 table1(NElements,9)
      real*8 curve1,curve2,SI
      real*8 NormCurve

      ! Parameters
      !
      real*8 zero,one,two,three,half,third,four,Pi,two_third
      parameter(zero=0.d0,one=1.d0,two=2.d0,three=3.d0,half=0.5d0,
     +     third=1.d0/3.d0,two_third=2.d0/3.d0,four=4.d0,Pi=3.1415926d0)

c      open(unit=80, file='C:\Abaqus_Working_Dir\FeFp\aaMSGS')



      ! pour initial coordinates into the global variable
      if (totaltime.lt. 10.0) then
          do km=1,nblock
             inicoord(nElement(km),1) = coordMp(km,1)
             inicoord(nElement(km),2) = coordMp(km,2)
             inicoord(nElement(km),3) = coordMp(km,3)
          enddo
      end if



          do km=1,nblock
             ! store all the CURRENT coordinates 
             ! at the reduced integration point location 
             StorageOld(nElement(km)-BaseElement,1) = coordMp(km,1)
             StorageOld(nElement(km)-BaseElement,2) = coordMp(km,2)
             StorageOld(nElement(km)-BaseElement,3) = coordMp(km,3)
          enddo




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


         ! calculate the current element 
         ! 
         CurrentElement  = nElement(km)-BaseElement


c         !  store the current coord for all at each element 
c         ! 
c         StorageOld(CurrentElement,1) = coordMp(km,1)
c         StorageOld(CurrentElement,2) = coordMp(km,2)
c         StorageOld(CurrentElement,3) = coordMp(km,3)








         ! take a snapshot of the coorf while loop is done 
         ! 
c         if (CurrentElement .eq. NElements) then
c
             StorageSnap = StorageOld

c         endif 



         
c         ! here I call the search algorithm to find my neibors 
c         ! 
c         call GoSearch(StorageSnap,CurrentElement,
c     +                    n1,n2,n3,n4,n5,n6,n7,n8)



         ! call search table 
         ! 

         call SearchTable1(table1)




c       write(*,*),'elem=',CurrentElement
c       write(*,*),'n1=',n1
c       write(*,*),'n2=',n2
c       write(*,*),'n3=',n3
c       write(*,*),'n4=',n4
c       write(*,*),'n5=',n5
c       write(*,*),'n6=',n6
c       write(*,*),'n7=',n7
c       write(*,*),'n8=',n8

       point(1,1) = StorageSnap(CurrentElement,1) 
       point(2,1) = StorageSnap(CurrentElement,2)
       point(3,1) = StorageSnap(CurrentElement,3)  

       neighbor1(1,1) = StorageSnap(table1(CurrentElement,2),1)
       neighbor1(2,1) = StorageSnap(table1(CurrentElement,2),2)
       neighbor1(3,1) = StorageSnap(table1(CurrentElement,2),3)

       neighbor2(1,1) = StorageSnap(table1(CurrentElement,3),1)
       neighbor2(2,1) = StorageSnap(table1(CurrentElement,3),2)
       neighbor2(3,1) = StorageSnap(table1(CurrentElement,3),3)

       neighbor3(1,1) = StorageSnap(table1(CurrentElement,4),1)
       neighbor3(2,1) = StorageSnap(table1(CurrentElement,4),2)
       neighbor3(3,1) = StorageSnap(table1(CurrentElement,4),3)


       neighbor4(1,1) = StorageSnap(table1(CurrentElement,5),1)
       neighbor4(2,1) = StorageSnap(table1(CurrentElement,5),2)
       neighbor4(3,1) = StorageSnap(table1(CurrentElement,5),3)

       neighbor5(1,1) = StorageSnap(table1(CurrentElement,6),1)
       neighbor5(2,1) = StorageSnap(table1(CurrentElement,6),2)
       neighbor5(3,1) = StorageSnap(table1(CurrentElement,6),3)

       neighbor6(1,1) = StorageSnap(table1(CurrentElement,7),1)
       neighbor6(2,1) = StorageSnap(table1(CurrentElement,7),2)
       neighbor6(3,1) = StorageSnap(table1(CurrentElement,7),3)

       neighbor7(1,1) = StorageSnap(table1(CurrentElement,8),1)
       neighbor7(2,1) = StorageSnap(table1(CurrentElement,8),2)
       neighbor7(3,1) = StorageSnap(table1(CurrentElement,8),3)

       neighbor8(1,1) = StorageSnap(table1(CurrentElement,9),1)
       neighbor8(2,1) = StorageSnap(table1(CurrentElement,9),2)
       neighbor8(3,1) = StorageSnap(table1(CurrentElement,9),3)



c       write(*,*),'time=', totaltime

c         write(*,*),'ele=', CurrentElement



         if (totaltime .lt. 0.1) then 

             curvature = 0.0 

         elseif (totaltime .ge. 0.1) then 

             call calc_normal(point,neighbor1,neighbor2,neighbor3,
     +                       neighbor4,neighbor5,neighbor6,
     +                       neighbor7,neighbor8,normal)

              call calc_curvature(point,normal,neighbor1,neighbor2,
     +                          neighbor3,neighbor4,neighbor5,
     +                          neighbor6,neighbor7,neighbor8,
     +                          curvature,curve1,curve2)

         endif 


        ! herer I calculate the Shape index (SI)
        ! 


         call SIFunc(curve1,curve2,SI)

 









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
     +                     coordx,coordy,coordz,curvature,NormCurve)

         else
            !
            ! Perform explicit time integration procedure
            !
            call integ_gray(matProps,nprops,F_tau,dt,T_tau,thetag_t,thetag_tau,
     +                     coordx,coordy,coordz,curvature,NormCurve)

         endif
         !---------------------------------------------------------------
         

         ! Update state variables
         !
         stateNew(km,1) = thetag_tau ! growth parameter at time tau

         stateNew(km,2) = coordx
         stateNew(km,3) = coordy
         stateNew(km,4) = coordz

         stateNew(km,5) = NormCurve
         stateNew(km,6) = SI


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

      end subroutine vumatXtrArg_gray
***********************************************************************
      subroutine vumatXtrArg_graz (
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

      integer i,j,l,i1,j1,ii,jj,kk,ll,km,ifail,CurrentElement

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
      real*8 n1,n2,n3,n4,n5,n6,n7,n8
      real*8 point(3,1),neighbor1(3,1),neighbor2(3,1),neighbor3(3,1)
      real*8 neighbor4(3,1),neighbor5(3,1),neighbor6(3,1),neighbor7(3,1)
      real*8 neighbor8(3,1),normal
      real*8 table2(NElements,9)
      real*8 curve1,curve2,SI
      real*8 NormCurve

      ! Parameters
      !
      real*8 zero,one,two,three,half,third,four,Pi,two_third
      parameter(zero=0.d0,one=1.d0,two=2.d0,three=3.d0,half=0.5d0,
     +     third=1.d0/3.d0,two_third=2.d0/3.d0,four=4.d0,Pi=3.1415926d0)

c      open(unit=80, file='C:\Abaqus_Working_Dir\FeFp\aaMSGS')



      ! pour initial coordinates into the global variable
      if (totaltime.lt. 10.0) then
          do km=1,nblock
             inicoord(nElement(km),1) = coordMp(km,1)
             inicoord(nElement(km),2) = coordMp(km,2)
             inicoord(nElement(km),3) = coordMp(km,3)
          enddo
      end if



          do km=1,nblock
             ! store all the CURRENT coordinates 
             ! at the reduced integration point location 
             StorageOld2(nElement(km)-BaseElement2,1) = coordMp(km,1)
             StorageOld2(nElement(km)-BaseElement2,2) = coordMp(km,2)
             StorageOld2(nElement(km)-BaseElement2,3) = coordMp(km,3)
          enddo




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


         ! calculate the current element 
         ! 
         CurrentElement  = nElement(km)-BaseElement2


c         !  store the current coord for all at each element 
c         ! 
c         StorageOld(CurrentElement,1) = coordMp(km,1)
c         StorageOld(CurrentElement,2) = coordMp(km,2)
c         StorageOld(CurrentElement,3) = coordMp(km,3)








         ! take a snapshot of the coorf while loop is done 
         ! 
c         if (CurrentElement .eq. NElements) then
c
             StorageSnap2 = StorageOld2

c         endif 



         
c         ! here I call the search algorithm to find my neibors 
c         ! 
c         call GoSearch(StorageSnap,CurrentElement,
c     +                    n1,n2,n3,n4,n5,n6,n7,n8)



         ! call search table 
         ! 

         call SearchTable2(table2)




c       write(*,*),'elem=',CurrentElement
c       write(*,*),'n1=',n1
c       write(*,*),'n2=',n2
c       write(*,*),'n3=',n3
c       write(*,*),'n4=',n4
c       write(*,*),'n5=',n5
c       write(*,*),'n6=',n6
c       write(*,*),'n7=',n7
c       write(*,*),'n8=',n8

       point(1,1) = StorageSnap2(CurrentElement,1) 
       point(2,1) = StorageSnap2(CurrentElement,2)
       point(3,1) = StorageSnap2(CurrentElement,3)  

       neighbor1(1,1) = StorageSnap2(table2(CurrentElement,2),1)
       neighbor1(2,1) = StorageSnap2(table2(CurrentElement,2),2)
       neighbor1(3,1) = StorageSnap2(table2(CurrentElement,2),3)

       neighbor2(1,1) = StorageSnap2(table2(CurrentElement,3),1)
       neighbor2(2,1) = StorageSnap2(table2(CurrentElement,3),2)
       neighbor2(3,1) = StorageSnap2(table2(CurrentElement,3),3)

       neighbor3(1,1) = StorageSnap2(table2(CurrentElement,4),1)
       neighbor3(2,1) = StorageSnap2(table2(CurrentElement,4),2)
       neighbor3(3,1) = StorageSnap2(table2(CurrentElement,4),3)


       neighbor4(1,1) = StorageSnap2(table2(CurrentElement,5),1)
       neighbor4(2,1) = StorageSnap2(table2(CurrentElement,5),2)
       neighbor4(3,1) = StorageSnap2(table2(CurrentElement,5),3)

       neighbor5(1,1) = StorageSnap2(table2(CurrentElement,6),1)
       neighbor5(2,1) = StorageSnap2(table2(CurrentElement,6),2)
       neighbor5(3,1) = StorageSnap2(table2(CurrentElement,6),3)

       neighbor6(1,1) = StorageSnap2(table2(CurrentElement,7),1)
       neighbor6(2,1) = StorageSnap2(table2(CurrentElement,7),2)
       neighbor6(3,1) = StorageSnap2(table2(CurrentElement,7),3)

       neighbor7(1,1) = StorageSnap2(table2(CurrentElement,8),1)
       neighbor7(2,1) = StorageSnap2(table2(CurrentElement,8),2)
       neighbor7(3,1) = StorageSnap2(table2(CurrentElement,8),3)

       neighbor8(1,1) = StorageSnap2(table2(CurrentElement,9),1)
       neighbor8(2,1) = StorageSnap2(table2(CurrentElement,9),2)
       neighbor8(3,1) = StorageSnap2(table2(CurrentElement,9),3)



c       write(*,*),'time=', totaltime

c         write(*,*),'ele=', CurrentElement



         if (totaltime .lt. 0.1) then 

             curvature = 0.0 

         elseif (totaltime .ge. 0.1) then 

             call calc_normal(point,neighbor1,neighbor2,neighbor3,
     +                       neighbor4,neighbor5,neighbor6,
     +                       neighbor7,neighbor8,normal)

              call calc_curvature(point,normal,neighbor1,neighbor2,
     +                          neighbor3,neighbor4,neighbor5,
     +                          neighbor6,neighbor7,neighbor8,
     +                          curvature,curve1,curve2)

         endif 


        ! herer I calculate the Shape index (SI)
        ! 


         call SIFunc(curve1,curve2,SI)

 









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
     +                     coordx,coordy,coordz,curvature,NormCurve)

         else
            !
            ! Perform explicit time integration procedure
            !
            call integ_gray(matProps,nprops,F_tau,dt,T_tau,thetag_t,thetag_tau,
     +                     coordx,coordy,coordz,curvature,NormCurve)

         endif
         !---------------------------------------------------------------
         

         ! Update state variables
         !
         stateNew(km,1) = thetag_tau ! growth parameter at time tau

         stateNew(km,2) = coordx
         stateNew(km,3) = coordy
         stateNew(km,4) = coordz

         stateNew(km,5) = NormCurve
         stateNew(km,6) = SI


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

      end subroutine vumatXtrArg_graz
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
      a0(2,1) = 0.0 
      a0(3,1) = 1.0 
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
     +                      curvature,NormCurve)

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


      ! unit norm 
      a0(1,1) = 0.0
      a0(2,1) = 0.0
      a0(3,1) = 1.0
      

      !!!!!!!!!!!!!!!!!!!!!!!!!!! dummy step !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(dtime.lt.zero) then

      
       thetag_tau = thetag_t


      ! update  kinematics 
      ! 
      ! iso 
c      Fg_tau  = thetag_tau*Iden
      ! area 
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
      ! iso 
c      Fg_tau  = thetag_tau*Iden
      ! area 
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
      subroutine integ_graz(Props,nprops,F_tau,dtime,T_tau,
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

      ! unit norm 
      a0(1,1) = 0.0
      a0(2,1) = 0.0
      a0(3,1) = 1.0      

      !!!!!!!!!!!!!!!!!!!!!!!!!!! dummy step !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(dtime.lt.zero) then

      
       thetag_tau = thetag_t


      ! update  kinematics 
      ! 
      ! iso 
c      Fg_tau  = thetag_tau*Iden
      ! area 
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
      ! iso 
c      Fg_tau  = thetag_tau*Iden
      ! area 
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


      end subroutine integ_graz
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
