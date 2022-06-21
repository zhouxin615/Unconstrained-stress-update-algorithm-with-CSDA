!**********************************************************************************
!* unconstrained stress update algorithm with FullCSDA                            *
!* Copyright (C)  Lu, D.C., Zhou, X., Zhang, Y.N.                                 *
!* All Rights Reserved.                                                           *
!* Unauthorized copying of this file, via any medium is strictly prohibited       *
!* Copyright (C) 2021 Zhou,X. (zhouxin615@126.com) Lu,D.C. (dechun@bjut.edu.cn)   *
!* The code can be used for academic research on the premise that the paper       *
!* "Lu, D.C., Zhou, X., Zhang, Y.N."                                              *
!* An unconstrained stress updating algorithm with FullCSDA for elastoplastic     *
!*  soil models.                                                                  *     
!**********************************************************************************

!**********************************************************************************
!*                              Very important!!!                                 * 
!* The subroutine K0e0 needs to be used when considering the change of initial    *
!* void ratio e0 along the depth direction, e.g., strip footing example and       *
!* pile foundation example. e0 is a constant in the Cylinder example.             *
!********************************************************************************** 
      
      
      
      subroutine umat(stress,statev,ddsdde,sse,spd,scd,
     1 rpl,ddsddt,drplde,drpldt,
     2 stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,
     3 ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     4 celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,jstep,kinc)
c
! include 'ABA_PARAM.INC'
      CHARACTER*80 CMNAME
      integer ntens, ndi, nshr, nstatv, nprops, noel, npt,
     & layer, kspt, kstep, kinc, inittension
c
      double precision stress(ntens), statev(15),
     &  ddsdde(ntens,ntens), ddsddt(ntens), drplde(ntens),jstep(4),
     &  stran(ntens), dstran(ntens), time(2), predef(1), dpred(1),
     &  props(nprops), coords(3), drot(3,3), dfgrd0(3,3), dfgrd1(3,3)
      double precision sse, spd, scd, rpl, drpldt, dtime, temp, 
     &  dtemp, pnewdt, celent
! This code is only available for 3D stress state, so NTENS=6
      double precision lambda,kappa,nu,e0,GK,norm_g,tol,e1,
     & p,q,pc,dphi,Mf,cp,ck,r,Beta,Ftol,cd,f,dev,BKw,GKw,
     & alpha,alphaj,F0,F1,rho,zeta,pc1,gi,hsize,u

      double precision s(6,1),dstra1(6,1),De(6,6),Jg(8,8),g(8,1),
     & x(8,1),xold(8,1),dx(4,1),II(6,6),IP(6,6),Ivol(6,6),Isym(6,6),
     & st(6,1),dum(1,1),paras(5,1),delta(6,1),IJg(8,8),hx(8,1),Dw(6,6),
     1 dstra(6,1),x4(4,1),x4old(4,1),sdold(6,1),nw(6,1),Jg4(4,4),
     2 IJg4(4,4),g4(4,1)

      integer K1,K2,Iteration,Iteration1
      DATA Maxit,Mmax /100, 10/
!   The Material parameters of MCC model --------------------------------------------------------------
!   props(1) - lambda   ! compression index
!   props(2) - kappa    ! swell index
!   props(3) - Mf       ! critical state stress ratio
!   props(4) - nu       ! Poisson's ratio
!   props(5) - e1       ! initial void ratio  
!   props(6) - u        ! fractional order
 
      Beta = 1.D-12/2.0d0        ! smoothing parameter
      Ftol = sqrt(2.0d0*Beta)    ! error tolerance for the NMTR method
      tol = 1.D-14               ! tolerance that is a threshold to determine when to use the L'hospital rule to avoid a zero denominator
      rho = 1.D-4
      zeta = 0.1d0
      hsize = 1D-16
      hx = 0.0d0

      delta = 0.0d0
      delta(1:3,1) = 1.0d0
!   initial values of stress, the strain increment and internal variables
      do K1=1,6
	  s(K1,1) = -stress(K1)
        dstra1(K1,1) = -dstran(K1)
      end do     
      
      do K1=1,3
        dstra(K1,1) = -dstran(K1)
      end do
      do K1=4,6
        dstra(K1,1)= -dstran(K1)/2.0d0
      end do

      lambda = props(1)
      kappa = props(2)
      Mf = props(3)
      nu = props(4)
      e1 = props(5)
      u = props(6)
      
      if (time(2).LT.1.D-12) then
         e1 = props(5) 
         call K0e0(e0,s,e1,lambda,kappa,Mf,coords,tol) 
         statev(6) = e0
      else    
         e0 = statev(6)
      end if   

      cp = (1.0D0+e0)/(lambda-kappa)
      ck = (1.0D0+e0)/kappa
      r = 3.0d0*(1.0d0-2.0d0*nu)/(2.0d0*(1.0d0+nu))
      paras(:,1) = (/Mf,cp,ck,r,u/)
      
      call TenMat1(II,IP,Ivol,Isym)     ! compute the fourth and second order unit tensors in matrix form
      Dw = 0.0d0
      Dw = (3.0d0 - 2.0d0*r)*Ivol+2.0d0*r*Isym  ! elastic stiffness matrix
      
      call pqs(s,p,q)              ! calculate hydrostatic stress, the generalized shear stress, and the deviatoric stress at the previous step.
      if (time(2).LT.1.D-7) then
         pc =(Mf*Mf*p*p+q*q)/(Mf*Mf*p)    
         pc1 = dexp((e1-e0-kappa*dlog(p))/(lambda-kappa))
         pc = max(pc,pc1)
      else    
         pc = statev(1)                        ! pre-consolidation pressure at the previous step
      end if
      
      dphi = 0.0d0                             ! the initial value of plastic multiplier
      xold(1:6,1) = s(1:6,1)
      xold(7,1) = pc
      xold(8,1) = 0.0d0
      x = xold 
      
      p = (s(1,1)+s(2,1)+s(3,1))/3.0D0
      q = sqrt(((s(1,1)-s(2,1))**2.0D0+(s(2,1)-s(3,1))**2.0D0
     1 +(s(3,1)-s(1,1))**2.0D0+6.0D0*(s(4,1)**2.0D0+s(5,1)**2.0D0
     2 +s(6,1)**2.0D0))/2.0D0)
      x4old(1:4,1)=(/p,q,pc,dphi/)
      
      sdold(1:3,1) = s(1:3,1) - p
      sdold(4:6,1) = s(4:6,1)
! four variable: the secant volume modulus(BKw), the secant shear modulus (GKw)
      dev = dstra1(1,1)+dstra1(2,1)+dstra1(3,1)   ! the total volume strain
      BKw = p*(exp(ck*dev)-1.0d0)/dev 
      if (abs(dev).LT.tol) then
         BKw = p*ck
      end if

      GKw = BKw*r
    
      De = 0.0d0
      Dw = (3.0d0 - 2.0d0*r)*Ivol+2.0d0*r*Isym 
      De = BKw*Dw                                  ! elastic stiffness matrix
      st = matmul(De,dstra1) + s                   ! calculate the trial stress   
      
      p = (st(1,1)+st(2,1)+st(3,1))/3.0D0
      q = sqrt(((st(1,1)-st(2,1))**2.0D0+(st(2,1)-st(3,1))**2.0D0
     1 +(st(3,1)-st(1,1))**2.0D0+6.0D0*(st(4,1)**2.0D0+st(5,1)**2.0D0
     2 +st(6,1)**2.0D0))/2.0D0)
      x4(1:4,1)=(/p,q,pc,dphi/)
      
      call norms(st,cd)
      cd = max(cd,1.0d0)**3.0d0   ! dimensional parameter in the smoothing function
      x(1:6,1) = st(1:6,1)
      
      Iteration = 0.0d0
       nw=0
       call gfunc4(gi,g4,Dw,paras,x4,x4old,cd,dstra1,sdold,nw,hx,5,tol)

      Dum = matmul(transpose(g4),g4)
      norm_g =  sqrt(Dum(1,1))

      do while ((norm_g.GT.Ftol).and.(Iteration .LT. Maxit))
         Iteration =Iteration + 1
         call Jacob4(Jg4,Dw,paras,x4,x4old,cd,dstra1,sdold,hsize,tol)

         call inverse(Jg4,IJg4,4)
         dx = -matmul(IJg4,g4)

         
         alpha = 1.0d0
         F0 = 0.5d0*norm_g
         Iteration1 = 0
         
         do while (Iteration1.LT.Mmax) 
            call gfunc4(gi,g4,Dw,paras,x4+alpha*dx,x4old,cd,dstra1,
     1                  sdold,nw,hx,5,tol)
           Dum = matmul(transpose(g4),g4)
           norm_g =  sqrt(Dum(1,1)) 
           F1 = 0.5d0*norm_g
              
            if (F1.LT.((1.0d0-2.0d0*rho*alpha)*F0)) then
                go to 192
            else
                alpha = max(zeta*alpha,F0/(F0+F1))
            end if
            
            Iteration1 = Iteration1 + 1
         end do   
         
192      continue   
         
        x4 = x4 + alpha*dx         
         
        call gfunc4(gi,g4,Dw,paras,x4,x4old,cd,dstra1,sdold,nw,hx,5,tol)
         
        Dum = matmul(transpose(g4),g4)
        norm_g =  sqrt(Dum(1,1)) 
      end do

! output the four independent variables in the residual function
      p=x4(1,1)
      q=x4(2,1)
      pc = x4(3,1)
      dphi = x4(4,1)
      s = p*delta + sqrt(2.0d0/3.0d0)*q*nw

! compute the smoothing consistent tangent operator
      x(1:6,1) = s(1:6,1)
      x(7,1)=pc
      x(8,1)=dphi
      call CTONum(DDSDDE,Dw,paras,x,xold,cd,dstra1,hsize,tol)

! The time step is required to decrease when the solution of NMTR is divergent   
      if ((Iteration .GE. 100).or.isnan(norm_g).or.
     2isnan(norm2(DDSDDE)).or.isnan(norm2(s))) then
         PNEWDT = 0.25d0
         DDSDDE = -1.D60*II
         s = -1.D10
         pc = -1.D10         
      end if
       
      do K1=1,6
         stress(K1) = -s(K1,1)
      end do
      
!   statev(1) - pc        ! pre-consolidation pressure 
!   statev(2) - Iteration ! the iteration number of NMTR
!   statev(3) - norm_g    ! 2-norm of residual function vector
!   statev(4) - f         ! the value of yield function 
!   statev(5) - statev(5) ! Counter 
      call pqs(s,p,q)              ! calculate hydrostatic stress, the generalized shear stress, and the deviatoric stress at the previous step.      
      f = q*q/(Mf*Mf) + p*(p-pc)
      statev(1) = pc
      statev(2) = Iteration
      statev(3) = norm_g
      statev(4) = f
      statev(5) = g(8,1)
      statev(7) = statev(7) + 1
      statev(8) = dphi
      
      end