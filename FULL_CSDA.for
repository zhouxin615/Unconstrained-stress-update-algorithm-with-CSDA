!**********************************************************************************
!* A robust stress update algorithm for elastoplastic models without analytical   *
!* derivation of the consistent tangent operator and loading/unloading estimation *
!* Copyright (C) Lu, D.C., Zhang, Y.N., Zhou, X., Su, C.C., Gao, Z.W., Du, X.L.   *
!* All Rights Reserved.                                                           *
!* Unauthorized copying of this file, via any medium is strictly prohibited       *
!* Copyright (C) 2023 Zhou,X. (zhouxin615@126.com) Lu,D.C. (dechun@bjut.edu.cn)   *
!*                    Zhang, Y.N. (1262477296@qq.com)                             *
!* The code can be used for academic research on the premise that the paper       *
!* "Lu, D.C., Zhang, Y.N., Zhou, X., Su, C.C., Gao, Z.W., Du, X.L.                *
!* A robust stress update algorithm for elastoplastic models without analytical   *
!* derivation of the consistent tangent operator and loading/unloading estimation.*
!* Int J Numer Anal Methods Geomech. 2023;47:1022–1050." is cited                *
!**********************************************************************************
      
      subroutine umat(stress,statev,ddsdde,sse,spd,scd,
     1 rpl,ddsddt,drplde,drpldt,
     2 stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,
     3 ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     4 celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,jstep,kinc)
c
!      include 'ABA_PARAM.INC'
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


      Beta = 1.D-12/2.0d0
      Ftol = sqrt(2.0d0*Beta)   
      tol = 1.D-14             ! tolerance that is a threshold to determine when to use the L'hospital rule to avoid a zero denominator
      rho = 1.D-4
      zeta = 0.1d0
      hsize = 1D-16
      hx = 0.0d0

      delta = 0.0d0
      delta(1:3,1) = 1.0d0

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
      e0 = props(5)
      Mf = sqrt((2.0-u)*Mf*Mf)
!      if (time(2).LT.1.D-12) then
!         e1 = props(5) 
!         call K0e0(e0,s,e1,lambda,kappa,Mf,coords,tol) 
!         statev(6) = e0
!      else    
!         e0 = statev(6)
!      end if   


      cp = (1.0D0+e0)/(lambda-kappa)
      ck = (1.0D0+e0)/kappa
      r = 3.0d0*(1.0d0-2.0d0*nu)/(2.0d0*(1.0d0+nu))
      paras(:,1) = (/Mf,cp,ck,r,u/)
      
      call TenMat1(II,IP,Ivol,Isym)     ! compute the fourth and second order unit tensors in matrix form
      Dw = 0.0d0
      Dw = (3.0d0 - 2.0d0*r)*Ivol+2.0d0*r*Isym  ! elastic stiffness matrix
      
      call pqs(s,p,q)              ! calculate hydrostatic stress, the generalized shear stress, and the deviatoric stress at the previous step.

      if (time(2).LT.1.D-7) then
!         pc =(Mf*Mf*p*p+q*q)/(Mf*Mf*p)    
!         pc1 = dexp((e1-e0-kappa*dlog(p))/(lambda-kappa))
!         pc = max(pc,pc1)
         pc = 110.d0
      else    
         pc = statev(1)                        ! pre-consolidation pressure at the previous step
      end if
      dphi = 0.0d0                             ! the initial value of plastic multiplier
! four variable: hydrostatic pressure,generalized shear stress,
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

      dev = dstra1(1,1)+dstra1(2,1)+dstra1(3,1)   ! the total volume strain
      BKw = p*(exp(ck*dev)-1.0d0)/dev 
      if (abs(dev).LT.tol) then
         BKw = p*ck
      end if

      GKw = BKw*r
    
      De = 0.0d0
      Dw = (3.0d0 - 2.0d0*r)*Ivol+2.0d0*r*Isym 
      De = BKw*Dw  ! elastic stiffness matrix
      st = matmul(De,dstra1) + s               ! calculate the trial stress   
      
      p = (st(1,1)+st(2,1)+st(3,1))/3.0D0
      q = sqrt(((st(1,1)-st(2,1))**2.0D0+(st(2,1)-st(3,1))**2.0D0
     1 +(st(3,1)-st(1,1))**2.0D0+6.0D0*(st(4,1)**2.0D0+st(5,1)**2.0D0
     2 +st(6,1)**2.0D0))/2.0D0)
      x4(1:4,1)=(/p,q,pc,dphi/)
      
      call norms(st,cd)
      cd = max(cd,1.0d0)**3.0d0
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
         F0 = 0.5d0*norm_g*norm_g
         Iteration1 = 0
         
         do while (Iteration1.LT.Mmax) 
            call gfunc4(gi,g4,Dw,paras,x4+alpha*dx,x4old,cd,dstra1,
     1                  sdold,nw,hx,5,tol)
           Dum = matmul(transpose(g4),g4)
           norm_g =  sqrt(Dum(1,1)) 
           F1 = 0.5d0*norm_g*norm_g
              
            if (F1.LT.((1.0d0-2.0d0*rho*alpha)*F0)) then
                go to 192
            else
                alpha = max(zeta*alpha,F0/(F0+F1))
            endif
            
            Iteration1 = Iteration1 + 1
        end do        
        
          
192     continue    

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

      
! The subroutine computing the fourth unit order tensors in matrix form     
	subroutine TenMat1(II,IP,Ivol,Isym)

      double precision II(6,6),IP(6,6),Ivol(6,6),Isym(6,6)
	integer K1,K2

      II = 0.0d0
      IP = 0.0d0
      Ivol = 0.0d0
      Isym = 0.0d0
      
      do  K1 = 1,3
        do  K2 = 1,3
           Ivol(K1,K2) = 1.0d0/3.0d0
           IP(K1,K2) = -1.0d0/3.0d0
        end do
      end do
      do  K1 = 1,3
        Isym(K1,K1) = 1.0d0
        Isym(K1+3,K1+3) = 1.0d0/2.0d0
        II(K1,K1) = 1.0d0
        II(K1+3,K1+3) = 1.0d0
        IP(K1,K1) = 2.0d0/3.0d0
        IP(K1+3,K1+3) = 1.0d0/2.0d0
      end do
      return
      
      end        
 
! The subroutine computing the hydrostatic pressure, the generalized shear stress, and the deviatoric stress
      subroutine pqs(s,p,q)
        double precision s(6,1)
        double precision p,q 
        integer K1
        
        p = (s(1,1)+s(2,1)+s(3,1))/3.0D0
        q = sqrt(((s(1,1)-s(2,1))**2.0D0+(s(2,1)-s(3,1))**2.0D0
     1 +(s(3,1)-s(1,1))**2.0D0+6.0D0*(s(4,1)**2.0D0+s(5,1)**2.0D0
     2 +s(6,1)**2.0D0))/2.0D0)
        return
      end      
      


       subroutine pqs16(s,p,q)
      
        double precision s(6,1)
        double precision p,q 
        integer K1
        p = (s(1,1)+s(2,1)+s(3,1))/3.0D0
        q = sqrt(((s(1,1)-s(2,1))**2.0D0+(s(2,1)-s(3,1))**2.0D0
     1 +(s(3,1)-s(1,1))**2.0D0+6.0D0*(s(4,1)**2.0D0+s(5,1)**2.0D0
     2 +s(6,1)**2.0D0))/2.0D0)
        return
      end    
! 求逆矩阵
      subroutine inverse(A,IA,N)
      implicit none
      integer :: i,j,N
      double precision :: A(N,N), IA(N,N),B(N,N),HIA(N),HA(N),Co_number
         ! 先把IA设定成单位矩阵
      forall(i=1:N,j=1:N,i==j) IA(i,j) = 1.0D0
      forall(i=1:N,j=1:N,i/=j) IA(i,j)=0.0D0
         ! 保存原先的矩阵A, 使用B来计算
         B=A
         ! 把B化成对角线矩阵(除了对角线外，都为0)
      call Upper(B,IA,N) !先把B化成上三角矩阵
      call Lower(B,IA,N) !再把B化成下三角矩阵
         ! 求解
      forall(i=1:N) IA(i,:) = IA(i,:)/B(i,i)
      forall(i=1:N) HIA(i) = sum(abs(IA(i,:)))
      forall(i=1:N) HA(i) = sum(abs(A(i,:)))
      Co_number = dlog((maxval(HA))*(maxval(HIA)))
      if (Co_number.GT.12) then
       !write(*,*) '条件数过大,Co_number=',Co_number    
      endif    

      return
      end subroutine

      ! 求上三角矩阵的子程序
      subroutine Upper(M,S,N)
      implicit none
      integer :: N
      double precision  :: M(N,N)
      double precision  :: S(N,N)
      integer :: i,j
      double precision :: E
      do i=1,N-1
      do j=i+1,N
         E = M(j,i)/M(i,i)
         M(j,i:N)=M(j,i:N)-M(i,i:N)*E
         S(j,:)=S(j,:)-S(i,:)*E
      enddo
      enddo
      return
      end subroutine Upper
      
      !求下三角矩阵的子程序
      subroutine Lower(M,S,N)
      implicit none
      integer :: N
      double precision  :: M(N,N)
      double precision  :: S(N,N)
      integer :: i,j
      double precision :: E
      do i=N,2,-1
      do j=i-1,1,-1
         E = M(j,i)/M(i,i)
         M(j,1:N)=M(j,1:N)-M(i,1:N)*E
         S(j,:)=S(j,:)-S(i,:)*E
      enddo
      enddo
      return
      end subroutine Lower   

	subroutine norms(A,B)
	double precision  A(6,1),B
      B = 0.0d0
      B = sqrt(A(1,1)**2.0d0+A(2,1)**2.0d0+A(3,1)**2.0d0+
     12.0d0*(A(4,1)**2.0d0+A(5,1)**2.0d0+A(6,1)**2.0d0))
      return
      end   

	subroutine normc(A,B)
	complex(kind=8)  A(6,1),B
      B = 0.0d0
      B = sqrt(A(1,1)**2.0d0+A(2,1)**2.0d0+A(3,1)**2.0d0+
     12.0d0*(A(4,1)**2.0d0+A(5,1)**2.0d0+A(6,1)**2.0d0))
      return
      end 
      
      
	subroutine K0e0(e0,s,e1,lambda,kappa,Mf,coords,tol)
	double precision  s(6,1),e0,p,q,p0,e1,lambda,kappa,Mf,tol,Z,
     1 coords(3),VSTRESS,HSTRESS

      call pqs(s,p,q)
      

      Z = coords(3)
      VSTRESS = -Z*6.0d0+50.d0 
	HSTRESS=0.609452015076834d0*VSTRESS  !foundation: K0 = 1-sin(φ)
	p=(VSTRESS+2.0d0*HSTRESS)/3.0d0
	q=VSTRESS-HSTRESS      
      
      e0=e1-lambda*dlog(q*q/Mf/Mf/p+p)+kappa*dlog(q*q/Mf/Mf/p/p+1.0)
      
      return
      end    
      
      
! The subroutine computing the smoothing continuum tangent operator  
	subroutine gfuni(gi,g,Dw,paras,x,xold,cd,dstra1,KK,tol)
	implicit none
	double precision p,q,pc,dphi,Mf,cp,ck,r,u,cd,BKw,dev,deve
      double precision gi,tol,Beta,pold,pcold,f,hpc
      double precision g(8,1),Dw(6,6),paras(5,1),x(8,1),
     & xold(8,1),dstra1(6,1),s(6,1),sold(6,1),gs(6,1),sd(6,1),Dum(1,1),
     & De(6,6),Deii(1,6),Dum6(6,1)
      integer K1,KK
    
      Beta = 1.D-12/2.0d0
      
      Mf = paras(1,1)
      cp = paras(2,1)
      ck = paras(3,1)
      r = paras(4,1)
      u = paras(5,1)
      
      s(1:6,1) = x(1:6,1)
      sold(1:6,1) = xold(1:6,1)
      pc = x(7,1)
      pcold = xold(7,1)
      pold = (xold(1,1)+xold(2,1)+xold(3,1))/3.0d0
      dphi = x(8,1)
      call pqs(s,p,q)              ! calculate hydrostatic stress, the generalized shear stress, and the deviatoric stress at the previous step.
      gs = 0.0d0
      call gsfun(gs,paras,x,tol)
      
      dev = dstra1(1,1)+dstra1(2,1)+dstra1(3,1)
      if ((q.LE.tol).or.(p.LE.tol).or.(u.EQ.1)) then
            hpc = (2.0d0*p-pc)
            deve = dev-dphi*hpc 
      else 
            hpc=q*q*p**(-u)/(Mf*Mf*gamma(1.0d0-u)) + 2.0d0*p**(2.0d0-u)/
     1 (gamma(3.0d0-u)) - pc*p**(1.0d0-u)/(gamma(2.0d0-u))          
            deve = dev-dphi*hpc
      endif      
      
      BKw = pold*(dexp(ck*deve)-1.0d0)/deve 
      if (abs(deve).LT.tol) then
         BKw = pold*ck
      end if
      
      De = 0.0d0
      De = BKw*Dw  ! elastic stiffness matrix

      gi = 0.0d0
      g = 0.0d0
      Deii = 0.0d0
      Dum  = 0.0d0
      Dum6 = 0.0d0
      
      if (KK.LE.6) then
         Deii(1,1:6) = De(KK,1:6)
         Dum = matmul(Deii,(dstra1-dphi*gs))
         gi = s(KK,1) - (sold(KK,1)+Dum(1,1))
      elseif (KK.EQ.7) then
         gi = pc - pcold*dexp(cp*dphi*hpc) 
      elseif (KK.EQ.8) then
         f = q*q/(Mf*Mf) + p*(p-pc)
         gi = sqrt(cd*cd*dphi*dphi + f*f + 2.0d0*Beta) - cd*dphi + f 
      elseif (KK.EQ.9) then
         Dum6 =  matmul(De,(dstra1-dphi*gs))
         g(1:6,1)=s(1:6,1)-(sold(1:6,1)+Dum6(1:6,1))
         g(7,1) = pc - pcold*dexp(cp*dphi*hpc) 
         f = q*q/(Mf*Mf) + p*(p-pc)
         g(8,1) = sqrt(cd*cd*dphi*dphi + f*f + 2.0d0*Beta)-cd*dphi+ f   
      else    
         write(*,*) "error&KK value, the allowable range of KK is [1,9]"
      end if      
      
      return
      
      end      



! The subroutine computing the smoothing continuum tangent operator  
	subroutine gfunc(gi,g,Dw,paras,x,xold,cd,dstra1,hx,KK,tol)
	implicit none
      integer, parameter :: rp = kind(1.d0)
      complex(kind=8) p,q,pc,dphi,BKw,deve,pold,pcold,f,hpc,fpu,fqu
      complex(kind=8) s(6,1),sold(6,1),gs(6,1),De(6,6),sd(6,1),gic,
     & Deii(1,6),gc(8,1),xc(8,1),Dum(1,1)
	double precision Mf,cp,ck,r,u,cd,dev,tol,Beta,gi,qtol,ptol,devetol
      double precision g(8,1),Dw(6,6),paras(5,1),x(8,1),
     & xold(8,1),dstra1(6,1),hx(8,1),Dum6(6,1)

      integer K1,KK    
       Beta = 1.D-12/2.0d0
      
      
      Mf = paras(1,1)
      cp = paras(2,1)
      ck = paras(3,1)
      r = paras(4,1)
      u = paras(5,1)
 
      xc(1:8,1) = cmplx(x(1:8,1), hx(1:8,1),kind=8)      
      s(1:6,1) = xc(1:6,1)
      pc = xc(7,1)
      dphi = xc(8,1)

     
      p = (s(1,1)+s(2,1)+s(3,1))/3.0D0
      q = sqrt(((s(1,1)-s(2,1))**2.0D0+(s(2,1)-s(3,1))**2.0D0
     1 +(s(3,1)-s(1,1))**2.0D0+6.0D0*(s(4,1)**2.0D0+s(5,1)**2.0D0
     2 +s(6,1)**2.0D0))/2.0D0)
      sd(1:3,1) = s(1:3,1) - p
      sd(4:6,1) = s(4:6,1)
      dev = dstra1(1,1)+dstra1(2,1)+dstra1(3,1)

      gs = 0.0d0
      
      qtol=sqrt(real(q)**2.D0+imag(q)**2.D0)
      ptol=sqrt(real(p)**2.D0+imag(p)**2.D0)
      
      if ((real(q).LE.tol).or.(real(p).LE.tol).or.(abs(u-1).LE.tol)
     1   .or.(u.EQ.1))then
        do  K1=1,3
            gs(K1,1) = 3.0d0*sd(K1,1)/(Mf*Mf)+(2.0d0*p-pc)/3.0d0
            gs(K1+3,1) = 6.0d0*sd(K1+3,1)/(Mf*Mf)
        end do        
        hpc = (2.0d0*p-pc)
      else
        fpu = q*q*p**(-u)/(Mf*Mf*gamma(1.0d0-u)) + 2.0d0*p**(2.0d0-u)/
     1 (gamma(3.0d0-u)) - pc*p**(1.0d0-u)/(gamma(2.0d0-u))
        fqu = 2.0d0*q**(2.0d0-u)/(Mf*Mf*gamma(3.0d0-u)) + 
     1  p*(p-pc)*q**(-u)/gamma(1.0d0-u)
        do  K1=1,3
            gs(K1,1) = fpu/3.0d0 + 3.0d0*fqu*sd(K1,1)/(2.0d0*q)
            gs(K1+3,1) = 6.0d0*fqu*sd(K1+3,1)/(2.0d0*q)
        end do
        hpc = fpu          
      endif      
      deve = dev-dphi*hpc
      sold(1:6,1) = cmplx(xold(1:6,1), 0.0d0,kind=8)
      pcold = cmplx(xold(7,1), 0.0d0,kind=8)
      pold = cmplx((xold(1,1)+xold(2,1)+xold(3,1))/3.0d0, 0.0d0,kind=8)

      BKw = pold*(exp(ck*deve)-1.0d0)/deve 
      
      devetol=sqrt((real(deve))**2.0d0+imag(deve)**2.0d0)
      if (abs(real(deve)).LT.tol) then
         BKw = pold*ck
      end if
      
      De = 0.0d0
      De = BKw*Dw  ! elastic stiffness matrix

      
      Deii = (0.0d0,0.0d0)  
      Dum  = (0.0d0,0.0d0)  
      
      gi  = 0.0d0  !(1,1)复数,输出g的虚部
      gic = (0.0d0,0.0d0) 
      gc = (0.0d0,0.0d0)  
      g = 0.0d0
      Dum6 = 0.0d0      

      
      if (KK.LE.6) then
         Deii(1,1:6) = De(KK,1:6)
         Dum = matmul(Deii,(dstra1-dphi*gs))
         gic = s(KK,1) - (sold(KK,1)+Dum(1,1))
         gi = imag(gic)                          
      elseif (KK.EQ.7) then
         gic = pc - pcold*exp(cp*dphi*hpc) 
         gi = imag(gic)
      elseif (KK.EQ.8) then
         f = q*q/(Mf*Mf) + p*(p-pc)
         gic = sqrt(cd*cd*dphi*dphi + f*f + 2.0d0*Beta) - cd*dphi + f 
         gi = imag(gic)
      elseif (KK.EQ.9) then
         Dum6 =  matmul(De,(dstra1-dphi*gs))
         gc(1:6,1)=s(1:6,1)-(sold(1:6,1)+Dum6(1:6,1))
         gc(7,1) = pc - pcold*exp(cp*dphi*hpc) 
         f = q*q/(Mf*Mf) + p*(p-pc)
         gc(8,1) = sqrt(cd*cd*dphi*dphi + f*f + 2.0d0*Beta)-cd*dphi+ f
         g = dble(gc)                     

      else    
         write(*,*) "error&KK value, the allowable range of KK is [1,9]"
      end if      
      
      return
      
      end  



! The subroutine computing the smoothing continuum tangent operator  
	subroutine gsfun(gs,paras,x,tol)

	implicit none
	double precision p,q,pc,Mf,u,tol,fpu,fqu
      double precision gs(6,1),paras(5,1),x(8,1),sd(6,1),s(6,1),
     & delta(6,1)
      integer K1
      Mf = paras(1,1)
      u = paras(5,1)
	s(1:6,1) = x(1:6,1)

      call pqs(s,p,q)
      
      pc = x(7,1)
      delta = 0.0d0
      delta(1:3,1) = 1.0d0
      sd = s - p*delta
      
      if ((real(q).LE.tol).or.(real(p).LE.tol).or.(abs(u-1.0d0).LE.tol)
     1   .or.(u.EQ.1)) then
        do  K1=1,3  
            gs(K1,1) = 3.0d0*sd(K1,1)/(Mf*Mf)+(2.0d0*p-pc)/3.0d0
            gs(K1+3,1) = 6.0d0*sd(K1+3,1)/(Mf*Mf)
        end do
      else
        fpu = q*q*p**(-u)/(Mf*Mf*gamma(1.0d0-u)) + 2.0d0*p**(2.0d0-u)/
     1 (gamma(3.0d0-u)) - pc*p**(1.0d0-u)/(gamma(2.0d0-u))
        fqu = 2.0d0*q**(2.0d0-u)/(Mf*Mf*gamma(3.0d0-u)) + 
     1  p*(p-pc)*q**(-u)/gamma(1.0d0-u)
        do  K1=1,3
            gs(K1,1) = fpu/3.0d0 + 3.0d0*fqu*sd(K1,1)/(2.0d0*q)
            gs(K1+3,1) = 6.0d0*fqu*sd(K1+3,1)/(2.0d0*q)
        end do
      endif      
      
      return
      
      end        
      
      
! The subroutine computing the smoothing continuum tangent operator  
	subroutine Jacob(Jg,Dw,paras,x,xold,cd,dstra1,hsize,tol)
	implicit none
	double precision hsize,gi,cd,tol
      double precision Jg(8,8),Dw(6,6),paras(5,1),x(8,1),
     1 xold(8,1),dstra1(6,1),hx(8,1),g(8,1)
      integer K1,K2
      
      
      Jg = 0.0d0
      gi = 0.0d0
      
      do  K1 = 1,8
        do  K2 = 1,8
           hx = 0.0d0
           hx(K2,1) = hsize
           call gfunc(gi,g,Dw,paras,x,xold,cd,dstra1,hx,K1,tol)
           Jg(K1,K2) = gi/hsize
        end do

      end do      

      return
      
      end
      
      
      
! The subroutine computing the smoothing continuum tangent operator  
	subroutine CTONum(DDSDDE,Dw,paras,x,xold,cd,dstra1,hsize,tol)
	implicit none
	double precision p,q,pc,dphi,Mf,cp,ck,r,u,cd,BKw,dev,deve
      double precision tol,Beta,pold,pcold,BK,dBK,hsize
      double precision paras(5,1),x(8,1),Dum6(6,1),delta(6,1),
     & xold(8,1),dstra1(6,1),s(6,1),sold(6,1),gs(6,1),sd(6,1),
     & Dw(6,6),deltaT(1,6),Ke(6,6),IJg(8,8),Jg(8,8),CTO(6,6),DDSDDE(6,6)
      integer K1,K2
 
      Beta = 1.D-12/2.0d0
      deltaT = 0.0d0
      deltaT(1,1:3) = 1.0d0
      
      
      Mf = paras(1,1)
      cp = paras(2,1)
      ck = paras(3,1)
      r = paras(4,1)
      u = paras(5,1)
      
      s(1:6,1) = x(1:6,1)
      sold(1:6,1) = xold(1:6,1)
      pc = x(7,1)
      pcold = xold(7,1)
      pold = (xold(1,1)+xold(2,1)+xold(3,1))/3.0d0
      dphi = x(8,1)
      call pqs16(s,p,q)            
      gs = 0.0d0
      call gsfun(gs,paras,x,tol)
      
      dev = dstra1(1,1)+dstra1(2,1)+dstra1(3,1)       
      
      if ((real(q).LE.tol).or.(real(p).LE.tol).or.(u.EQ.1)) then
            deve = dev-dphi*(2.0d0*p-pc) 
      else    
            deve = dev-dphi*(q*q/(Mf*Mf)/p**u/gamma(1.0d0-u)+2.0d0*p**
     1   (2.0d0-u)/gamma(3.0d0-u)-pc*p**(1.0d0-u)/gamma(2.0d0-u))
      endif      
      
      
      BKw = pold*(exp(ck*deve)-1.0d0)/deve
      BK = ck*pold*exp(ck*deve)
      dBK = (BK-BKw)/deve
      if (abs(deve).LT.tol) then
         BKw = pold*ck
         dBK = pold*ck*ck
      end if
      
      Ke =0.0d0
      Dum6 = 0.0d0
      Dum6 = matmul(Dw,(dstra1-dphi*gs))
      Ke = matmul(Dum6,deltaT)
      Ke = -(BKw*Dw + dBK*Ke)             
 
      call Jacob(Jg,Dw,paras,x,xold,cd,dstra1,hsize,tol)

      call inverse(Jg,IJg,8)

      CTO = 0.0d0
      CTO(1:6,1:6) = IJg(1:6,1:6)
      DDSDDE = -matmul(CTO,Ke)
      return
      
      end

! The subroutine computing the smoothing continuum tangent operator  
	subroutine gfunc4(gi,g,Dw,paras,x4,x4old,cd,dstra1,sdold,nw,hx,KK,tol)
	implicit none

      complex(kind=8) p,q,pc,dphi,BKw,deve,f,hpc,fpu,eta,
     &                gic,dl,norm_sdt,norm_sd,t,t1,t2,Gkw
      complex(kind=8) sd(6,1),gc(4,1),x4c(4,1),sdt(6,1),nwc(6,1) 
	double precision Mf,cp,ck,r,u,cd,dev,tol,Beta,gi,pold,pcold,qtol,ptol,
     &                 devetol
      double precision g(4,1),Dw(6,6),paras(5,1),x4(4,1),delta(6,1),
     & dstra1(6,1),hx(4,1),sdold(6,1),dgamma(6,1),x4old(4,1),nw(6,1)
      integer K1,KK
      Beta = 1.D-12/2.0d0
      Mf = paras(1,1)
      cp = paras(2,1)
      ck = paras(3,1)
      r = paras(4,1)
      u = paras(5,1)
      delta = 0.0d0
      delta(1:3,1) = 1.0d0
      pold = x4old(1,1)
      pcold = x4old(3,1)

      x4c(1:4,1) = cmplx(x4(1:4,1), hx(1:4,1),kind=8)
      
      p = x4c(1,1)
      q = x4c(2,1)
      pc = x4c(3,1)
      dphi = x4c(4,1)

      qtol=sqrt(real(q)**2.D0+imag(q)**2.D0)
      ptol=sqrt(real(p)**2.D0+imag(p)**2.D0)
      
      dev = dstra1(1,1)+dstra1(2,1)+dstra1(3,1)
      if ((real(q).LE.tol).or.(real(p).LE.tol).or.(abs(u-1).LE.tol)
     1    .or.(u.EQ.1)) then
          hpc = (2.0d0*p-pc)

      else
       fpu = q*q*p**(-u)/(Mf*Mf*gamma(1.0d0-u)) + 2.0d0*p**(2.0d0-u)/
     1 (gamma(3.0d0-u)) - pc*p**(1.0d0-u)/(gamma(2.0d0-u))

       hpc = fpu  
      end if
      
      deve = dev-dphi*hpc

      devetol=sqrt((real(deve))**2.0d0+imag(deve)**2.0d0)
      
      BKw = pold*(exp(ck*deve)-1.0d0)/deve 
      if (abs(real(deve)).LT.tol) then
         BKw = pold*ck
        
      end if
        Gkw=BKw*r
      
      gi  = 0.0d0  
      gic = (0.0d0,0.0d0) 
      gc = (0.0d0,0.0d0)  
      g = 0.0d0  

      
      if ((real(q).LE.tol).or.(real(p).LE.tol).or.(abs(u-1).LE.tol)
     1   .or.(u.EQ.1))then
            eta = 1.0d0/(1.0d0+6.0d0*Gkw*dphi/(Mf*Mf))
      else    
            dl = q**(1.0d0-u)/Mf**2.0d0/gamma(3.0d0-u)+
     1   p*(p-pc)*q**(-u-1.0d0)/2.0d0/gamma(1.0d0-u)
            eta = 1.0d0/(1.0d0+6.0d0*Gkw*dphi*dl)
      end if 
      
      dgamma(1:3,1) = dstra1(1:3,1) - delta(1:3,1)*dev/3.0d0
      dgamma(4:6,1) = dstra1(4:6,1)/2.0D0
      sdt = sdold+2.0d0*Gkw*dgamma

      call normc(sdt,norm_sdt)


      sd = eta*(sdold+2.0d0*Gkw*dgamma)   
      call normc(sd,norm_sd)
         nwc = 0.0d0
         nwc = sd/norm_sd                   ! the direction of deviatoric stress tensor 
         if (abs(real(norm_sd).LT.tol)) then   
            nwc = delta/sqrt(3.0d0)
         end if 
         nw=dble(nwc)

         f = q*q/(Mf*Mf) + p*(p-pc)

      if ((real(q).LE.tol).or.(real(p).LE.tol).or.(abs(u-1).LE.tol)
     1     .or.(u.EQ.1)) then
         if (KK.EQ.1) then
            gic = (p - pold*exp(ck*deve))
            gi = imag(gic)                          
         elseif (KK.EQ.2) then
            gic = (q - sqrt(3.0d0/2.0d0)*eta*norm_sdt)   
            gi = imag(gic)
         elseif (KK.EQ.3) then
            gic = (pc - pcold*exp(cp*dphi*(2.0d0*p-pc))) 
            gi = imag(gic)
         elseif (KK.EQ.4) then 
            gic = sqrt(cd*cd*dphi*dphi + f*f + 2.0d0*Beta)-cd*dphi+ f
            gi = imag(gic)
         elseif (KK.EQ.5) then 
            gc(1,1) = (p - pold*exp(ck*deve))       
            gc(2,1) = (q - sqrt(3.0d0/2.0d0)*eta*norm_sdt)    
            gc(3,1) = (pc - pcold*exp(cp*dphi*(2.0d0*p-pc)))  
            gc(4,1) = sqrt(cd*cd*dphi*dphi + f*f + 2.0d0*Beta)-cd*dphi+f
            g = dble(gc)
         else     
             write(*,*) "error"
         end if
      else
          if (KK.EQ.1) then
            gic = (p - pold*exp(ck*deve))
            gi = imag(gic)                          
          elseif (KK.EQ.2) then
            gic = (q - sqrt(3.0d0/2.0d0)*eta*norm_sdt)
            gi = imag(gic)

         elseif (KK.EQ.3) then
            gic =(pc-pcold*exp(cp*dphi*(q*q/(Mf*Mf)/p**u/gamma(1.0d0-u)
     1           +2.0d0*p**(2.0d0-u)/gamma(3.0d0-u)
     2            -pc*p**(1.0d0-u)/gamma(2.0d0-u))))
            gi = imag(gic)
         elseif (KK.EQ.4) then 
            gic = sqrt(cd*cd*dphi*dphi + f*f + 2.0d0*Beta) - cd*dphi + f
            gi = imag(gic)
         elseif (KK.EQ.5) then 
            gc(1,1) = (p - pold*exp(ck*deve))
            gc(2,1) = (q - sqrt(3.0d0/2.0d0)*eta*norm_sdt) 
            gc(3,1) = (pc-pcold*exp(cp*dphi*(q*q/(Mf*Mf)/p**u
     1               /gamma(1.0d0-u)+2.0d0*p**(2.0d0-u)/gamma(3.0d0-u)
     2               -pc*p**(1.0d0-u)/gamma(2.0d0-u))))
            gc(4,1) = sqrt(cd*cd*dphi*dphi + f*f + 2.0d0*Beta)-cd*dphi+f
            g = dble(gc)
         else     
             write(*,*) "error"
         end if
      end  if
         
      
      return
      
      end  


! The subroutine computing the smoothing continuum tangent operator  
	subroutine Jacob4(Jg,Dw,paras,x4,x4old,cd,dstra1,sdold,hsize,tol)
	implicit none
	double precision hsize,gi,cd,tol
      double precision Jg(4,4),Dw(6,6),paras(5,1),x4(4,1),
     1 x4old(4,1),dstra1(6,1),hx(4,1),g(4,1),sdold(6,1),nw(6,1)
      integer K1,K2
      
      
      Jg = 0.0d0
      gi = 0.0d0
      
      do  K1 = 1,4
        do  K2 = 1,4
           hx = 0.0d0
           hx(K2,1) = hsize
           call gfunc4(gi,g,Dw,paras,x4,x4old,cd,dstra1,sdold,nw,hx,K1,
     1                 tol)
           Jg(K1,K2) = gi/hsize

        end do

      end do      

      return
      
      end