      subroutine BelPnmdt(pnm,dpt1,dpt2,maxn,t)
      !The normalized associated Legendre functions and thier derivative to θ.
      !计算规格化Pnm及其对sita一、二阶导数t=cos_sita
      !Improved Belikov recursion algorithm for Pnm.规格化Pnm采用改进的Belikov递推算法
      !Non-singular recursive algorithm for derivative of Pnm to θ. 一、二阶导数非奇异递推算法
	implicit none
	integer::maxn,n,m,k,kk,k1,k2,k3,astat(6)
	real*8::pnm((maxn+2)**2),dpt1((maxn+2)**2),dpt2((maxn+2)**2)
	real*8::t,u,a,b,c,d,e
	real*8::anm,bnm,cnm,da
!---------------------------------------------------------------------------
      pnm=0.d0;dpt1=0.d0;dpt2=0.d0
      u=dsqrt(1-t**2);pnm(1)=1.d0
      pnm(2)=dsqrt(3.d0)*t;pnm(3)=dsqrt(3.d0)*u
	do n=1,maxn
        a=dsqrt(2.d0*n+1.d0)/dsqrt(2.d0*n-1.d0)
        b=dsqrt(2.d0*(n-1.d0)*(2.d0*n+1.d0))/dsqrt((2.d0*n-1.d0)*n)
        kk=n*(n+1)/2+1;k1=n*(n-1)/2+1;k2=n*(n-1)/2+2
        pnm(kk)=a*t*pnm(k1)-b*u/2.d0*pnm(k2)
	  do m=1,n
          kk=n*(n+1)/2+m+1;k1=n*(n-1)/2+m+1
          k2=n*(n-1)/2+m+2;k3=n*(n-1)/2+m
          c=a/n*dsqrt(dble(n**2-m**2));d=a/n/2.d0*dsqrt(dble((n-m)*(n-m-1)))
          e=a/n/2.d0*dsqrt(dble((n+m)*(n+m-1)))
          if(m==1)e=e*dsqrt(2.d0)
          pnm(kk)=c*t*pnm(k1)-d*u*pnm(k2)+e*u*pnm(k3)
	  enddo
      enddo
      !derivative of Pnm to θ. 计算一、二阶导数
      dpt1=0.d0;dpt2=0.d0
      dpt1(1)=0.d0;dpt1(2)=-dsqrt(3.d0)*u;dpt1(3)=dsqrt(3.d0)*t
      dpt2(1)=0.d0;dpt2(2)=-dsqrt(3.d0)*t;dpt2(3)=-dsqrt(3.d0)*u
	do n=2,maxn
        kk=n*(n+1)/2+1;k1=n*(n+1)/2+2 !m=0
        dpt1(kk)=-dsqrt(dble(n*(n+1))/2.d0)*pnm(k1)
        kk=n*(n+1)/2+2;k1=n*(n+1)/2+1;k2=n*(n+1)/2+3 !m=1
        anm=dsqrt(dble(2*n))*dsqrt(dble(n+1))/2.d0
        bnm=-dsqrt(dble(n-1))*dsqrt(dble(n+2))/2.d0
	  dpt1(kk)=anm*pnm(k1)+bnm*pnm(k2)
        kk=n*(n+1)/2+1;k1=n*(n+1)/2+3 !m=0
        anm=-dble(n*(n+1))/2.d0
        bnm=dsqrt(dble(n*(n-1)))*dsqrt(dble((n+1)*(n+2)))/dsqrt(8.d0)
	  dpt2(kk)=anm*pnm(kk)+bnm*pnm(k1)
        kk=n*(n+1)/2+2;k1=n*(n+1)/2+4 !m=1
        anm=-dble(2*n*(n+1)+(n-1)*(n+2))/4.d0
        bnm=dsqrt(dble((n-2)*(n-1)))*dsqrt(dble((n+2)*(n+3)))/4.d0
	  dpt2(kk)=anm*pnm(kk)+bnm*pnm(k1)
	  do m=2,n
          kk=n*(n+1)/2+m+1;k1=n*(n+1)/2+m;k2=n*(n+1)/2+m+2
          anm=dsqrt(dble(n+m))*dsqrt(dble(n-m+1))/2.d0
          bnm=-dsqrt(dble(n-m))*dsqrt(dble(n+m+1))/2.d0
	    dpt1(kk)=anm*pnm(k1)+bnm*pnm(k2)
          k1=n*(n+1)/2+m-1;k2=n*(n+1)/2+m+3
          anm=dsqrt(dble((n-m+1)*(n-m+2)))*dsqrt(dble((n+m-1)*(n+m)))/4.d0
          bnm=-dble((n-m+1)*(n+m)+(n-m)*(n+m+1))/4.d0
          cnm=dsqrt(dble((n-m-1)*(n-m)))*dsqrt(dble((n+m+1)*(n+m+2)))/4.d0
	    dpt2(kk)=anm*pnm(k1)+bnm*pnm(kk)+cnm*pnm(k2)
	  enddo
      enddo
904   return
      end
