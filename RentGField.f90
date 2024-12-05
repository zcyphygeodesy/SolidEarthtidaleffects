      subroutine RentGField(maxn,rln,cnm,snm,gvm,GRS)
!输入maxn－最大阶数
!输入cnm、snm－位系数
!输出grav(4)-剩余扰动位、剩余扰动重力东北天分量
!位系数总数=(maxn+1)**2-3
!cnm/snm最大个数NMAX=(maxn+1)*(maxn+2)/2
!2006年10月5日,章传银
      implicit none
	integer::maxn
	real*8::cnm((maxn+2)**2),snm((maxn+2)**2)
	real*8::cosml,sinml,tn(5),gr,pi,RAD,XYZ(3),ZCFD(3)
	integer n,m,kk,i
	real*8 rln(3),rr,rlat,rlon,gm,ae,gvm(9)
	real*8 tn0,p01,p02,p11,p12,GRS(6)
	real*8::BLH(3),NFD(5),t,u,sp
	real*8,allocatable::pnm(:),dpt1(:),dpt2(:)
!---------------------------------------------------------------------------
   	gm=GRS(1);ae=GRS(2);pi=datan(1.d0)*4.d0; RAD=pi/180.d0
	rr=rln(1);rlat=rln(2);rlon=rln(3)
 	allocate(pnm((maxn+2)**2),dpt1((maxn+2)**2),dpt2((maxn+2)**2))
      t=dsin(rlat*RAD);u=dcos(rlat*RAD)
      call RLAT_BLH(GRS,rln,BLH);call GNormalfd(BLH,NFD,GRS);gr=NFD(2)
      call BelPnmdt(pnm,dpt1,dpt2,maxn,t)
      gvm=0.d0
	do n=2,maxn
        tn=0.d0
	  do m=0,n
          kk=n*(n+1)/2+m
          cosml=dcos(dble(m)*rlon*RAD);sinml=dsin(dble(m)*rlon*RAD)
	    tn(1)=tn(1)+(cnm(kk)*cosml+snm(kk)*sinml)*pnm(kk+1)
	    tn(2)=tn(2)+dble(m)*(cnm(kk)*sinml-snm(kk)*cosml)*pnm(kk+1)
	    tn(3)=tn(3)+(cnm(kk)*cosml+snm(kk)*sinml)*dpt1(kk+1)
	    tn(4)=tn(4)+(cnm(kk)*cosml+snm(kk)*sinml)*dpt2(kk+1)
	    tn(5)=tn(5)+dble(m**2)*(cnm(kk)*cosml+snm(kk)*sinml)*pnm(kk+1)
	  enddo
        sp=dexp(dble(n)*dlog(ae/rr))
	  gvm(1)=gvm(1)+tn(1)*sp
	  gvm(2)=gvm(2)+tn(1)*sp
	  gvm(3)=gvm(3)+tn(1)*sp*dble(n-1)
	  gvm(4)=gvm(4)+tn(1)*sp*dble(n+1)
	  gvm(5)=gvm(5)-tn(3)*sp
	  gvm(6)=gvm(6)-tn(2)*sp
	  gvm(7)=gvm(7)+tn(1)*sp*dble(n+1)*dble(n+2)     
	  gvm(8)=gvm(8)+tn(4)*sp                         
	  gvm(9)=gvm(9)+tn(5)*sp                         
	enddo
	gvm(1)=gm/rr*gvm(1)
	gvm(2)=(gm/rr/gr)*gvm(2)*1.d3
	gvm(3)=(gm/rr**2)*gvm(3)*1.0e8
	gvm(4)=(gm/rr**2)*gvm(4)*1.0e8
	gvm(5)=gm/rr**2/gr*gvm(5)/RAD*36.d5
	gvm(6)=gm/rr**2/gr*gvm(6)/u/RAD*36.d5
	gvm(7)=-gm/rr**3*gvm(7)*1.d14 !0.01mE
	gvm(8)=gm/rr**3*gvm(8)*1.d14 !0.01mE
	gvm(9)=gm/rr**3*gvm(9)/u**2*1.d14 !0.01mE
      deallocate(pnm,dpt1,dpt2)
	return
      end
