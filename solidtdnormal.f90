      subroutine solidtdnormal(rln,maxn,cnm,snm,tdn,GRS)
!由位系数变化计算固体潮直接影响和间接影响之和。标称勒夫数，顾及位移勒夫数纬度相关性
!tdn(13)高程异常,地面重力,扰动重力,地倾斜,垂线偏差,水平,径向,“11”,梯度3分量
	integer::maxn,mode,nn,n,m,kk,i,j
	real*8::cnm(80),snm(80)
	real*8::rln(3),gm,tdn(24),tn(12),gr,GRS(6),cosml,sinml
	real*8::BLH(3),NFD(5),t,u,knm(80),hnm(80),lnm(80)
 	real*8::rr,rlat,rlon,ae,pi,RAD,dk,dh,dl,sinfi2
	real*8,allocatable::pnm(:),dpt1(:),dpt2(:)
!---------------------------------------------------------------------------
    	tdn(1:14)=0.d0;gm=GRS(1);ae=GRS(2)
   	pi=datan(1.d0)*4.d0;RAD=pi/180.d0
 	allocate(pnm((maxn+2)**2),dpt1((maxn+2)**2),dpt2((maxn+2)**2))
      rr=rln(1);rlat=rln(2);rlon=rln(3)
      t=dsin(rlat*RAD);u=dcos(rlat*RAD)
      call RLAT_BLH(GRS,rln,BLH);call GNormalfd(BLH,NFD,GRS);gr=NFD(2)
      call BelPnmdt(pnm,dpt1,dpt2,maxn,t)
      knm=0.d0; hnm=0.d0; lnm=0.d0
	knm(1)=0.30190d0;knm(2)=0.29830d0;knm(3)=0.30102d0
	knm(4)=0.093d0;knm(5)=0.093d0;knm(6)=0.093;knm(7)=0.094d0
	hnm(1)=0.6078d0;hnm(2)=0.6078d0;hnm(3)=0.6078d0
	hnm(4)=0.292d0;hnm(5)=0.292d0;hnm(6)=0.292d0;hnm(7)=0.292d0
	lnm(1)=0.0847d0;lnm(2)=0.0847d0;lnm(3)=0.0847d0
	lnm(4)=0.015d0;lnm(5)=0.015d0;lnm(6)=0.015d0;lnm(7)=0.015d0
      sinfi2=(3.d0*t**2-1.d0)/2.d0
      hnm(1:3)=hnm(1:3)-0.0006d0*sinfi2;lnm(1:3)=lnm(1:3)-0.0006d0*sinfi2
      kk=0
      do m=0,2
         j=m+1
         kk=kk+1;tdn(14+kk)=cnm(j)*(1.d0+knm(j))
         if(m>0)then
           kk=kk+1;tdn(14+kk)=snm(j)*(1.d0+knm(j))
         endif
      enddo
      do n=2,maxn
         tn=0.d0
         do m=0,n
            j=n*(n+1)/2+m-2
            cosml=dcos(rlon*RAD*m);sinml=dsin(rlon*RAD*m)
	      tn(1)=tn(1)+(cnm(j)*cosml+snm(j)*sinml)*pnm(j+3)*(1.d0+knm(j))
	      tn(2)=tn(2)+(cnm(j)*cosml+snm(j)*sinml)*pnm(j+3)*(1.d0+(2.d0*hnm(j)-dble(n+1)*knm(j))/dble(n))
	      tn(3)=tn(3)+(cnm(j)*cosml+snm(j)*sinml)*pnm(j+3)*(1.d0+knm(j))
	      tn(4)=tn(4)+(cnm(j)*cosml+snm(j)*sinml)*dpt1(j+3)*(1.d0+knm(j)-hnm(j))
	      tn(5)=tn(5)+dble(m)*(cnm(j)*sinml-snm(j)*cosml)*pnm(j+3)*(1.d0+knm(j)-hnm(j))
	      tn(6)=tn(6)+(cnm(j)*cosml+snm(j)*sinml)*dpt1(j+3)*(1.d0+knm(j))
	      tn(7)=tn(7)+dble(m)*(cnm(j)*sinml-snm(j)*cosml)*pnm(j+3)*(1.d0+knm(j))
	      tn(8)=tn(8)+dble(m)*(cnm(j)*sinml-snm(j)*cosml)*pnm(j+3)*lnm(j)
	      tn(9)=tn(9)+(cnm(j)*cosml+snm(j)*sinml)*dpt1(j+3)*lnm(j)
	      tn(10)=tn(10)+(cnm(j)*cosml+snm(j)*sinml)*pnm(j+3)*hnm(j)
	      tn(11)=tn(11)+(cnm(j)*cosml+snm(j)*sinml)*dpt2(j+3)*(1.d0+knm(j))
	      tn(12)=tn(12)+dble(m**2)*(cnm(j)*cosml+snm(j)*sinml)*pnm(j+3)*(1.d0+knm(j))
         enddo
         tn=tn*dexp(dble(n)*dlog(ae/rr))
 	   tdn(1)=tdn(1)+tn(1)
	   tdn(2)=tdn(2)+tn(2)*dble(n+1)
	   tdn(3)=tdn(3)+tn(3)*dble(n+1)
	   tdn(4)=tdn(4)+tn(4)
	   tdn(5)=tdn(5)+tn(5)
	   tdn(6)=tdn(6)+tn(6)
	   tdn(7)=tdn(7)+tn(7)
	   tdn(8)=tdn(8)+tn(8)
	   tdn(9)=tdn(9)+tn(9)
 	   tdn(10)=tdn(10)+tn(10)
 	   tdn(12)=tdn(12)+tn(1)*dble(n+1)*dble(n+2)
	   tdn(13)=tdn(13)+tn(11)
 	   tdn(14)=tdn(14)+tn(12)
      enddo
	tdn(1)=tdn(1)*gm/rr/gr*1.d3
	tdn(2)=tdn(2)*gm/rr**2*1.0e8
	tdn(3)=tdn(3)*gm/rr**2*1.0e8
	tdn(4)=tdn(4)*gm/rr**2/gr/RAD*36.d5
	tdn(5)=tdn(5)*gm/rr**2/gr/u/RAD*36.d5
	tdn(6)=tdn(6)*gm/rr**2/gr/RAD*36.d5
	tdn(7)=tdn(7)*gm/rr**2/gr/u/RAD*36.d5
	tdn(8)=-tdn(8)*gm/rr/gr/u*1.d3
 	tdn(9)=-tdn(9)*gm/rr/gr*1.d3
	tdn(10)=tdn(10)*gm/rr/gr*1.d3
      tdn(11)=tdn(10)-tdn(1)
	tdn(12)=-gm/rr**3*tdn(12)*1.d14 !0.01mE
	tdn(13)=gm/rr**3*tdn(13)*1.d14 !0.01mE
	tdn(14)=gm/rr**3*tdn(14)/u**2*1.d14 !0.01mE
      deallocate(pnm,dpt1,dpt2)
 	return
	end
