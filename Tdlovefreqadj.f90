      subroutine Tdlovefreqadj(rln,cnm,snm,tdn,GRS)
!计算勒夫数频率相关性校正
!tdn(24)高程异常,地面重力,扰动重力,地倾斜,垂线偏差,水平,径向,“11”,梯度3分量
	integer::maxn,mode,nn,n,m,kk,i,j
	real*8::cnm(80),snm(80),ch(3),sh(3),cl(3),sl(3)
	real*8::rln(3),gm,tdn(24),tn(12),gr,GRS(6),cosml,sinml
	real*8::BLH(3),NFD(5),t,u,knm(80),hnm(80),lnm(80)
 	real*8::rr,rlat,rlon,ae,pi,RAD,dk,dh,dl,pnm(40),dpt1(40),dpt2(40)
!---------------------------------------------------------------------------
    	tdn(1:14)=0.d0;gm=GRS(1);ae=GRS(2);pi=datan(1.d0)*4.d0;RAD=pi/180.d0
      ch(1:3)=cnm(4:6);sh(1:3)=snm(4:6);cl(1:3)=cnm(7:9);sl(1:3)=snm(7:9)
      rr=rln(1);rlat=rln(2);rlon=rln(3)
      t=dsin(rlat*RAD);u=dcos(rlat*RAD)
      call RLAT_BLH(GRS,rln,BLH);call GNormalfd(BLH,NFD,GRS);gr=NFD(2)
      call BelPnmdt(pnm,dpt1,dpt2,3,t)
         tn=0.d0
         do m=0,2
            j=m+1
            cosml=dcos(rlon*RAD*m);sinml=dsin(rlon*RAD*m)
	      tn(1)=tn(1)+(cnm(j)*cosml+snm(j)*sinml)*pnm(j+3)!高程异常
	      tn(2)=tn(2)-1.5d0*(cnm(j)*cosml+snm(j)*sinml)*pnm(j+3)!地面重力
	      tn(2)=tn(2)+(ch(j)*cosml+sh(j)*sinml)*pnm(j+3)
	      tn(3)=tn(3)+(cnm(j)*cosml+snm(j)*sinml)*pnm(j+3)!扰动重力
	      tn(4)=tn(4)+(cnm(j)*cosml+snm(j)*sinml)*dpt1(j+3)!地倾斜
	      tn(4)=tn(4)-(ch(j)*cosml+sh(j)*sinml)*dpt1(j+3)
	      tn(5)=tn(5)+dble(m)*(cnm(j)*sinml-snm(j)*cosml)*pnm(j+3)
	      tn(5)=tn(5)-dble(m)*(ch(j)*sinml-sh(j)*cosml)*pnm(j+3)
	      tn(6)=tn(6)+(cnm(j)*cosml+snm(j)*sinml)*dpt1(j+3)!垂线偏差
	      tn(7)=tn(7)+dble(m)*(cnm(j)*sinml-snm(j)*cosml)*pnm(j+3)
	      tn(8)=tn(8)+dble(m)*(cl(j)*sinml-sl(j)*cosml)*pnm(j+3)!水平位移
	      tn(9)=tn(9)+(cl(j)*cosml+sl(j)*sinml)*dpt1(j+3)
	      tn(10)=tn(10)+(ch(j)*cosml+sh(j)*sinml)*pnm(j+3)!径向位移
	      tn(11)=tn(11)+(cnm(j)*cosml+snm(j)*sinml)*dpt2(j+3)
	      tn(12)=tn(12)+dble(m**2)*(cnm(j)*cosml+snm(j)*sinml)*pnm(j+3)
         enddo
         tn=tn*(ae/rr)**2
 	   tdn(1)=tdn(1)+tn(1)
	   tdn(2)=tdn(2)+tn(2)*dble(3)
	   tdn(3)=tdn(3)+tn(3)*dble(3)
	   tdn(4)=tdn(4)+tn(4)
	   tdn(5)=tdn(5)+tn(5)
	   tdn(6)=tdn(6)+tn(6)
	   tdn(7)=tdn(7)+tn(7)
	   tdn(8)=tdn(8)+tn(8)
	   tdn(9)=tdn(9)+tn(9)
 	   tdn(10)=tdn(10)+tn(10)
 	   tdn(12)=tdn(12)+tn(1)*dble(12)
	   tdn(13)=tdn(13)+tn(11)
 	   tdn(14)=tdn(14)+tn(12)
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
 	return
	end
