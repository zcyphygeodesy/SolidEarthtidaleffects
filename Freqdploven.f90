      subroutine Freqdploven(tdt,cnm,snm,frd,eop)
!eop(6)-x",y",UT1_JD,LODs,dx",dy"
!fr(71*9)-频率相关项IERS2010,1~48-Tab6.5a,49~69Tab6.5b,70~71-ab6.5c
!2006年10月5日,章传银;2015年7月修改
      implicit none
	integer::n,m,kk,i,j,dod(6),bita(6),bs
	real*8::tdt,t,cnm(80),snm(80),frd(200,8),eop(6),th0
	real*8::pi,RAD,th,v0,sigma,f(5),a0,a1,a2,cc(3),ss(3)
	real*8::frk(71,8),frh(30,8),frl(18,8),GMST2000
!---------------------------------------------------------------------------
   	pi=datan(1.d0)*4.d0;RAD=pi/180.d0;cnm=0.d0;snm=0.d0
	t=(tdt-2451545.d0)/36525.d0 !J2000.0 世纪
      frk(1:71,1:8)=frd(1:71,1:8)
      frh(1:30,1:8)=frd(72:101,1:8)
      frl(1:18,1:8)=frd(102:119,1:8)
      !five-vector of fundamental arguments Fj(the Delaunay variables 
      !l, l', F, D, omga) of nutation theory (IERS2010, 5.43)
	f(1)=134.96340251d0+(1717915923.2178d0*t+31.8792d0*t**2
     >	     +0.051635d0*t**3-0.0002447d0*t**4)/36.d2
	f(2)=357.52910918d0+(129596581.0481d0*t-0.5532d0*t**2
     >	     +0.000136d0*t**3-0.00001149d0*t**4)/36.d2
	f(3)=93.27209062d0+(1739527262.8478d0*t-12.7512d0*t**2
     >	     +0.001037d0*t**3-0.00000417d0*t**4)/36.d2
	f(4)=297.85019547d0+(1602961601.209d0*t-6.3706d0*t**2
     >	     +0.006593d0*t**3-0.00003196d0*t**4)/36.d2
	f(5)=125.04455501d0+(-6962890.5431d0*t+7.4722d0*t**2
     >     +0.007702d0*t**3-0.00005939d0*t**4)/36.d2
	f=dmod(f*RAD,2.d0*pi)
	th0=GMST2000(tdt,eop(3),tdt,0.d0)
      th0=th0+pi;a0=4.4228d-8;a1=-3.1274d-8;a2=3.1274d-8
      !位勒夫数频率相关性校正k
      do i=1,2!m=2半日
 	  th=0.d0
	  do j=1,5
	    th=th+f(j)*frk(i,j)
	  enddo
	  th=2.d0*th0-th
	  cnm(3)=cnm(3)+a2*((frk(i,6)-144.d0)*dcos(th)-frk(i,7)*dsin(th))*frk(i,8)*1.d-10
	  snm(3)=snm(3)-a2*((frk(i,6)-144.d0)*dsin(th)+frk(i,7)*dcos(th))*frk(i,8)*1.d-10
      enddo
	do i=3,50!m=1周日48
 	  th=0.d0
	  do j=1,5
	    th=th+f(j)*frk(i,j)
	  enddo
	  th=th0-th
	  cnm(2)=cnm(2)+a1*((frk(i,6)-130.d0)*dsin(th)+frk(i,7)*dcos(th))*frk(i,8)*1.d-10
	  snm(2)=snm(2)+a1*((frk(i,6)-130.d0)*dcos(th)-frk(i,7)*dsin(th))*frk(i,8)*1.d-10
	enddo
	do i=51,71!m=0长周期21
	  th=0.d0
	  do j=1,5
	    th=th+f(j)*frk(i,j)
        enddo
        th=-th
	  cnm(1)=cnm(1)+a0*(frk(i,6)*dcos(th)-frk(i,7)*dsin(th))*frk(i,8)*1.d-10
	enddo
      cc=0.d0;ss=0.d0!径向位移球谐系数频率相关性校正h
      i=1!m=2半日
 	  th=0.d0
	  do j=1,5
	    th=th+f(j)*frh(i,j)
	  enddo
	  th=2.d0*th0-th
	  cc(3)=cc(3)+a2*(frh(i,6)*dcos(th)-frh(i,7)*dsin(th))*frh(i,8)*1.d-9
	  ss(3)=ss(3)-a2*(frh(i,6)*dsin(th)+frh(i,7)*dcos(th))*frh(i,8)*1.d-9
	do i=2,25!m=1周日24
 	  th=0.d0
	  do j=1,5
	    th=th+f(j)*frh(i,j)
	  enddo
	  th=th0-th
	  cc(2)=cc(2)+a1*(frh(i,6)*dsin(th)+frh(i,7)*dcos(th))*frh(i,8)*1.d-9
	  ss(2)=ss(2)+a1*(frh(i,6)*dcos(th)-frh(i,7)*dsin(th))*frh(i,8)*1.d-9
	enddo
	do i=26,30!m=0长周期5
	  th=0.d0
	  do j=1,5
	    th=th+f(j)*frh(i,j)
        enddo
        th=-th
	  cc(1)=cc(1)+a0*(frh(i,6)*dcos(th)-frh(i,7)*dsin(th))*frh(i,8)*1.d-9
      enddo
      cnm(4:6)=cc(1:3);snm(4:6)=ss(1:3)
      cc=0.d0;ss=0.d0!水平位移球谐系数频率相关性校正l
      i=1!m=2半日
 	  th=0.d0
	  do j=1,5
	    th=th+f(j)*frl(i,j)
	  enddo
	  th=2.d0*th0-th
	  cc(3)=cc(3)+a2*(frl(i,6)*dcos(th)-frl(i,7)*dsin(th))*frl(i,8)*1.d-9
	  ss(3)=ss(3)-a2*(frl(i,6)*dsin(th)+frl(i,7)*dcos(th))*frl(i,8)*1.d-9
	do i=2,13!m=1周日13
 	  th=0.d0
	  do j=1,5
	    th=th+f(j)*frl(i,j)
	  enddo
	  th=th0-th
	  cc(2)=cc(2)+a1*(frl(i,6)*dsin(th)+frl(i,7)*dcos(th))*frl(i,8)*1.d-9
	  ss(2)=ss(2)+a1*(frl(i,6)*dcos(th)-frl(i,7)*dsin(th))*frl(i,8)*1.d-9
	enddo
	do i=14,18!m=0长周期5
	  th=0.d0
	  do j=1,5
	    th=th+f(j)*frl(i,j)
        enddo
        th=-th
	  cc(1)=cc(1)+a0*(frl(i,6)*dcos(th)-frl(i,7)*dsin(th))*frl(i,8)*1.d-9
      enddo
      cnm(7:9)=cc(1:3);snm(7:9)=ss(1:3)
	return
	end
