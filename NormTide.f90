      subroutine NormTide(tm,val,nn,cnm,snm,rln,GRS,TRStCRS)
!������Ԫtd(UT1_JD)��N��(������)���峱λϵ��-ֱ��Ӱ��
!val(nn)��DE405���ĳ���
!2015��7��21��,�´���
      implicit none
	integer::nn,n,m,kk,i,j,sp
	real*8::tm,val(nn),rln(3),XYZ(3),RLATLON(3)
	real*8::cnm(80),snm(80),cosml,sinml,tmp
	real*8::pnm(80),pnm1(80),din(4),GRS(6),RRD(6),gm(11),grav(6)
	real*8::rr,rlat,rlon,ae,pi,RAD,auday,TRStCRS(3,3)
!---------------------------------------------------------------------------
	!auday=(val(7)*1.d3)**3/(8.64d4*8.64d4);!DE405
	auday=(val(10)*1.d3)**3/(8.64d4*8.64d4);!DE440-val(10),DE405-val(7)
	!gm(1:9)=val(9:17); gm(11)=val(18);gm=gm*auday!DE405
	gm(1:9)=val(12:20); gm(11)=val(21);gm=gm*auday!DE440
	gm(3)=GRS(1);gm(10)=gm(3)/val(11) !DE400-val(11),DE405-val(8)
   	ae=GRS(2);pi=datan(1.d0)*4.d0;RAD=pi/180.d0
      cnm=0.d0; snm=0.d0
	do i=1,11
	  if(abs(i-3)>0)then!��������i=10,̫��11,����2,����5,����4,ˮ��1,����6,������7,������8,ڤ����9
     	    call PLEPH(tm,i,3,RRD)
	    !XYZ=RRD(1:3)*val(7)*1.d3!DE405
	    XYZ=RRD(1:3)*val(10)*1.d3!DE240
          do j=1,3
            tmp=RRD(j)
            tmp=XYZ(j)
          enddo
          XYZ=matmul(transpose(TRStCRS),XYZ)
          call XYZ_RLAT(GRS,XYZ,RLATLON)
	    rr=RLATLON(1);rlat=RLATLON(2);rlon=RLATLON(3)
!	    ����ȫ��Legendre����pnm,����pnm���±��cnm��snm��3
          pnm=0.d0;pnm1=0.d0
          call PlmBar_d(pnm, pnm1, 6, rlat)
          do n=2,6
             do m=0,n
                j=n*(n+1)/2-2+m
                cosml=dcos(dble(m)*rlon*RAD);sinml=dsin(dble(m)*rlon*RAD)
 	          cnm(j)=cnm(j)+pnm(j+3)*cosml*gm(i)/gm(3) 
     >	              *dexp(dble(n+1)*dlog(ae/rr))/dble(2*n+1)
	          snm(j)=snm(j)+pnm(j+3)*sinml*gm(i)/gm(3)
     >	              *dexp(dble(n+1)*dlog(ae/rr))/dble(2*n+1)
            enddo
          enddo
	  endif
	enddo
	end
