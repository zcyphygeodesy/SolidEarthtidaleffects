      subroutine GravFdSldTd(td,BLH,tdn,eop,GRS,val,nn,frd,TRStCRS)
         !Calculate the solid tidal effects on all-element geodetic variations.
         !�����ز���ȫҪ��(10��)���峱ЧӦ
!eop(6)-x",y",UT1_JD,LODs,dx",dy"
!val(nn)-DE405/DE440 the astronomical constants !���ĳ���
!frd(71*9)-The corrections for frequency dependence  !Ƶ�������
!IERS2010,1~48-Tab6.5a,49~69Tab6.5b,70~71-ab6.5c
!-------------------------------------------------------------
      implicit none
	integer::nn,i,n,m,j,maxn,kk
	real*8::td,rln(3),BLH(3),grav(9),eop(6),GRS(6),val(nn),frd(200,8)
	real*8::ae,TRStCRS(3,3),tdn(24),pp,ftd(24)
      real*8::k4(3),cnm(80),snm(80),cc(80),ss(80)
!---------------------------------------------------------------------
	ae=GRS(2); cnm=0.d0; snm=0.d0;tdn(1:14)=0.d0
      call BLH_RLAT(GRS,BLH,rln)
      !Calculate the geopotential coefficient variations of the Earth's tidal generating potential (TGP) from celestial bodies
	!����N��(������)���۹��峱λϵ��
      call NormTide(td,val,nn,cnm,snm,rln,GRS,TRStCRS)
      !Calculate the degree-4 geopotential coefficient variations
      !����4��λϵ���仯
      cc=0.d0;ss=0.d0;grav=0.d0
      k4(1)=-.89d-3; k4(2)=-.8d-3; k4(3)=-.57d-3
      do i=1,3
         cc(7+i)=cnm(i)*k4(i); ss(7+i)=snm(i)*k4(i)
      enddo
      call RentGField(4,rln,cc,ss,grav,GRS)
      ! solid tidal effects on all-element geodetic variations using the nominal Love numbers.
      !�ñ���շ�������ȫҪ�ع��峱�α�
      call solidtdnormal(rln,6,cnm,snm,tdn,GRS)
      tdn(1)=tdn(1)+grav(2);tdn(2:3)=tdn(2:3)+grav(4)
      tdn(4:5)=tdn(4:5)+grav(5:6);tdn(6:7)=tdn(6:7)+grav(5:6)
      tdn(12:14)=tdn(12:14)+grav(7:9)
      !Calculate the geopoential coefficient adjustments caused by Love number frequency dependence
      !����λ�շ���Ƶ����ص��µ�λϵ���仯��
	call Freqdploven(td,cnm,snm,frd,eop)
      !Calculate the solid tidal effect adjustments caused by Love number frequency dependence
      !ȫҪ�ع��峱ЧӦ�շ���Ƶ�������У����
	ftd=0.d0;call Tdlovefreqadj(rln,cnm,snm,ftd,GRS)
      tdn(1:14)=tdn(1:14)+ftd(1:14)
9002  return
      end
