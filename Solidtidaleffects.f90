!  Solidtidaleffects.f90 
!
!  FUNCTIONS:
!  Solidtidaleffects - Entry point of console application.
!

      program Solidtidaleffects
      implicit none
	character out1(800),flEtd(800)
	character*800::outfl,fini
	character*6 NAM(900)
	character*800::line,str,astr,stm
	integer fln(2),i,j,k,NUM,nn,IY,IM,ID,ih,imn,kln,kk,sn,kpar(10)
      real*8 TAI_UTC,td1,td2,TRStCRS(3,3),mjd,mjd01,mjd02,tm,sec,rec(800)
	real*8::GRS(6),BLH(3),tmp,td,VAL(900),SSS(3),dt(2,6),cnm(80),snm(80)
	real*8 eop(80000,6),eoput1(6),frd(200,8),tdn(24),ENR(3),mjd1,mjd2
	real*8 dcin(30),dout(30),tm1,tm2
	real*8 plon, plat, phgt, bgntm, endtm, tmdlt,ltm
	integer::status=0
	character*800::jeph405fl
      common /JEPH/ jeph405fl
!---------------------------------------------------------------------
      GRS(1)= 3.986004415d14; GRS(2)=6378137.d0; GRS(3)=1.0826359d-3
      GRS(4) = 7.292115d-5; GRS(5)=1.d0/298.25641153d0
      !Input longitude (degree decimal), latitude (degree decimal), ellipsoidal height (m) 
      !输入经纬度大地高
	plon=121.24d0; plat=29.4281; phgt=17.830
      !Input starting time (long integer time), ending time, time interval (minute) 输入起止时间与时间间隔
      bgntm=20160701;endtm=20160704; tmdlt=10.d0/1440.d0
 	!Read IERS EOP file IERSeopc04.dat
      !MJD(day);EOP(i,1:6)= x(") y(") UT1-UT(s) LOD(s) dX(") dY(")
      open(unit=8,file="IERSeopc04.dat",status="old",iostat=status)
      if(status/=0)then
          dout(12)=1.d0; goto 902 
      endif
      do i=1,14
         read(8,'(a)') line
      enddo
	i=1
      do while(.not.eof(8))  
         read(8,'(a)') line
         call PickRecord(line,kln,rec,sn)
         if(sn<8)goto 407
         tm=rec(4);eop(i,1:6)=rec(5:10)
	   if(i==1)tm1=tm
         i=i+1
407      continue
	enddo
903   close(8)
      tm2=tm
      !Read the Love number correction file frqadjlovekhl.txt for frequency dependence
      !读取勒夫数频率相关系数文件frqadjlovekhl.txt
      open(unit=8,file="frqadjlovekhl.txt",status="old",iostat=status)
      if(status/=0)then!
          dout(13)=1.d0; goto 902 
      endif
      j=0;read(8,'(a)') line
      do i=1,119
        read(8,'(a)') line
        call PickRecord(line,kln,rec,sn)
	  j=j+1; frd(j,1:8)=rec(3:10)
	enddo
904   close(8)
      !Read the JPL Planetary ephemeris DE405/D440 file header to obtain the astronomical constants val(NUM).
      !读DE405/D440星历头文件,获取天文常数val(NUM)
      write(jeph405fl,*)"JPLEPH.440"
      NUM=645!DE405-156,DE440-645
	call CONST(NAM,VAL,SSS,NUM) !from ReadDE405.f90
      !Transform the long integer time (date) agreed by ETideLoad to year, month, day, hour, minute and second
      !ETideLoad格式日期tm转年月日时分秒。
      call tmcnt(bgntm,IY,IM,ID,ih,imn,sec)
      !Gregorian Calendar to Julian Date.
      call CAL2JD (IY,IM,ID,mjd,j)
      mjd01=mjd+dble(ih)/24.d0+dble(imn)/1440.d0+dble(sec)/864.d2 !GPS_MJD
      call tmcnt(endtm,IY,IM,ID,ih,imn,sec)
      call CAL2JD (IY,IM,ID,mjd,j)
      mjd02=mjd+dble(ih)/24.d0+dble(imn)/1440.d0+dble(sec)/864.d2 !GPS_MJD
      if(mjd02<mjd01)then
          dout(14)=1.d0; goto 902
      endif
      open(unit=10,file='reslt.txt',status="replace") !Output file 输出文件
      BLH(1)=plat;BLH(2)=plon;BLH(3)=phgt
      write(10,'(a8,2F12.6,F10.3,F15.6)')'Forcast',plon,plat,phgt,mjd01
      mjd=mjd01;mjd02=mjd02+1.d-6
      do while(mjd<mjd02)
         call EOP_UT1_UTC(mjd,tm1,eop,eoput1)
         if(mjd>tm2.or.mjd<tm1)eoput1=0.d0  !EOP=0.d0
         !Calculate the conversion matrix from ITRS to GCRS at mjd epoch time.
         !计算mjd历元时刻地球参考框架到地心惯性框架的转换矩阵TRStCRS
	   call TRStoCRS(mjd,eoput1,TRStCRS) 
	   td=mjd+2400000.5d0+32.184d0/864.d2
         tdn(1:14)=0.d0
         !Calculate the solid tidal effects on all-element geodetic variations.
         !计算大地测量全要素(10种)固体潮效应
	   call GravFdSldTd(td,BLH,tdn,eoput1,GRS,val,NUM,frd,TRStCRS)
         call mjdtotm(mjd,ltm);call tmtostr(ltm,stm)
         write(10,'(a15,F12.6,20F12.4)')adjustl(stm),mjd-mjd01,(tdn(i),i=1,14)!tdh(14): solid tidal effects. 固体潮效应（14列）
         mjd=mjd+tmdlt
906      continue
	enddo
902	continue
      close(10)
      write (*,*)'    Complete the computation! The results are saved in the file reslt.txt.'
      pause
      end
!
!************************************************************************
!
      subroutine tmcnt(tm,iyr,imo,idy,ihr,imn,sec)
      !Transform the long integer time (date) agreed by ETideLoad to year, month, day, hour, minute and second
      !ETideLoad格式日期tm转年月日时分秒。
      implicit none
	integer::iyr,imo,idy,ihr,imn,kln
	real*8::tm,sec,tmp,dj0,fd
	character*40::tmstr,astr,aymd,ahms
!-------------------------------------------------------------------------
      write(astr,*)tm
      astr=trim(adjustl(astr));kln=len(astr)
      read(astr(1:4),*)iyr
      read(astr(5:6),*)imo
      read(astr(7:8),*)idy
      ihr=0;imn=0;sec=0.d0
      if(kln>9)read(astr(9:10),*)ihr
      if(kln>11)read(astr(11:12),*)imn
      if(kln>13)read(astr(13:kln),*)sec
      continue
      end
!
!************************************************************************
!
      subroutine mjdtotm(mjd0,ltm)
      implicit none
	integer::IY,IM,ID,ihr,imn,k
	real*8::ltm,mjd0,mjd,sec,tmp
      mjd=mjd0
      call JD2CAL(2400000.5d0,mjd,IY,IM,ID,tmp,k)
      mjd=dble(nint(mjd*1.d6))*1.d-6;tmp=tmp+1.d-8
      ihr=floor(tmp*24.d0);tmp=tmp*24.d0-dble(ihr)
      imn=floor(tmp*6.d1);tmp=tmp*6.d1-dble(imn)
      sec=tmp*6.d1;
      if(nint(sec)-59.5>0)then
        sec=0.d0;imn=imn+1
      endif
      if(imn-59.5>0)then
        imn=0;ihr=ihr+1
      endif
      if(sec>5.d-2)ltm=IY*1.d10+IM*1.d8+ID*1.d6+ihr*1.d4+imn*1.d2+sec
      if(sec<5.d-2.and.imn>0.1)ltm=IY*1.d8+IM*1.d6+ID*1.d4+ihr*1.d2+imn
      if(sec<5.d-2.and.imn<0.1.and.ihr>0.1)ltm=IY*1.d6+IM*1.d4+ID*1.d2+ihr
      if(sec<5.d-2.and.imn<0.1.and.ihr<0.1)ltm=IY*1.d4+IM*1.d2+ID
      if(ltm<1901010131.d0)ltm=ltm*1.d2
      end
!
!************************************************************************
!
      subroutine tmtostr(tm,tmstr)
      implicit none
	character*800::tmstr,astr
	real*8::tm
!-------------------------------------------------------------------------
      write(astr,'(F16.0)')tm
      astr=trim(adjustl(astr))
      write(tmstr,*)astr(1:len_trim(astr)-1)
      return
      end


