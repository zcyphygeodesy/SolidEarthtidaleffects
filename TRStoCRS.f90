      SUBROUTINE TRStoCRS(tm,eoput1,TRStCRS)
!-------------------------------------------------------------
      implicit none
	real*8::tm,TRStCRS(3,3)
	real*8::tm0,dt(2,6),eoput1(6)
	real*8 pdt,td,td1,td2,XP,YP,SP,ERA2000
      real*8 RPOM(3,3),THETA,RBPN(3,3)
	integer*4 k
!-----------------------------------------------------------------------
	td=tm+19.d0/864.d2+eoput1(3)/864.d2 !UT1_MJD
	call XYS2000A(2400000.5d0,td,XP,YP,SP)
	call POM2000(XP,YP,SP,RPOM)
	THETA=ERA2000(td,2400000.5d0)
	call RBPN2000(XP,YP,SP,RBPN)
	call T2C2000(RPOM,THETA,RBPN,TRStCRS)
      end
!
!*************************************************
!
      SUBROUTINE EOP_UT1_UTC(tm,tm0,eop,eoput1)
!-------------------------------------------------------------
      implicit none
	real*8::tm
	real*8::eop(80000,6),tm0,dt(2,6),eoput1(6)
	real*8 pdt,td,td1,td2
	integer*4 k
!-----------------------------------------------------------------------
	td=tm+19.d0/864.d2
	td1=dble(nint(td-0.5d0));td2=dble(nint(td+0.5d0))
	pdt=td-td1;k=nint(td1)-tm0+1;dt(1:2,1:6)=eop(k:k+1,1:6)
	eoput1(1:6)=dt(1,1:6)+pdt*(dt(2,1:6)-dt(1,1:6))
      end
