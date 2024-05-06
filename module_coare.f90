!========================================================================
! The module to calculate heat fluxes via COARE (2.6 or 4.0)
! 
! components:
!   ---COARE26Z (UWIND,VWIND,ZU,TA,ZT,QV,ZQ,PA,TS,DLW,DSW,HSB,HLB,NET,USTRESS,VSTRESS,LAT)
!   ---COARE40VN(UWIND,VWIND,ZU,TA,ZT,QV,ZQ,PA,TS,DLW,DSW,HSB,HLB,NET,USTRESS,VSTRESS,LAT)
!   ---COARE26Z_2d (UWIND,VWIND,ZU,TA,ZT,QV,ZQ,PA,TS,DLW,DSW,HSB,HLB,NET,USTRESS,VSTRESS,LAT)
!   ---COARE40VN_2d(UWIND,VWIND,ZU,TA,ZT,QV,ZQ,PA,TS,DLW,DSW,HSB,HLB,NET,USTRESS,VSTRESS,LAT)

!
!   Variable  |    Standard Name         | Units
!   ----INPUT----
!    UWIND    | wind on x-direction      | m/s
!    VWIND    | wind on y-direction      | m/s
!    ZU       | height of WIND           | m
!    TA       | temperature of air       | K
!    ZT       | height of TA             | m
!    QV       | water vapor mixing ratio | kg kg-1
!    PA       | surface pressure         | Pa 
!    TS       | temperature of sea       | K
!    DLW*     | longwave radiation       | W m-2 (positive: ocean gains energy)
!    DSW      | shortwave radiation      | W m-2 (positive: ocean gains energy)
!    LAT      | latitude                 | degree (positive: N) 
!   ----OUTPUT----
!    DLW*     | longwave radiation       | W m-2 (positive: ocean gains energy)
!    HSB      | sensible heat flux       | W m-2 (positive: ocean gains energy)
!    HLB      | latent heat flux         | W m-2 (positive: ocean gains energy)
!    NET      | net heat flux            | W m-2 (positive: ocean gains energy)
!    USTRESS  | wind stress on x         | Pa
!    VSTRESS  | wind stress on y         | Pa
!========================================================================

MODULE module_coare


CONTAINS


SUBROUTINE COARE26Z_2d (UWIND,VWIND,ZU,TA,ZT,QV,ZQ,PA,TS,DLW,DSW,HSB,HLB,NET,USTRESS,VSTRESS,LAT)
  !
  implicit none
  !
  real, dimension(:,:), intent(in   ) :: UWIND, VWIND, TA, QV, PA, TS, DSW, LAT
  real, intent(in   )                 :: ZU, ZT, ZQ
  real, dimension(:,:), intent(inout) :: DLW
  real, dimension(:,:), intent(  out) :: HSB, HLB, NET, USTRESS, VSTRESS

  integer :: i, j, nx, ny

  nx = size(TA, 1)
  ny = size(TA, 2)

  do i = 1, nx
    do j = 1, ny
      call COARE26Z (UWIND(i,j),VWIND(i,j),ZU,TA(i,j),ZT,QV(i,j),ZQ,PA(i,j),TS(i,j), &
                       DLW(i,j),DSW(i,j),HSB(i,j),HLB(i,j),NET(i,j),                   &
                       USTRESS(i,j),VSTRESS(i,j),LAT(i,j))
    end do
  end do
  !
END SUBROUTINE 

SUBROUTINE COARE40VN_2d(UWIND,VWIND,ZU,TA,ZT,QV,ZQ,PA,TS,DLW,DSW,HSB,HLB,NET,USTRESS,VSTRESS,LAT)
  !
  implicit none
  !
  real, dimension(:,:), intent(in   ) :: UWIND, VWIND, TA, QV, PA, TS, DSW, LAT
  real, intent(in   )                 :: ZU, ZT, ZQ
  real, dimension(:,:), intent(inout) :: DLW
  real, dimension(:,:), intent(  out) :: HSB, HLB, NET, USTRESS, VSTRESS

  integer :: i, j, nx, ny

  nx = size(TA, 1)
  ny = size(TA, 2)

  do i = 1, nx
    do j = 1, ny
      call COARE40VN(UWIND(i,j),VWIND(i,j),ZU,TA(i,j),ZT,QV(i,j),ZQ,PA(i,j),TS(i,j), &
                       DLW(i,j),DSW(i,j),HSB(i,j),HLB(i,j),NET(i,j),                   &
                       USTRESS(i,j),VSTRESS(i,j),LAT(i,j))
    end do
  end do
  !
END SUBROUTINE 

!   ---COARE40VN(UWIND,VWIND,ZU,TA,ZT,QV,ZQ,PA,TS,DLW,DSW,HSB,HLB,NET,USTRESS,VSTRESS,LAT)
!---------------------------------------------------------------------------
! COARE26Z
! The original version is written by Dr. Song Hu. 
! I made FOUR changes:
! 1) Move the latitude into the arguments.
! 2) Temperatures are all in unit of Kevin.
! 3) Use u10, v10, q2 as input, rather than wspd10 and RH(relative humidity);
! 4) Recalculate the Longwave, and calculate the net heat flux
! 
!   Variable  |    Standard Name         | Units
!    UWIND    | wind on x-direction      | m/s
!    VWIND    | wind on y-direction      | m/s
!    ZU       | height of WIND           | m
!    TA       | temperature of air       | K
!    ZT       | height of TA             | m
!    QV       | water vapor mixing ratio | kg kg-1
!    PA       | surface pressure         | Pa 
!    TS       | temperature of sea       | K
!    DLW      | longwave radiation       | W m-2 (positive: downward --->
!                                           ocean gain energy, sign changed)
!    DSW      | shortwave radiation      | W m-2 (positive: downward --->
!                                           ocean gain energy, same sign)
!    HSB      | sensible heat flux       | W m-2 (positive: ocean gain energy)
!    HLB      | latent heat flux         | W m-2 (positive: ocean gain energy)
!    NET      | net heat flux            | W m-2 (positive: ocean gain energy)
!    USTRESS  | wind stress on x         | Pa
!    VSTRESS  | wind stress on y         | Pa
!    LAT      | latitude                 | degree (positive: N) 
!
! Siqi Li
!---------------------------------------------------------------------------
SUBROUTINE COARE26Z(UWIND,VWIND,ZU,TA,ZT,QV,ZQ,PA,TS,DLW,DSW,HSB,HLB,NET,USTRESS,VSTRESS,LAT)
!SUBROUTINE COARE26Z(UR,ZU,TA,ZT,RH,ZQ,PA,TS,DLW,DSW,tau,hsb,hlb)
!  A=coare26(u,zu,Ta,zt,rh,zq,Pa,Ts,dlw,dsw)
! Simplified non-vectorized version of coare2.6 code
! with cool skin option retained but warm layer and 
! surface wave options removed, and rain set to zero. 
! Assumes input are single scalars and that u is the
! magitude of the difference between wind and surface 
! current vectors, if latter available. Output:
! A = [tau qsen qlat Cd Ch Ce Cdn_10 Chn_10 Cen_10].
        implicit none

!---> Siqi Li
! I rewrite this part, just to make sure the input variables do not change
! during the calculation. Values of TA and TS do not change, though their 
! intent is 'inout'.
  real, intent(in)    :: UWIND, VWIND, ZU, TA, ZT, QV, ZQ, PA, TS, DSW, LAT
  real, intent(inout) :: DLW
  real, intent(out)   :: HSB, HLB, NET, USTRESS, VSTRESS
!<--- Siqi Li
        real  U
        real  UR
        real  US
!        real  TS
        real  TSea  ! Siqi
        real  T
!        real  TA
        real  RH
        real  P
!        real  PA
        real  RL
!        real  DLW
        real  RS
!        real  DSW
        real  QSAT26S
        real  QS
        real  Q
        real  ZI
        real  BETA
        real  VON
        real  FDG
        real  TDK
        real  GRVS
        real  GRAV
        real  RGAS
        real  CPA
        real  CPV
        real  RHOA
        real  VISA
        real  AL
        real  BE
        real  CPW
        real  RHOW
        real  VISW
        real  TCW
        real  BIGC
        real  WETC
        real  RNS
        real  RNL
        real  DU
        real  DT
!        real  ZT
        real  DQ
        real  UG
        real  DTER
        real  DQER
        real  UT
        real  SQRT
        real  U10
!        real  ZU
        real  USR
        real  ZO10
        real  CD10
        real  CH10
        real  CT10
        real  ZOT10
        real  CD
        real  CT
        real  CC
        real  RIBCU
        real  RIBU
        real  NITS
        real  ZETU
        real  PSIU_26S
        real  PSIT_26S
        real  TSR
!        real  ZQ
        real  QSR
        real  TKT
        real  CHARN
        real  ZET
        real  ZO
        real  ZOQ
        real  ZOT
        real  LOG
        real  BF
!        real  HSB
!        real  HLB
        real  QOUT
        real  DELS
        real  QCOL
        real  ALQ
        real  XLAMX
        real  TAU
        real  CH
        real  CE
        integer I
        real CDN_10,CHN_10,CEN_10
           REAL LE,L,Jcool,L10!,LAT
           real rr
        real AG, theta

! Set jcool=0 if Ts is surface, =1 if Ts is bulk.
! rcb checked 6/9/04
! set jcool=1 if Ts is bulk, 0 if Ts is true skin jcool=1;
           JCOOL=1.
! rename variables from fairall et al coare3 code
 
!---> Siqi Li
! There are three things done before the COARE26Z in wrf_to_fvcom.F90
! I moved them into this subroutine.
! 1. Calculate the wind speed
! 2. Calculate the relative humidity
! 3. Convert the unit of temperature (Sea and Air) from C to K
   UR=sqrt(UWIND**2+VWIND**2)  ! Add by Siqi Li
   
   RH=TA-273.16
   RH=6.112*exp(17.67*RH/(RH+243.5))
   RH=0.622*RH/(PA-RH)
   RH=QV*100/RH
   if (RH > 100.) RH=100.

!   TA=TA-273.16
!   TS=TS-273.16
!<--- Siqi Li

! wind speed (m/s) at height zr (m)
           U=UR
!  surface current speed in the wind direction(m/s)
           US=0*UR
!  water temperature (deg C)
           TSea=TS - 273.16
!  BULK AIR TEMPERATURE (C) AT HEIGHT ZT(m)
            T=TA - 273.16
!  RELATIVE HUMIDITY (%) AT HEIGHT zq(M)
           RH=RH
!  SURFACE PRESSURE (mb)
!           P=PA
           P=PA/100. 
! Siqi Li, do not use PA(Pa), or the latent will have large error
!  DOWNWARD LONGWAVE RADIATION (W/m2)
           RL=DLW
!  DOWNWARD SHORTWAVE RADIATION (W/m2)
          RS=DSW
!  CONVERT RH TO SPECIFIC HUMIDITY (G/KG)
            CALL QSAT26(TSea,P,QSAT26S)
            QS=0.98*QSAT26S/1000.
!  SPECIFIC HUMIDITY OF AIR (G/KG)  
            CALL QSAT26(T,P,QSAT26S)
            Q=(0.01*RH)*QSAT26S/1000.
!   SET RAIN TO ZERO RAIN=0*U
!   SET RAIN RATE (MM/HR) - KEEP AS OPTION

!  ***********SET LOCAL CONSTANTS *********
!    PBL HEIGHT (M)
        ZI=600.
!     LATITUDE (DEG,N=+)- GEORGES BANK
!---> Siqi Li, Set this as an input
!        LAT=42.
!<--- Siqi Li
! ************SET CONSTANTS **************
        BETA=1.2
        VON=0.4
        FDG=1.00
        TDK=273.16
        CALL GRV(LAT,GRVS)
        GRAV=GRVS
! ************AIR CONSTANTS **************
        RGAS=287.1
        LE=(2.501-0.00237*TSea)*1000000.
        CPA=1004.67
        CPV=CPA*(1.+0.84*Q)
        rhoa=P*100./(Rgas*(t+tdk)*(1+0.61*q));
        VISA=1.326*0.00001*(1+6.542*0.001*T+8.301*0.000001*T*T-4.84*0.000000001*T*T*T)
! ***********COOL SKIN CONSTANTS******************
        AL=2.1*0.00001*((TSea+3.2)**0.79)
        BE=0.026
        CPW=4000.
        RHOW=1022.
        VISW=0.000001
        TCW=0.6
        BIGC=16.*GRAV*CPW*((RHOW*VISW)**3)/(TCW*TCW*RHOA*RHOA)
        WETC=0.622*LE*QS/(RGAS*((TSea+TDK)**2))
! ***************COMPUTE AUX STUFF *********
         RNS=RS*0.945
         RNL=0.97*(5.67*0.00000001*((TSea-0.3*JCOOL+TDK)**4)-RL)


! **************BEGIN BULK LOOP ***********

! **************FIRST GUESS ***************
         DU=U-US
         DT=TSea-T-0.0098*ZT
         DQ=QS-Q
!         TA=T+TDK
         UG=0.5
         DTER=0.3
         DQER=WETC*DTER
         UT=SQRT(DU*DU+UG*UG) 
         U10=UT*LOG(10./1e-4)/LOG(ZU/1e-4)
         USR=0.035*U10

         zo10=0.011*usr*usr/grav+0.11*visa/usr
         Cd10=(von/log(10./zo10))**2
         Ch10=0.00115
         Ct10=Ch10/sqrt(Cd10)
         zot10=10./exp(von/Ct10)
         Cd=(von/log(zu/zo10))**2
         Ct=von/log(zt/zot10)
         CC=von*Ct/Cd
         Ribcu=-zu/(zi*0.004*(Beta**3))
         Ribu=-grav*zu*((dt-dter*jcool)+.61*ta*dq)/(ta*(ut**2))
!        same as edson
         nits=6.   
         if(Ribu.lt.0)then
          zetu=CC*Ribu/(1.+Ribu/Ribcu)
         else
          zetu=CC*Ribu/(1.+27./(9*Ribu*cc))
          endif
         L10=zu/zetu
         if(zetu.gt.50)then
          nits=1
          endif
         CALL psiu_26(zu/L10,psiu_26s)
          usr=ut*von/(log(zu/zo10)-psiu_26s)
         CALL psit_26(zt/L10,psit_26s)
          tsr=-(dt-dter*jcool)*von*fdg/(log(zt/zot10)-psit_26s)
         CALL psit_26(zq/L10,psit_26s)
          qsr=-(dq-wetc*dter*jcool)*von*fdg/(log(zq/zot10)-psit_26s)
          tkt=.001
          charn=0.011
          if(ut.gt.10)then 
          charn=0.011+(ut-10)/(18-10)*(0.018-0.011)
          endif
          if(ut.gt.18)then
          charn=0.018
          endif

!*************** bulk loop ************
!            do I=1,nits
             do I=1,1
           zet=von*grav*zu/ta*(tsr*(1.+0.61*Q)+0.61*ta*qsr)/(usr*usr)/(1.+0.61*Q)
           zo=charn*usr*usr/grav+0.11*visa/usr
           rr=zo*usr/visa
           L=zu/zet
           zoq=amin1(1.15e-4,5.5e-5/(rr**0.6)) 
           zot=zoq
           CALL psiu_26(zu/L,psiu_26s)
           usr=ut*von/(log(zu/zo)-psiu_26s)
           CALL psit_26(zt/L,psit_26s)
           tsr=-(dt-dter*jcool)*von*fdg/(log(zt/zot)-psit_26s)
           CALL psit_26(zq/L,psit_26s)
           qsr=-(dq-wetc*dter*jcool)*von*fdg/(log(zq/zoq)-psit_26s) 
           Bf=-grav/ta*usr*(tsr+.61*ta*qsr)               
           if(Bf.gt.0)then
            ug=Beta*((Bf*zi)**0.333)
           else
            ug=0.2
           endif
           ut=sqrt(du*du+ug*ug) 
            Rnl=0.97*(5.67*0.00000001*((ts-dter*jcool+tdk)**4)-Rl)
            hsb=-rhoa*cpa*usr*tsr
            hlb=-rhoa*Le*usr*qsr
            qout=Rnl+hsb+hlb
 

            dels=Rns*(.065+11*tkt-  6.6*0.00001/(tkt*(1-exp(-tkt/8.0*0.0001))))
            qcol=qout-dels
            alq=Al*qcol+be*hlb*cpw/Le
            if (alq.gt.0)then
            xlamx=6./(1.+(bigc*alq/usr**4)**0.75)**0.333
            tkt=xlamx*visw/(sqrt(rhoa/rhow)*usr)
            else
            xlamx=6.0
           tkt=amin1(.01,xlamx*visw/(sqrt(rhoa/rhow)*usr))
            endif
!           cool skin
            dter=qcol*tkt/tcw
            dqer=wetc*dter
            enddo
!  of  end bulk iter loop

!****** compute fluxes ******************************************
!             wind stress
            tau=rhoa*usr*usr*du/ut
!              sensible heat flux
            hsb=rhoa*cpa*usr*tsr
!              latent heat flux
            hlb=rhoa*Le*usr*qsr
!****** compute transfer coeffs relative to ut @ meas. ht ********
            Cd=tau/rhoa/ut/amax1(.1,du)
            Ch=-usr*tsr/ut/(dt-dter*jcool)
            Ce=-usr*qsr/(dq-dqer*jcool)/ut

!****** compute 10-m neutral coeff relative to ut ****************
            Cdn_10=von*von/log(10./zo)/log(10./zo)
            Chn_10=von*von*fdg/log(10./zo)/log(10./zot)
            Cen_10=von*von*fdg/log(10./zo)/log(10./zoq)

!******** rain heat flux (save to use if desired) *************
! dwat=2.11e-5*((t+tdk)/tdk)^1.94; %! water vapour diffusivity
! dtmp=(1.+3.309e-3*t-1.44e-6*t*t)*0.02411/(rhoa*cpa); %!heat diffusivity
! alfac= 1/(1+(wetc*Le*dwat)/(cpa*dtmp)); %! wet bulb factor
! RF= rain*alfac*cpw*((ts-t-dter*jcool)+(Qs-Q-dqer*jcool)*Le/cpa)/3600;
!**************************************************************
!----------------------------------------------------------

!---> Siqi Li

!! The TS was converted back again to K. Change it into C
!! TS = TS + 273.15

! Below are added according to wrf_to_fvcom_26z.F90.
  ! Stefan-Boltzmann Law
  DLW = DLW- 0.98*5.6697*((TS*0.01)**4)
  NET = DSW + DLW + HSB + HLB

! Calculate the wind stress
  USTRESS = UWIND/UR * tau
  VSTRESS = VWIND/UR * tau

!<--- Siqi Li

            end subroutine COARE26Z
            subroutine psit_26(zet,psi)
! computes temperature structure function
        implicit none
        real  ZET
        real  X
        real  PSIK
        real  PSIC
        real  F
        real  PSI
        real  C

            if (zet.lt.0)then
            x=(1.-15.*zet)**.5
            psik=2.*log((1.+x)/2.)
            x=(1.-34.15*zet)**.3333
            psic=1.5*log((1.+x+x*x)/3.)-sqrt(3.)*atan((1.+2.*x)/sqrt(3.))+4.*atan(1.)/sqrt(3.)
            f=zet*zet/(1.+zet*zet)
            psi=(1.-f)*psik+f*psic
            else
            c=amin1(50.,.35*zet)
            psi=-((1.+2./3.*zet)**1.5+.6667*(zet-14.28)/exp(c)+8.525)
            endif
            end subroutine psit_26
  
!----------------------------------------------------------
            subroutine psiu_26(zet,psi)
! computes velocity structure function
        implicit none
        real  ZET
        real  X
        real  PSIK
        real  PSIC
        real  F
        real  PSI
        real  C

            if(zet.lt.0)then
             x=(1.-15.*zet)**.25
             psik=2.*log((1.+x)/2.)+ log((1.+x*x)/2.)-2.*atan(x)+2.*atan(1.)
             x=(1.-10.15*zet)**.3333
             psic=1.5*log((1.+x+x*x)/3.)-sqrt(3.)*atan((1.+2.*x)/sqrt(3.))+4.*atan(1.)/sqrt(3.)
             f=zet*zet/(1.+zet*zet)
             psi=(1.-f)*psik+f*psic
            else
            c=amin1(50.,.35*zet)
            psi=-((1.+1.0*zet)**1.0+.667*(zet-14.28)/exp(c)+8.525)
            endif
            end subroutine psiu_26
!-----------------------------------------------------------
            subroutine qsat26(T,P,qs)
! computes saturation specific humidity
        implicit none
        real es,T,P,qs

         es=6.112*exp(17.502*T/(T+241.0))*(1.0007+3.46*0.000001*P)
            qs=es*622/(P-.378*es)
            end subroutine qsat26

!----------------------------------------- function g=grv(lat)
           subroutine grv(lat,g)
        implicit none
        real gamma,c1,c2,c3,c4,phi,lat,pi,x,g
           pi=3.1415926
! computes g given lat in deg
           gamma=9.7803267715
           c1=0.0052790414
           c2=0.0000232718
           c3=0.0000001262
           c4=0.0000000007
           phi=lat*pi/180.
           x=sin(phi)
           g=gamma*(1.+c1*x**2+c2*x**4+c3*x**6+c4*x**8)
           end subroutine grv


!---------------------------------------------------------------------------
! COARE40VN
! The original version is coare40vn.m:
! Vectorized version of COARE 3 code (Fairall et al, 2003) with 
! modification based on the CLIMODE, MBL and CBLAST experiments 
! (Edson et al., JPO, 43, 2013).
!
! Dr. Zhongxiang Wu changed it into Fortran.
! Siqi Li made some changes and tested it.
!
! I made FOUR changes:
! 1) Delete some unused variables 'zi,rain,cp,sigH'
! 2) Temperatures are all in unit of Kevin.
! 3) Use u10, v10, q2 as input, rather than wspd10 and RH(relative humidity);
! 4) Recalculate the Longwave, and calculate the net heat flux
! 
!   Variable  |    Standard Name         | Units
!    UWIND    | wind on x-direction      | m/s
!    VWIND    | wind on y-direction      | m/s
!    ZU       | height of WIND           | m
!    TA       | temperature of air       | K
!    ZT       | height of TA             | m
!    QV       | water vapor mixing ratio | kg kg-1
!    PA       | surface pressure         | hPa (mb)
!    TS       | temperature of sea       | K
!    DLW      | longwave radiation       | W m-2 (positive: downward --->
!                                           ocean gain energy, sign changed)
!    DSW      | shortwave radiation      | W m-2 (positive: downward --->
!                                           ocean gain energy, same sign)
!    HSB      | sensible heat flux       | W m-2 (positive: ocean gain energy)
!    HLB      | latent heat flux         | W m-2 (positive: ocean gain energy)
!    NET      | net heat flux            | W m-2 (positive: ocean gain energy)
!    LAT      | latitude                 | degree (positive: N) 
!
! Siqi Li
!---------------------------------------------------------------------------

subroutine COARE40VN(UWIND,VWIND,ZU,TA,ZT,QV,ZQ,PA,TS,DLW,DSW,HSB,HLB,NET,USTRESS,VSTRESS,LAT)

!      subroutine
!      coare40vn(u,zu,t,zt,rh,zq,P,ts,Rl,Rs,tau,hsb,hlb,lat,zi,rain,cp,sigH,fmiss)
! No-vectorized verion - Revised from the vectorized of coare40vn.m
! Zhongxiang Wu
! 4/8/2016
        
! Vectorized version of COARE 3 code (Fairall et al, 2003) with 
! modification based on the CLIMODE, MBL and CBLAST experiments 
! (Edson et al., JPO, 43, 2013). 
!
! The current version of the code includes the wind-speed, wave-age and
! sea-state dependent parameterizations of the Charnock variabile as
! described in Edson et al. (2013).  The parameterization is chosen by the 
! inputed values of
!    cp = phase speed of dominant waves (m/s)  
!  sigH = significant wave height (m)
!
! An important component of this code is whether the inputed ts 
! represents the skin temperature of a near surface temperature.  
! How this variable is treated is determined by the jcool parameter:
! set jcool=1 if Ts is bulk ocean temperature (default),
!     jcool=0 if Ts is true ocean skin temperature. 
!********************************************************************
!
! The code assumes u,t,rh,ts are vectors; 
! sensor heights zu,zt,zl, latitude lat, and PBL height zi are constants;
! air pressure P and radiation Rs,Rl may be vectors or constants. 
! Default values are assigned for P,Rs,Rl,lat,and zi if these data are not 
! available.  Input NaNs to indicate no data. Defaults should be set to 
! representative regional values if possible.
!
! Inputs:  
!
!     u = relative wind speed (m/s) at heigth zu 
!    zu = height of wind speed measurement (m)
!     t = bulk air temperature (degC) at height zt
!    zt = height of temperature measurement (m)
!    rh = relative humidity (%) at height zq
!    zq = height of relative humidity measurement (m)
!     P = surface air pressure (mb) (default = 1015)
!    ts = water temperature (degC) see jcool below
!    Rs = downward shortwave radiation (W/m^2) (default = 150) 
!    Rl = downward longwave radiation (W/m^2) (default = 370)
!   lat = latitude (default = +45 N)
!    zi = PBL height (m) (default = 600m)
!  rain = rain rate (mm/hr) - not required for turbulent flux estimates,
!         set to NaN if unavailable.
!    cp = phase speed of dominant waves (m/s)  
!  sigH = significant wave height (m)
!

! The user controls the output.  This is currently set as:
! 
! A=[usr tau hsb hlb hbb hsbb wbar  tsr qsr zot zoq Cd Ch Ce  L 
!      1   2   3   4   5    6    7    8   9  10  11 12 13 14 15
!    zet dter dqer tkt Urf Trf Qrf RHrf UrfN Rnl Le rhoa UN U10 U10N
!     16   17   18  19  20  21  22   23   24  25 26   27 28  29   30
! Cdn_10 Chn_10 Cen_10 RF Qs Evap T10 Q10 RH10];
!     31     32     33 34 35   36  37  38   39        
!  where
!
!   usr = friction velocity that includes gustiness (m/s)
!   tau = wind stress (N/m^2)
!   hsb = sensible heat flux into ocean (W/m^2)
!   hlb = latent heat flux into ocean (W/m^2)
!   hbb = buoyany flux into ocean (W/m^2)
!   hsbb = "sonic" buoyancy flux measured directly by sonic anemometer 
!   wbar = the vertical velocity required to "Webb Correct" direct
!          covariance fluxes
!   tsr = temperature scaling parameter (K)
!   qsr = specific humidity scaling parameter (g/Kg)
!   zot = thermal roughness length (m)
!   zoq = moisture roughness length (m)
!   Cd = wind stress transfer (drag) coefficient at height zu   
!   Ch = sensible heat transfer coefficient (Stanton number) at height zu   
!   Ce = latent heat transfer coefficient (Dalton number) at height zu
!    L = Obukhov length scale (m) 
!  zet = Monin-Obukhov stability parameter zu/L 
! dter = cool-skin temperature depression (degC)
! dqer = cool-skin humidity depression (degC)
!  tkt = cool-skin thickness (m)
!  Urf = wind speed at reference height (user can select height below)
!  Trf = temperature at reference height
!  Qrf = specific humidity at reference height
! RHrf = relative humidity at reference height
! UrfN = neutral value of wind speed at reference height
!  Rnl = Upwelling IR radiation computed by COARE
!   Le = latent heat of vaporization
! rhoa = density of air
!   UN = neutral value of wind speed at zu
!  U10 = wind speed adjusted to 10 m
! UN10 = neutral value of wind speed at 10m
!Cdn_10 = neutral value of drag coefficient at 10m    
!Chn_10 = neutral value of Stanton number at 10m    
!Cen_10 = neutral value of Dalton number at 10m
!    RF = sensible heat flux due to rain (W/m^2)
!    Qs = surface value of specific humidity (g/kg) without cool skin
!         correction
!  Evap = evaporation rate from surface in mm/hr
!   T10 = temperarure at 10 m 
!   Q10 = specific humidity at 10 m
!  RH10 = relative humidity at 10 m
!

! Notes: 1) u is the relative wind speed, i.e., the magnitude of the
!           difference between the wind (at zu) and ocean surface current 
!           vectors.
!        2) Set jcool=0 in code if ts is true surface skin temperature,
!           otherwise ts is assumed the bulk temperature and jcool=1.
!        3) Set P=NaN to assign default value if no air pressure data 
!           available. 
!        4) Set Rs=NaN, Rl=NaN if no radiation data available.  This assigns 
!           default values to Rs, Rl so that cool skin option can be applied. 
!        5) Set lat=NaN and/or zi=NaN to assign default values if latitude
!           and/or PBL height not given. 
!        6) The code to compute the heat flux caused by precipitation is 
!           included if rain data is available (default is no rain).
!        7) Code updates the cool-skin temperature depression dter and thickness
!           tkt during iteration loop for consistency.
!        8) Number of iterations set to nits = 6.

! Reference:
!
!  Fairall, C.W., E.F. Bradley, J.E. Hare, A.A. Grachev, and J.B. Edson (2003),
!  Bulk parameterization of air sea fluxes: updates and verification for the 
!  COARE algorithm, J. Climate, 16, 571-590.

!  ----------------------------------------------------------

!---> Siqi Li
! I rewrite this part, just to make sure the input variables do not change
! during the calculation. Values of TA and TS do not change, though their 
! intent is 'inout'.
  real, intent(in)    :: UWIND, VWIND, ZU, TA, ZT, QV, ZQ, PA, TS, DSW, LAT
  real, intent(inout) :: DLW
  real, intent(out)   :: HSB, HLB, NET, USTRESS, VSTRESS
  real                :: tsea
!<--- Siqi Li

        real u,rh,zi,rain,cp,sigH
        real Pdef,Rsdef,Rldef,Latdef,Zidef
        real lapse,L10
        data Pdef,Rsdef,Rldef,Latdef,Zidef/ &
             1030,150,  370,  45,    600/
       
! convert input to column vectors
!u=u(:);t=t(:);rh=rh(:);P=P(:);ts=ts(:);
!Rs=Rs(:);Rl=Rl(:);lat=lat(:);zi=zi(:);
!zu=zu(:);zt=zt(:);zq=zq(:);
!rain=rain(:);
!N=length(u);

! set local variables to default values if input is NaN
!        if(p.eq.fmiss) p=pdef
!        if(Rs.eq.fmiss) Rs=Rsdef
!        if(Rl.eq.fmiss) Rl=Rldef
!        if(Lat.eq.fmiss) Lat=Latdef
!        if(Zi.eq.fmiss)Zi=Zidef
         Zi=Zidef
        
        waveage=0
        seastate=0
        if(cp.ne.0) then
           waveage=1
           if(sigh.ne.0) seastate=1
        endif

!************************************************************************
! Check on which parameterization you are using.  You can remove the 
! pause once you are familiar with how the code selects the 
! appropriate parameterization based on the input variables.
!************************************************************************

        if (waveage.eq.1 .and. seastate.eq.1)  &
           print *,' Use the seastate dependent parameterization.'
        if (waveage.eq.1 .and. seastate.eq.0)  &
           print *, ' Use the waveage dependent parameterization.'

! --->Siqi Li
        U=sqrt(UWIND**2+VWIND**2) 

   RH=TA-273.16
   RH=6.112*exp(17.67*RH/(RH+243.5))
   RH=0.622*RH/(P-RH)
   RH=QV*100/RH
   if (RH > 100.) RH=100.

        t=ta-273.16
        tsea=ts-273.16
   RS=dsw
   RL=dlw
!   P=PA
   P=PA/100.      ! Siqi Li, 2022-03-03
! <---Siqi Li

! input variable u is assumed relative wind speed (magnitude of difference
! between wind and surface current vectors). to follow orginal Fairall code, set
! surface current speed us=0. if us data are available, construct u prior to
! using this code.
        us = 0
! convert rh to specific humidity
        Qs = qsat26sea(tsea,P)/1000    ! surface water specific humidity (g/kg)
        CALL qsat26air(t,P,rh,Q,Pv)   ! specific humidity of air (g/kg)
        Q=Q/1000                      ! specific humidity of air (kg/kg)

!***********  set constants **********************************************
        zref=10
        Beta = 1.2
        von  = 0.4
        fdg  = 1.00 ! Turbulent Prandtl number
        tdk  = 273.16
        grav = grvf(lat)

!***********  air constants **********************************************
        Rgas = 287.05
        Le   = (2.501-.00237*tsea)*1e6
        cpa  = 1004.67
        cpv  = cpa*(1+0.84*Q)
        rhoa = P*100./(Rgas*(t+tdk)*(1+0.61*Q))
        rhodry = (P-Pv)*100./(Rgas*(t+tdk))
        visa = 1.326e-5*(1+6.542e-3*t+8.301e-6*t**2-4.84e-9*t**3)
        lapse=grav/cpa

!***********  cool skin constants  ***************************************
        Al   = 2.1e-5*(tsea+3.2)**0.79
        be   = 0.026
        cpw  = 4000
        rhow = 1022
        visw = 1e-6
        tcw  = 0.6
        bigc = 16*grav*cpw*(rhow*visw)**3/(tcw**2*rhoa**2)
        wetc = 0.622*Le*Qs/(Rgas*(ts+tdk)**2)

!***********  net radiation fluxes ***************************************
        Rns = 0.945*Rs ! albedo correction
        Rnl = 0.97*(5.67e-8*(tsea-0.3*jcool+tdk)**4-Rl) ! initial value

!****************  begin bulk loop ********************************************

!***********  first guess ************************************************
        rovcp=Rgas/cpa
        lapse=grav/cpa

        du = u-us
        dt = tsea-t-lapse*zt
        dq = Qs-Q     
!        ta = t+tdk   ! Siqi
        tv = ta*(1+0.61*Q)
        ug = 0.5
        dter  = 0.3
        dqer = dter*wetc 
        ut    = sqrt(du**2+ug**2)
        u10   = ut*log(10/1e-4)/log(zu/1e-4)
        usr   = 0.035*u10
        zo10  = 0.011*usr**2/grav + 0.11*visa/usr
        Cd10  = (von/log(10/zo10))**2
        Ch10  = 0.00115
        Ct10  = Ch10/sqrt(Cd10)
        zot10 = 10/exp(von/Ct10)
        Cd    = (von/log(zu/zo10))**2
        Ct    = von/log(zt/zot10)
        CC    = von*Ct/Cd
        Ribcu = -zu/(zi*0.004*Beta**3)
        Ribu  = -grav*zu/ta*((dt-dter*jcool)+.61*ta*dq)/ut**2 !(ta*(tu**2)) in 26z
        if(Ribu.lt.0) then
           zetu=CC*Ribu/(1+Ribu/Ribcu)
        else
           zetu = CC*Ribu*(1.0+27/9*Ribu/CC)
!          zetu=CC*Ribu/(1.+27/(9*Ribu*cc))   ! in 26z          
        endif
        
!!! Late ---- what for this        
!        k50=find(zetu>50) ! stable with very thin M-O length relative to zu
        
        L10 = zu/zetu
        gf=ut/du
        usr = ut*von/(log(zu/zo10)-psiu_26f(zu/L10))
        
!---> Changed by Siqi Li
!        tsr = -(dt-dter*jcool)*von*fdg/(log(zt/zot10)-psit_26(zt/L10))
        tsr = -(dt-dter*jcool)*von*fdg/(log(zt/zot10)-psit_26f(zt/L10))

!        qsr = -(dq-dqer*jcool)*von*fdg/(log(zq/zot10)-psit_26(zq/L10))
        qsr = -(dq-dqer*jcool)*von*fdg/(log(zq/zot10)-psit_26f(zq/L10))
!<--- Changed by Siqi Li
        tkt = 0.001
      
!**********************************************************
!  The Charnock variable for COARE 3.0
!**********************************************************
        if(ut.gt.10) then 
           charn=0.011+(ut-10)/(18-10)*(0.018-0.011)
        elseif(ut.gt.18) then
           charn=0.018
        else
           charn = 0.011
        endif
!**********************************************************
!  The following gives the new formulation for the
!  Charnock variable in COARE 3.5
!**********************************************************
        charnC=0.11
        umax=19
        a1=0.0017
        a2=-0.0050
        if(u10.gt.umax) then  !wind-speed dependent coefficients
           charnC=a1*umax+a2
        else
           charnC=a1*u10+a2  
        endif
        
        A=0.114  !wave-age dependent coefficients
        B=0.622
        charnW=A*(usr/cp)**B

        Ad=0.091  !Sea-state/wave-age dependent coefficients
        Bd=2.0
        zoS=sigH*Ad*(usr/cp)**Bd
        charnS=zoS*grav/usr/usr

        nits=10 ! number of iterations
                ! Note: nits=1 in 26z_v1.0
!**************  bulk loop **************************************************

        do i=1,nits
           zet=von*grav*zu/tv*(tsr +0.61*ta*qsr)/(usr**2)
           if (waveage.eq.1) then
              if (seastate.eq.1) then
                 charn=charnS
              else
                 charn=charnW
              endif
           else
              charn=charnC
           endif
           L=zu/zet
           zo=charn*usr*usr/grav+0.11*visa/usr ! surface roughness
           rr=zo*usr/visa
           zoq=min(1.6e-4,6e-3/rr**1.6) !These thermal roughness lengths give
           zot=zoq                   !Stanton and Dalton numbers for COARE 4.0
           cdhf=von/(log(zu/zo)-psiu_26f(zu/L))
!---> Changed by Siqi Li
!           cqhf=von*fdg/(log(zq/zoq)-psit_26(zq/L))
           cqhf=von*fdg/(log(zq/zoq)-psit_26f(zq/L))
!           cthf=von*fdg/(log(zt/zot)-psit_26(zt/L))
           cthf=von*fdg/(log(zt/zot)-psit_26f(zt/L))
!<--- Changed by Siqi Li
           usr=ut*cdhf
           qsr=-(dq-dqer*jcool)*cqhf
           tsr=-(dt-dter*jcool)*cthf
           tvsr=tsr+0.61*ta*qsr
           tssr=tsr+0.51*ta*qsr
           Bf=-grav/tv*usr*tvsr
           if(Bf.gt.0) then
              ug=max(0.2,Beta*(Bf*zi)**0.333)
           else
              ug=0.2
           endif
           ut=sqrt(du**2+ug**2)
           gf=ut/du
           hsb=-rhoa*cpa*usr*tsr
           hlb=-rhoa*Le*usr*qsr
           qout=Rnl+hsb+hlb
           dels=Rns*(0.065+11*tkt-6.6e-5/tkt*(1-exp(-tkt/8.0e-4)))
           qcol=qout-dels
           alq=Al*qcol+be*hlb*cpw/Le
           xlamx=6.0
           tkt=min(0.01, xlamx*visw/(sqrt(rhoa/rhow)*usr))
           if(alq.gt.0) then
              xlamx=6/(1+(bigc*alq/usr**4)**0.75)**0.333
              tkt=xlamx*visw/(sqrt(rhoa/rhow)*usr)
           endif
           dter=qcol*tkt/tcw
           sst=tsea-dter  
           dqer=Qs-qsat26sea(sst,P)/1000 !wetc.*dter
           Rnl=0.97*(5.67e-8*(tsea-dter*jcool+tdk)**4-Rl) ! update dter
           if (i.eq.1) then ! save first iteration solution for case of zetu>50
              if(zetu.gt.50) then ! stable with very thin M-O length relative to zu
                 usr50=usr
                 tsr50=tsr
                 qsr50=qsr
                 L50=L
                 zet50=zet
                 dter50=dter
                 dqer50=dqer
                 tkt50=tkt
              endif
           end if
           u10N = usr/von/gf*log(10/zo)
           if (waveage.eq.1) then
              if (seastate.eq.1) then
                 zoS=sigH*Ad*(usr/cp)**Bd-0.11*visa/usr
                 charnS=zoS*grav/usr/usr
              else
                 charnW=A*(usr/cp)**B
              end if
           else
              charnC=a1*u10N+a2
              if(u10N.gt.umax)then
                 charnC=a1*umax+a2
              endif
           end if
        end do
        ! insert first iteration solution for case with zetu>50
        if(zetu.gt.50) then
           usr=usr50
           tsr=tsr50
           qsr=qsr50
           L=L50
           zet=zet50
           dter=dter50
           dqer=dqer50
           tkt=tkt50
        endif
!****************  compute fluxes  ********************************************
        tau=rhoa*usr*usr/gf      ! wind stress
        hsb=-rhoa*cpa*usr*tsr     ! sensible heat flux
        hlb=-rhoa*Le*usr*qsr      ! latent heat flux
        hbb=-rhoa*cpa*usr*tvsr    ! buoyancy flux
        hsbb=-rhoa*cpa*usr*tssr   ! sonic heat flux
        wbar=1.61*hlb/Le/(1+1.61*Q)/rhoa+hsb/rhoa/cpa/ta !Useful for Webb Correction
        !hlwebb=rhoa*wbar*Q*Le
        Evap=1000*hlb/Le/1000*3600   !mm/hour

!*****  compute transfer coeffs relative to ut @ meas. ht  ********************
        Cd=tau/rhoa/ut/max(.1,du)
        Ch=-usr*tsr/ut/(dt-dter*jcool)
        Ce=-usr*qsr/(dq-dqer*jcool)/ut
!***  compute 10-m neutral coeff relative to ut (output if needed) ************
        Cdn_10=1000*von**2./log(10./zo)**2
        Chn_10=1000*von**2*fdg/log(10./zo)/log(10./zot)
        Cen_10=1000*von**2*fdg/log(10./zo)/log(10./zoq)
        
!***  compute 10-m neutral coeff relative to ut (output if needed) ************
!  Find the stability functions
!*********************************
        zrf_u=2             !User defined reference heights to
        zrf_t=2             !compute values at zrf.
        zrf_q=2
        psi=psiu_26f(zu/L)
        psi10=psiu_26f(10./L)
        psirf=psiu_26f(zrf_u/L)
        psiT=psit_26f(zt/L)
        psi10T=psit_26f(10./L)
        psirfT=psit_26f(zrf_t/L)
        psirfQ=psit_26f(zrf_q/L)
        gf=ut/du

!*********************************************************
!  Determine the wind speeds relative to ocean surface
!  Note that usr is the friction velocity that includes 
!  gustiness usr = sqrt(Cd) S, which is equation (18) in
!  Fairall et al. (1996)
!*********************************************************
        S = ut
        U = du
        S10 = S + usr/von*(log(10./zu)-psi10+psi)
        U10 = S10/gf
        ! or U10 = U + usr/von/gf*(log(10/zu)-psi10+psi)
        Urf = U + usr/von/gf*(log(zrf_u/zu)-psirf+psi)
        UN = U + psi*usr/von/gf
        U10N = U10 + psi10*usr/von/gf
        UrfN = Urf + psirf*usr/von/gf
        
        UN2 = usr/von/gf*log(zu/zo)
        U10N2 = usr/von/gf*log(10./zo)
        UrfN2  = usr/von/gf*log(zrf_u/zo)
        
!******** rain heat flux (save to use if desired) *****************************
        if(rain.gt.0) then
           dwat=2.11e-5*((t+tdk)/tdk)**1.94 !! water vapour diffusivity
           dtmp=(1. + 3.309e-3*t - 1.44e-6*t*t)*0.02411/(rhoa*cpa) !! heat diffusivity
           dqs_dt=Q*Le/(Rgas*(t+tdk)**2) !! Clausius-Clapeyron
           alfac= 1./(1+0.622*(dqs_dt*Le*dwat)/(cpa*dtmp)) !! wet bulb factor
           RF= rain*alfac*cpw*((tsea-t-dter*jcool)+ &
                (Qs-Q-dqer*jcool)*Le/cpa)/3600
        else
           RF=0
        end if

        lapse=grav/cpa
        SST=tsea-dter*jcool

        T = t
        T10 = T + tsr/von*(log(10./zt)-psi10T+psiT) + lapse*(zt-10)
        Trf = T + tsr/von*(log(zrf_t/zt)-psirfT+psiT) + lapse*(zt-zrf_t)
        TN = T + psiT*tsr/von
        T10N = T10 + psi10T*tsr/von
        TrfN = Trf + psirfT*tsr/von

        TN2 = SST + tsr/von*log(zt/zot)-lapse*zt
        T10N2 = SST + tsr/von*log(10./zot)-lapse*10
        TrfN2 = SST + tsr/von*log(zrf_t/zot)-lapse*zrf_t

        dqer=(Qs-qsat26sea(SST,P)/1000)*jcool !wetc*dter*jcool
        SSQ=Qs-dqer
        SSQ=SSQ*1000
        Q=Q*1000
        qsr=qsr*1000
        Q10 = Q + qsr/von*(log(10./zq)-psi10T+psiT)
        Qrf = Q + qsr/von*(log(zrf_q/zq)-psirfQ+psiT)
        QN = Q + psiT*qsr/von/sqrt(gf)
        Q10N = Q10 + psi10T*qsr/von
        QrfN = Qrf + psirfQ*qsr/von
        
        QN2 = SSQ + qsr/von*log(zq/zoq)
        Q10N2 = SSQ + qsr/von*log(10./zoq)
        QrfN2 = SSQ + qsr/von*log(zrf_q/zoq)
        RHrf=RHcalc(Trf,P,Qrf/1000)
        RH10=RHcalc(T10,P,Q10/1000)

!--->Siqi Li
        ! Converted the temperature from K into C
!        TA=T+tdk
!        TS=TS+tdk
        ! Heat flux
        hsb=-hsb
        hlb=-hlb
        rl=rl-0.98*5.6697*((TS*0.01)**4)
        dlw=rl
        net=dsw+dlw+hsb+hlb
        ! Calculate the wind stress
        USTRESS = UWIND/U * tau       ! Siqi Li, 2022-03-03
        VSTRESS = VWIND/U * tau
!        USTRESS = UWIND/UR * tau
!        VSTRESS = VWIND/UR * tau
!<---Siqi Li

!****************  output  ****************************************************

!A=[usr tau hsb hlb hbb hsbb wbar  tsr qsr zot zoq Cd Ch Ce  L zet dter dqer tkt
!Urf Trf Qrf RHrf UrfN Rnl Le rhoa UN U10 U10N Cdn_10 Chn_10 Cen_10 RF Qs Evap
!T10 Q10 RH10 gf]
!   1   2   3   4   5   6    7      8   9  10  11  12 13 14 15 16   17   18  19
!   20  21  22  23   24   25 26  27  28  29  30    31     32     33   34 35  36
!   37  38   39  40]
END SUBROUTINE COARE40VN

!------------------------------------------------------------------------------
      function psit_26f(zet)
! computes temperature structure function
        dzet=min(50.,0.35*zet) ! stable
        psit_26f=-((1+0.6667*zet)**1.5+0.6667*(zet-14.28)*exp(-dzet)+8.525)

!        k=find(zet<0) ! unstable
        if(zet.lt.0) then ! unstable
           x=(1-16*zet)**0.5
           psik=2*log((1+x)/2)
           x=(1-34.15*zet)**0.3333
           psic=1.5*log((1+x+x**2)/3.)-sqrt(3.) &
                *atan((1+2*x)/sqrt(3.))+4*atan(1.)/sqrt(3.)
           f=zet**2./(1+zet**2)
           psit_26f=(1-f)*psik+f*psic
        endif
      end function psit_26f
!------------------------------------------------------------------------------
      function psiu_26f(zet)
! computes velocity structure function
        dzet=min(50.,0.35*zet) ! stable
        a=0.7
        b=3./4.
        c=5.
        d=0.35
        psiu_26f=-(a*zet+b*(zet-c/d)*exp(-dzet)+b*c/d)
!        k=find(zet<0) ! unstable
        if(zet.lt.0) then ! unstable
           x=(1-16*zet)**0.25
           psik=2*log((1+x)/2)+log((1+x*x)/2)-2*atan(x)+2*atan(1.)
           x=(1-10.15*zet)**0.3333
           psic=1.5*log((1+x+x**2)/3)-sqrt(3.) &
                *atan((1+2*x)/sqrt(3.))+4*atan(1.)/sqrt(3.)
           f=zet**2./(1+zet**2)
           psiu_26f=(1-f)*psik+f*psic
        endif
      end function psiu_26f

!------------------------------------------------------------------------------
      function psiu_40(zet)                             
! computes velocity structure function
        dzet=min(50.,0.35*zet) ! stable
        a=1
        b=3./4.
        c=5.
        d=0.35
        psiu_40=-(a*zet+b*(zet-c/d)*exp(-dzet)+b*c/d)
        !k=find(zet<0) ! unstable
        if(zet.lt.0) then
           x=(1-18*zet)**0.25
           psik=2*log((1+x)/2)+log((1+x*x)/2)-2*atan(x)+2*atan(1.)
           x=(1-10*zet)**0.3333
           psic=1.5*log((1+x+x**2)/3.)-sqrt(3.0)  &
                *atan((1+2*x)/sqrt(3.))+4*atan(1.)/sqrt(3.)
           f=zet**2./(1+zet**2)
           psiu_40=(1-f)*psik+f*psic
        endif
      end function psiu_40
!------------------------------------------------------------------------------
      function bucksat(T,P)
! computes saturation vapor pressure [mb]
! given T [degC] and P [mb]
        bucksat=6.1121*exp(17.502*T/(T+240.97))*(1.0007+3.46e-6*P)
      end function bucksat
!------------------------------------------------------------------------------
      function qsat26sea(T,P)
! computes surface saturation specific humidity [g/kg]
! given T [degC] and P [mb]
        ex=bucksat(T,P)
        es=0.98*ex ! reduction at sea surface
        qsat26sea=622*es/(P-0.378*es)
      end function qsat26sea
!------------------------------------------------------------------------------
      subroutine qsat26air(T,P,rh,q,em)

! computes saturation specific humidity [g/kg]
! given T [degC] and P [mb]
        es=bucksat(T,P)
        em=0.01*rh*es
        q=622*em/(P-0.378*em)
      end subroutine qsat26air
!------------------------------------------------------------------------------
      function grvf(lat)
! computes g [m/sec**2] given lat in deg
        real lat,pi
        parameter (pi=3.1415926)

        gamma=9.7803267715
        c1=0.0052790414
        c2=0.0000232718
        c3=0.0000001262
        c4=0.0000000007
        phi=lat*pi/180
        x=sin(phi)
        grvf=gamma*(1+c1*x**2+c2*x**4+c3*x**6+c4*x**8)
      end function grvf

!------------------------------------------------------------------------------
      function RHcalc(T,P,Q)
! computes relative humidity given T,P, & Q

        es=6.1121*exp(17.502*T/(T+240.97))*(1.0007+3.46e-6*P)
        em=Q*P/(0.378*Q+0.622)
        RHcalc=100*em/es
      end function RHcalc

!------------------------------------------------------------------------------
! Calculate evaporation from latent heat flux and sst
! (Yu, L., 2007. Global variations in oceanic evaporation (1958–2005): The role 
!  of the changing wind speed. Journal of climate, 20(21), pp.5376-5390.)
! 
! latent      --- latent heat flux           W/m2
! sst         --- sea surface temperature    K
! evaporation --- evaporation                m/s
!
! Siqi Li, SMAST
! 2022-03-29
!------------------------------------------------------------------------------
  SUBROUTINE CALC_EVAPORATION(latent, sst, evaporation)

    implicit none

    real, intent(in   )  :: latent, sst
    real, intent(  out)  :: evaporation
    real, parameter      :: rho0 = 1023  ! (ocean density,kg/m3)

    evaporation = latent / rho0 / ((2.501-0.00237*(sst-273.15))*1e6)
    
  END SUBROUTINE CALC_EVAPORATION

!------------------------------------------------------------------------------
! Calculate precipitation from rainc and rainnc
! (WRF manual)
! 
! rainc1        --- cumulus precipitation at time 1       mm
! rainnc1       --- grid-scale precipitation at time 1    mm
! rainc2        --- cumulus precipitation at time 2       mm
! rainnc2       --- grid-scale precipitation at time 2    mm
! dt            --- time difference                       s
! precipitation --- precipitaion                          m/s
!
! Siqi Li, SMAST
! 2022-03-29
!------------------------------------------------------------------------------
  SUBROUTINE CALC_PRECIPITATION(rainc1, rainnc1, rainc2, rainnc2, dt, precipitation)

    implicit none

    real, intent(in   )  :: rainc1, rainnc1, rainc2, rainnc2, dt
    real, intent(  out)  :: precipitation

    precipitation = ( (rainc2+rainnc2) - (rainc1+rainnc1) ) / dt / 1000.
    
  END SUBROUTINE CALC_PRECIPITATION
  !
END MODULE module_coare
