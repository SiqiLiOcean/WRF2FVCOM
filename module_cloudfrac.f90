!========================================================================
! Calculate the could fraction                                          !
!                                                                       !
!                                                                       !
! Siqi Li, SMAST                                                        !
! 2024-03-01                                                            !
!========================================================================
MODULE module_cloudfrac

CONTAINS

SUBROUTINE CALC_cloudfrac(P, PB, QVAPOR, T, cloudfrac)
  !
  implicit none
  !
  real, dimension(:), intent(in   ) :: P, PB, QVAPOR, T
  real,               intent(  out) :: cloudfrac
  integer                           :: nz
  real, parameter                   :: T_BASE  = 300.0
  real, parameter                   :: RD      = 287.0
  real, parameter                   :: CP      = 1004.5
  real, parameter                   :: P1000MB = 100000.0
  real, dimension(:), allocatable   :: pres, tc, es, ws, rh
  real                              :: lowc, midc, highc
  integer                           :: k, kc_high, kc_mid, kc_low
  

  nz = size(P, 1)
  allocate(pres(nz), tc(nz), es(nz), ws(nz), rh(nz))

  ! Calculate the pressure
  pres = P + PB
  
  ! Calculate the true temperature
  tc = (pres/P1000MB)**(RD/CP) * (T+T_BASE) - 273.15

  ! Calucate the relative humidity
  es = 6.1094 * exp(17.625*tc/(tc+243.04))
  ws = 0.622 * es / (pres/100.-(1-0.622)*es)
  rh = QVAPOR / ws * 100.
  
  ! Calucate the clouds fraction of different layers
  kc_low = 0
  kc_mid = 0
  kc_high = 0
  lowc = 0.0
  midc = 0.0
  highc = 0.0
  DO k = 1, nz-1
    IF ( pres(k) > 97000.) kc_low = k
    IF ( pres(k) > 80000.) kc_mid = k
    IF ( pres(k) > 45000.) kc_high = k
  END DO
  DO k = 1, nz-1
    IF (k>=kc_low .AND. k<kc_mid ) THEN
      lowc = max(rh(k), lowc)
    ELSEIF (k>=kc_mid .AND. k<kc_high) THEN
      midc = max(rh(k), midc)
    ELSEIF (k>=kc_high) THEN
      highc = max(rh(k), highc)
    ENDIF
  END DO
  lowc  = 4.0*lowc/100. - 3.0
  midc  = 4.0*midc/100. - 3.0
  highc = 2.5*highc/100. - 1.5
  lowc  = MIN(lowc, 1.0)
  lowc  = MAX(lowc, 0.0)
  midc  = MIN(midc, 1.0)
  midc  = MAX(midc, 0.0)
  highc = MIN(highc, 1.0)
  highc = MAX(highc, 0.0)

  cloudfrac = 1.0 - (1.0-lowc) * (1.0-midc) * (1.0-highc)

  deallocate(pres, tc, es, ws, rh)
  !
  !
END SUBROUTINE CALC_cloudfrac

END MODULE module_cloudfrac
