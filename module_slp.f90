!========================================================================
! Calculate the sea level pressure                                      !
!                                                                       !
! The original code is from ARWpost                                     !
!                                                                       !
! Siqi Li, SMAST                                                        !
! 2024-03-01                                                            !
!========================================================================
MODULE MODULE_SLP


CONTAINS

!---------------------------------------------------------------------------
! CALC_SLP
! The original version is in module_calc_slp.f90 in ARWPost_3.1
! Download from
! http://www2.mmm.ucar.edu/wrf/users/download/get_sources.html#post_processing
!
! 1) There are five inputs of the column at one grid.
! 2) The vertical dimension lengths of P, PB, T, QVAPOR are bottom_top.
! 3) The vertical dimension lengths of PH, PHB are bottom_top_stag (bottom_top+1).
!
!   Variable  |    Standard Name                                 | Units
!    P(:)     | Perturbation pressure                            | Pa
!    PB(:)    | Base state pressure                              | Pa
!    T(:)     | perturbation potential temperature (theta-t0)    | K
!    QVAPOR(:)| Water vapor mixing ratio                         | kg kg-1
!    PH(:)    | Perturbation geopotential                        | m2 s-2
!    PHB(:)   | Base-state geopotentia                           | m2 s-2
!    SLP(:)   | Sea level pressure                               | Pa
!
! Siqi Li
!---------------------------------------------------------------------------
  SUBROUTINE CALC_SLP(P, PB, T, QVAPOR, PH, PHB, SLP)
    !
    implicit none
    !
    real, intent(in)   :: P(:), PB(:), T(:), QVAPOR(:), PH(:), PHB(:)
    real, intent(out)  :: SLP
    !
    real, parameter    :: G=9.81
    real, parameter    :: PCONST=10000.
    real, parameter    :: GAMMA=0.0065
    real, parameter    :: Rd=287.04
    real, parameter    :: Cp=7.0*Rd/2.
    real, parameter    :: p0=100000.
    real, parameter    :: RCP=Rd/Cp
    logical, parameter :: traditional_comp=.true.
    real, parameter    :: TC=273.16+17.5
    !
    integer            :: nz, level, klo, khi
    real, allocatable  :: PRES(:), GEOPT(:), TK(:), z(:)
    real               :: plo, phi, tlo, thi, zlo, zhi,             &
                          p_at_pconst, t_at_pconst, z_at_pconst,    &
                          t_surf, t_sea_level, z_half_lowest
    logical            :: found, l1, l2, l3
    !
    integer            :: i, k
    !
    !
    nz=size(P)
    !
    allocate(PRES(nz), GEOPT(nz), z(nz), TK(nz))
    PRES=P+PB
    GEOPT=(PHB(1:nz)+PHB(2:nz+1))*0.5 + (PH(1:nz)+PH(2:nz+1))*0.5
    z=GEOPT/G
    do i=1,nz
      TK(i)=(T(i)+300.)*(PRES(i)/p0)**RCP
    end do
    !
    ! 1) Find least zeta level that is PCONST Pa above the surface.
    ! We later use this level to extrapolate a surface pressure and
    ! temperature, which is supposed to reduce the effect of the
    ! diurnal heating cycle in the pressure field.
    level=-1
    k=1
    found=.false.
    do while (.not. found .and. k<=nz)
      if (PRES(k)<PRES(1)-PCONST) then      ! Ki-hwan Kim found the typo
        level=k
        found=.true.
      end if
      k=k+1
    end do

    if (level==-1) then
      write(*,*) "Error in CALC_SLP, error 1."
    end if
    !
    ! 2) Get temperature PCONST Pa above surface. Use this to
    ! extrapolate the temperature at the surface and down to sea level.
    klo=max(level-1,1)
    khi=min(klo+1,nz-1)

    if (klo==khi) then
      write(*,*) "Error in CALC_SLP, error 2."
    end if

    plo=PRES(klo)
    phi=PRES(khi)
    tlo=TK(klo)*(1.+0.608*QVAPOR(klo))
    thi=TK(khi)*(1.+0.608*QVAPOR(khi))
    zlo=z(klo)
    zhi=z(khi)

    p_at_pconst=PRES(1)-PCONST
    t_at_pconst=thi-(thi-tlo)*log(p_at_pconst/phi)*log(plo/phi)
    z_at_pconst=zhi-(zhi-zlo)*log(p_at_pconst/phi)*log(plo/phi)

    t_surf=t_at_pconst*(PRES(1)/p_at_pconst)**(GAMMA*Rd/G)
    t_sea_level=t_at_pconst+GAMMA*z_at_pconst

    ! 3) If we follow a traditional computation, there is a correction
    ! to the sea level temperature if both the surface and sea level
    ! temperatures are *too* hot.
    if (traditional_comp) then
      l1=t_sea_level<TC
      l2=t_surf<=TC
      l3=.not. l1
      if (l2 .and. l3) then
        t_sea_level=TC
      else
        t_sea_level=TC-0.005*(t_surf-TC)**2.
      end if
    end if
    !
    ! 4) Calculate the sea level pressure.
    z_half_lowest=z(1)
    slp=PRES(1)*exp( (2.*G*z_half_lowest) / (Rd*(t_sea_level+t_surf)) )
    !slp=slp/100. ! FVCOM need the pressure in Pa, so DO NOT divide 100.
    !
    !
    deallocate(PRES, GEOPT, z, TK)
    !
  END SUBROUTINE CALC_SLP
  
END MODULE MODULE_SLP
