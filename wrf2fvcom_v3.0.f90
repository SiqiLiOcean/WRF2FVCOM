!=======================================================================!
! Adjust the WRF output to FVCOM input via COARE                        !
!                                                                       !
! Compilation :                                                         !
!ifort module_nc.f90                        \                           !
!      module_slp.f90                       \                           !
!      module_cloudfrac.f90                 \                           !
!      module_coare.f90                     \                           !
!      wrf2fvcom_v3.0.f90                   \                           !
!      -L${nc_path}/lib -lnetcdff -lnetcdf  \                           !
!      -I${nc_path}/include                 \                           !
!      -o wrf2fvcom                                                     !
!                                                                       !
! Usage :                                                               !
!   ./wrf2fvcom -i input.nc -o output.nc                                !
!   ./wrf2fvcom -i list.dat -o output.nc -l -s                          !
!   ./wrf2fvcom -i list.dat -o output.nc -l -s -ice                     !
!      -i   input file name / filelist                                  !
!      -o   output file name                                            !
!      -l   flag of filelist                      (optional)            !
!      -s   flag of successive files              (optional)            !
!      -x1  start x index                         (optional)            !
!      -nx  x length                              (optional)            !
!      -y1  start y index                         (optional)            !
!      -ny  y length                              (optional)            !
!      -t1  start time index                      (optional)            !
!      -nt  time length                           (optional)            !
!      -ice output variables for the ice module   (optional)            !
!      -slp calculate SLP instead of PSFC         (optional)            !
!      -h   help information                                            !
!                                                                       !
! Siqi Li, SMAST                                                        !
! 2022-03-29                                                            !
!                                                                       !
! Version 3.0                                                           !
!                                                                       !
! Updated:                                                              !
! 2022-09-11  Siqi Li       Added x1,nx,y1,ny                           !
! 2022-09-19  Siqi Li       Changed longitude (LON) to the ranging 0-360!
! 2022-10-07  Chenyu Zhang  Changed NC global attribute "source'        !
! 2022-12-01  Siqi Li       Changed the SLP calculating method          !
! 2023-01-31  Siqi Li       Added variables required by the ice module  !
! 2023-10-26  Siqi Li       SLP is optional to be calculated            !
! 2023-10-30  Siqi Li       Fixed the bug of BUCKET_MM                  !
!=======================================================================!
PROGRAM wrf2fvcom
  !
  use module_nc
  use module_coare
  use module_slp
  use module_cloudfrac
  !
  implicit none
  !
  ! Read the command line
  character(len=200) :: files, fout
  logical            :: flag_list, flag_successive, flag_ice, flag_slp
  integer            :: x1, nx, y1, ny, t1, nt, nz
  ! Variables
!  integer            :: nx, ny
  real, allocatable, dimension(:,:) :: U10, V10, T2, Q2, PSFC, SST, SLP,    &
                                       LONG, SHORT, SENSIBLE, LATENT, NET,  &
                                       USTRESS, VSTRESS,                    &
                                       LON, LAT,                            &
                                       RAINC, RAINNC, RAINC_p, RAINNC_p,    &
                                       precipitation, evaporation,          &
                                       SPQ, SAT, cloud_cover
  real, allocatable, dimension(:,:,:) :: P, PB, T, QVAPOR, PH, PHB
                               integer, allocatable, dimension(:,:) :: I_RAINC, I_RAINNC, I_RAINC_p, I_RAINNC_p
  character(len=19)  :: Times                                     
  real               :: ZU, ZT, ZQ, BUCKET_MM, dt
  
  character(len=200), allocatable :: fin(:)
  type(type_NC)      :: nc
  integer            :: ncid, dimid
  integer            :: nf
  integer, dimension(3):: start, count, start_out, count_out
  integer            :: i, j, k, it, iout
  

  ! Parameters
  ZU = 10.
  ZT = 2.
  ZQ = 2.
  dt = 3600.  ! output time step (s)
  !
  ! Read the command line
  CALL read_args(files, fout, flag_list, flag_successive, flag_ice, flag_slp, x1, nx, y1, ny, t1, nt)
  
  ! Read the input file name(s)
  if (flag_list) then
    nf = 0
    open(11, file=files)
    do while (.not. eof(11))
      read(11, '(A200)')
      nf = nf + 1
    end do
    rewind(11)
    allocate(fin(nf))
    do k = 1, nf
      read(11, '(A200)') fin(k)
    end do
    close(11)

  else
    nf = 1
    allocate(fin(1))
    fin(1) = files
    
  end if

  ! Read the nx and ny dimensions
  if (nx==0) then
    CALL nc_read_dim(fin(1), 'west_east', nx)
    nx = nx - x1 + 1
  end if
  if (ny==0) then 
    CALL nc_read_dim(fin(1), 'south_north', ny)
    ny = ny - y1 + 1
  end if
  CALL nc_read_dim(fin(1), 'bottom_top', nz)

  allocate(lon(nx,ny), lat(nx,ny))
  allocate(U10(nx,ny), V10(nx,ny), T2(nx,ny), Q2(nx,ny), PSFC(nx,ny), SST(nx,ny))
  allocate(LONG(nx,ny), SHORT(nx,ny), SENSIBLE(nx,ny), LATENT(nx,ny), NET(nx,ny))
  allocate(USTRESS(nx,ny), VSTRESS(nx,ny))
  allocate(RAINC(nx,ny), RAINNC(nx,ny), I_RAINC(nx,ny), I_RAINNC(nx,ny))
  allocate(RAINC_p(nx,ny), RAINNC_p(nx,ny), I_RAINC_p(nx,ny), I_RAINNC_p(nx,ny))
  allocate(precipitation(nx,ny), evaporation(nx,ny))
  allocate(P(nx,ny,nz), PB(nx,ny,nz), T(nx,ny,nz), QVAPOR(nx,ny,nz))
  allocate(PH(nx,ny,nz+1), PHB(nx,ny,nz+1))
  allocate(SLP(nx,ny))
  if (flag_ice) allocate(SPQ(nx,ny), SAT(nx,ny), cloud_cover(nx,ny))

  ! Create the output NetCDF file
  nc%name = fout
  ! Dimensions
  nc%dims(1)%name   = 'south_north'
  nc%dims(1)%length = ny
  nc%dims(2)%name   = 'west_east'
  nc%dims(2)%length = nx
  nc%dims(3)%name   = 'DateStrLen'
  nc%dims(3)%length = 19
  nc%dims(4)%name   = 'Time'
  nc%dims(4)%length = -1
  ! Variables
  nc%vars(1)%name       = 'XLONG'
  nc%vars(1)%xtype      = 'float'
  nc%vars(1)%dims(1:2)  = (/'west_east', 'south_north'/)
  nc%vars(2)%name       = 'XLAT'
  nc%vars(2)%xtype      = 'float'
  nc%vars(2)%dims(1:2)  = (/'west_east', 'south_north'/)  
  nc%vars(3)%name       = 'U10'
  nc%vars(3)%xtype      = 'float'
  nc%vars(3)%dims(1:3)  = (/'west_east', 'south_north', 'Time'/)  
  nc%vars(4)%name       = 'V10'
  nc%vars(4)%xtype      = 'float'
  nc%vars(4)%dims(1:3)  = (/'west_east', 'south_north', 'Time'/)  
  nc%vars(5)%name       = 'Stress_U'
  nc%vars(5)%xtype      = 'float'
  nc%vars(5)%dims(1:3)  = (/'west_east', 'south_north', 'Time'/)  
  nc%vars(6)%name       = 'Stress_V'
  nc%vars(6)%xtype      = 'float'
  nc%vars(6)%dims(1:3)  = (/'west_east', 'south_north', 'Time'/)  
  nc%vars(7)%name       = 'Net_Heat'
  nc%vars(7)%xtype      = 'float'
  nc%vars(7)%dims(1:3)  = (/'west_east', 'south_north', 'Time'/)  
  nc%vars(8)%name       = 'Shortwave'
  nc%vars(8)%xtype      = 'float'
  nc%vars(8)%dims(1:3)  = (/'west_east', 'south_north', 'Time'/)  
  nc%vars(9)%name       = 'Longwave'
  nc%vars(9)%xtype      = 'float'
  nc%vars(9)%dims(1:3)  = (/'west_east', 'south_north', 'Time'/)  
  nc%vars(10)%name      = 'Sensible'
  nc%vars(10)%xtype     = 'float'
  nc%vars(10)%dims(1:3) = (/'west_east', 'south_north', 'Time'/)  
  nc%vars(11)%name      = 'Latent'
  nc%vars(11)%xtype     = 'float'
  nc%vars(11)%dims(1:3) = (/'west_east', 'south_north', 'Time'/)  
  nc%vars(12)%name      = 'Evaporation'
  nc%vars(12)%xtype     = 'float'
  nc%vars(12)%dims(1:3) = (/'west_east', 'south_north', 'Time'/)  
  nc%vars(13)%name      = 'Precipitation'
  nc%vars(13)%xtype     = 'float'
  nc%vars(13)%dims(1:3) = (/'west_east', 'south_north', 'Time'/)  
  nc%vars(14)%name      = 'SLP'
  nc%vars(14)%xtype     = 'float'
  nc%vars(14)%dims(1:3) = (/'west_east', 'south_north', 'Time'/)  
  nc%vars(15)%name      = 'Times'
  nc%vars(15)%xtype     = 'char'
  nc%vars(15)%dims(1:2) = (/'DateStrLen', 'Time'/)  
  if (flag_ice) then 
    nc%vars(16)%name      = 'SPQ'
    nc%vars(16)%xtype     = 'float'
    nc%vars(16)%dims(1:3) = (/'west_east', 'south_north', 'Time'/)  
    nc%vars(17)%name      = 'SAT'
    nc%vars(17)%xtype     = 'float'
    nc%vars(17)%dims(1:3) = (/'west_east', 'south_north', 'Time'/)
    nc%vars(18)%name      = 'cloud_cover'
    nc%vars(18)%xtype     = 'float'
    nc%vars(18)%dims(1:3) = (/'west_east', 'south_north', 'Time'/)  
  end if
!  nc%vars(19)%name      = 'T2'
!  nc%vars(19)%xtype     = 'float'
!  nc%vars(19)%dims(1:3) = (/'west_east', 'south_north', 'Time'/)

  ! Define the header
  CALL nc_def_header(nc)
  ! Add attributes
  CALL nc_put_att(fout, 'GLOBAL', 'source', 'wrf2fvcom version 3.0')
  if (flag_slp) then
    CALL nc_put_att(fout, 'GLOBAL', 'SLP', 'calculated')
    
  else
    CALL nc_put_att(fout, 'GLOBAL', 'SLP', 'read from PSFC')
  endif
  CALL nc_put_att(fout, 'XLONG', 'description', 'Longitude')
  CALL nc_put_att(fout, 'XLONG', 'units', 'degree_east')
  CALL nc_put_att(fout, 'XLAT', 'description', 'Latitude')
  CALL nc_put_att(fout, 'XLAT', 'units', 'degree_north')
  CALL nc_put_att(fout, 'U10', 'description', 'U-wind at 10 m')
  CALL nc_put_att(fout, 'U10', 'units', 'm/s')
  CALL nc_put_att(fout, 'V10', 'description', 'V-wind at 10 m')
  CALL nc_put_att(fout, 'V10', 'units', 'm/s')
  CALL nc_put_att(fout, 'Stress_U', 'description', 'U-wind stress at sea surface')
  CALL nc_put_att(fout, 'Stress_U', 'units', 'Pa')
  CALL nc_put_att(fout, 'Stress_V', 'description', 'V-wind stress at sea surface')
  CALL nc_put_att(fout, 'Stress_V', 'units', 'Pa')
  CALL nc_put_att(fout, 'Net_Heat', 'description', 'Net heat flux, positive for ocean gaining energy')
  CALL nc_put_att(fout, 'Net_Heat', 'units', 'W/m2')
  CALL nc_put_att(fout, 'Shortwave', 'description', 'Shortwave radiation, positive for ocean gaining energy')
  CALL nc_put_att(fout, 'Shortwave', 'units', 'W/m2')
  CALL nc_put_att(fout, 'Longwave', 'description', 'Longwave radiation, positive for ocean gaining energy')
  CALL nc_put_att(fout, 'Longwave', 'units', 'W/m2')
  CALL nc_put_att(fout, 'Sensible', 'description', 'Sensible heat flux, positive for ocean gaining energy')
  CALL nc_put_att(fout, 'Sensible', 'units', 'W/m2')
  CALL nc_put_att(fout, 'Latent', 'description', 'Latent heat flux, positive for ocean gaining energy')
  CALL nc_put_att(fout, 'Latent', 'units', 'W/m2')
  CALL nc_put_att(fout, 'Evaporation', 'description', 'Evaporation, positive for ocean gaining water')
  CALL nc_put_att(fout, 'Evaporation', 'units', 'm/s')
  CALL nc_put_att(fout, 'Precipitation', 'description', 'Precipitation, positive for ocean gaining water')
  CALL nc_put_att(fout, 'Precipitation', 'units', 'm/s')
  CALL nc_put_att(fout, 'SLP', 'description', 'Sea level pressure')
  CALL nc_put_att(fout, 'SLP', 'units', 'PA')
!  CALL nc_put_att(fout, 'T2', 'description', 'Temp at 2m')
!  CALL nc_put_att(fout, 'T2', 'units', 'K')
  CALL nc_put_att(fout, 'Times', 'description', 'GMT time')
  CALL nc_put_att(fout, 'Times', 'units', 'yyyy-mm-dd_HH:MM:SS')
  if (flag_ice) then
    CALL nc_put_att(fout, 'SPQ', 'description', 'Specific humidity at 2m')
    CALL nc_put_att(fout, 'SPQ', 'units', 'kg/kg')
    CALL nc_put_att(fout, 'SAT', 'description', 'Aire temperature at 2m')
    CALL nc_put_att(fout, 'SAT', 'units', 'degree C')
    CALL nc_put_att(fout, 'cloud_cover', 'description', 'Total cloud cover')
    CALL nc_put_att(fout, 'cloud_cover', 'units', '1')
  end if

  ! Write longitude and latitude
  CALL nc_read_var(fin(1), 'XLONG', LON, (/x1, y1/), (/nx, ny/))
  WHERE (LON<0.0) LON = LON + 360.0
  CALL nc_read_var(fin(1), 'XLAT', LAT, (/x1, y1/), (/nx, ny/))

  CALL nc_put_var(fout, 'XLONG', LON)
  CALL nc_put_var(fout, 'XLAT', LAT)


  ! Read the Variables
  iout = 0
  do k = 1, nf
    
    write(*,'(A4, A)') '---', trim(fin(k))

    if (nt==0) then
      CALL nc_read_dim(fin(k), 'Time', nt)
      nt = nt - t1 + 1
    end if

 
    CALL nc_read_att(fin(k), 'GLOBAL', 'BUCKET_MM', BUCKET_MM)

    if (t1>1) then
      start = (/x1, y1, t1-1/)
      count = (/nx, ny, 1/)
    else
      start = (/x1, y1, t1/)
      count = (/nx, ny, 1/)
    end if
    if (flag_successive .and. k>1) then
      start = (/x1, y1, nt/)
      count = (/nx, ny, 1/)
      CALL nc_read_var(fin(k-1), 'RAINC',    RAINC_p,    start, count)
      CALL nc_read_var(fin(k-1), 'RAINNC',   RAINNC_p,   start, count)
      if (BUCKET_MM>0) then
        CALL nc_read_var(fin(k-1), 'I_RAINC',  I_RAINC_p,  start, count)
        CALL nc_read_var(fin(k-1), 'I_RAINNC', I_RAINNC_p, start, count)      
      else
        I_RAINC_p = 0
        I_RAINNC_p = 0
      end if
    else
      CALL nc_read_var(fin(k), 'RAINC',    RAINC_p,    start, count)
      CALL nc_read_var(fin(k), 'RAINNC',   RAINNC_p,   start, count)
      if (BUCKET_MM>0) then
        CALL nc_read_var(fin(k), 'I_RAINC',  I_RAINC_p,  start, count)
        CALL nc_read_var(fin(k), 'I_RAINNC', I_RAINNC_p, start, count)
      else
        I_RAINC_p = 0
        I_RAINNC_p = 0
      end if
    end if


    do it = t1, t1+nt-1

      CALL nc_read_var(fin(k), 'Times',    Times,    (/1,it/), (/19,1/))
      print*, '  -', TRIM(Times(1:19))
      start = (/x1, y1, it/)
      count = (/nx, ny, 1/)
      CALL nc_read_var(fin(k), 'U10',      U10,      start, count)
      CALL nc_read_var(fin(k), 'V10',      V10,      start, count)
      CALL nc_read_var(fin(k), 'T2',       T2,       start, count)
      CALL nc_read_var(fin(k), 'Q2',       Q2,       start, count)
      CALL nc_read_var(fin(k), 'PSFC',     PSFC,     start, count)
      CALL nc_read_var(fin(k), 'P',        P,       (/x1, y1, 1, it/),   (/nx, ny, nz, 1/))
      CALL nc_read_var(fin(k), 'PB',       PB,      (/x1, y1, 1, it/),   (/nx, ny, nz, 1/))
      CALL nc_read_var(fin(k), 'T',        T,       (/x1, y1, 1, it/),   (/nx, ny, nz, 1/))
      CALL nc_read_var(fin(k), 'QVAPOR',   QVAPOR,  (/x1, y1, 1, it/),   (/nx, ny, nz, 1/))
      CALL nc_read_var(fin(k), 'PH',       PH,      (/x1, y1, 1, it/), (/nx, ny, nz+1, 1/))
      CALL nc_read_var(fin(k), 'PHB',      PHB,     (/x1, y1, 1, it/), (/nx, ny, nz+1, 1/))
      CALL nc_read_var(fin(k), 'SST',      SST,      start, count)
      CALL nc_read_var(fin(k), 'GLW',      LONG,     start, count)
      CALL nc_read_var(fin(k), 'SWDOWN',   SHORT,    start, count)
      CALL nc_read_var(fin(k), 'RAINC',    RAINC,    start, count)
      CALL nc_read_var(fin(k), 'RAINNC',   RAINNC,   start, count)
      if (BUCKET_MM>0) then
        CALL nc_read_var(fin(k), 'I_RAINC',  I_RAINC,  start, count)
        CALL nc_read_var(fin(k), 'I_RAINNC', I_RAINNC, start, count)
      else
        I_RAINC = 0
        I_RAINNC = 0
      end if

      do i = 1, nx
        do j = 1, ny
          ! Sea-level pressure
          if (flag_slp) then
            CALL CALC_SLP(P(i,j,:), PB(i,j,:), T(i,j,:), QVAPOR(i,j,:), PH(i,j,:), PHB(i,j,:), SLP(i,j))
          else
            SLP = PSFC
          endif
          ! Heat flux and wind stress
          CALL COARE26Z (U10(i,j),V10(i,j),ZU, T2(i,j),ZT, Q2(i,j),ZQ, SLP(i,j), SST(i,j), &
                         LONG(i,j),SHORT(i,j),SENSIBLE(i,j),LATENT(i,j),NET(i,j),    &
                         USTRESS(i,j),VSTRESS(i,j),lat(i,j))
          ! Evaporation
          CALL CALC_evaporation(LATENT(i,j), SST(i,j), evaporation(i,j))
          ! Cloud cover
          if (flag_ice) CALL CALC_cloudfrac(P(i,j,:), PB(i,j,:), QVAPOR(i,j,:), T(i,j,:), cloud_cover(i,j))
        end do
      end do

      ! Precipitation  
      precipitation = ((RAINC-RAINC_p) + BUCKET_MM*(I_RAINC-I_RAINC_p) +  &
                      (RAINNC-RAINNC_p) + BUCKET_MM*(I_RAINNC-I_RAINNC_p))&
                      /dt/1000
      RAINC_p = RAINC
      I_RAINC_p = I_RAINC
      RAINNC_p = RAINNC
      I_RAINNC_p = I_RAINNC

      ! Specific humidity
      if (flag_ice) SPQ = Q2 / (Q2+1)

      ! Aire temperature in degree C
      if (flag_ice) SAT = T2 - 273.15

      iout = iout + 1
      start_out = (/1, 1, iout/)
      count_out = (/nx, ny, 1/)
      CALL nc_put_var(fout, 'U10', U10, start_out, count_out)
      CALL nc_put_var(fout, 'V10', V10, start_out, count_out)
      CALL nc_put_var(fout, 'Stress_U', USTRESS, start_out, count_out)
      CALL nc_put_var(fout, 'Stress_V', VSTRESS, start_out, count_out)
      CALL nc_put_var(fout, 'Net_Heat', NET, start_out, count_out)
      CALL nc_put_var(fout, 'Shortwave', SHORT, start_out, count_out)
      CALL nc_put_var(fout, 'Longwave', LONG, start_out, count_out)
      CALL nc_put_var(fout, 'Sensible', SENSIBLE, start_out, count_out)
      CALL nc_put_var(fout, 'Latent', LATENT, start_out, count_out)
      CALL nc_put_var(fout, 'Evaporation', evaporation, start_out, count_out)
      CALL nc_put_var(fout, 'Precipitation', precipitation, start_out, count_out)
      CALL nc_put_var(fout, 'SLP', SLP, start_out, count_out)
!      CALL nc_put_var(fout, 'T2', T2, start_out, count_out)
      if (flag_ice) then
        CALL nc_put_var(fout, 'SPQ', SPQ, start_out, count_out)
        CALL nc_put_var(fout, 'SAT', SAT, start_out, count_out)
        CALL nc_put_var(fout, 'cloud_cover', cloud_cover, start_out, count_out)
      end if
      CALL nc_put_var(fout, 'Times', Times, (/1,iout/), (/19,1/))

    end do

  end do

  deallocate(fin)
  deallocate(U10, V10, T2, Q2, PSFC, SST)
  deallocate(LONG, SHORT, SENSIBLE, LATENT, NET)
  deallocate(USTRESS, VSTRESS)
  deallocate(lon, lat)
  deallocate(RAINC, RAINNC, I_RAINC, I_RAINNC)
  deallocate(RAINC_p, RAINNC_p, I_RAINC_p, I_RAINNC_p)
  deallocate(precipitation, evaporation)
  deallocate(P, PB, T, QVAPOR, PH, PHB, SLP)
  if (flag_ice) then
    deallocate(SPQ, SAT, cloud_cover)
  endif
  !

CONTAINS

  SUBROUTINE read_args(files, fout, flag_list, flag_successive, flag_ice, flag_slp, x1, nx, y1, ny, t1, nt)
    implicit none

    character(len=200), intent(out) :: files, fout
    logical, intent(out)            :: flag_list, flag_successive, flag_ice, flag_slp
    integer, intent(out)            :: x1, nx, y1, ny, t1, nt
    integer                         :: numarg, i, iargc
    character(len=200)              :: str, dummy

    numarg = iargc()

    if (numarg == 0) then
      CALL help_info
    end if

    i = 1
    fin = ''
    fout = ''
    flag_list = .false.
    flag_successive = .false.
    flag_ice = .false.
    flag_slp = .false.
    x1 = 1
    nx = 0
    y1 = 1
    ny = 0
    t1 = 1
    nt = 0
    do while (i<=numarg)
      call getarg(i, dummy)

      if (dummy(1:1) == '-') then
        SELECT CASE (trim(dummy))
          CASE ('-h')
            CALL help_info
          CASE ('-i')
            i = i + 1
            CALL getarg(i, files)
          CASE ('-o')
            i = i + 1
            CALL getarg(i, fout)
          CASE ('-l')
            flag_list = .true.
          CASE ('-s')
            flag_successive = .true.
          CASE ('-ice')
            flag_ice = .true.
          CASE ('-slp')
            flag_slp = .true.
          CASE ('-x1')
            i = i + 1
            CALL getarg(i, str)
            read(str, *) x1
          CASE ('-nx')
            i = i + 1
            CALL getarg(i, str)
            read(str, *) nx
          CASE ('-y1')
            i = i + 1
            CALL getarg(i, str)
            read(str, *) y1
          CASE ('-ny')
            i = i + 1
            CALL getarg(i, str)
            read(str, *) ny
          CASE ('-t1')
            i = i + 1
            CALL getarg(i, str)
            read(str, *) t1
          CASE ('-nt')
            i = i + 1
            CALL getarg(i, str)
            read(str, *) nt   
          CASE DEFAULT
            CALL help_info
        END SELECT

      else
        CALL help_info
      end if

    i = i + 1
    end do

    if (len_trim(files)==0) stop 'Input file name not specified'
    if (len_trim(fout)==0) stop 'Output file name not specified'

  END SUBROUTINE
      
    

  SUBROUTINE help_info
      print*, 'Usage:'
      print*, '  -i    : input file name / filelist'
      print*, '  -o    : output file name'
      print*, '  -l    : flag of filelist'
      print*, '  -s    : flag of successive files'
      print*, '  -x1   : start x index (optional)'
      print*, '  -nx   : x length (optional)'
      print*, '  -y1   : start y index (optional)'
      print*, '  -ny   : y length (optional)'
      print*, '  -t1   : start time index (optional)'
      print*, '  -nt   : time length (optional)'
      print*, '  -ice  : output variables for ice module'
      print*, '  -slp  : calculate slp'
      print*, '  -h    : help information'
      stop
  END SUBROUTINE
END PROGRAM wrf2fvcom
