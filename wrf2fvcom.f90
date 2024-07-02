!=======================================================================!
! Adjust the WRF output to FVCOM input via COARE                        !
!                                                                       !
! Compilation :                                                         !
! ifort module_nc.f90                        \                           !
!       module_slp.f90                       \                           !
!       module_cloudfrac.f90                 \                           !
!       module_coare.f90                     \                           !
!       wrf2fvcom_v3.0.f90                   \                           !
!       -L${nc_path}/lib -lnetcdff -lnetcdf  \                           !
!       -I${nc_path}/include                 \                           !
!       -o wrf2fvcom                                                     !
!                                                                       !
! Usage :                                                               !
!   ./wrf2fvcom -i input.nc -o output.nc                                !
!   ./wrf2fvcom -i list.dat -o output.nc -l -s                          !
!   ./wrf2fvcom -i list.dat -o output.nc -l -s -ice                     !
!      -i          input file name / filelist                           !
!      -o          output file name                                     !
!      -l          flag of filelist                      (optional)     !
!      -s          flag of successive files              (optional)     !
!      -x1         start x index                         (optional)     !
!      -nx         x length                              (optional)     !
!      -y1         start y index                         (optional)     !
!      -ny         y length                              (optional)     !
!      -t1         start time index                      (optional)     !
!      -nt         time length                           (optional)     !
!      -ice        output variables for the ice module   (optional)     !
!      -slp        calculate SLP instead of PSFC         (optional)     !
!      -proj       options in PROJ (use quote)           (optional)     !
!      -v          COARE version, use 2.6 or 4.0(def)    (optional)     !
!      -land_wind  factor of the wind speed on land      (optional)     !
!      -h          help information                                     !
!                                                                       !
! Siqi Li, SMAST                                                        !
! 2022-03-29                                                            !
!                                                                       !
! Version 3.4                                                           !
!                                                                       !
! Updated:                                                              !
! 2022-09-11  Siqi Li       Added x1,nx,y1,ny                           !
! 2022-09-19  Siqi Li       Changed longitude (LON) to the ranging 0-360!
! 2022-10-07  Chenyu Zhang  Changed NC global attribute "source'        !
! 2022-12-01  Siqi Li       Changed the SLP calculating method          !
! 2023-01-31  Siqi Li       Added variables required by the ice module  !
! 2023-10-26  Siqi Li       SLP is optional to be calculated            !
! 2023-10-30  Siqi Li       Fixed the bug of BUCKET_MM                  !
! 2024-02-29  Siqi Li       Added the option of xx and yy               !
! 2024-03-01  Siqi Li       Added the option to select COARE version    !
! 2024-03-15  Siqi Li       Corrected shortwave radiation with albedo   !
! 2024-06-28  Siqi Li       Added the factor of land-wind               !
!=======================================================================!
PROGRAM wrf2fvcom
  !
  USE module_nc
  USE module_coare
  USE module_slp
  USE module_cloudfrac
  !
  IMPLICIT NONE
  !
  ! Read the command line
  CHARACTER(LEN=200)                   :: files
  CHARACTER(LEN=200)                   :: fout
  LOGICAL                              :: flag_list
  LOGICAL                              :: flag_successive
  LOGICAL                              :: flag_ice
  LOGICAL                              :: flag_slp
  CHARACTER(LEN=100)                   :: proj_ref
  CHARACTER(LEN=10)                    :: version
  INTEGER                              :: x1
  INTEGER                              :: nx
  INTEGER                              :: y1
  INTEGER                              :: ny
  INTEGER                              :: t1
  INTEGER                              :: nt
  INTEGER                              :: nz
  REAL                                 :: land_wind
  ! Variables
  REAL, ALLOCATABLE, DIMENSION(:,:)    :: XLAND
  REAL, ALLOCATABLE, DIMENSION(:,:)    :: U10
  REAL, ALLOCATABLE, DIMENSION(:,:)    :: V10
  REAL, ALLOCATABLE, DIMENSION(:,:)    :: T2
  REAL, ALLOCATABLE, DIMENSION(:,:)    :: Q2
  REAL, ALLOCATABLE, DIMENSION(:,:)    :: PSFC
  REAL, ALLOCATABLE, DIMENSION(:,:)    :: SST
  REAL, ALLOCATABLE, DIMENSION(:,:)    :: SLP
  REAL, ALLOCATABLE, DIMENSION(:,:)    :: LONG
  REAL, ALLOCATABLE, DIMENSION(:,:)    :: SHORT
  REAL, ALLOCATABLE, DIMENSION(:,:)    :: SENSIBLE
  REAL, ALLOCATABLE, DIMENSION(:,:)    :: LATENT
  REAL, ALLOCATABLE, DIMENSION(:,:)    :: NET
  REAL, ALLOCATABLE, DIMENSION(:,:)    :: ALBEDO
  REAL, ALLOCATABLE, DIMENSION(:,:)    :: USTRESS
  REAL, ALLOCATABLE, DIMENSION(:,:)    :: VSTRESS
  REAL, ALLOCATABLE, DIMENSION(:,:)    :: LON
  REAL, ALLOCATABLE, DIMENSION(:,:)    :: LAT
  REAL, ALLOCATABLE, DIMENSION(:,:)    :: XX
  REAL, ALLOCATABLE, DIMENSION(:,:)    :: YY
  REAL, ALLOCATABLE, DIMENSION(:,:)    :: RAINC
  REAL, ALLOCATABLE, DIMENSION(:,:)    :: RAINNC
  REAL, ALLOCATABLE, DIMENSION(:,:)    :: RAINC_p
  REAL, ALLOCATABLE, DIMENSION(:,:)    :: RAINNC_p
  REAL, ALLOCATABLE, DIMENSION(:,:)    :: precipitation
  REAL, ALLOCATABLE, DIMENSION(:,:)    :: evaporation
  REAL, ALLOCATABLE, DIMENSION(:,:)    :: SPQ
  REAL, ALLOCATABLE, DIMENSION(:,:)    :: SAT
  REAL, ALLOCATABLE, DIMENSION(:,:)    :: cloud_cover
  REAL, ALLOCATABLE, DIMENSION(:,:,:)  :: P
  REAL, ALLOCATABLE, DIMENSION(:,:,:)  :: PB
  REAL, ALLOCATABLE, DIMENSION(:,:,:)  :: T
  REAL, ALLOCATABLE, DIMENSION(:,:,:)  :: QVAPOR
  REAL, ALLOCATABLE, DIMENSION(:,:,:)  :: PH
  REAL, ALLOCATABLE, DIMENSION(:,:,:)  :: PHB
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: I_RAINC
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: I_RAINNC
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: I_RAINC_p
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: I_RAINNC_p
  CHARACTER(LEN=19)                    :: Times                                     
  REAL                                 :: ZU
  REAL                                 :: ZT
  REAL                                 :: ZQ
  REAL                                 :: BUCKET_MM
  REAL                                 :: dt
  ! Read files
  CHARACTER(LEN=200), ALLOCATABLE      :: fin(:)
  TYPE(type_NC)                        :: nc
  INTEGER                              :: ncid
  INTEGER                              :: dimid
  INTEGER                              :: nf
  INTEGER, DIMENSION(3)                :: start
  INTEGER, DIMENSION(3)                :: count
  INTEGER, DIMENSION(3)                :: start_out
  INTEGER, DIMENSION(3)                :: count_out
  INTEGER                              :: i, j, k, it, iout, iv, id


  ! Parameters
  ZU = 10.
  ZT = 2.
  ZQ = 2.
  dt = 3600.  ! output time step (s)
  
  ! Read the command line
  CALL read_args(files, fout, flag_list, flag_successive, flag_ice, flag_slp, proj_ref, version, x1, nx, y1, ny, t1, nt, land_wind)
  
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
  allocate(XLAND(nx,ny))
  allocate(U10(nx,ny), V10(nx,ny), T2(nx,ny), Q2(nx,ny), PSFC(nx,ny), SST(nx,ny))
  allocate(LONG(nx,ny), SHORT(nx,ny), SENSIBLE(nx,ny), LATENT(nx,ny), NET(nx,ny))
  allocate(ALBEDO(nx,ny))
  allocate(USTRESS(nx,ny), VSTRESS(nx,ny))
  allocate(RAINC(nx,ny), RAINNC(nx,ny), I_RAINC(nx,ny), I_RAINNC(nx,ny))
  allocate(RAINC_p(nx,ny), RAINNC_p(nx,ny), I_RAINC_p(nx,ny), I_RAINNC_p(nx,ny))
  allocate(precipitation(nx,ny), evaporation(nx,ny))
  allocate(P(nx,ny,nz), PB(nx,ny,nz), T(nx,ny,nz), QVAPOR(nx,ny,nz))
  allocate(PH(nx,ny,nz+1), PHB(nx,ny,nz+1))
  allocate(SLP(nx,ny))
  if (flag_ice) allocate(SPQ(nx,ny), SAT(nx,ny), cloud_cover(nx,ny))
  if (LEN_TRIM(proj_ref) /= 0) allocate(xx(nx,ny), yy(nx,ny))

  ! Create the output NetCDF file
  nc%name = fout
  ! Dimensions
  id = 1
  nc%dims(id)%name   = 'south_north'
  nc%dims(id)%length = ny
  id = id + 1
  nc%dims(id)%name   = 'west_east'
  nc%dims(id)%length = nx
  id = id + 1
  nc%dims(id)%name   = 'DateStrLen'
  nc%dims(id)%length = 19
  id = id + 1
  nc%dims(id)%name   = 'Time'
  nc%dims(id)%length = -1
  ! Variables
  iv = 1
  nc%vars(iv)%name      = 'XLONG'
  nc%vars(iv)%xtype     = 'float'
  nc%vars(iv)%dims(1:2) = (/'west_east', 'south_north'/)
  iv = iv + 1
  nc%vars(iv)%name      = 'XLAT'
  nc%vars(iv)%xtype     = 'float'
  nc%vars(iv)%dims(1:2) = (/'west_east', 'south_north'/)  
  iv = iv + 1
  nc%vars(iv)%name      = 'U10'
  nc%vars(iv)%xtype     = 'float'
  nc%vars(iv)%dims(1:3) = (/'west_east', 'south_north', 'Time'/)  
  iv = iv + 1  
  nc%vars(iv)%name      = 'V10'
  nc%vars(iv)%xtype     = 'float'
  nc%vars(iv)%dims(1:3) = (/'west_east', 'south_north', 'Time'/)  
  iv = iv + 1  
  nc%vars(iv)%name      = 'Stress_U'
  nc%vars(iv)%xtype     = 'float'
  nc%vars(iv)%dims(1:3) = (/'west_east', 'south_north', 'Time'/)  
  iv = iv + 1  
  nc%vars(iv)%name      = 'Stress_V'
  nc%vars(iv)%xtype     = 'float'
  nc%vars(iv)%dims(1:3) = (/'west_east', 'south_north', 'Time'/)  
  iv = iv + 1  
  nc%vars(iv)%name      = 'Net_Heat'
  nc%vars(iv)%xtype     = 'float'
  nc%vars(iv)%dims(1:3) = (/'west_east', 'south_north', 'Time'/)  
  iv = iv + 1  
  nc%vars(iv)%name      = 'Shortwave'
  nc%vars(iv)%xtype     = 'float'
  nc%vars(iv)%dims(1:3) = (/'west_east', 'south_north', 'Time'/)  
  iv = iv + 1  
  nc%vars(iv)%name      = 'Longwave'
  nc%vars(iv)%xtype     = 'float'
  nc%vars(iv)%dims(1:3) = (/'west_east', 'south_north', 'Time'/)  
  iv = iv + 1  
  nc%vars(iv)%name      = 'Sensible'
  nc%vars(iv)%xtype     = 'float'
  nc%vars(iv)%dims(1:3) = (/'west_east', 'south_north', 'Time'/)  
  iv = iv + 1  
  nc%vars(iv)%name      = 'Latent'
  nc%vars(iv)%xtype     = 'float'
  nc%vars(iv)%dims(1:3) = (/'west_east', 'south_north', 'Time'/)  
  iv = iv + 1  
  nc%vars(iv)%name      = 'Evaporation'
  nc%vars(iv)%xtype     = 'float'
  nc%vars(iv)%dims(1:3) = (/'west_east', 'south_north', 'Time'/)  
  iv = iv + 1  
  nc%vars(iv)%name      = 'Precipitation'
  nc%vars(iv)%xtype     = 'float'
  nc%vars(iv)%dims(1:3) = (/'west_east', 'south_north', 'Time'/)  
  iv = iv + 1  
  nc%vars(iv)%name      = 'SLP'
  nc%vars(iv)%xtype     = 'float'
  nc%vars(iv)%dims(1:3) = (/'west_east', 'south_north', 'Time'/)  
  iv = iv + 1  
  nc%vars(iv)%name      = 'Times'
  nc%vars(iv)%xtype     = 'char'
  nc%vars(iv)%dims(1:2) = (/'DateStrLen', 'Time'/)  
  if (flag_ice) then  
    iv = iv + 1 
    nc%vars(iv)%name      = 'SPQ'
    nc%vars(iv)%xtype     = 'float'
    nc%vars(iv)%dims(1:3) = (/'west_east', 'south_north', 'Time'/)  
    iv = iv + 1  
    nc%vars(iv)%name      = 'SAT'
    nc%vars(iv)%xtype     = 'float'
    nc%vars(iv)%dims(1:3) = (/'west_east', 'south_north', 'Time'/)  
    iv = iv + 1
    nc%vars(iv)%name      = 'cloud_cover'
    nc%vars(iv)%xtype     = 'float'
    nc%vars(iv)%dims(1:3) = (/'west_east', 'south_north', 'Time'/)  
  end if
  if (LEN_TRIM(proj_ref) /= 0) then  
    iv = iv + 1
    nc%vars(iv)%name      = 'xx'
    nc%vars(iv)%xtype     = 'float'
    nc%vars(iv)%dims(1:2) = (/'west_east', 'south_north'/)  
    iv = iv + 1  
    nc%vars(iv)%name      = 'yy'
    nc%vars(iv)%xtype     = 'float'
    nc%vars(iv)%dims(1:2) = (/'west_east', 'south_north'/)  
  end if 
!  iv = iv + 1
!  nc%vars(iv)%name      = 'T2'
!  nc%vars(iv)%xtype     = 'float'
!  nc%vars(iv)%dims(1:3) = (/'west_east', 'south_north', 'Time'/)

  ! Define the header
  CALL nc_def_header(nc)
  ! Add attributes
  CALL nc_put_att(fout, 'GLOBAL', 'source', 'wrf2fvcom version 3.0')
  if (flag_slp) then
    CALL nc_put_att(fout, 'GLOBAL', 'SLP', 'calculated')
  else
    CALL nc_put_att(fout, 'GLOBAL', 'SLP', 'read from PSFC')
  endif
  SELECT CASE (TRIM(VERSION))
  CASE ('2.6')
    CALL nc_put_att(fout, 'GLOBAL', 'COARE_version', '2.6')
  CASE ('4.0')
    CALL nc_put_att(fout, 'GLOBAL', 'COARE_version', '4.0')
  CASE DEFAULT
    STOP 'UNKNOWN COARE VERSION. USE 2.6 OR 4.0'
  END SELECT
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
  if (LEN_TRIM(proj_ref) /= 0) then
    CALL nc_put_att(fout, 'xx', 'description', 'Cartesian Coordinate X')
    CALL nc_put_att(fout, 'xx', 'units', 'm')
    CALL nc_put_att(fout, 'yy', 'description', 'Cartesian Coordinate Y')
    CALL nc_put_att(fout, 'yy', 'units', 'm')
  end if

  ! Write longitude and latitude
  CALL nc_read_var(fin(1), 'XLONG', LON, (/x1, y1/), (/nx, ny/))
  WHERE (LON<0.0) LON = LON + 360.0
  CALL nc_read_var(fin(1), 'XLAT', LAT, (/x1, y1/), (/nx, ny/))

  CALL nc_put_var(fout, 'XLONG', LON)
  CALL nc_put_var(fout, 'XLAT', LAT)

  ! PROJ
  if (LEN_TRIM(proj_ref) /= 0) then
    CALL GEO2XY(fout, proj_ref, LON, LAT, XX, YY)
    CALL nc_put_var(fout, 'xx', XX)
    CALL nc_put_var(fout, 'yy', YY)
  end if


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
      CALL nc_read_var(fin(k), 'XLAND',    XLAND,    start, count)
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
      CALL nc_read_var(fin(k), 'ALBEDO',   ALBEDO,   start, count)
      CALL nc_read_var(fin(k), 'RAINC',    RAINC,    start, count)
      CALL nc_read_var(fin(k), 'RAINNC',   RAINNC,   start, count)
      if (BUCKET_MM>0) then
        CALL nc_read_var(fin(k), 'I_RAINC',  I_RAINC,  start, count)
        CALL nc_read_var(fin(k), 'I_RAINNC', I_RAINNC, start, count)
      else
        I_RAINC = 0
        I_RAINNC = 0
      end if

      ! Adjust the wind on land 
      WHERE (XLAND==1) U10 = U10 * land_wind
      WHERE (XLAND==1) V10 = V10 * land_wind
          
      do i = 1, nx
        do j = 1, ny
          ! Calculate the net shortwave radiation at surface
          SHORT(i,j) = SHORT(i,j) * (1-ALBEDO(i,j))
          
          ! Sea-level pressure
          if (flag_slp) then
            CALL CALC_SLP(P(i,j,:), PB(i,j,:), T(i,j,:), QVAPOR(i,j,:), PH(i,j,:), PHB(i,j,:), SLP(i,j))
          else
            SLP = PSFC
          endif

          ! Heat flux and wind stress
          SELECT CASE (TRIM(VERSION))
          CASE ('2.6')
            CALL COARE26Z (U10(i,j),V10(i,j),ZU, T2(i,j),ZT, Q2(i,j),ZQ, SLP(i,j), SST(i,j), &
                           LONG(i,j),SHORT(i,j),SENSIBLE(i,j),LATENT(i,j),NET(i,j),          &
                           USTRESS(i,j),VSTRESS(i,j),lat(i,j))
          CASE ('4.0')
            CALL COARE40VN(U10(i,j),V10(i,j),ZU, T2(i,j),ZT, Q2(i,j),ZQ, SLP(i,j), SST(i,j), &
                           LONG(i,j),SHORT(i,j),SENSIBLE(i,j),LATENT(i,j),NET(i,j),          &
                           USTRESS(i,j),VSTRESS(i,j),lat(i,j))
          END SELECT
          
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
  if (flag_ice) deallocate(SPQ, SAT, cloud_cover)
  if (LEN_TRIM(proj_ref) /= 0) deallocate(xx, yy)

CONTAINS

  SUBROUTINE GEO2XY(fout, proj_ref, LON, LAT, XX, YY)

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)      :: fout
    CHARACTER(LEN=*), INTENT(IN)      :: proj_ref
    REAL, DIMENSION(:,:), INTENT(IN)  :: LON
    REAL, DIMENSION(:,:), INTENT(IN)  :: LAT
    REAL, DIMENSION(:,:), INTENT(OUT) :: XX
    REAL, DIMENSION(:,:), INTENT(OUT) :: YY

    CHARACTER(LEN=100)                :: fgeo
    CHARACTER(LEN=100)                :: fxy
    INTEGER                           :: nx
    INTEGER                           :: ny
    CHARACTER(LEN=300)                :: cmd

    INTEGER                           :: i
    INTEGER                           :: j

    nx = size(LON, 1)
    ny = size(LAT, 2)

    fgeo = TRIM(fout)//'_geo'
    fxy = TRIM(fout)//'_xy'
    OPEN(14, FILE=fgeo)
    DO j = 1, ny
      DO i = 1, nx
        WRITE(14, '(2F20.6)') LON(i,j), LAT(i,j)
      END DO
    END DO
    CLOSE(14)

    cmd = 'proj '//TRIM(proj_ref)//' +unit=m -f "%.2f" '//TRIM(fgeo)//' > '//TRIM(fxy)
    CALL SYSTEM(cmd)

    OPEN(15, FILE=fxy)
    DO j = 1, ny
      DO i = 1, nx
        READ(15, *) XX(i,j), YY(i,j)
      END DO
    END DO
    CLOSE(15)

    cmd = 'rm '//TRIM(fgeo)//' '//TRIM(fxy)
    CALL SYSTEM(cmd)

  END SUBROUTINE GEO2XY


  SUBROUTINE read_args(files, fout, flag_list, flag_successive, flag_ice, flag_slp, proj_ref, version, x1, nx, y1, ny, t1, nt, land_wind)
    IMPLICIT NONE

    CHARACTER(len=200), INTENT(out) :: files
    CHARACTER(len=200), INTENT(out) :: fout
    LOGICAL, INTENT(out)            :: flag_list
    LOGICAL, INTENT(out)            :: flag_successive
    LOGICAL, INTENT(out)            :: flag_ice
    LOGICAL, INTENT(out)            :: flag_slp
    CHARACTER(len=*), INTENT(out)   :: proj_ref
    CHARACTER(len=*), INTENT(out)   :: version
    INTEGER, INTENT(out)            :: x1
    INTEGER, INTENT(out)            :: nx
    INTEGER, INTENT(out)            :: y1
    INTEGER, INTENT(out)            :: ny
    INTEGER, INTENT(out)            :: t1
    INTEGER, INTENT(out)            :: nt
    REAL, INTENT(OUT)               :: land_wind
    INTEGER                         :: numarg
    INTEGER                         :: i
    INTEGER                         :: iargc
    CHARACTER(len=200)              :: str
    CHARACTER(len=200)              :: dummy

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
    proj_ref = ''
    version = '4.0'
    x1 = 1
    nx = 0
    y1 = 1
    ny = 0
    t1 = 1
    nt = 0
    land_wind = 1.0
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
          CASE ('-proj')
            i = i + 1
            CALL getarg(i, proj_ref)
          CASE ('-v')
            i = i + 1
            CALL getarg(i, version)
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
          CASE ('-land_wind')
            i = i + 1
            CALL getarg(i, str)
            read(str, *) land_wind           
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
      print*, '  -i         : input file name / filelist'
      print*, '  -o         : output file name'
      print*, '  -l         : flag of filelist'
      print*, '  -s         : flag of successive files'
      print*, '  -x1        : start x index (optional)'
      print*, '  -nx        : x length (optional)'
      print*, '  -y1        : start y index (optional)'
      print*, '  -ny        : y length (optional)'
      print*, '  -t1        : start time index (optional)'
      print*, '  -nt        : time length (optional)'
      print*, '  -ice       : output variables for ice module'
      print*, '  -slp       : calculate slp'
      print*, '  -proj      : options in PROJ (use quote)'
      print*, '  -land_wind : factor of land-wind'
      print*, '  -h         : help information'
      stop
  END SUBROUTINE
END PROGRAM wrf2fvcom
