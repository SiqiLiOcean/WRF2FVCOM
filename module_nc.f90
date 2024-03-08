!========================================================================
! The module to deal with NetCDF format files
! 
! components:
!   File operation
!    --nc_create(file, ncid)
!    --nc_open(file, ncid, flag_write)
!    --nc_close(ncid)
!    --nc_def_header(nc)
!   Mode
!    --nc_redef(ncid)
!    --nc_enddef(ncid)
!   Dimension
!    --nc_def_dim(file, dim_name, dim_length)
!    --nc_read_dim(file, dim_name, dim_length)
!   Varaibel
!    --nc_def_var(file, var_name, var_type, var_dims)
!    --nc_put_var(file, var_name, value, start, count)
!    --nc_read_var(file, var_name, value, start, count)
!   Attribute
!    --nc_put_att(file, var_name, att_name, att_value)
!    --nc_read_att(file, var_name, att_name, att_value)
!   Others     
!    --nc_check(status)
!
!========================================================================
MODULE module_nc 
  !
  USE NETCDF
  !
  IMPLICIT NONE
  !
  INTERFACE nc_put_att
    MODULE PROCEDURE nc_put_att_char
    MODULE PROCEDURE nc_put_att_int
    MODULE PROCEDURE nc_put_att_real
  END INTERFACE

  INTERFACE nc_read_att
    MODULE PROCEDURE nc_read_att_char
    MODULE PROCEDURE nc_read_att_int
    MODULE PROCEDURE nc_read_att_real
  END INTERFACE

  INTERFACE nc_put_var
    MODULE PROCEDURE nc_put_var_char
    MODULE PROCEDURE nc_put_var_int_1d
    MODULE PROCEDURE nc_put_var_int_2d
    MODULE PROCEDURE nc_put_var_int_3d
    MODULE PROCEDURE nc_put_var_int_4d
    MODULE PROCEDURE nc_put_var_real_1d
    MODULE PROCEDURE nc_put_var_real_2d
    MODULE PROCEDURE nc_put_var_real_3d
    MODULE PROCEDURE nc_put_var_real_4d
  END INTERFACE

  INTERFACE nc_read_var
    MODULE PROCEDURE nc_read_var_char
    MODULE PROCEDURE nc_read_var_int_1d
    MODULE PROCEDURE nc_read_var_int_2d
    MODULE PROCEDURE nc_read_var_int_3d
    MODULE PROCEDURE nc_read_var_int_4d
    MODULE PROCEDURE nc_read_var_real_1d
    MODULE PROCEDURE nc_read_var_real_2d
    MODULE PROCEDURE nc_read_var_real_3d
    MODULE PROCEDURE nc_read_var_real_4d
  END INTERFACE

  TYPE type_ATT
    character(len=20)              :: name
    character(len=6)               :: xtype
    integer                        :: value_int
    real                           :: value_real
    character(len=100)             :: value_char
  END TYPE type_ATT

  TYPE type_DIM
    character(len=20)              :: name
    integer                        :: length
  END TYPE type_DIM

  TYPE type_VAR
    character(len=20)              :: name
    character(len=6)               :: xtype
    character(len=20)              :: dims(5)
    type(type_ATT)                 :: atts(300)
  END TYPE type_VAR

  TYPE type_NC
    character(len=200)             :: name
    type(type_DIM)    :: dims(20)
    type(type_VAR)    :: vars(300)
    type(type_ATT)    :: atts(300)
  END TYPE type_NC
 
CONTAINS

  SUBROUTINE nc_create(file, ncid)
    !
    character(len=*), intent(in) :: file
    integer, intent(out)         :: ncid
    !
    call nc_check( nf90_create(file, nf90_clobber, ncid) )
    call nc_check( nf90_close(ncid) )
    !
  END SUBROUTINE

  SUBROUTINE nc_open(file, ncid, flag_write)
    !
    character(len=*), intent(in) :: file
    integer, intent(out)         :: ncid
    integer, intent(in), optional:: flag_write
    integer                      :: mode
    !
    if (present(flag_write) .and. flag_write==1) then
      mode = nf90_write 
    else
      mode = nf90_nowrite
    end if
    call nc_check( nf90_open(file, mode, ncid) )
    !
  END SUBROUTINE

  SUBROUTINE nc_close(ncid)
    !
    integer, intent(in)          :: ncid
    !
    call nc_check( nf90_close(ncid) )
    !
  END SUBROUTINE
  !
  SUBROUTINE nc_def_header(nc)
    !
    type(type_nc), intent(in) :: nc
    integer                   :: ncid
    integer, allocatable      :: dimids(:)
    integer                   :: n_dim, n_var, n_att
    integer                   :: i, j

    CALL nc_create(nc%name, ncid)
    ! CALL nc_close(ncid)

 !   CALL system('chmod u+w '//trim(nc%name))

 
    ! Set global att
    n_att = count(len_trim(nc%atts(:)%name)<20)
    do i = 1, n_att
      SELECT CASE (trim(nc%atts(i)%xtype))
      CASE ('CHAR', 'char')
        CALL nc_put_att(nc%name, 'GLOBAL',  &
              trim(nc%atts(i)%name), trim(nc%atts(i)%value_char))
      CASE ('INT', 'int')
        CALL nc_put_att(nc%name, 'GLOBAL',  &
              trim(nc%atts(i)%name), nc%atts(i)%value_int)
      CASE ('FLOAT', 'float')
        CALL nc_put_att(nc%name, 'GLOBAL',  &
              trim(nc%atts(i)%name), nc%atts(i)%value_real)
      CASE DEFAULT
        print*, 'Unknown type'//trim(nc%atts(i)%xtype)
        STOP 
      END SELECT
    end do

    ! Set dimensions
    n_dim = count(len_trim(nc%dims(:)%name)<20)
    do i = 1, n_dim
      CALL nc_def_dim(nc%name, nc%dims(i)%name, nc%dims(i)%length)
    end do

    ! Set variables
    n_var = count(len_trim(nc%vars(:)%name)<20) 
    do i = 1, n_var
      n_dim = count(len_trim(nc%vars(i)%dims)<20)
      allocate(dimids(n_dim))
      CALL nc_open(nc%name, ncid)
      do j = 1, n_dim
        CALL nc_check( nf90_inq_dimid(ncid, trim(nc%vars(i)%dims(j)), dimids(j)) )
      end do
      CALL nc_close(ncid)

      CALL nc_def_var(nc%name, nc%vars(i)%name, nc%vars(i)%xtype, dimids)
      deallocate(dimids)

      ! Set variable att
      n_att = count(len_trim(nc%vars(i)%atts%name)<20)
      do j = 1, n_att
        SELECT CASE (trim(nc%vars(i)%atts(j)%xtype))
        CASE ('CHAR', 'char')
          CALL nc_put_att(nc%name, nc%vars(i)%name,  &
                trim(nc%vars(i)%atts(j)%name), trim(nc%vars(i)%atts(j)%value_char))
        CASE ('INT', 'int')
          CALL nc_put_att(nc%name, nc%vars(i)%name,  &
                trim(nc%vars(i)%atts(j)%name), nc%vars(i)%atts(j)%value_int)
        CASE ('FLOAT', 'float')
          CALL nc_put_att(nc%name, nc%vars(i)%name,  &
                trim(nc%vars(i)%atts(j)%name), nc%vars(i)%atts(j)%value_real)
        CASE DEFAULT
          print*, 'Unknown type'//trim(nc%vars(i)%atts(j)%xtype)//'of var '//trim(nc%vars(i)%name)
          STOP 
        END SELECT
      end do

    end do
      
    !
  END SUBROUTINE
  !
  SUBROUTINE nc_redef(ncid)
    !
    integer, intent(in)         :: ncid
    !
    call nc_check( nf90_redef(ncid) )
    !
  END SUBROUTINE

  SUBROUTINE nc_enddef(ncid)
    !
    integer, intent(in)         :: ncid
    !
    call nc_check( nf90_enddef(ncid) )
    !
  END SUBROUTINE
  !
  SUBROUTINE nc_def_dim(file, dim_name, dim_length)
    !
    character(len=*), intent(in)  :: file
    character(len=*), intent(in)  :: dim_name
    integer, intent(in)           :: dim_length
    integer                       :: ncid, dimid
    !
    CALL nc_open(file, ncid, 1)

    CALL nc_redef(ncid)

    if (dim_length>0) then
      call nc_check( nf90_def_dim(ncid, dim_name, dim_length, dimid) )
    else
      call nc_check( nf90_def_dim(ncid, dim_name, nf90_unlimited, dimid) )
    endif

    CALL nc_enddef(ncid)

    CALL nc_close(ncid)
    !
  END SUBROUTINE
  !
  SUBROUTINE nc_read_dim(file, dim_name, dim_length)
    !
    character(len=*), intent(in) :: file
    character(len=*), intent(in) :: dim_name
    integer, intent(out)         :: dim_length
    integer                      :: ncid, dimid
    !
    call nc_open(file, ncid)

    call nc_check( nf90_inq_dimid(ncid, dim_name, dimid) )

    call nc_check( nf90_inquire_dimension(ncid, dimid, len=dim_length) )

    call nc_close(ncid)
    !
  END SUBROUTINE
  !
  SUBROUTINE nc_def_var(file, var_name, var_type, var_dims)
    !
    character(len=*), intent(in)     :: file
    character(len=*), intent(in)     :: var_name
    character(len=*), intent(in)     :: var_type
    integer, intent(in)              :: var_dims(:)          
    integer                          :: ncid, varid
    !
    CALL nc_open(file, ncid, 1)

    CALL nc_redef(ncid)

    SELECT CASE(trim(var_type))
    CASE ('CHAR', 'char')
      call nc_check( nf90_def_var(ncid, var_name, nf90_char, var_dims, varid) )
    CASE ('INT', 'int')
      call nc_check( nf90_def_var(ncid, var_name, nf90_int, var_dims, varid) )
    CASE ('FLOAT', 'float')
      call nc_check( nf90_def_var(ncid, var_name, nf90_float, var_dims, varid) )
    CASE ('DOUBLE', 'double')
      call nc_check( nf90_def_var(ncid, var_name, nf90_double, var_dims, varid) )
    CASE DEFAULT
      print*, '|', var_type ,'|'
      print*, 'Wrong variable type. Use CHAR, INT, FLOAT, or DOUBLE.'
      stop
    END SELECT

    CALL nc_enddef(ncid)

    CALL nc_close(ncid)
    !
  END SUBROUTINE
  !
  SUBROUTINE nc_put_att_char(file, var_name, att_name, att_value)
    !
    character(len=*), intent(in)      :: file
    character(len=*), intent(in)      :: var_name
    character(len=*), intent(in)      :: att_name
    character(len=*), intent(in)      :: att_value
    integer                           :: ncid, varid
    !
    CALL nc_open(file, ncid, 1)
    CALL nc_redef(ncid)

    if (trim(var_name) == 'GLOBAL') then
      varid = nf90_global
    else
      CALL nc_check( nf90_inq_varid(ncid, var_name, varid) )
    end if

    call nc_check( nf90_put_att(ncid, varid, att_name, att_value) )

    CALL nc_enddef(ncid)

    CALL nc_close(ncid)
    !
  END SUBROUTINE 
  !
  SUBROUTINE nc_put_att_int(file, var_name, att_name, att_value)
    !
    character(len=*), intent(in)      :: file
    character(len=*), intent(in)      :: var_name
    character(len=*), intent(in)      :: att_name
    integer, intent(in)               :: att_value
    integer                           :: ncid, varid
    !
    CALL nc_open(file, ncid, 1)

    CALL nc_redef(ncid)

    if (trim(var_name) == 'GLOBAL') then
      varid = nf90_global
    else
      CALL nc_check( nf90_inq_varid(ncid, var_name, varid) )
    end if

    call nc_check( nf90_put_att(ncid, varid, att_name, att_value) )

    CALL nc_enddef(ncid)

    CALL nc_close(ncid)
    !
  END SUBROUTINE
  !
  SUBROUTINE nc_put_att_real(file, var_name, att_name, att_value)
    !
    character(len=*), intent(in)      :: file
    character(len=*), intent(in)      :: var_name
    character(len=*), intent(in)      :: att_name
    real, intent(in)                  :: att_value
    integer                           :: ncid, varid
    !
    CALL nc_open(file, ncid, 1)
   
    CALL nc_redef(ncid)
    
    if (trim(var_name) == 'GLOBAL') then
      varid = nf90_global
    else
      CALL nc_check( nf90_inq_varid(ncid, var_name, varid) )
    end if

    call nc_check( nf90_put_att(ncid, varid, att_name, att_value) )

    CALL nc_enddef(ncid)

    CALL nc_close(ncid)
    !
  END SUBROUTINE
  !
  SUBROUTINE nc_read_att_char(file, var_name, att_name, att_value)
    !
    character(len=*), intent(in)      :: file
    character(len=*), intent(in)      :: var_name
    character(len=*), intent(in)      :: att_name
    character(len=*), intent(out)     :: att_value
    integer  :: ncid, varid
    !
    CALL nc_open(file, ncid)

    if (trim(var_name) == 'GLOBAL') then
      varid = nf90_global
    else
      CALL nc_check( nf90_inq_varid(ncid, var_name, varid) )
    end if

    CALL nc_check( nf90_get_att(ncid, varid, att_name, att_value) )

    CALL nc_close(ncid)
    !
  END SUBROUTINE 
  !
  SUBROUTINE nc_read_att_int(file, var_name, att_name, att_value)
    !
    character(len=*), intent(in)      :: file
    character(len=*), intent(in)      :: var_name
    character(len=*), intent(in)      :: att_name
    integer, intent(out)              :: att_value
    integer  :: ncid, varid
    !
    CALL nc_open(file, ncid)

    if (trim(var_name) == 'GLOBAL') then
      varid = nf90_global
    else
      CALL nc_check( nf90_inq_varid(ncid, var_name, varid) )
    end if

    CALL nc_check( nf90_get_att(ncid, varid, att_name, att_value) )

    CALL nc_close(ncid)
    !
  END SUBROUTINE 
  !
  SUBROUTINE nc_read_att_real(file, var_name, att_name, att_value)
    !
    character(len=*), intent(in)      :: file
    character(len=*), intent(in)      :: var_name
    character(len=*), intent(in)      :: att_name
    real, intent(out)                 :: att_value
    integer  :: ncid, varid
    !
    CALL nc_open(file, ncid)

    if (trim(var_name) == 'GLOBAL') then
      varid = nf90_global
    else
      CALL nc_check( nf90_inq_varid(ncid, var_name, varid) )
    end if

    CALL nc_check( nf90_get_att(ncid, varid, att_name, att_value) )

    CALL nc_close(ncid)
    !
  END SUBROUTINE 
  !
  SUBROUTINE nc_put_var_char(file, var_name, value, start, count)
    !
    character(len=*), intent(in)      :: file
    character(len=*), intent(in)      :: var_name
    character(len=*), intent(in)      :: value
    integer, optional, intent(in)     :: start(:)
    integer, optional, intent(in)     :: count(:)
    integer                           :: ncid, varid
    !
    call nc_open(file, ncid, 1)

    call nc_check( nf90_inq_varid(ncid, var_name, varid) )

    if (present(start)) then
      if (present(count)) then
        call nc_check( nf90_put_var(ncid, varid, value, start=start, count=count) )
      else
        call nc_check( nf90_put_var(ncid, varid, value, start=start) )
      endif
    else
      if (present(count)) then
        call nc_check( nf90_put_var(ncid, varid, value, count=count) )
      else
        call nc_check( nf90_put_var(ncid, varid, value) )
      endif
    endif

    CALL nc_close(ncid)
    !
  END SUBROUTINE 
  !
  SUBROUTINE nc_put_var_int_1d(file, var_name, value, start, count)
    !
    character(len=*), intent(in)      :: file
    character(len=*), intent(in)      :: var_name
    integer, intent(in)               :: value(:)
    integer, optional, intent(in)     :: start(:)
    integer, optional, intent(in)     :: count(:)
    integer                           :: ncid, varid
    !
    call nc_open(file, ncid, 1)

    call nc_check( nf90_inq_varid(ncid, var_name, varid) )

    if (present(start)) then
      if (present(count)) then
        call nc_check( nf90_put_var(ncid, varid, value, start=start, count=count) )
      else
        call nc_check( nf90_put_var(ncid, varid, value, start=start) )
      endif
    else
      if (present(count)) then
        call nc_check( nf90_put_var(ncid, varid, value, count=count) )
      else
        call nc_check( nf90_put_var(ncid, varid, value) )
      endif
    endif

    CALL nc_close(ncid)
    !
  END SUBROUTINE 
  !
  SUBROUTINE nc_put_var_int_2d(file, var_name, value, start, count)
    !
    character(len=*), intent(in)      :: file
    character(len=*), intent(in)      :: var_name
    integer, intent(in)               :: value(:,:)
    integer, optional, intent(in)     :: start(:)
    integer, optional, intent(in)     :: count(:)
    integer                           :: ncid, varid
    !
    call nc_open(file, ncid, 1)

    call nc_check( nf90_inq_varid(ncid, var_name, varid) )

    if (present(start)) then
      if (present(count)) then
        call nc_check( nf90_put_var(ncid, varid, value, start=start, count=count) )
      else
        call nc_check( nf90_put_var(ncid, varid, value, start=start) )
      endif
    else
      if (present(count)) then
        call nc_check( nf90_put_var(ncid, varid, value, count=count) )
      else
        call nc_check( nf90_put_var(ncid, varid, value) )
      endif
    endif

    CALL nc_close(ncid)
    !
  END SUBROUTINE 
  !
  SUBROUTINE nc_put_var_int_3d(file, var_name, value, start, count)
    !
    character(len=*), intent(in)      :: file
    character(len=*), intent(in)      :: var_name
    integer, intent(in)               :: value(:,:,:)
    integer, optional, intent(in)     :: start(:)
    integer, optional, intent(in)     :: count(:)
    integer                           :: ncid, varid
    !
    call nc_open(file, ncid, 1)

    call nc_check( nf90_inq_varid(ncid, var_name, varid) )

    if (present(start)) then
      if (present(count)) then
        call nc_check( nf90_put_var(ncid, varid, value, start=start, count=count) )
      else
        call nc_check( nf90_put_var(ncid, varid, value, start=start) )
      endif
    else
      if (present(count)) then
        call nc_check( nf90_put_var(ncid, varid, value, count=count) )
      else
        call nc_check( nf90_put_var(ncid, varid, value) )
      endif
    endif

    CALL nc_close(ncid)
    !
  END SUBROUTINE 
  !
  SUBROUTINE nc_put_var_int_4d(file, var_name, value, start, count)
    !
    character(len=*), intent(in)      :: file
    character(len=*), intent(in)      :: var_name
    integer, intent(in)               :: value(:,:,:,:)
    integer, optional, intent(in)     :: start(:)
    integer, optional, intent(in)     :: count(:)
    integer                           :: ncid, varid
    !
    call nc_open(file, ncid, 1)

    call nc_check( nf90_inq_varid(ncid, var_name, varid) )

    if (present(start)) then
      if (present(count)) then
        call nc_check( nf90_put_var(ncid, varid, value, start=start, count=count) )
      else
        call nc_check( nf90_put_var(ncid, varid, value, start=start) )
      endif
    else
      if (present(count)) then
        call nc_check( nf90_put_var(ncid, varid, value, count=count) )
      else
        call nc_check( nf90_put_var(ncid, varid, value) )
      endif
    endif

    CALL nc_close(ncid)
    !
  END SUBROUTINE 
  !
  SUBROUTINE nc_put_var_real_1d(file, var_name, value, start, count)
    !
    character(len=*), intent(in)      :: file
    character(len=*), intent(in)      :: var_name
    real, intent(in)                  :: value(:)
    integer, optional, intent(in)     :: start(:)
    integer, optional, intent(in)     :: count(:)
    integer                           :: ncid, varid
    !
    call nc_open(file, ncid, 1)

    call nc_check( nf90_inq_varid(ncid, var_name, varid) )

    if (present(start)) then
      if (present(count)) then
        call nc_check( nf90_put_var(ncid, varid, value, start=start, count=count) )
      else
        call nc_check( nf90_put_var(ncid, varid, value, start=start) )
      endif
    else
      if (present(count)) then
        call nc_check( nf90_put_var(ncid, varid, value, count=count) )
      else
        call nc_check( nf90_put_var(ncid, varid, value) )
      endif
    endif

    CALL nc_close(ncid)
    !
  END SUBROUTINE 
  !
  SUBROUTINE nc_put_var_real_2d(file, var_name, value, start, count)
    !
    character(len=*), intent(in)      :: file
    character(len=*), intent(in)      :: var_name
    real, intent(in)                  :: value(:,:)
    integer, optional, intent(in)     :: start(:)
    integer, optional, intent(in)     :: count(:)
    integer                           :: ncid, varid
    !
    call nc_open(file, ncid, 1)

    call nc_check( nf90_inq_varid(ncid, var_name, varid) )

    if (present(start)) then
      if (present(count)) then
        call nc_check( nf90_put_var(ncid, varid, value, start=start, count=count) )
      else
        call nc_check( nf90_put_var(ncid, varid, value, start=start) )
      endif
    else
      if (present(count)) then
        call nc_check( nf90_put_var(ncid, varid, value, count=count) )
      else
        call nc_check( nf90_put_var(ncid, varid, value) )
      endif
    endif

    CALL nc_close(ncid)
    !
  END SUBROUTINE 
  !
  SUBROUTINE nc_put_var_real_3d(file, var_name, value, start, count)
    !
    character(len=*), intent(in)      :: file
    character(len=*), intent(in)      :: var_name
    real, intent(in)                  :: value(:,:,:)
    integer, optional, intent(in)     :: start(:)
    integer, optional, intent(in)     :: count(:)
    integer                           :: ncid, varid
    !
    call nc_open(file, ncid, 1)

    call nc_check( nf90_inq_varid(ncid, var_name, varid) )

    if (present(start)) then
      if (present(count)) then
        call nc_check( nf90_put_var(ncid, varid, value, start=start, count=count) )
      else
        call nc_check( nf90_put_var(ncid, varid, value, start=start) )
      endif
    else
      if (present(count)) then
        call nc_check( nf90_put_var(ncid, varid, value, count=count) )
      else
        call nc_check( nf90_put_var(ncid, varid, value) )
      endif
    endif

    CALL nc_close(ncid)
    !
  END SUBROUTINE 
  !
  SUBROUTINE nc_put_var_real_4d(file, var_name, value, start, count)
    !
    character(len=*), intent(in)      :: file
    character(len=*), intent(in)      :: var_name
    real, intent(in)                  :: value(:,:,:,:)
    integer, optional, intent(in)     :: start(:)
    integer, optional, intent(in)     :: count(:)
    integer                           :: ncid, varid
    !
    call nc_open(file, ncid, 1)

    call nc_check( nf90_inq_varid(ncid, var_name, varid) )

    if (present(start)) then
      if (present(count)) then
        call nc_check( nf90_put_var(ncid, varid, value, start=start, count=count) )
      else
        call nc_check( nf90_put_var(ncid, varid, value, start=start) )
      endif
    else
      if (present(count)) then
        call nc_check( nf90_put_var(ncid, varid, value, count=count) )
      else
        call nc_check( nf90_put_var(ncid, varid, value) )
      endif
    endif

    CALL nc_close(ncid)
    !
  END SUBROUTINE 
  !
  SUBROUTINE nc_read_var_char(file, var_name, value, start, count)
    !
    character(len=*), intent(in)      :: file
    character(len=*), intent(in)      :: var_name
    integer, optional, intent(in)     :: start(:)
    integer, optional, intent(in)     :: count(:)
    character(len=*), intent(out)     :: value
    integer                           :: ncid, varid
    !
    call nc_open(file, ncid)

    call nc_check( nf90_inq_varid(ncid, var_name, varid) )

    if (present(start)) then
      if (present(count)) then
        call nc_check( nf90_get_var(ncid, varid, value, start=start, count=count) )
      else
        call nc_check( nf90_get_var(ncid, varid, value, start=start) )
      endif
    else
      if (present(count)) then
        call nc_check( nf90_get_var(ncid, varid, value, count=count) )
      else
        call nc_check( nf90_get_var(ncid, varid, value) )
      endif
    endif

    call nc_close(ncid)
    !
  END SUBROUTINE 
  !
  SUBROUTINE nc_read_var_int_1d(file, var_name, value, start, count)
    !
    character(len=*), intent(in)      :: file
    character(len=*), intent(in)      :: var_name
    integer, optional, intent(in)     :: start(:)
    integer, optional, intent(in)     :: count(:)
    integer, intent(out)              :: value(:)
    integer                           :: ncid, varid
    !
    call nc_open(file, ncid)

    call nc_check( nf90_inq_varid(ncid, var_name, varid) )

    if (present(start)) then
      if (present(count)) then
        call nc_check( nf90_get_var(ncid, varid, value, start=start, count=count) )
      else
        call nc_check( nf90_get_var(ncid, varid, value, start=start) )
      endif
    else
      if (present(count)) then
        call nc_check( nf90_get_var(ncid, varid, value, count=count) )
      else
        call nc_check( nf90_get_var(ncid, varid, value) )
      endif
    endif

    call nc_close(ncid)
    !
  END SUBROUTINE 
  !
  SUBROUTINE nc_read_var_int_2d(file, var_name, value, start, count)
    !
    character(len=*), intent(in)      :: file
    character(len=*), intent(in)      :: var_name
    integer, optional, intent(in)     :: start(:)
    integer, optional, intent(in)     :: count(:)
    integer, intent(out)              :: value(:,:)
    integer                           :: ncid, varid
    !
    call nc_open(file, ncid)

    call nc_check( nf90_inq_varid(ncid, var_name, varid) )

    if (present(start)) then
      if (present(count)) then
        call nc_check( nf90_get_var(ncid, varid, value, start=start, count=count) )
      else
        call nc_check( nf90_get_var(ncid, varid, value, start=start) )
      endif
    else
      if (present(count)) then
        call nc_check( nf90_get_var(ncid, varid, value, count=count) )
      else
        call nc_check( nf90_get_var(ncid, varid, value) )
      endif
    endif

    call nc_close(ncid)
    !
  END SUBROUTINE 
  !
  SUBROUTINE nc_read_var_int_3d(file, var_name, value, start, count)
    !
    character(len=*), intent(in)      :: file
    character(len=*), intent(in)      :: var_name
    integer, optional, intent(in)     :: start(:)
    integer, optional, intent(in)     :: count(:)
    integer, intent(out)              :: value(:,:,:)
    integer                           :: ncid, varid
    !
    call nc_open(file, ncid)

    call nc_check( nf90_inq_varid(ncid, var_name, varid) )

    if (present(start)) then
      if (present(count)) then
        call nc_check( nf90_get_var(ncid, varid, value, start=start, count=count) )
      else
        call nc_check( nf90_get_var(ncid, varid, value, start=start) )
      endif
    else
      if (present(count)) then
        call nc_check( nf90_get_var(ncid, varid, value, count=count) )
      else
        call nc_check( nf90_get_var(ncid, varid, value) )
      endif
    endif

    call nc_close(ncid)
    !
  END SUBROUTINE 
  !
  SUBROUTINE nc_read_var_int_4d(file, var_name, value, start, count)
    !
    character(len=*), intent(in)      :: file
    character(len=*), intent(in)      :: var_name
    integer, optional, intent(in)     :: start(:)
    integer, optional, intent(in)     :: count(:)
    integer, intent(out)              :: value(:,:,:,:)
    integer                           :: ncid, varid
    !
    call nc_open(file, ncid)

    call nc_check( nf90_inq_varid(ncid, var_name, varid) )

    if (present(start)) then
      if (present(count)) then
        call nc_check( nf90_get_var(ncid, varid, value, start=start, count=count) )
      else
        call nc_check( nf90_get_var(ncid, varid, value, start=start) )
      endif
    else
      if (present(count)) then
        call nc_check( nf90_get_var(ncid, varid, value, count=count) )
      else
        call nc_check( nf90_get_var(ncid, varid, value) )
      endif
    endif

    call nc_close(ncid)
    !
  END SUBROUTINE 
  !
  SUBROUTINE nc_read_var_real_1d(file, var_name, value, start, count)
    !
    character(len=*), intent(in)      :: file
    character(len=*), intent(in)      :: var_name
    integer, optional, intent(in)     :: start(:)
    integer, optional, intent(in)     :: count(:)
    real, intent(out)                 :: value(:)
    integer                           :: ncid, varid
    !
    call nc_open(file, ncid)

    call nc_check( nf90_inq_varid(ncid, var_name, varid) )

    if (present(start)) then
      if (present(count)) then
        call nc_check( nf90_get_var(ncid, varid, value, start=start, count=count) )
      else
        call nc_check( nf90_get_var(ncid, varid, value, start=start) )
      endif
    else
      if (present(count)) then
        call nc_check( nf90_get_var(ncid, varid, value, count=count) )
      else
        call nc_check( nf90_get_var(ncid, varid, value) )
      endif
    endif

    call nc_close(ncid)
    !
  END SUBROUTINE 
  !
  SUBROUTINE nc_read_var_real_2d(file, var_name, value, start, count)
    !
    character(len=*), intent(in)      :: file
    character(len=*), intent(in)      :: var_name
    integer, optional, intent(in)     :: start(:)
    integer, optional, intent(in)     :: count(:)
    real, intent(out)                 :: value(:,:)
    integer                           :: ncid, varid
    !
    call nc_open(file, ncid)

    call nc_check( nf90_inq_varid(ncid, var_name, varid) )

    if (present(start)) then
      if (present(count)) then
        call nc_check( nf90_get_var(ncid, varid, value, start=start, count=count) )
      else
        call nc_check( nf90_get_var(ncid, varid, value, start=start) )
      endif
    else
      if (present(count)) then
        call nc_check( nf90_get_var(ncid, varid, value, count=count) )
      else
        call nc_check( nf90_get_var(ncid, varid, value) )
      endif
    endif

    call nc_close(ncid)
    !
  END SUBROUTINE 
  !
  SUBROUTINE nc_read_var_real_3d(file, var_name, value, start, count)
    !
    character(len=*), intent(in)      :: file
    character(len=*), intent(in)      :: var_name
    integer, optional, intent(in)     :: start(:)
    integer, optional, intent(in)     :: count(:)
    real, intent(out)                 :: value(:,:,:)
    integer                           :: ncid, varid
    !
    call nc_open(file, ncid)

    call nc_check( nf90_inq_varid(ncid, var_name, varid) )

    if (present(start)) then
      if (present(count)) then
        call nc_check( nf90_get_var(ncid, varid, value, start=start, count=count) )
      else
        call nc_check( nf90_get_var(ncid, varid, value, start=start) )
      endif
    else
      if (present(count)) then
        call nc_check( nf90_get_var(ncid, varid, value, count=count) )
      else
        call nc_check( nf90_get_var(ncid, varid, value) )
      endif
    endif

    call nc_close(ncid)
    !
  END SUBROUTINE 
  !
  SUBROUTINE nc_read_var_real_4d(file, var_name, value, start, count)
    !
    character(len=*), intent(in)      :: file
    character(len=*), intent(in)      :: var_name
    integer, optional, intent(in)     :: start(:)
    integer, optional, intent(in)     :: count(:)
    real, intent(out)                 :: value(:,:,:,:)
    integer                           :: ncid, varid
    !
    call nc_open(file, ncid)

    call nc_check( nf90_inq_varid(ncid, var_name, varid) )

    if (present(start)) then
      if (present(count)) then
        call nc_check( nf90_get_var(ncid, varid, value, start=start, count=count) )
      else
        call nc_check( nf90_get_var(ncid, varid, value, start=start) )
      endif
    else
      if (present(count)) then
        call nc_check( nf90_get_var(ncid, varid, value, count=count) )
      else
        call nc_check( nf90_get_var(ncid, varid, value) )
      endif
    endif

    call nc_close(ncid)
    !
  END SUBROUTINE 
  !
  !
  SUBROUTINE nc_check(status)
    !
    integer, intent(in)          :: status
    !
    if (status /= nf90_noerr) then
      print*, trim(nf90_strerror(status))
      stop "Stopped"
    endif
  END SUBROUTINE


END MODULE module_nc
