# Compiler and flags
FC = ifort
FFLAGS = -O3
LDFLAGS = -lnetcdff -lnetcdf

# NetCDF include and library paths
#NETCDF = /usr/medm/install/APPS/INTEL/22/with_mpich-4.0.2/netcdf-4.9.0
NETCDF = /usr/local/share/x86_64/netcdf/4.9.2/intel/2023.2.0

# List of source files
SRCS = module_cloudfrac.f90 module_coare.f90 module_nc.f90 module_slp.f90 wrf2fvcom.f90

# List of object files
OBJS = $(SRCS:.f90=.o)
MODS = $(SRCS:.f90=.mod)

# Executable name
EXEC = wrf2fvcom

all: $(EXEC)

$(EXEC): $(OBJS)
	$(FC) $(FFLAGS) -L$(NETCDF)/lib $(OBJS) -o $(EXEC) $(LDFLAGS)

%.o: %.f90
	$(FC) $(FFLAGS) -I$(NETCDF)/include -c $< -o $@

clean:
	rm -f $(OBJS) $(MODS) $(EXEC)



