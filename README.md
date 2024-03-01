# WRF2FVCOM
WRF2FVCOM is a collection of Fortran modules and a main program used for converting data from the Weather Research and Forecasting (WRF) model to the Finite Volume Coastal Ocean Model (FVCOM) format. This tool is useful for researchers and scientists working in the fields of atmospheric and oceanic modeling.

## Requirements
- Fortran compiler (e.g., gfortran)
- NetCDF library

## Compilation
To compile the WRF2FVCOM tool, use the provided Makefile. Ensure that the NetCDF library is installed on your system and update the paths in the Makefile accordingly.

```bash
make
```

## Usage
After compiling, you can run the WRF2FVCOM tool as follows.
- Basic usage:
```bash
./wrf2fvcom -i input.nc -o output.nc
```

- Include forcing variables for FVCOM-ICE module:
```bash
./wrf2fvcom -i input.nc -o output.nc -ice
```

- Calculate sea-level pressure for air pressure (rather than surface pressure):
```bash
./wrf2fvcom -i input.nc -o output.nc -slp
```

- Add Cartesian Coordinates (take the Gulf of Maine for an example):
```bash
./wrf2fvcom -i input.nc -o output.nc -proj '+init=nad83:1802'
```

- Use COAREv4.0 algorithum (rather than COAREv2.6):
```bash
./wrf2fvcom -i input.nc -o output.nc -v 4.0
```

- Extract one region and certain times, for example, extracting x-grid 3-12, y-grid 5-9, time 13 to the end :
```bash
./wrf2fvcom -i input.nc -o output.nc -x1 3 -nx 10 -y1 5 -ny 5 -t1 13
```

- Convert several wrf files using a list (list.dat):
```bash
./wrf2fvcom -i list.dat -o output.nc -l
```

- Convert several wrf files using a list (list.dat), where the listed files are successively archived from one WRF case:
```bash
./wrf2fvcom -i list.dat -o output.nc -l -s
```


## Files
- `module_cloudfrac.f90`: Module for cloud fraction calculations.
- `module_coare.f90`: Module for bulk flux calculations.
- `module_nc.f90`: Module for NetCDF file operations.
- `module_slp.f90`: Module for sea level pressure calculations.
- `wrf2fvcom_v3.1.f90`: Main program for converting WRF output to FVCOM format.

## Contributing
Contributions to WRF2FVCOM are welcome! If you find any issues or have suggestions for improvements, please submit a pull request or open an issue on GitHub.

## License
This project is licensed under the MIT License - see the LICENSE file for details.
