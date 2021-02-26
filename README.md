# NetCDF in Python

This repository contains tutorials and utilities for accessing and manipulating NetCDF files in Python.

## Introduction to NetCDF
NetCDF (or Network Common Data Form) files are commonly used to store multidimensional geographic data. This
multi-dimensionality of the files is often the cause of confusion. NetCDF files can include the following dimensions:

- Longitude (or x-axis)
- Latitude (or y-axis)
- Depth/Height (or z-axis)
- Time (or t-axis)

### X/Y dimensions
Since a NetCDF file is intended for geographic data, NetCDF files will likely always contain Longitude/Latitude (or
x/y) dimensions. While each cell can be laid out in a regular grid (like graph paper), NetCDF does not mandate this.
Other layout formats, such as curvilinear or skewed grids are also common.

### Depth/Height dimension
This is an optional dimension, depending on the data being represented. The Z dimension can be either "depth" (+ve =
below the surface, -ve = above the surface) or "height" ( +ve = above the surface, -ve = below the surface).

### Time dimension
This is an optional dimension, depending on the data being represented. Normally the time dimension contains either
absolute values (eg: 1-July-2020 15:37 UTC) or relative values (relative to a specific date/time, such as number of
milliseconds since start of 1-Jan-1990 UTC).

## Executing

This repository includes a number of Jupyter Notebooks-based tutorials designed to run on 
[Google Colab](https://colab.research.google.com/notebooks/intro.ipynb).

### Visualising NetCDF
Links: [View on Github](./notebooks/visualising-netcdf.ipynb) | [Open in Colab](https://colab.research.google.com/github/eatlas/netcdf-python/blob/master/notebooks/visualising-netcdf.ipynb) 

A notebook that introduces the reader to accessing a NetCDF file, retrieving metadata, and displaying data on a map.

### Regridding NetCDF
Links: [View on Github](./notebooks/regrid-netcdf.ipynb) | [Open in Colab](https://colab.research.google.com/github/eatlas/netcdf-python/blob/master/notebooks/regrid-netcdf.ipynb)

A notebook that introduces non-rectilinear grids, and uses the [Regridding module](./regrid) from this repository to
regrid a curvilinear grid to a rectilinear grid, displaying before and after data on maps.

## Modules

This repository includes one or more Python modules.

### Regridding
Links: [View on Github](./regrid)

A Python module to regrid a non-rectilinear grid to a rectilinear grid.

## Development environment
To perform further development on this repository, it is helpful to have a standardised development environment. To that
end, a [Docker](https://www.docker.com/) image has been created to provide a standardised Python 3 environment with the 
necessary NetCDF libraries installed.

__Note:__ this development environment is only intended for further development within this repository, NOT for running
the tutorial Jupyter notebooks. The Jupyter notebooks are executed via Google Colab.

The use this development environment, your computer must have Docker installed, which is beyond the scope of this 
article.

To start the Docker virtual environment on a Linux host machine:
```shell
./start-environment.sh 
```

## Datasets

While there are many NetCDF-based datasets available on the internet, the tutorials in this repository reference the
following datasets:

- [eReefs Hydrodynamic and BioGeoChemical models of the Great Barrier Reef (GBR)](https://ereefs.aims.gov.au/ereefs-aims)
  (downloads: [raw](https://dapds00.nci.org.au/thredds/catalogs/fx3/catalog.html)
  or
  [regridded and aggregated](http://thredds.ereefs.aims.gov.au/thredds/s3catalogue/aims-ereefs-public-prod/derived/ncaggregate/ereefs/catalog.html))
- [GBR photosynthetically active radiation (PAR) at 8 m depth](https://eatlas.org.au/data/uuid/eebd1438-2d4e-4f60-9055-27e6b9e58c3a)
  ([download](https://maps.eatlas.org.au/thredds/catalog/NESP-TWQ-5-3_Benthic-light/xr_par8/orig/xr_par8_daily/catalog.html))

### Acknowledgements
These tutorials build on ideas and content from the following sources:
- [Read NetCDF Data with Python](https://towardsdatascience.com/read-netcdf-data-with-python-901f7ff61648) 
  (10-July-2020) by [Konrad Hafen](https://khafen.medium.com/)
