"""
Author: Farzad
Created on: 2024-10-24
"""

# %% [CELL 0] import libraries
import numpy as np

# import netCDF4 as nc
import xarray as xr

# import matplotlib.pyplot as plt
import rasterio as rio
from affine import Affine

# %% [CELL 0] Functions


def save_nc_as_raster(data_array, output_file, crs, transform, nodata=-9999):
    """
    Save a xarray DataArray as a raster file
    """
    profile = {
        "driver": "GTiff",
        "height": data_array.shape[0],  # Number of rows
        "width": data_array.shape[1],  # Number of columns
        "count": 1,  # Number of bands
        "dtype": data_array.dtype,  # Data type (e.g., 'float32')
        "crs": crs,  # Coordinate reference system
        "transform": transform,  # Affine transformation
        "nodata": nodata,  # Nodata value (if applicable)
    }

    with rio.open(output_file, "w", **profile) as dst:
        dst.write(data_array.values[0].astype(profile["dtype"]), 1)
    print(f"Saved raster to {output_file}")


# %% [CELL 1] Load the data


# Load the data
hist_data_folder = "/share/home/farzad/datacube/Data/CHELSA/Historical/V1.1"
data_xr_ds = xr.open_dataset(hist_data_folder + "/CHELSA_EUR11_hd_1987_V1.1.nc")

# Print the dataset summary to see all details including variable names
print(data_xr_ds.attrs)

# Check dataset details to understand attributes and coordinates
print(data_xr_ds.attrs)  # Attributes (may contain CRS or transform info)
print(data_xr_ds.encoding)  # Encoding (might contain transform information)

# let us see variable name
var_name = list(data_xr_ds.data_vars.keys())[1]

# let us see the data
data_xr_ds
# let us see the variables
data_xr_ds.variables

# let us see the coordinates
data_xr_ds.coords
# let us see the dimensions
data_xr_ds.dims
# let us see the data
data_xr_ds.data_vars
# let us see the crs
data_xr_ds.crs
# let us see time values
data_xr_ds.time.values
transform = data_xr_ds.crs.GeoTransform



# give 4 corner of intrested area
left = 19  # lon_min
right = 29  # lon_max
bottom = 34  # lat_min
top = 42  # lat_max

# clip the data on the region of interest
clipped_data = data_xr_ds.sel(lon=slice(left, right), lat=slice(bottom, top))

# ley us plot the data
clipped_data["heating_degree_days_per_time_period"].isel(time=0).plot(cmap="viridis")
transform_clipped = Affine(
    float(0.0083333333),  # a
    float(0),  # b
    float(left),  # c (x-coordinate of the upper-left corner)
    float(0),  # d
    float(-0.0083333333),  # e
    float(top),  # f (y-coordinate of the upper-left corner)
)
# visualize the data
hdd_xr_data_arra_clipped = clipped_data["heating_degree_days_per_time_period"]
hdd_xr_data_arra_clipped.isel(time=0).plot(cmap="viridis")


# invert to raster and save the data
path_to_save = (
    "/share/home/farzad/datacube/Data/CHELSA/Historical/V1.1/hdd_italy_1987.tif"
)
hdd_xr_data_arra_clipped.rio.to_raster(
    path_to_save, transform=transform_clipped, nodata=-9999
)

# %% [CELL 2] Load the data
arr = rio.open("/share/home/farzad/climax/FIRE/provaclip.tiff").profile


profile = {
    "driver": "GTiff",
    "height": hdd_xr_data_arra_clipped.shape[1],  # Number of rows
    "width": hdd_xr_data_arra_clipped.shape[2],  # Number of columns
    "count": 1,  # Number of bands
    "dtype": "float32",  # Data type (e.g., 'float32')
    "crs": "EPSG:4326",  # Coordinate reference system
    "transform": transform_italy,  # Affine transformation
    "nodata": -999,  # Nodata value (if applicable)
}

{
    "driver": "GTiff",
    "dtype": "float32",
    "nodata": -9999.0,
    "width": 2810,
    "height": 2554,
    "count": 1,
    "crs": CRS.from_wkt(
        'LOCAL_CS["ETRS89-extended / LAEA Europe",UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'
    ),
    "transform": Affine(100.0, 0.0, 3488731.355109199, 0.0, -100.0, 2241986.6503754416),
    "blockxsize": 2810,
    "blockysize": 1,
    "tiled": False,
    "interleave": "band",
}


with rio.open(
    "/share/home/farzad/datacube/Data/CHELSA/Historical/V1.1/hdd_italy_1987.tif",
    "w",
    **profile,
) as dst:
    dst.write(hdd_array, 1)
