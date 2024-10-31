"""
Author: Farzad
Created on: 2024-10-22
"""

# This is to downlaod data from the web CHELSA data

# %% [CELL 0 ] Import libraries
import os
import urllib.request
import zipfile
import shutil
import wget
import cdsapi
import pyproj
from affine import Affine
import xarray as xr
import rasterio as rio
import matplotlib.pyplot as plt
from rasterio import CRS
from pyproj import Proj, transform, Transformer
from pyproj.transformer import TransformerGroup


# %% functions
def download_data(urls: list, output_folder: str, overwrite=False):
    # Create the output directory if it does not exist
    os.makedirs(output_folder, exist_ok=True)

    for url in urls:
        # Get the filename from the URL
        name_file = url.split("/")[-1]
        file_path = os.path.join(output_folder, name_file)

        # Check if the file exists and overwrite is not set
        if os.path.exists(file_path) and not overwrite:
            print(f"{file_path} already exists. Skipping download.")
            continue

        # Try to download the file
        try:
            wget.download(url, out=output_folder)
            print(f"\nDownloaded data from {url} to {file_path}")
        except NameError as e:
            print(f"Error: {url} not found - {e}")
        except Exception as e:
            print(f"An unexpected error occurred while downloading {url}: {e}")


def download_data_get_clipped_raster(
    urls: list,
    output_folder: str,
    reference_dem_path,
    dem_crs="EPSG:3035",
    overwrite=False,
):
    # Create the output directory if it does not exist
    os.makedirs(output_folder, exist_ok=True)
    tg = TransformerGroup(dem_crs, "EPSG:4326", always_xy=True)
    tg.transformers[0].definition
    bounds_3035 = rio.open(reference_dem_path).bounds
    left, bottom = tg.transformers[0].transform(bounds_3035.left, bounds_3035.bottom)
    right, top = tg.transformers[0].transform(bounds_3035.right, bounds_3035.top)
    # Transform the bounds to EPSG:4326
    print("The bounding box (in EPSG:4326) for the DEM is: ", left, bottom, right, top)

    transform_clipped = Affine(
        float(0.0083333333),  # a
        float(0),  # b
        float(left),  # c (x-coordinate of the upper-left corner)
        float(0),  # d
        float(-0.0083333333),  # e
        float(top),  # f (y-coordinate of the upper-left corner)
    )
    for url in urls:
        # Get the filename from the URL
        name_file = url.split("/")[-1]
        file_path = os.path.join(output_folder, name_file)

        # Check if the file exists and overwrite is not set
        if os.path.exists(file_path) and not overwrite:
            print(f"{file_path} already exists. Skipping download.")
            continue

        # Try to download the file
        try:
            wget.download(url, out=output_folder)
            print(f"\nDownloaded data from {url} to {file_path}")
            data_xr_ds = xr.open_dataset(file_path)
            var_name = list(data_xr_ds.data_vars.keys())[-1]

            clipped_data = data_xr_ds.sel(
                lon=slice(left, right), lat=slice(bottom, top)
            )

            clipped_data = clipped_data[var_name]
            clipped_data.rio.to_raster(
                file_path.replace(".nc", ".tif"),
                transform=transform_clipped,
                nodata=-999,
            )

            # let us delet nc file
            os.remove(file_path)
        except NameError as e:
            print(f"Error: {url} not found - {e}")
        except Exception as e:
            print(f"An unexpected error occurred while downloading {url}: {e}")


# %% [CELL 1] Download data CHELSA ----------------------------------

# download Meta data ------------------------------------------------
name_pdfs = [
    "RCA4_bioclims",
    "RACMO22E_bioclims",
    "HIRHAM5_bioclims",
    "CSC-REMO2009_bioclims",
]
pdfs_links = [
    f"https://os.zhdk.cloud.switch.ch/chelsav2/EUR11/{name_pdf}.pdf"
    for name_pdf in name_pdfs
] + [
    "https://os.zhdk.cloud.switch.ch/chelsav2/EUR11/documentation/CHELSA_EUR11_technical_documentation.pdf"
]
output_folder = "/share/home/farzad/datacube/Data/Meta_data"
download_data(pdfs_links, output_folder)

# let us download the annual data bio01 that mean observation (It has only one variable)----------------------------------
path_data = [
    f"https://os.zhdk.cloud.switch.ch/chelsav2/EUR11/obs/annual/V2.1/bio01d/CHELSA_EUR11_obs_bio01d_{year}_V.2.1.nc"
    for year in range(1991, 2021)
]
output_folder = "/share/home/farzad/datacube/Data/CHELSA/Historical/V2.1"
download_data(path_data, output_folder)


variables_v1 = ["cdd", "cwd", "npp", "pr_timmean", "pr_timstd"]
years = list(range(1981, 2006))
path_data = [
    f"https://os.zhdk.cloud.switch.ch/chelsav2/EUR11/obs/annual/V1.0/CHELSA_EUR11_{variable_v1}_{year}.nc"
    for variable_v1 in variables_v1
    for year in years
]
output_folder = "/share/home/farzad/datacube/Data/CHELSA/Historical/V1.0"
download_data(path_data, output_folder)

variables_v11 = [
    "cdd",
    "cfd",
    "csu",
    "cwd",
    "etr",
    "fd",
    "hd",
    "npp",
    "tas_timmean",
    "tas_timstd",
    "tasmax_timmax",
    "tasmax_timstd",
    "tasmin_timmin",
    "tasmin_timstd",
]
years = list(range(1981, 2005))
path_data = [
    f"https://os.zhdk.cloud.switch.ch/chelsav2/EUR11/obs/annual/V1.1/CHELSA_EUR11_{variable_v11}_{year}_V1.1.nc"
    for variable_v11 in variables_v11
    for year in years
]
output_folder = "/share/home/farzad/datacube/Data/CHELSA/Historical/V1.1"
download_data(path_data, output_folder, overwrite=False)

# %% [CELL 2] Download data ERA5 (HISTORICAL) ----------------------------------

variables = [
    "2m_dewpoint_temperature",
    "2m_temperature",
    #  "skin_temperature",
    #  "soil_temperature_level_1",
    #  "soil_temperature_level_2",
    #  "soil_temperature_level_3",
    #  "soil_temperature_level_4",
    #  "lake_bottom_temperature",
    #  "lake_ice_depth",
    #  "lake_ice_temperature",
    #  "lake_mix_layer_depth",
    #  "lake_mix_layer_temperature",
    #  "lake_shape_factor",
    #  "lake_total_layer_temperature",
    #  "snow_albedo",
    #  "snow_cover",
    #  "snow_density",
    #  "snow_depth",
    #  "snow_depth_water_equivalent",
    #  "snowfall",
    #  "snowmelt",
    #  "temperature_of_snow_layer",
    #  "skin_reservoir_content",
    #  "volumetric_soil_water_layer_1",
    #  "volumetric_soil_water_layer_3",
    #  "volumetric_soil_water_layer_4",
    #  "forecast_albedo",
    #  "surface_latent_heat_flux",
    #  "surface_net_solar_radiation",
    #  "surface_net_thermal_radiation",
    #  "surface_sensible_heat_flux",
    #  "surface_solar_radiation_downwards",
    #  "surface_thermal_radiation_downwards",
    #  "evaporation_from_bare_soil",
    #  "evaporation_from_open_water_surfaces_excluding_oceans",
    #  "evaporation_from_the_top_of_canopy",
    #  "evaporation_from_vegetation_transpiration",
    #  "potential_evaporation",
    #  "runoff",
    #  "snow_evaporation",
    #  "sub_surface_runoff",
    #  "surface_runoff",
    #  "total_evaporation",
    "10m_u_component_of_wind",
    "10m_v_component_of_wind",
    #  "surface_pressure",
    "total_precipitation",
    "leaf_area_index_high_vegetation",
    "leaf_area_index_low_vegetation",
]

years = list(range(1950, 2021))
years_str = [str(year) for year in years]

leap_years = [
    1952,
    1956,
    1960,
    1964,
    1968,
    1972,
    1976,
    1980,
    1984,
    1988,
    1992,
    1996,
    2000,
    2004,
    2008,
    2012,
    2016,
    2020,
    2024,
]

months = list(range(1, 13))
months_str = [str(month).zfill(2) for month in months]

days = list(range(1, 32))
days_str = [str(day).zfill(2) for day in days]

times = [
    "00:00",
    "01:00",
    "02:00",
    "03:00",
    "04:00",
    "05:00",
    "06:00",
    "07:00",
    "08:00",
    "09:00",
    "10:00",
    "11:00",
    "12:00",
    "13:00",
    "14:00",
    "15:00",
    "16:00",
    "17:00",
    "18:00",
    "19:00",
    "20:00",
    "21:00",
    "22:00",
    "23:00",
]

dataformat = "netcdf"  # 'grib'
downloadformat = "unarchived"  # 'zip'
area = [50, 12, 42, 50]  # Italy


# test
variables = ["2m_temperature"]
years_str = "2020"
months_str = "01"
days = "01"
times = "00:00"
dataformat = "netcdf"  # 'grib'
downloadformat = "unarchived"  # 'zip'
area = [50, 12, 42, 50]  # Italy

dataset = "reanalysis-era5-land"
request = {
    "variable": variables,
    "year": years_str,
    "month": months_str,
    "day": days,
    "time": times,
    "data_format": dataformat,
    "download_format": downloadformat,
    "area": area,
}

client = client = cdsapi.Client(
    url="https://cds.climate.copernicus.eu/api",
    key="594b0504-516a-491d-b89d-f5a14b80c313",  # Replace with your username and API key
)
output_folder = "/share/home/farzad/datacube/Data/ERA5/1_historical"
os.makedirs(output_folder, exist_ok=True)
output_name = output_folder + f"/{variable}_{year}_{month}_{day}_{time}.nc"
client.retrieve(dataset, request).download()


# %% [CELL 3] Test

path_data = [
    f"https://os.zhdk.cloud.switch.ch/chelsav2/EUR11/obs/annual/V2.1/bio01d/CHELSA_EUR11_obs_bio01d_{year}_V.2.1.nc"
    for year in range(1991, 2021)
]
output_folder = "/share/home/farzad/datacube/Data/CHELSA/Historical/V2.1"
reference_dem_path = "/share/home/farzad/World_bank/Europe/Italy/DEM/dem_3035.tif"
download_data_get_clipped_raster(
    path_data, output_folder, reference_dem_path, overwrite=True
)


# open

# %%
