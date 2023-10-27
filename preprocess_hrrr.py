import pylab as plt
import pandas as pd
import cartopy.crs as ccrs
import numpy as np
import xarray as xr
import warnings
import pyproj

warnings.filterwarnings("ignore")

"""
# This is another way to do the rotation.
# It should be equivalent to Kyle's way.
def uv_to_earth_new(lon, lat, proj, u, v):
    m = proj.get_factors(lon, lat)

    dy_dphi = m.dy_dphi
    dy_dlam = m.dy_dlam

    # alpha = -np.arctan2(-np.cos(np.deg2rad(lat)) * (1.0 / dy_dlam), (1.0 / dy_dphi))
    alpha = np.arctan2(dy_dlam/np.cos(np.deg2rad(lat)), dy_dphi)

    return (
        u * np.cos(alpha) + v * np.sin(alpha),
        v * np.cos(alpha) - u * np.sin(alpha),
        alpha,
    )


def uv_to_grid_new(lon, lat, proj, u, v):
    m = proj.get_factors(lon, lat)

    dy_dphi = m.dy_dphi
    dy_dlam = m.dy_dlam

    alpha = - np.arctan2(dy_dlam/np.cos(np.deg2rad(lat)), dy_dphi)

    return (
        u * np.cos(alpha) + v * np.sin(alpha),
        v * np.cos(alpha) - u * np.sin(alpha),
        alpha,
    )
"""

"""
# This is the WPS way of rotating the winds from HRRR grid-relative to earth-relative.
# This specific subroutine is a translation of the Fortran code found at
# This is from the HRRR FAQ webpage at https://rapidrefresh.noaa.gov/faq/HRRR.faq.html.
def uv_to_earth_noaa(lon, u, v):
    rlon = -97.5
    rlat = 38.5
    rotcon_p = np.sin(np.deg2rad(rlat))

    alpha = rotcon_p * np.deg2rad(lon - rlon)

    return (
        np.cos(alpha) * u + np.sin(alpha) * v,
        -np.sin(alpha) * u + np.cos(alpha) * v,
        alpha,
    )
"""


def uv_to_earth_kyle(lon, lat, proj, u, v):
    m = proj.get_factors(lon, lat)

    dy_dlam = m.dy_dlam
    dx_dlam = m.dx_dlam

    theta = -np.arctan2(dy_dlam, dx_dlam)

    return (
        u * np.cos(theta) - v * np.sin(theta),
        u * np.sin(theta) + v * np.cos(theta),
    )


def uv_to_grid_kyle(lon, lat, proj, u, v):
    m = proj.get_factors(lon, lat)

    dy_dlam = m.dy_dlam
    dx_dlam = m.dx_dlam

    theta = np.arctan2(dy_dlam, dx_dlam)

    return (
        u * np.cos(theta) - v * np.sin(theta),
        u * np.sin(theta) + v * np.cos(theta),
    )


def get_nearest_point(projection, longitude, latitude):
    # based on the 1799 x 1059 grid size and dx = 3000 m
    xa = np.arange(-2.697e06, 2.697e06 + 3000.0, 3000.0)
    ya = np.arange(-1.587e06, 1.587e06 + 3000.0, 3000.0)

    pc = ccrs.PlateCarree()
    # x, y = projection.transform_point(longitude, latitude, ccrs.PlateCarree())

    proj_2 = pyproj.Transformer.from_crs(pc, projection)

    x, y = proj_2.transform(longitude, latitude)

    xi = np.argmin(np.abs(x - xa))
    yi = np.argmin(np.abs(y - ya))

    return xi, yi


# Time range to be processed
# DATES = pd.date_range("2021-03-08 00:00", "2021-03-16 00:00", freq="1H")
# of_name = (
#     DATES[0].strftime("%y-%m-%d-%Hz") + "_to_" + DATES[-1].strftime("%y-%m-%d-%Hz")
# )
DATES = pd.date_range("2021-03-15 00:00", "2021-03-17 00:00", freq="1H")
of_name = (
    DATES[0].strftime("%y-%m-%d-%Hz") + "_to_" + DATES[-1].strftime("%y-%m-%d-%Hz")
)

# input zarr files
# produced by convert_hrrr_nat_to_zarr.py and convert_hrrr_sfc_to_zarr.ipynb
sfc_zarr_file = "/home/xiao169/lidarbuoy/zarr/" + of_name + "_surface.zarr"
levels_zarr_file = "/home/xiao169/lidarbuoy/zarr/" + of_name + "_levels.zarr"

# output file location
out_file_header = "/home/xiao169/scratch/pinacles_in.new/pinacles_in_" + of_name

globe = ccrs.Globe(ellipse="sphere", semimajor_axis=6370000, semiminor_axis=6370000)
# map projections for HRRR
crsproj_hrrr = ccrs.LambertConformal(
    central_longitude=-97.5,
    central_latitude=38.5,
    standard_parallels=[38.5],
    globe=globe,
)
pyproj_hrrr = pyproj.Proj(crsproj_hrrr)

# map projections for PINACLES
lat_center = 35.00
lon_center = -123.71
proj_lon_center = -68.0
crsproj_pinacles = ccrs.LambertConformal(
    central_longitude=proj_lon_center,
    central_latitude=lat_center,
    standard_parallels=(lat_center, lat_center),
    globe=globe,
)
pyproj_pinacles = pyproj.Proj(crsproj_pinacles)

# Reading in surface variables
ds = xr.open_zarr(sfc_zarr_file)
# Convert longitude from 0_360 to -180_180 and rename vairables
ds = ds.assign_coords(longitude=(((ds.longitude + 180) % 360) - 180))
ds = ds.rename_vars({"t": "SST"})
ds = ds.rename_vars({"t2m": "T2m"})
ds = ds.rename_vars({"sh2": "QV2m"})
ds = ds.rename_vars({"u10": "U10m"})
ds = ds.rename_vars({"v10": "V10m"})
ds = ds.rename_vars({"sp": "PSFC"})

# Reading in sigma-level variables
ds_p = xr.open_zarr(levels_zarr_file)
ds_p = ds_p.assign_coords(longitude=(((ds_p.longitude + 180) % 360) - 180))
ds_p = ds_p.rename_vars({"pres": "P"})
ds_p = ds_p.rename_vars({"t": "T"})
ds_p = ds_p.rename_vars({"gh": "Z"})
ds_p = ds_p.rename_vars({"q": "QV"})
ds_p = ds_p.rename_vars({"clwmr": "QC"})
ds_p = ds_p.rename_vars({"unknown": "QI"})
ds_p = ds_p.rename_vars({"u": "U"})
ds_p = ds_p.rename_vars({"v": "V"})

# slice of hrrr domain to extract
# width is just a rough estimate
width = 100
buoy_lat = 35.71074
buoy_lon = -121.84606
xi_buoy, yi_buoy = get_nearest_point(crsproj_hrrr, buoy_lon, buoy_lat)
xi, yi = get_nearest_point(crsproj_hrrr, lon_center, lat_center)
xslice = slice(xi - width, xi + width)
yslice = slice(yi - width, yi + width)
print(
    f"PINACLES domain center, lon = {ds.longitude.values[yi, xi]}, lat = {ds.latitude.values[yi, xi]}"
)
print("Slices of HRRR domain to extract:")
print(xslice, yslice)

nt = ds.time.shape[0]
print(f"Total time steps to output: {nt}")

for ti in np.arange(nt):
    print(f"Processing time step {ti} ...")
    tslice = slice(ti, ti + 1)
    out_file = f"{out_file_header}_{ti:04d}.nc"

    # First prepare the SST
    sst = ds.SST.isel(time=tslice, x=xslice, y=yslice)
    sst.to_netcdf(out_file, mode="w", unlimited_dims="time")

    # land = ds.LAND.isel(x=xslice, y=yslice)
    # land.to_netcdf(out_file, 'a')

    t2m = ds.T2m.isel(time=tslice, x=xslice, y=yslice)
    t2m.to_netcdf(out_file, mode="a", unlimited_dims="time")

    qv2m = ds.QV2m.isel(time=tslice, x=xslice, y=yslice)
    qv2m.to_netcdf(out_file, mode="a", unlimited_dims="time")

    u10m = ds.U10m.isel(time=tslice, x=xslice, y=yslice)
    v10m = ds.V10m.isel(time=tslice, x=xslice, y=yslice)
    lon = v10m.longitude.values
    lat = u10m.latitude.values

    # Rotate wind into PINACLES Projection
    for tii in range(u10m.values.shape[0]):
        u_earth, v_earth = uv_to_earth_kyle(
            lon, lat, pyproj_hrrr, u10m.values[tii, :, :], v10m.values[tii, :, :]
        )
        u10m.data[tii, :, :], v10m.data[tii, :, :] = uv_to_grid_kyle(
            lon, lat, pyproj_pinacles, u_earth, v_earth
        )

    u10m.to_netcdf(out_file, mode="a", unlimited_dims="time")
    v10m.to_netcdf(out_file, mode="a", unlimited_dims="time")

    pres = ds.PSFC.isel(time=tslice, x=xslice, y=yslice)
    pres.to_netcdf(out_file, mode="a", unlimited_dims="time")

    p = ds_p.P.isel(time=tslice, x=xslice, y=yslice)
    p.to_netcdf(out_file, mode="a", unlimited_dims="time")

    t = ds_p.T.isel(time=tslice, x=xslice, y=yslice)
    t.to_netcdf(out_file, mode="a", unlimited_dims="time")

    z = ds_p.Z.isel(time=tslice, x=xslice, y=yslice)
    z.to_netcdf(out_file, mode="a", unlimited_dims="time")

    qv = ds_p.QV.isel(time=tslice, x=xslice, y=yslice)
    qv.to_netcdf(out_file, mode="a", unlimited_dims="time")

    qc = ds_p.QC.isel(time=tslice, x=xslice, y=yslice)
    qc.to_netcdf(out_file, mode="a", unlimited_dims="time")

    qi = ds_p.QI.isel(time=tslice, x=xslice, y=yslice)
    qi.to_netcdf(out_file, mode="a", unlimited_dims="time")

    u = ds_p.U.isel(time=tslice, x=xslice, y=yslice)
    v = ds_p.V.isel(time=tslice, x=xslice, y=yslice)

    lon = u.longitude.values
    lat = u.latitude.values

    # Rotate winds to PINACLES Projeciton
    for tii in range(u.values.shape[0]):
        for ki in range(u.values.shape[1]):
            u_earth, v_earth = uv_to_earth_kyle(
                lon, lat, pyproj_hrrr, u.values[tii, ki, :, :], v.values[tii, ki, :, :]
            )
            u.data[tii, ki, :, :], v.data[tii, ki, :, :] = uv_to_grid_kyle(
                lon, lat, pyproj_pinacles, u_earth, v_earth
            )

    u.to_netcdf(out_file, mode="a", unlimited_dims="time")
    v.to_netcdf(out_file, mode="a", unlimited_dims="time")
