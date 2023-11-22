import numpy as np
import xarray as xr
import pandas as pd
import cartopy.crs as ccrs
import pyproj
import warnings
import os

os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"


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
    # Bill's WRF domain
    lon_c, lat_c = -129.0, 37.0
    nx, ny = 699, 699
    half_n = nx // 2
    dx = 3000.0  # m
    xa = np.arange(-half_n * dx, half_n * dx + dx, dx)
    ya = np.arange(-half_n * dx, half_n * dx + dx, dx)

    pc = ccrs.PlateCarree()

    proj_2 = pyproj.Transformer.from_crs(pc, projection)

    x, y = proj_2.transform(longitude, latitude)
    x_c, y_c = proj_2.transform(lon_c, lat_c)

    xi = np.argmin(np.abs(x - x_c - xa))
    yi = np.argmin(np.abs(y - y_c - ya))

    return xi, yi


globe = ccrs.Globe(ellipse="sphere", semimajor_axis=6370000, semiminor_axis=6370000)

# map projections for WRF
""" Copied from WRF output file
                :CEN_LAT = 37.00002f ;
                :CEN_LON = -129.f ;
                :TRUELAT1 = 30.f ;
                :TRUELAT2 = 42.f ;
                :MOAD_CEN_LAT = 37.00002f ;
                :STAND_LON = -100.f ;
"""
crsproj_wrf = ccrs.LambertConformal(
    central_longitude=-100.0,
    central_latitude=37.0,
    standard_parallels=[30.0, 42.0],
    globe=globe,
)
pyproj_wrf = pyproj.Proj(crsproj_wrf)

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

start_time = pd.to_datetime("2021-03-09 12:00")
start_tstring = start_time.strftime("%Y-%m-%d_%H:%M:%S")
# DATES = pd.date_range("2021-03-10 00:00", "2021-03-11 00:00", freq="30T")
DATES = pd.date_range("2021-03-09 12:00", "2021-03-11 00:00", freq="30T")

xr.backends.file_manager.FILE_CACHE.clear()

for DATE in DATES:
    datestring = DATE.strftime("%Y-%m-%d_%H_%M_%S")
    in_file = f"/lustre/eaglefs/projects/lidarbuoy/d3m088/wrf/era5/run_wrf_20210309/output_morr/wrf2pin_d01_{datestring}"
    out_file = (
        f"/lustre/eaglefs/scratch/xiao169/pinacles_in.wrf/pinacles_in_{datestring}.nc"
    )
    print(f"Input: {in_file}")
    print(f"Output: {out_file}")
    ds = xr.open_dataset(in_file, decode_times=False)
    new_unit = f"minutes since {start_tstring}"
    ds.XTIME.attrs["description"] = new_unit
    ds.XTIME.attrs["units"] = new_unit
    ds = xr.decode_cf(ds)
    print(f"Time: {ds.XTIME.values[0]}")
    if DATE == DATES[0]:
        width = 100
        xi, yi = get_nearest_point(crsproj_wrf, lon_center, lat_center)
        xslice = slice(xi - width, xi + width)
        yslice = slice(yi - width, yi + width)
        print(
            f"PINACLES domain center on WRF grid, lon = {ds.XLONG.values[0, yi, xi]}, lat = {ds.XLAT.values[0, yi, xi]}"
        )
        print("Slices of HRRR domain to extract:")
        print(xslice, yslice)
    ds = ds.rename_dims(
        {"west_east": "x", "south_north": "y", "Time": "t", "bottom_top": "z"}
    )
    ds = ds.rename_vars(
        {
            "XTIME": "time",
            "XLAT": "latitude",
            "XLONG": "longitude",
            "Q2": "QV2m",
            "T2": "T2m",
            "QVAPOR": "QV",
            "QCLOUD": "QC",
            "QICE": "QI",
            "U10": "U10m",
            "V10": "V10m",
        }
    )
    ds["latitude"] = xr.DataArray(
        data=ds.latitude[0, :, :],
        dims=["y", "x"],
        attrs={"units": "degree_north", "description": "LATITUDE, SOUTH IS NEGATIVE"},
    )
    ds["longitude"] = xr.DataArray(
        data=ds.longitude[0, :, :],
        dims=["y", "x"],
        attrs={"units": "degree_east", "description": "LONGITUDE, WEST IS NEGATIVE"},
    )
    dout = xr.Dataset()
    # 2D: SST, PSFC, QV2m, T2m, U10m, V10m
    dout["SST"] = ds.SST.isel(x=xslice, y=yslice)
    dout["PSFC"] = ds.PSFC.isel(x=xslice, y=yslice)
    dout["QV2m"] = ds.QV2m.isel(x=xslice, y=yslice)
    dout["T2m"] = ds.T2m.isel(x=xslice, y=yslice)
    u10m, v10m = ds.U10m.isel(x=xslice, y=yslice), ds.V10m.isel(x=xslice, y=yslice)
    lon = u10m.longitude.values
    lat = u10m.latitude.values
    u_earth, v_earth = uv_to_earth_kyle(lon, lat, pyproj_wrf, u10m.values, v10m.values)
    u10m.values, v10m.values = uv_to_grid_kyle(
        lon, lat, pyproj_pinacles, u_earth, v_earth
    )
    dout["U10m"] = u10m
    dout["V10m"] = v10m
    # 3D: QV, QC, QI, U, V, P, T, Z
    dout["QV"] = ds.QV.isel(x=xslice, y=yslice)
    dout["QC"] = ds.QC.isel(x=xslice, y=yslice)
    dout["QI"] = ds.QI.isel(x=xslice, y=yslice)
    ui, vi = ds.U.values, ds.V.values
    u = (ui[:, :, :, 1:] + ui[:, :, :, :-1]) * 0.5
    v = (vi[:, :, 1:, :] + vi[:, :, :-1, :]) * 0.5
    u_earth, v_earth = uv_to_earth_kyle(
        lon, lat, pyproj_wrf, u[:, :, yslice, xslice], v[:, :, yslice, xslice]
    )
    u[:, :, yslice, xslice], v[:, :, yslice, xslice] = uv_to_grid_kyle(
        lon, lat, pyproj_pinacles, u_earth, v_earth
    )
    u_dressed = xr.DataArray(data=u, dims=["t", "z", "y", "x"], attrs={"units": "m/s"})
    v_dressed = xr.DataArray(data=v, dims=["t", "z", "y", "x"], attrs={"units": "m/s"})
    dout["U"] = u_dressed.isel(x=xslice, y=yslice)
    dout["V"] = v_dressed.isel(x=xslice, y=yslice)
    p_dressed = xr.DataArray(
        data=ds.P.values + ds.PB.values,
        dims=["t", "z", "y", "x"],
        attrs={"units": "Pa"},
    )
    dout["P"] = p_dressed.isel(x=xslice, y=yslice)
    g = 9.81
    z = (ds.PH.values[:, :-1, :, :] + ds.PHB.values[:, 1:, :, :]) * 0.5 / g
    z_dressed = xr.DataArray(data=z, dims=["t", "z", "y", "x"], attrs={"units": "m"})
    dout["Z"] = z_dressed.isel(x=xslice, y=yslice)
    rd = 287.0 # From module_model_constants.F in WRF code
    rv = 461.6
    t0 = 300.0  # K
    # cp = 7.*rd/2. 
    t = (ds.THM.values + t0) * ((p_dressed.values/1.0e5)**(2./7.)) / (1.0 + rv * ds.QV.values / rd)
    t_dressed = xr.DataArray(data=t, dims=["t", "z", "y", "x"], attrs={"units": "K"})
    dout["T"] = t_dressed.isel(x=xslice, y=yslice)
    try:
        dout.to_netcdf(out_file, unlimited_dims=["t"])
    except PermissionError:
        print("AGAIN?")
    except:
        print("What Happened?")
    else:
        print("Good!")
    finally:
        dout.close()
        ds.close()
        del ds, dout
        xr.backends.file_manager.FILE_CACHE.clear()
