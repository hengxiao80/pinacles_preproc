# %%
import pandas as pd
from herbie import FastHerbie
import numpy as np
import xarray as xr
import warnings
warnings.filterwarnings('ignore')

# %%
#Specify the time range over which to generate zarr files
DATES = pd.date_range('2021-03-16 01:00', '2021-03-17 00:00', freq='1H')

# %%
# Create zarr filename prefix
of_name = DATES[0].strftime("%y-%m-%d-%Hz") + '_to_' + DATES[-1].strftime("%y-%m-%d-%Hz")

# %%
# Define native model level fields to be downloaded
level_name = {} 
level_name["U"] = ":(?:U)GRD:[0-9]+ hybrid"
level_name["V"] = ":(?:V)GRD:[0-9]+ hybrid"
level_name["P"] = ":(?:PRES):[0-9]+ hybrid"
level_name["T"] = ":(?:TMP):[0-9]+ hybrid"
level_name["QV"] = ":(?:SPFH):[0-9]+ hybrid"
level_name["QC"] = ":(?:CLMR):[0-9]+ hybrid"
level_name["QI"] = ":(?:CIMIXR):[0-9]+ hybrid"
level_name["Z"] =":(?:HGT:)[0-9]+ hybrid"

# %%
# for DATE in tqdm.tqdm(DATES):
    # with io.capture_output() as captured:
new_data = True
for i in range(len(DATES)):
    tslice = slice(i,i+1)
    Hnat = FastHerbie(DATES[tslice], model="hrrr", product="nat", save_dir='/home/xiao169/scratch/', fxx=[1], verbose=True)
    searchstring = f"({'|'.join(list(level_name.values()))})"
    m = Hnat.xarray(searchstring, max_threads=4)
    m.drop_vars('gribfile_projection')
    del m.attrs['local_grib'] 
    del m.attrs['remote_grib'] 
    print(m)
    if new_data: 
        m.to_zarr('/home/xiao169/scratch/' + of_name + '_levels.zarr', mode="a")
        new_data = False
    else:
        m.to_zarr('/home/xiao169/scratch/' + of_name + '_levels.zarr',append_dim="time")
    m.close() 


