import subprocess
import glob

file_list = sorted(glob.glob("pinacles_in_21-03-15-00z_to_21-03-17-00z_*.nc"))
print(file_list[25:])

header = "pinacles_in_21-03-16-01z_to_21-03-17-00z"
for i, f in enumerate(file_list[25:]):
    subprocess.call(["ncatted", "-a", "units,valid_time,o,c,hours since 2021-03-08T01:00:00", f, f"{header}_{i:04d}_in.nc"])
    subprocess.call(["ncap2", "-s", "valid_time+=168", f"{header}_{i:04d}_in.nc", f"{header}_{i:04d}.nc"])
