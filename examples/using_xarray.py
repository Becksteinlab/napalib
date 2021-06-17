# xarray is used to store most data in napalib due to its ability
# to handle dimensions in a clean way. I highly recommend learning
# this library for use outside of napalib. I will give a brief 
# overview of how to use the data generated.

# let's say you ran the collection method for RMSD on a trajectory
# 

import xarray as xr

da = xr.load_dataarray("anton_358_OF_0.nc")

print(da.shape)
# (3, 3, 2, 4, 62500)

print(da)
#Coordinates:
#* time       (time) float64 0.0 240.0 480.0 720.0 ... 1.5e+07 1.5e+07 1.5e+07
#* crystal    (crystal) object 'OF' 'IF' 'OCC' 'EIF'
#* alignment  (alignment) object 'c' 'd' 'f'
#* selection  (selection) object 'c' 'd' 'f'
#* protomer   (protomer) object 'A' 'B'

# note that crystal should be thought of as reference and is only called
# crystal because we were originally only interested relative to the
# crystal strctures. Changing this would break existing analysis and plotting
# code. 
# OF -> outward facing crystal
# IF -> inward facing crystal
# OCC -> cluster center for occluded 
# EIF -> cluster center for the observed inward facing structure

# here we can see that the dataarray has a shape matching that of the
# coordinates. This allows for fast selection of very large datasets
# without needing to be concerned with management of indices. To get
# the RMSD trace of the core domain when the whole protein was used as
# the alignment relative to the OF crystal structure

core_on_full_A_rOF = da.sel(protomer="A", alignment="f", selection="c", crystal="OF")

# take the values out as numpy arrays if you wish
times = core_on_full_A_rOF.time.values
timeseries = core_on_full_A_rOF.values
