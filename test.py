import h5pyd

f = h5pyd.File("/nrel/nsrdb/v3/nsrdb_2022.h5", 'r')
print(list(f.keys()))
