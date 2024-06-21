#Process to load different files
import h5py

def load_3D_temperature_data(file):
    data_dict = {}
    with h5py.File(file, "r") as g:    
        data_dict["temperatures"]=g["temperature"][()] # this is (healpix x distance_bin)
        data_dict["distance_slices"] = g["distance_slices"][()]
        data_dict["dEBV"] = g["dEBV_smooth"][()]
        data_dict["rhos"] = g["rhos"][()] 
        data_dict["betas"]= g["betas"][()]
        data_dict["nside"] = g.attrs["nside"]
        data_dict["nr_of_distance_bins"] = g.attrs["nr_of_distance_bins"]
        data_dict["healpix_ordering"] = g.attrs["healpix_ordering"] #nest
        data_dict["run_name"] = g.attrs["original_run_name"]
        g.close()
    return data_dict

def load_just_temperature_data(file):
    '''
    A function that loads the 3D dust temperature map and returns a dictionary of the data. This dictionary contains the temperature, distance slices, the number of distance bins, the HEALPix resolution Nside, and the healpix ordering (nesting). 
    Parameters: 
        file------------------------------hdf5 dust temperature file
    Returns:
        data_dict-------------------------dictionary

    '''
    data_dict = {}
    with h5py.File(file, "r") as g:    
        data_dict["temperatures"]=g["temperature"][()] # this is (healpix x distance_bin)
        data_dict["distance_slices"] = g["distance_slices"][()]
        data_dict["nside"] = g.attrs["nside"]
        data_dict["nr_of_distance_bins"] = g.attrs["nr_of_distance_bins"]
        data_dict["healpix_ordering"] = g.attrs["healpix_ordering"] #nest
        g.close()
    return data_dict