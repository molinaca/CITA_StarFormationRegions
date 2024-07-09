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

def load_sftracer(tracer_data):
    '''
    A function that loads the longitude, latitude and distance of a star formation tracer. It converts them from degrees to radians, and converts the latitude into colatitude. If no distance is found, the function only returns theta and phi. 

    Parameters: 
        tracer_data-------------------------------pandas dataframe, must have a component 'l' and a component 'b', optionally can have 'D'

    Returns:
        long--------------------------------------array-like, longitude in degrees
        lat---------------------------------------array-like, latitude in degrees
        distance----------------------------------array-like, distance in kpc
    '''
    long = tracer_data['l']
    lat = tracer_data['b']

    if 'D' in tracer_data: #If the file has distance create a distance variable and return it
        distance = tracer_data['D'] 
        
        return long, lat, distance
    else: #If there is no distance only return that and phi
        print("No distance measurement")
        return long, lat