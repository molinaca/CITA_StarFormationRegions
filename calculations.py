#This script contains any functions needed to perform calculations with arrays in the main script
import numpy as np
import ephem
import healpy as hp
from constants import h, c, k
import matplotlib.pyplot as plt



## 1: Scientific Functions

#These are all functions that are used to portray scientific calculations. 

#Wien's Law 
def WiensLaw(T):
    '''
    Function that uses Wien's law to calculate the maximum frequency of a black body radiation curve at a given temperature

    Parameters:
    T: temperature in Kelvin
    h: planck constant from constants.py
    k: boltzmann constant from constants.py

    Output:
    nu_max: maximum frequency in Hz
    '''
    nu_max = 2.824*k*T/h
    return nu_max

## 2: Manipulating Arrays 

def get_temptracers_at_freq(Tmap, nu=None, normalize=True, method='planck', limits=None):
    '''
    Function to find the temperatures at a frequency nu that will act as tracers for the RGB channels.
    Includes different methods to do this. 

    Parameters:
    nu (optional): frequency in Hz used in planck method
    Tmap: temperature map in Kelvin
    normalize (optional): boolean to normalize the Planck function
    method (optional): method to use to calculate the temperature tracers. Options are 'planck', 'custom'
    limts (optional): limits for the custom method, of form [R_limit, G_limit]

    h: planck constant from constants.py
    c: speed of light from constants.py
    k: boltzmann constant from constants.py
    x_d: variable used to simplify planck function

    Output:
    rad : planck function at temperature in Tmap in units W/m^2/Hz/sr
    '''
    if method == 'planck':
        x_d = h*nu/(k*Tmap)
        rad = 2*h*nu**3/c**2/(np.exp(x_d)-1)
        if normalize:
            print("Normalizing by max")
            freq_max = WiensLaw(Tmap)
            x_d_max = h*freq_max/(k*Tmap)
            rad_max = 2*h*freq_max**3/c**2/(np.exp(x_d_max)-1)
            rad = rad/rad_max
        return rad 
    
    if method == 'custom':
        R_limit, G_limit = limits

        channel_1 = np.where(Tmap <R_limit, Tmap, 0.01)
        channel_2 = np.where((Tmap >= R_limit) & (Tmap < G_limit), Tmap, 0.01)
        channel_3 = np.where(Tmap >= G_limit, Tmap, 0.01)

        channel_array = np.array([channel_1, channel_2, channel_3])
        return channel_array

def increase_temp_res(data_dict, nside_new):
    '''
    Function used to increase the resolution of the temperature map

    Parameters:
    data_dict: dictionary containing the data
    nside_new: new nside to be used

    Output:
    Ts_new: new temperature map with the new nside
    '''
    Ts = data_dict["temperatures"]
    model_nslices = data_dict["nr_of_distance_bins"]

    Ts_new = np.zeros((model_nslices, hp.nside2npix(nside_new))) #Create new temperature with the new nside, note that have to switch axis because of the way ud_grade works
    
    #Increase resolution at each distance slice
    for ds_index in range(model_nslices):
        Ts_new[ds_index] = hp.ud_grade(Ts[:,ds_index], nside_out=nside_new, order_in='NESTED', order_out='NESTED') 

    return Ts_new

def custom_window_func(Tmap, R_limit, G_limit):
    '''
    Function used to just create custom windows for RGB in the temperature map. Creates a very stark contrast between the three channels.

    Parameters: 
    Tmap: temperature map, in form (distance_bin x pixel)
    R_limit: limit for the red channel, in Kelvin
    G_limit: limit for the green channel, in Kelvin

    Output:
    channel_array: array of shape (3 x distance_bin x pixel) containing the three channels
    '''

    channel_1 = np.where(Tmap <R_limit, Tmap, 0.01)
    channel_2 = np.where((Tmap >= R_limit) & (Tmap < G_limit), Tmap, 0.01)
    channel_3 = np.where(Tmap >= G_limit, Tmap, 0.01)

    channel_array = np.array([channel_1, channel_2, channel_3])
    return channel_array

def multiply_dEBVandTtracer(data_dict, dEBVmap, tracermap, frequency):
    '''
    Function that multiplied the dust density (dEBV) by the temperature tracer map at each distance slice and frequency. 

    Parameters:
    data_dict: dictionary containing the data, used to obtain number of distance slices
    dEBVmap: dust density map with shape (distance_bin x pixel)
    tracermap: temperature tracer map with shape (frequency x distance_bin x pixel)
    frequency: array containing the frequencies

    Output:
    dens_temp: map of shape (frequency x distance_bin x pixel) containing the multiplication of dEBV and B 
    '''
    model_nslices = data_dict["nr_of_distance_bins"]
    nfreq = len(frequency) #number of frequencies

    dens_temp = np.zeros((nfreq, dEBVmap.shape[0], dEBVmap.shape[1])) #shape (freq x distance x pixel)

    #Calculate multiplication at each distance slice and frequency
    for f_index in range(nfreq):
    
        for ds_index in range(model_nslices):
            dens_temp[f_index, ds_index] = dEBVmap[ds_index]*tracermap[f_index, ds_index]
    
    return dens_temp

def normalize_multiplied_array(data_dict, dens_temp, frequency):
    '''
    Function that normalizes the newly multiplied map of dust density and temperature emission across the whole sky.  

    Parameters:
    data_dict: dictionary containing the data, used to obtain number of distance slices
    dens_temp: map of shape (frequency x distance_bin x pixel) containing the multiplication of dEBV and B
    frequency: array containing the frequencies

    Output:
    dens_temp_norm: normalized map of shape (frequency x distance_bin x pixel) containing the multiplication of dEBV and B
    '''
    model_nslices = data_dict["nr_of_distance_bins"]
    nfreq = len(frequency) #number of frequencies

    dens_temp_norm = np.zeros((nfreq, dens_temp.shape[1], dens_temp.shape[2])) #shape (freq x distance x pixel)

    max_distance = np.zeros(model_nslices) #array containing max emission at each distance slice
    max_freq = np.zeros(nfreq) #array contaning max emission over all distance slices for each frequency

    for f_index in range(nfreq):
        for ds_index in range(model_nslices):
            max_distance[ds_index] = np.max(dens_temp[f_index, ds_index]) #calculate max emission at each distance

        max_freq[f_index] = np.max(max_distance) #choose overall max emission for each frequency

    max_overall = np.max(max_freq) #choose overall max emission over all frequencies

    #Normalize
    for f_index in range(nfreq):
        for ds_index in range(model_nslices):
            dens_temp_norm[f_index, ds_index] = dens_temp[f_index, ds_index]/max_overall #divide by max so is now between 0 and 1

    return dens_temp_norm

### 2.1: Functions for sf tracers
def remove_tracers_in_mask(long, lat, distance_bins):
    '''
    Function that removes any tracers that are within the declination mask (any dec < -30 degrees). It will output new, filtered arrays. 

    Parameters:
    long: list within list of longitudes for each distance slice
    lat: list within list of latitudes for each distance slice
    distance_bins: number of distance bins

    Output:
    filtered_long: list within list of longitudes for each distance slice after filtering
    filtered_lat: list within list of latitudes for each distance slice after filtering
    '''
    
    #Initialize new, filtered arrays in the same form as long, lat
    filtered_long = [[] for i in range(distance_bins)]
    filtered_lat = [[] for i in range(distance_bins)]

    for ds_index in range(distance_bins):
        long_atdist = long[ds_index] #To make it easier, get list of long and lat at each distance slice
        lat_atdist = lat[ds_index]

        ntracers = len(long_atdist) #Only want number at each distance slice

        #Initialize ra anddec
        ra = np.zeros(ntracers)
        dec = np.zeros(ntracers)

        for index in range(ntracers):
            l = long_atdist[index] #get each long and lat
            b = lat_atdist[index]
            
            galactic = ephem.Galactic(l/180.0*np.pi,b/180.0*np.pi)
            equatorial = ephem.Equatorial(galactic, epoch=ephem.J2000)
            ra[index]=equatorial.ra/np.pi*180.0
            dec[index]=equatorial.dec/np.pi*180.0

            #If not within mask, add to new arrays
            if dec[index] > -30:
                filtered_long[ds_index].append(l)
                filtered_lat[ds_index].append(b)
        
    return filtered_long, filtered_lat

### 2.2: Density variation functions

def flag_regions(nside, dist_slices, map, primary_threshold, secondary_threshold, radius, lowtohigh=False):
    '''
    Function that based on two thresholds, identifies which areas of the dust map have high dEBV (dust density). It saves these areas on their
    own, which includes details about the region maximum, as well as a map with all the regions. 

    Parameters:
    nside: int, nside of dEBV map
    dist_slices: int, number of distance slices
    map: np.array, dust map of form (distance_slices, pixel)
    primary_threshold: float, threshold for high density regions
    secondary_threshold: float, threshold for parts of the map that could still be in high density regions
    radius: float, radius of region around pixel to look at, in degrees

    Output:
    region_info: list of lists of dictionaries, contains the maximum, pixel and pixel values for each region at each distance slice
    high_density_map: np.array, map of all high density regions in form (distance_slices, pixel)
    '''

    #Make region_info list and high_density map    

    region_info = [[] for i in range(dist_slices)]

    high_density_map = np.zeros_like(map)

    if lowtohigh == True:
        reverse = False
    else:
        reverse = True

    for ds_index in range(dist_slices):
        # Create a binary mask for high-density regions (True if high dEBV, false otherwise)

        current_map = map[ds_index] #Define to make it easier

        if lowtohigh == True:
            binary_mask = current_map < primary_threshold
        else:
            binary_mask = current_map > primary_threshold #Pick highest density regions

        # Only look at pixels where binary mask is True
        pixels_to_check = np.where(binary_mask)[0]
        checked_pixels = set() #Don't want to repeat same pixels

        pixels_to_check = sorted(pixels_to_check, key=lambda p: current_map[p], reverse=reverse)

        for pixel in pixels_to_check:
            if pixel not in checked_pixels: # If pixel has not been checked yet

                #Make vector for pixel we're looking at
                pixel_vec = hp.pix2vec(nside, pixel, nest=True)
                
                region = hp.query_disc(nside, pixel_vec, np.radians(radius), nest=True, inclusive=True) #pick region around pixel
                
                if lowtohigh == True:
                    region = [p for p in region if current_map[p] < secondary_threshold]
                else:
                    region = [p for p in region if current_map[p] > secondary_threshold] #make sure to only keep pixels within secondary threshold

                #Makes sure to only do this if region is not empty
                if region:
                    
                    # Mark the region in the map
                    high_density_map[ds_index, region] = current_map[region]
                    checked_pixels.update(region) #update pixels that have been checked
                    
                    # Calculate the max because that is where regions will be "flagged" at
                    max_dEBV_pixel = region[np.argmax(current_map[region])]
                    
                    theta, phi = hp.pix2ang(nside, pixel, nest=True)#mark this pixel in the map
                    
                    #Add info to dictionary
                    region_info[ds_index].append({
                        'center': (theta, phi),
                        'region_pixels': region,
                        'region_values': current_map[region],
                    })

    return region_info, high_density_map

def angular_distance(theta1, phi1, theta2, phi2):
    """
    Function to calculate the angular distance between two points on the sphere
    
    Parameters:
    theta1, phi1: float, spherical coordinates of point 1 in radians
    theta2, phi2: float, spherical coordinates of point 2 in radians

    Output:
    float, angular distance between the two points
    """
    dtheta = theta1 - theta2
    dphi = np.abs(phi1 - phi2)
    if dphi > np.pi:
        dphi = 2 * np.pi - dphi
    return np.sqrt(dtheta**2 + dphi**2)

def get_region_maps(region_info, nside, ds_index, filter = False, rot = None, radius = None, combine=False):

    '''
    Function that gets individual maps of each region, and has to option to filter regions around a point, and to combine all region maps
    into one large map

    Parameters:
    region_info: list of lists of dictionaries, contains the center, pixel and pixel values for each region at each distance slice
    nside: int, nside of dEBV map
    ds_index: int, distance slice index
    filter: bool, whether to filter regions around a point
    rot (only if filter = True): tuple, (lon, lat) of point to filter around in degrees
    radius (only if filter = True): float, radius around point to filter in radians
    combine: bool, whether to combine all region maps into one large map

    Output:
    region_maps: list of np.arrays, contains individual region maps
    combined_map: np.array, map of all regions combined, done if combine = True

    '''

    region_maps = [] #List of maps for each individual region
    

    infilter_count = 0 #Keep track of how many regions in filter

    if filter == True:

        print(rot) #center of region 

        theta_filter, phi_filter = np.radians(90. - rot[1]), np.radians(rot[0]) #convert to theta phi

    for region in region_info[ds_index]:

        #get variables of region 
        theta, phi = np.array(region['center'])
        region_pixels = region['region_pixels']
        region_dEBV = region['region_values']

        individ_map = np.zeros(hp.nside2npix(nside)) #Map for individual regions
        individ_map[region_pixels] = region_dEBV

        if filter == True:

            dist = angular_distance(theta, phi, theta_filter, phi_filter) #get angular distance between region and point

            if dist < radius:

                region_maps.append(individ_map.copy()) #add to list if in filter
                infilter_count += 1 #add to count

            else:
                continue

        else:
            print('Not filtering')
            region_maps.append(individ_map.copy()) #if just want individual regions without filter

    if filter == True:
        print(f"Number of regions in filter: {infilter_count}")

    #If want to make combined map
    if combine == True:
        combined_map = np.zeros(hp.nside2npix(nside))
        for map in region_maps:
            combined_map += map
        return combined_map
    
    return region_maps    

## 3: Functions related to RGB
def get_RGB(dens_temp):
    '''
    Function that gets the R, G and B color channels from the normalized map of dust density and temperature emission. It also scales the values
    so that they're in between 0 and 255 (uint8).

    Parameters:
    dens_temp: normalized map of shape (frequency x distance_bin x pixel) containing the multiplication of dEBV and B

    Output:
    R, G, B: arrays of type uint8 containing the R, G and B color channels.
    '''
    nslices = dens_temp.shape[1]
    pixels = dens_temp.shape[2]
    max_pixel = 255 #maximum value allowed in an image
    

    R_array = np.zeros((nslices, pixels))
    G_array = np.zeros((nslices, pixels))
    B_array = np.zeros((nslices, pixels))

    for ds_index in range(nslices):
        R_array[ds_index] = dens_temp[0, ds_index] #Choose channel based on frequency of dens_temp
        G_array[ds_index] = dens_temp[1, ds_index]
        B_array[ds_index] = dens_temp[2, ds_index]

    #Make it so that its in between 0 and 255 and convert to uint8
    R = (R_array*max_pixel).astype(np.uint8) 
    G = (G_array*max_pixel).astype(np.uint8)
    B = (B_array*max_pixel).astype(np.uint8)

    return R, G, B


def get_rgb_ratios(R, G, B, title, image_name):

    R_flat = R.flatten()
    G_flat = G.flatten()
    B_flat = B.flatten()

    ratio_RG = np.divide(R_flat, G_flat, where=G_flat!=0)
    ratio_RB = np.divide(R_flat, B_flat, where=B_flat!=0)
    ratio_GB = np.divide(G_flat, B_flat, where=B_flat!=0)

    plt.hist(ratio_RG, bins=30, color='olive', alpha=0.5, label='R/G')
    plt.hist(ratio_RB, bins=30, color='indigo', alpha=0.5, label='R/B')
    plt.hist(ratio_GB, bins=30, color='turquoise', alpha=0.5, label='G/B')
    plt.yscale('log')
    plt.xlabel('Color Depth Ratio')
    plt.ylabel('Number of Pixels')
    plt.title(title)
    plt.legend(loc='upper right')
    plt.savefig(image_name)
    plt.show()

    return ratio_RG, ratio_RB, ratio_GB

def get_rgb_roi(R, G, B, x, y, w, h):
    '''
    Function that gets the RGB values of a region of interest (roi) in an image

    Parameters:
    R, G, B : 2d numpy arrays that are from 0 to 255
    x, y : coordinates of the top left corner of the roi
    w, h : width and height of the roi

    Returns:
    R_roi, G_roi, B_roi : 2d numpy arrays of the roi
    '''
    R_roi = R[y:y+h, x:x+w]
    G_roi = G[y:y+h, x:x+w]
    B_roi = B[y:y+h, x:x+w]

    return R_roi, G_roi, B_roi
    
           
    
