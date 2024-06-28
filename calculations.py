#This script contains any functions needed to perform calculations with arrays in the main script
import numpy as np
import ephem
import healpy as hp
from constants import h, c, k
import matplotlib.pyplot as plt



## 1: Scientific Functions

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
    
           
    
