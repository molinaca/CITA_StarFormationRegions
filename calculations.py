#This script contains any functions needed to perform calculations with arrays in the main script
import numpy as np
import healpy as hp
from constants import h, c, k
import matplotlib.pyplot as plt



# 1: Scientific Functions

def Planck (nu, Tmap):
    '''
    Function to calculate the Planck function at a frequency nu for a given temperature map

    Parameters:
    nu: frequency in Hz
    Tmap: temperature map in Kelvin
    h: planck constant from constants.py
    c: speed of light from constants.py
    k: boltzmann constant from constants.py
    x_d: variable used to simplify planck function

    Output:
    B = 2*h*nu**3/c**2/(np.exp(x_d)-1) : planck function in units W/m^2/Hz/sr
    '''
    x_d = h*nu/(k*Tmap)
    return 2*h*nu**3/c**2/(np.exp(x_d)-1)

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

# 2: Manipulating Arrays 
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

def multiply_dEBVandB(data_dict, dEBVmap, Bmap, frequency):
    '''
    Function that multiplied the dust density (dEBV) by the temperature emission (B) at each distance slice and frequency. 

    Parameters:
    data_dict: dictionary containing the data, used to obtain number of distance slices
    dEBVmap: dust density map with shape (distance_bin x pixel)
    Bmap: temperature emission map with shape (frequency x distance_bin x pixel)
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
            dens_temp[f_index, ds_index] = dEBVmap[ds_index]*Bmap[f_index, ds_index]
    
    return dens_temp

def normalize_dEBVandB(data_dict, dens_temp, frequency):
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

# 3: Functions related to RGB
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
    R = (R_array*255).astype(np.uint8) 
    G = (G_array*255).astype(np.uint8)
    B = (B_array*255).astype(np.uint8)

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
           
    