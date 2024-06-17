#This script contains any functions required to make plots in the 2D_RGBmap.ipynb jupyter notebook
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os
from PIL import Image, ImageOps

# 1: Healpy plotting functions
def plot_3D_temperature_slice_maps(data_dict):

    '''
    Function to ONLY plot the original temperature map of shape (pixel x distance_bin) at each distance slice
    '''
    Ts = data_dict['temperatures']
    model_nslices = data_dict["nr_of_distance_bins"]
    model_dist_slices = data_dict["distance_slices"]
    for ds_index in range(model_nslices):                                 
        hp.mollview(Ts[:,ds_index],title=r"$T$ at distance slice "+str(ds_index) +\
                                   " at "+'{:.2f}'.format(model_dist_slices[ds_index])+" kpc",nest=True,min=10,max=25, unit='K')
        #plt.savefig("T_128_"+str(ds_index)+".png")

def plot_dEBV(data_dict):
    """ Plot the reddening in each distance bin
    """ 
    dEBV = data_dict["dEBV"]
    model_nslices = data_dict["nr_of_distance_bins"]
    model_dist_slices = data_dict["distance_slices"]
    for ds_index in range(model_nslices):
        hp.mollview(dEBV[ds_index], title="Differential E(B-V) at distance slice "+str(ds_index) +\
                     " at "+'{:.2f}'.format(model_dist_slices[ds_index])+" kpc", nest=True, max=1)
        
def plot_map(data_dict, map, min_map, max_map, title_map, unit_map):
    """ 
    Function to plot any map of the whole sky as long as it has the form distance bin x pixel

    Parameters:
    data_dict: dictionary containing the data
    map: map to be plotted with shape (distance_bin x pixel) so if it has freq have to select freq before
    min_map: float, minimum value of the map
    max_map: float, maximum value of the map
    title: string, title of the plot
    unit_map: string, unit of the map

    Output:
    Plots the map at each distance slice
    """ 
    model_nslices = data_dict["nr_of_distance_bins"]
    model_dist_slices = data_dict["distance_slices"]

    for ds_index in range(model_nslices):
        
        hp.mollview(map[ds_index],title=f"{title_map} at distance slice "+str(ds_index) +\
                                        " at "+'{:.2f}'.format(model_dist_slices[ds_index])+" kpc",nest=True,min=min_map, max=max_map, unit=unit_map)

def plot_map_region(map, distance, longitude, latitude, x, y, min_map, max_map, title_map, unit_map):
    '''
    Function used to plot a specific region of the map which requires hp.gnomview instead of hp.mollview. 

    Parameters:
    map: map to be plotted with shape (distance_bin x pixel) so if it has freq have to select freq before
    distance: int, distance slice to be plotted
    longitude: float, longitude (degrees) of the centre of the region to be plotted
    latitude: float, latitude (degrees) of the centre of the region to be plotted
    x: int, size of the x axis of plot
    y: int, size of the y axis of plot
    min_map: float, minimum color value of the map
    max_map: float, maximum color value of the map
    title_map: string, title of the plot
    unit_map: string, unit of the map

    Output:
    A gnomview plot of the map at the specified region
    '''
    hp.gnomview(map[distance], rot=(longitude,latitude), title=title_map, nest=True, xsize=x, ysize=y, min=min_map, max=max_map, unit=unit_map, notext=True)

# 2: Matplotlib plotting functions
def plot_RGB_histogram(R, G, B, title, image_name):

    '''
    Function to create of plot of histograms of the color depth of the RGB channels of an image. 

    Parameters:
    R, G, B: numpy arrays, the RGB channels of the image
    title: string, title of the plot
    image_name: string, name of the image to be saved that should include the path

    Output:
    A single histogram with the color depth of the RGB channels of the image
    '''
        
    #Plot R
    plt.hist(R.flatten(), bins=50, color='red', alpha=0.5, label='R') #arrays are 2d so we need to flatten them
    #Plot G
    plt.hist(G.flatten(), bins=50, color='green', alpha=0.5, label='G')
    #Plot B
    plt.hist(B.flatten(), bins=50, color='blue', alpha=0.5, label='B')

    plt.legend(loc='upper right')
    plt.yscale('log')
    plt.xlabel('Color Depth')
    plt.ylabel('Number of Pixels')
    plt.title(title)
    plt.savefig(image_name)
    plt.show()

# 3: Getting Images
def create_image(R, G, B):
    '''
    Function that when given R, G, B arrays that are already 2d and normalized, will create an image. 

    Parameters:
    R, G, B : 2d numpy arrays that are from 0 to 255

    Returns:
    RGB_image : RGB image object
    '''
    
    R_uint8 = R.astype(np.uint8)
    G_uint8 = G.astype(np.uint8)
    B_uint8 = B.astype(np.uint8)

    R_image = Image.fromarray(R_uint8)
    G_image = Image.fromarray(G_uint8)
    B_image = Image.fromarray(B_uint8)

    RGB_image = Image.merge('RGB', (R_image, G_image, B_image))

    return RGB_image

def get_sky_image(data_dict, R, G, B, scale=True):
    '''
    Function that creates an RGB image of the whole sky using the R, G and B color channels. It produces an image at each distance slice

    Parameters:
    data_dict: dictionary containing the data, used to obtain number of distance slices
    R, G, B: arrays of type uint8 containing the R, G and B color channels
    scale: boolean, if True, scales the RGB values so that they expand over the whole range

    Output:
    Saves the RGB image of the whole sky
    '''
    #Trying to create images for whole sky

    dist_nslice = data_dict["nr_of_distance_bins"]

    directory = 'RGB_images' #directory to save images

    target_size = (2000, 1000) #have to make smaller size or else will not load

    # Create the directory if it does not exist
    if not os.path.exists(directory):
        os.makedirs(directory)

    for ds_index in range(dist_nslice):

        #Now make 2d for image purposes
        R_array_2d = hp.mollview(R[ds_index], nest=True, return_projected_map=True, title="R Channel", cbar=False) #return_projected_map=True is what returns the 2d "array" needed
        plt.close() #Close because actual image is not important
        G_array_2d = hp.mollview(G[ds_index], nest=True, return_projected_map=True, title="G Channel", cbar=False)
        plt.close()  
        B_array_2d = hp.mollview(B[ds_index], nest=True, return_projected_map=True, title="B Channel", cbar=False)
        plt.close()  

        R_array = np.array(R_array_2d) #Convert to numpy array
        G_array = np.array(G_array_2d)
        B_array = np.array(B_array_2d)

        #If desired, scale the values so that they span entire range
        if scale==True:
            max_pixel = 255
            max_intensity = np.max([np.max(R_array), np.max(G_array), np.max(B_array)])
            scale_factor = max_pixel/max_intensity

            R_array = R_array*scale_factor
            G_array = G_array*scale_factor
            B_array = B_array*scale_factor

        #gnomview returns inverted array, so have to flip it to get correct image
        R_flipped = np.flipud(R_array)
        G_flipped = np.flipud(G_array)
        B_flipped = np.flipud(B_array)

        #Convert to uint8 for image
        R_uint8 = R_flipped.astype(np.uint8)
        G_uint8 = G_flipped.astype(np.uint8)
        B_uint8 = B_flipped.astype(np.uint8)

        #Convert into image
        R_image = Image.fromarray(R_uint8)
        G_image = Image.fromarray(G_uint8)
        B_image = Image.fromarray(B_uint8)

        RGB_image = Image.merge("RGB", (R_image, G_image, B_image))

        #Because image is of the whole sky, have to resize it so that it requires less computing power
        RGB_image = ImageOps.fit(RGB_image, target_size, method=Image.Resampling.LANCZOS)

        RGB_image.save(f'{directory}/allsky_rgb_{ds_index}.png')

def get_region_image(R, G, B, dist, longitude, latitude, x, y, scale=False):
    '''
    Function that gets the RGB image of a region in the sky. Differs from the whole sky due to the use of hp.gnomview. Also provides an option to scale
    the RGB values so that they expand over whole range and in turn makes the image brghter. 

    Parameters:
    R, G, B: arrays of type uint8 containing the R, G and B color channels *Note, does not have to be in uint8 format, code converts it into this anyways
    dist: int, distance slice to be plotted
    longitude: float, longitude (degrees) of the centre of the region to be plotted 
    latitude: float, latitude (degrees) of the centre of the region to be plotted
    x: int, size of the x axis of the plot
    y: int, size of the y axis of the plot
    scale: boolean, if True, scales the RGB values so that they expand over the whole range

    Output:
    R_uint, G_uint, B_uint: arrays of type uint8 containing the R, G and B color channels for that region
    RGB_img: Image object containing the RGB image of the region
    '''

    #Use hp.gnomview to get a 2d array of the region that will be used to create the image

    R_2darray = hp.gnomview(R[dist], rot=(longitude,latitude), nest=True, xsize=x, ysize=y, return_projected_map=True) #return_projected_map=True is what returns the 2d "array" needed 
    plt.close() #Close because actual image is not important
    G_2darray = hp.gnomview(G[dist], rot=(longitude,latitude), nest=True, xsize=x, ysize=y, return_projected_map=True)
    plt.close()
    B_2darray = hp.gnomview(B[dist], rot=(longitude,latitude), nest=True, xsize=x, ysize=y, return_projected_map=True)
    plt.close()

    #Output is not in numpy array so have to convert it
    R_array = np.array(R_2darray)
    G_array = np.array(G_2darray)
    B_array = np.array(B_2darray)

    #Code to scale values so that they span entire range
    if scale==True:
        max_pixel = 255
        max_intensity = np.max([np.max(R_array), np.max(G_array), np.max(B_array)])
        scale_factor = max_pixel/max_intensity

        R_array = R_array*scale_factor
        G_array = G_array*scale_factor
        B_array = B_array*scale_factor

    #gnomview returns inverted array, so have to flip it to get correct image
    R_flipped = np.flipud(R_array)
    G_flipped = np.flipud(G_array)
    B_flipped = np.flipud(B_array)

    #Convert to uint8 for image
    R_uint = R_flipped.astype(np.uint8)
    G_uint = G_flipped.astype(np.uint8)
    B_uint = B_flipped.astype(np.uint8)

    #Convert into image
    R_img = Image.fromarray(R_uint)
    G_img = Image.fromarray(G_uint)
    B_img = Image.fromarray(B_uint)

    #Make inmage
    RGB_img = Image.merge("RGB", (R_img, G_img, B_img))

    return R_uint, G_uint, B_uint, RGB_img

def create_panel(maps_dict, frequency, dist, longitude, latitude, plot_title, image_path):
    '''
    Function to create a 5 x 3 panel of images of the Cepheus LMC region. The function first saves the images and then calls them when creating the panel.
    
    The panel will contain the following images:
    - Original Temperature Map at nside 32
    - New Temperature Map at nside 1024
    - E(B-V) Map
    - Planck Map at 588, 1106, and 1765 GHz
    - Normalized Density x Temperature Map at 588, 1106, and 1765 GHz
    - RGB histogram and ratios
    - RGB scaled image

    Parameters:
    maps_dict: dictionary containing the maps
    frequency: numpy array, array of frequencies
    dist: int, distance slice to be plotted
    longitude: float, longitude (degrees) of the centre of the region to be plotted
    latitude: float, latitude (degrees) of the centre of the region to be plotted
    plot_title: string, title of the region
    image_path: string, path to save the images, also used to call histograms so make sure it is the same as the path for those

    Output:
    A 5 x 3 panel of images of the desired region
    '''
    Ts = maps_dict["Temp_og"]
    Ts_new = maps_dict["Temp_new"]
    dEBV = maps_dict["Density"]
    planck = maps_dict["Planck"]
    densxtemp = maps_dict["Normalized_denstemp"]

    ##First save images needed for the panel
    hp.gnomview(Ts[:,dist], rot=[longitude, latitude], xsize=2000, ysize=2000,title=f'$T$ of {plot_title} with nside 32', nest=True,min=10,max=25, unit='K', notext=True)
    plt.savefig(image_path + "T_og.png")
    plt.close()
    hp.gnomview(Ts_new[dist], rot=[longitude, latitude], xsize=2000, ysize=2000,title=f'$T$ of {plot_title} with nside 1024', nest=True,min=10,max=25, unit='K', notext=True)
    plt.savefig(image_path + "T_new.png")
    plt.close()
    hp.gnomview(dEBV[dist], rot=[longitude, latitude], xsize=2000, ysize=2000,title=f'$E(B-V)$ of {plot_title}', nest=True, max=1, notext=True)
    plt.savefig(image_path + "dEBV.png")
    plt.close()


    for f_index in range(3):
    
        hp.gnomview(planck[f_index, dist], rot=[longitude, latitude], xsize=2000, ysize=2000,title=f'Temperature Emission $B$ of {plot_title} at {int(frequency[f_index])} GHz', nest=True,min=1e-17, max=2e-15, unit=r'$\mathrm{W/M^{2}\cdot sr \cdot Hz}$', notext=True)
        plt.savefig(image_path + f"planck_{int(frequency[f_index])}.png")
        plt.close()
        hp.gnomview(densxtemp[f_index, dist], rot=[longitude, latitude], xsize=2000, ysize=2000,title=f'Normalized $B$ and $E(B-V)$ of {plot_title} at {int(frequency[f_index])} GHz', nest=True, min = 0, max = 0.2, unit=None, notext=True)
        plt.savefig(image_path + f"dxT_norm_{int(frequency[f_index])}.png")
        plt.close()

    ##Now create the panel
    

    #Creating panel of images of Cepheus LMC
    panel = plt.figure(figsize=(18, 30))
    gs = gridspec.GridSpec(5, 3, panel)

    #Get image paths and names to not repeat code

    img_names = ['T_og', 'T_new', 'dEBV', 'planck_588', 'planck_1106', 'planck_1765', 'dxT_norm_588', 'dxT_norm_1106', 'dxT_norm_1765', 'R', 'G', 'B', 'rgb_hist_scaled', 'rgb_hist_ratios', 'rgb_scaled'] #array of names

    for i in range(5):  # Assuming 5 rows
        for j in range(3):  # Assuming 3 columns
            panel.add_subplot(gs[i, j]) #add subplots
            img_index = i * 3 + j  #Can avoid iterating over the image index too
            img = Image.open(image_path + img_names[img_index] + '.png') #Have to string them all together

            if img_index < len(img_names): #Making a check 
                plt.imshow(img)
            plt.axis('off')  

            if img_names[img_index] == 'rgb_scaled': #This image has no title so have to manually add title
                plt.title(f'RGB Image of {plot_title}')

    plt.subplots_adjust(wspace=0, hspace=0, top=0.97)
    plt.suptitle(f'{plot_title} Panel')

    plt.savefig(image_path + 'panel.png', facecolor='white', edgecolor='none')

    plt.show()
