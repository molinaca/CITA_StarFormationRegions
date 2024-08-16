#This script contains any functions required to make plots in the 2D_RGBmap.ipynb jupyter notebook
import numpy as np
import ephem
import healpy as hp
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import storage_management as sm
import os
from PIL import Image, ImageOps
import calculations as calc

## 1: General Plotting Functions:

### 1.1: Healpy 
def plot_map(data_dict, map, min_map, max_map, title_map, unit_map, nside, declination_mask = False):
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
    npix = hp.nside2npix(nside)
    pixel_index_array = np.arange(npix)

    for ds_index in range(model_nslices):

        map_title = f"{title_map} at distance slice "+str(ds_index) +\
                                        " at "+'{:.2f}'.format(model_dist_slices[ds_index])+" kpc"
        
        plot_array = np.zeros(npix)

        if declination_mask ==True:
            dec = create_declination_mask(nside)
            data_masked = map[ds_index].copy()
            data_masked[dec<-30]=hp.UNSEEN
            plot_array[pixel_index_array]=data_masked

        else:
            plot_array[pixel_index_array]=map[ds_index]
            
        hp.mollview(plot_array, title = map_title,nest=True, min=min_map, max=max_map, unit=unit_map)
        plt.title(map_title, fontsize = 16)
        cbar = plt.gcf().axes[-1]
        cbar.tick_params(labelsize=15) 

def create_declination_mask(nside, nested=True):
        ### mask the map at np.abs(declinations) larger than 30
        ### we used nested data sets in this analysis
        ### inspired from https://stackoverflow.com/questions/29007648/pyephem-coordinate-transformation-galactic-to-equatorial
        ### also inspired from this https://healpy.readthedocs.io/en/latest/tutorial.html
        ### get the pixel indices
        npix = hp.nside2npix(nside)
        pixel_array = np.array(range(npix))
        ### get the lat long coordinates to then mask over them
        ls, bs=hp.pixelfunc.pix2ang(nside=nside, ipix=pixel_array, nest=nested, lonlat=True)
        #### converting galactic coordinates to equatorial coordinates
        ra = np.zeros(npix)
        dec = np.zeros(npix)
        for index in range(npix):
            l = ls[index]
            b = bs[index]
            galactic = ephem.Galactic(l/180.0*np.pi,b/180.0*np.pi)
            equatorial = ephem.Equatorial(galactic, epoch=ephem.J2000)
            ra[index]=equatorial.ra/np.pi*180.0
            dec[index]=equatorial.dec/np.pi*180.0
            #print('%.13f %.13f' % (equatorial.ra/np.pi*180.0, equatorial.dec/np.pi*180.0))
        return dec

        
def plot_map_region(map, distance, longitude, latitude, title_map, x=None, y=None, min_map=None, max_map=None, unit_map=None):
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
    hp.gnomview(map[distance], rot=(longitude,latitude), title=title_map, nest=True, xsize=x, ysize=y, min=min_map, max=max_map, 
                cbar=True, unit=unit_map, notext=True)
    plt.title(title_map, fontsize = 16)
    cbar = plt.gcf().axes[-1]
    cbar.tick_params(labelsize=15) 

def make_label_dict(title, cbar_label = None, x_label = None, y_label = None, sup_title = None):
    '''
    Make a dictionary containing labels for a plot. 

    Parameters:
    title: string, title of the plot

        Optional:
        cbar_label: string, label of the colorbar
        x_label: string, label of the x-axis
        y_label: string, label of the y-axis
        sup_title: string, super title of the plot

    Output:
    label_dict: dictionary, contains the labels from parameters

    '''
    label_dict = {'title': title}
    if cbar_label:
        label_dict['cbar_label'] = cbar_label
    if x_label:
        label_dict['x_label'] = x_label
    if y_label:
        label_dict['y_label'] = y_label
    if sup_title:
        label_dict['sup_title'] = sup_title
    return label_dict


def make_scatterplot_dict(long, lat, lonlat, color, marker, size, alpha=None, label=None):

    '''
    Makes a dictionary that contains the parameters used for hp.projscatter scatter plots. Note that long, lat must be a 1D array 
    (so this function must be run for each distance slice).

    Parameters:
    long, lat: longitude and latitude of the object in degrees
    lonlat: bool, if True, long and lat are in lonlat format, if False, they are in theta, phi format
    color: color of the marker
    marker: marker style
    size: size of the marker
    alpha (optional): transparency of the marker
    label (optional): label of the marker

    Output:
    scatter_plot_dict: dictionary containing the parameters for the scatter plot
    '''
    scatter_plot_dict = {
        'coords' : [long, lat],
        'lonlat' : lonlat,
        'color' :  color,
        'marker' : marker,
        #Optional arguments
        'alpha' : alpha,
        'size': size,
        'label' : label
    }

    return scatter_plot_dict

def make_scatter_plot(scatter_plot_dict, legend = True, scatter_type = 'healpy', subplot = False, ax = None):
    '''
    Uses the dictionary made by make_scatterplot_dict to create a scatter plot on a healpy map (could be gnonview or mollview).
    Provides the option to add a legend.

    Parameters:
    scatter_plot_dict: dictionary containing the parameters for the scatter plot
    legend: bool, if True, adds a legend to the plot

    Output:
    A scatter plot on a healpy map
    '''
    long, lat = scatter_plot_dict['coords']
    lonlat_bool = scatter_plot_dict['lonlat']
    color = scatter_plot_dict['color']
    marker = scatter_plot_dict['marker']
    alpha = scatter_plot_dict['alpha']
    size = scatter_plot_dict['size']
    label = scatter_plot_dict['label']

    if scatter_type == 'healpy':
        hp.projscatter(long, lat, lonlat=lonlat_bool, c=color, marker=marker, alpha=alpha, s=size, label=label)
    elif scatter_type == 'matplotlib':
        if subplot == True:
            ax.scatter(long, lat, c=color, marker=marker, alpha=alpha, s=size, label=label)
        else:
            plt.scatter(long, lat, c=color, marker=marker, alpha=alpha, s=size, label=label)
    else:
        print('Invalid scatter type, please choose healpy or matplotlib')
    
    if legend == True:
        plt.legend(fontsize=14)
        
def make_matplotlib_scatter_plot(scatter_plot_dict, proj_map, subplot = False, ax = None, boundaries = None):
    '''
    Performs conversions and calculations to simulate the effect of hp.projscatter but with matplotlib.scattter instead. Provides the option to include this in
    a subplot to as a standalone figure. 

    Parameters:
    scatter_plot_dict: dictionary, contains the following keys:
        'coords': list of longitudes and latitudes
        'lonlat': bool, if True, then the coordinates are in lonlat, if False, then in radians
        'color': string, color of the points
        'marker': string, marker of the points
        'alpha': float, transparency of the points
        'size': float, size of the points
        'label': string, label of the points
    proj_map: array, gnomonoic projection from hp.projector.GnomonicProj
    subplot: bool, if True, then the scatter plot will be included in a subplot, if False, then it will be a standalone figure
    ax: axis, if subplot is True, then the axis must be provided
    boundaries: list, boundaries of the plot, if provided, then points outside the boundaries will be masked

    Output:
    scatter plots of the dictionaries provided
    '''

    #Extract valuable information from dictionary
    coords = np.array(scatter_plot_dict['coords'])
    long, lat = coords[0], coords[1]
    lonlat_bool = scatter_plot_dict['lonlat']
    color = scatter_plot_dict['color']
    marker = scatter_plot_dict['marker']
    alpha = scatter_plot_dict['alpha']
    size = scatter_plot_dict['size']
    label = scatter_plot_dict['label']

    #Convert to lonlat if necessary
    if lonlat_bool == False:
        long, lat = calc.convert_to_lonlat(long, lat)

    #Use hp.projector to positions of x and y on projection
    x, y = proj_map.ang2xy(long, lat, lonlat=True)

    #Mask points outside the boundaries
    if boundaries:
        #Define boundaries
        x_lower, x_upper = boundaries[0], boundaries[1]
        y_lower, y_upper = boundaries[2], boundaries[3]

        #Mask points 
        x_mask = (x >= x_lower) & (x <= x_upper)
        y_mask = (y >= y_lower) & (y <= y_upper)

        x = x[x_mask & y_mask]
        y = y[x_mask & y_mask]

    #Now plot them depending on if subplot or not
    if subplot == True:
        ax.scatter(x, y, c=color, marker=marker, alpha=alpha, s=size, label=label)
        ax.legend(fontsize=14)
    elif subplot == False:
        plt.scatter(x, y, c=color, marker=marker, alpha=alpha, s=size, label=label)
        plt.legend(fontsize=14)

    else:
        print('Error: subplot must be True or False')


def get_fov(xsize, reso, degree=True):
    '''
    Get the field of view of a gnomview plot using the xsize and resolution, can return in degrees or radians. Only works of xsize == ysize. 

    Parameters:
    xsize: int, size of x axis in gnomview plot
    reso: float, resolution of plot
    degree: bool, true if want fov returned in degrees

    Output:
    fov: float, field of view of gnomview plot in degrees or radians
    '''
    #Get fov using xsize and reso
    fov = xsize*(reso/60)

    #If want fov in radians instead of degrees make conversion
    if degree==False:
        fov = np.radians(fov)

    return fov

def get_axis_skycoords(l, b, spacing, xsize, reso, degree=True, label = True, prec = None):
    '''
    Get x and y axis ticks and labels as sky coordinates from the fov of the gnomview plot. 

    Parameters:
    l, b: float, longitude and latitude in degrees, would also work for theta, phi or RA and DEC. 
    spacing: float, spacing between ticks
    xsize: int, xsize of gnomview plot
    reso: float, resolution of gnomview plot
    degree: bool, units of fov, default==True
    label: bool, will create custom x and y labels
    prec: int, optional, number of decimal points for labels
    '''
    #Get fov
    fov = get_fov(xsize, reso, degree=degree)

    #Want half of fov because will add and subtract from centre
    half_fov = fov/2

    #Get min and max for each axis
    x_lower = l - half_fov 
    x_upper = l + half_fov
    y_lower = b - half_fov
    y_upper = b + half_fov
    
    #Get x and y ticks
    x_ticks = np.linspace(x_lower, x_upper, spacing)
    y_ticks = np.linspace(y_lower, y_upper, spacing)

    #If also want to get labels make them here and return all
    if label == True:
        x_label = ['{:.{}f}'.format(x, prec) for x in x_ticks]
        y_label = ['{:.{}f}'.format(y, prec) for y in y_ticks]

        return x_ticks, y_ticks, x_label, y_label

    else:
        return x_ticks, y_ticks
        


## 2: Specific Plotting Functions

### 2.1: Temperature
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

### 2.2: dEBV
def plot_dEBV(data_dict):
    """ Plot the reddening in each distance bin
    """ 
    dEBV = data_dict["dEBV"]
    model_nslices = data_dict["nr_of_distance_bins"]
    model_dist_slices = data_dict["distance_slices"]
    for ds_index in range(model_nslices):
        hp.mollview(dEBV[ds_index], title="Differential E(B-V) at distance slice "+str(ds_index) +\
                     " at "+'{:.2f}'.format(model_dist_slices[ds_index])+" kpc", nest=True, max=1)

### 2.3: Flagging Regions with Certain Properties

def overplot_regions_mollview(region_info, map, dist_slices):

    '''
    Function that overplots the maximum/center of the high_density regions over an hp.mollview map

    Parameters:
    region_info: list of lists of dictionaries, contains the center, pixel and pixel values for each region at each distance slice
    map: np.array, map that regions will be plotted over in form (distance_slices, pixel)
    dist_slices: int, number of distance slices

    Output:
    plots map and a scatter plot of the maximum of the regions at each distance slice
    '''

    #Get centers of regions
    region_centers = calc.get_region_centers(region_info, dist_slices)

    for ds_index in range(dist_slices):
        hp.mollview(map[ds_index], title=f'High density regions slice {ds_index}', nest=True, cbar=True)
        plt.title(f'High density regions slice {ds_index}', fontsize = 16)

        #Make sure region_info exists
        if region_centers[ds_index].size > 0:
            hp.projscatter(region_centers[ds_index][0], region_centers[ds_index][1], lonlat=True, s=8, marker='o', color='red')
            plt.show()

        else:
            print(f"No high density regions at distance slice {dist_slices}")
            plt.close()

def overplot_region_gnomview(region_centers, map, ds_index, rot, title, xsize=None, ysize=None, unit=None, save=True, filename=None, show=False ):

    '''
    Function that overplots region centers on a gnomview map. 

    Parameters:
    region_info (list of lists): list of lists of dictionaries containing information about regions.
    map : np.array, map to plot
    ds_index : int, index of distance slice
    rot : tuple, center of map in (theta, phi) coordinates
    title : str, title of plot
    xsize (optional): int, size of x-axis in gnomview map
    ysize (optional): int, size of y-axis in gnomview map
    unit (optional): str, unit of map
    save (optional): bool, default = True, if True, saves the plot
    filename (optional): str, name of file to save plot, only use if save=True
    show (optional): bool, default = False, if True, shows the plot

    Returns:
    plot of gnomview map with region centers overplotted
    '''
    #Get centers of regions
    regions_atdist = region_centers[ds_index]

    if regions_atdist.size > 0:

        hp.gnomview(map[ds_index], rot=rot, title=title, xsize = xsize, ysize = ysize,unit=unit, nest=True, cbar=True, notext=True)
        hp.projscatter(regions_atdist[0], regions_atdist[1], s=15, marker='o', color='red')
        plt.legend(['Region Centers'], fontsize=12)
        plt.title(title, fontsize = 16)

        #Increase size of colorbar
        cbar = plt.gcf().axes[-1]
        cbar.tick_params(labelsize=15)

        if save == True:
            plt.savefig(filename, bbox_inches='tight', pad_inches=0.1)

        if show == True:
            plt.show()
        else:
            plt.close()

def make_gnomview_matplotlib(map, rot, size, labels, nside, ax, fig, grid = True, scatter = True, scatter_list = None):
    '''
    Function that simulated the same output as overplot_region_gnomview but with matplotlib. It uses hp.projector.GnomonicProj to project the map
    and provides the option to overplot a scatter plot on top of the map. Only works if a subplot, should change this in the future. 

    Parameters: 
    map: array, map to be plotted
    rot: list, rotation of the map in form [long, lat] in degrees
    size: int list, size of the map in form [xsize, ysize]
    labels: dictionary, contains the following keys:
        'title': string, title of the plot
        'cbar_label': string, label of the colorbar (optional)
        'x_label': string, label of the x-axis (optional)
        'y_label': string, label of the y-axis (optional)
        'sup_title': string, super title of the plot (optional)

    nside: int, nside of the map
    ax: axis, axis of the plot
    fig: figure, figure of the plot
    grid: bool, if True, then grid will be plotted, if False, then not
    scatter: bool, if True, then scatter plot will be plotted, if False, then not
    scatter_list: list, list of dictionaries containing the scatter plot information, only used if scatter = True

    Output:
    plot of the map with the scatter plot if scatter = True
    '''

    #Define variables
    l = rot[0]
    b = rot[1]
    xsize, ysize = size[0], size[1]

    #Get Gnomonic Projection
    proj_map = hp.projector.GnomonicProj(rot=[l, b], xsize=xsize, ysize=ysize)

    #Get field of view in degrees
    fov = proj_map.get_fov()
    fov_deg = np.degrees(fov)

    #In order to get image, need to define a function that GnomonicProj.projmap can use. Needs nside so defined within function
    def vec2pix_func(x, y, z):
        return hp.vec2pix(nside, x, y, z, nest=True)
    
    img = proj_map.projmap(map, vec2pix_func)

    #Get bounds in sky coordinates
    x_lower = l - fov_deg/2
    x_upper = l + fov_deg/2
    y_lower = b - fov_deg/2
    y_upper = b + fov_deg/2

    #Get bounds in from of extent
    edges_extent = proj_map.get_extent()

    #Plot image
    im = ax.imshow(img, origin='lower', extent=edges_extent)

    #Plot grid
    if grid:
        ax.grid(True, color='dimgrey')

    #Plot scatter plot
    if scatter:
        for plot in scatter_list:
            make_matplotlib_scatter_plot(plot, proj_map, subplot=True, ax=ax, boundaries=edges_extent)
    
    
    #x and y ticks will be extents
    x_extent = np.linspace(edges_extent[0], edges_extent[1], 5)
    y_extent = np.linspace(edges_extent[2], edges_extent[3], 5)

    ax.set_xticks(x_extent)  # Example for x ticks
    ax.set_yticks(y_extent)    # Example for y ticks

    #Set the labels using bounds
    x_ticks_label = np.linspace(x_lower, x_upper, 5)
    y_ticks_label = np.linspace(y_lower, y_upper, 5)

    ax.set_xticklabels(['{:.2f}'.format(x) for x in x_ticks_label[::-1]], fontsize=14)
    ax.set_yticklabels(['{:.2f}'.format(y) for y in y_ticks_label], fontsize=14)
    
    ax.set_xlabel('l (deg)', fontsize=15)
    ax.set_ylabel('b (deg)', fontsize=15)

    #Get labels from dictionary
    cbar_label = labels['cbar_label']
    title = labels['title']
    
    #Get color bar
    cbar = fig.colorbar(im, ax=ax, orientation='vertical', shrink=0.6)
    cbar.ax.tick_params(labelsize=14)
    cbar.set_label(cbar_label, fontsize=14)

    #Set title
    ax.set_title(title, fontsize=18)

## 3: RGB and Imaging
### 3.1: Getting RGB Images
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

def brighten_image(image, factor):

    '''
    Function to brighten the images by a factor

    Parameters:
    factor: factor to brighten the images by

    Output:
    R, G, B: brightened RGB arrays
    '''

    image_float = image.astype(np.float32) / 255.
    brightened_image = image_float * factor
    brightened_image = np.clip(brightened_image, 0.0, 1.0)
    image_float = np.clip(image_float, 0, 1)
    brightened_image = (brightened_image * 255).astype(np.uint8)
    return brightened_image

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
        R_array_2d = hp.mollview(R[ds_index], nest=True, return_projected_map=True, title="R Channel", cbar=False) 
        #return_projected_map=True is what returns the 2d "array" needed

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

def get_region_image(R, G, B, dist, longitude, latitude, x, y, scale=False, flip=True):
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

    R_2darray = hp.gnomview(R[dist], rot=(longitude,latitude), nest=True, xsize=x, ysize=y, return_projected_map=True) 
    #return_projected_map=True is what returns the 2d "array" needed 

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
    if flip ==True:
        R_flipped = np.flipud(R_array)
        G_flipped = np.flipud(G_array)
        B_flipped = np.flipud(B_array)
    else:
        R_flipped = R_array #Just don't want to have to change the name right now
        G_flipped = G_array
        B_flipped = B_array

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

### 3.2: Analysis      
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
    plt.xlabel('Color Depth', fontsize=14)
    plt.ylabel('Number of Pixels', fontsize=14)
    plt.title(title, fontsize = 16)
    plt.savefig(image_name, bbox_inches='tight', pad_inches=0.1)
    plt.show()

## 4: Panels
def create_panel(size, figsize, path, images, title, filename = None):
    panel, axs = plt.subplots(size[0], size[1], figsize=figsize)
    axs = np.atleast_2d(axs) 

    for i in range(size[0]):
        for j in range(size[1]):
            img = Image.open(path + f'/{images[i*size[1]+j]}.png')
            axs[i,j].imshow(img)
            axs[i,j].axis('off')
    
    plt.tight_layout(rect=[0, 0, 1, 0.96]) 
    panel.suptitle(title, fontsize=22)    
    if filename:
        plt.savefig(filename, bbox_inches='tight', pad_inches=0.1)
    plt.show() 

def create_rgb_panel(maps_dict, frequency, dist, longitude, latitude, plot_title, image_path):
    '''
    Function to create a 5 x 3 panel of images of the Cepheus LMC region. The function first saves the images and then calls them 
    when creating the panel.
    
    The panel will contain the following images:
    - Original Temperature Map at nside 32
    - New Temperature Map at nside 1024
    - E(B-V) Map
    - Temperature Tracer Map at chosen frequencies
    - Normalized Density x Temperature Map at chosen frequencies
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
    temptracer = maps_dict["Temperature Tracer"]
    densxtemp = maps_dict["Normalized_denstemp"]

    ##First save images needed for the panel

    #Can use plot_map_region for all of them except original temp map because format is pixel x distance 

    #Parameters for plots
    xsize = 2000
    ysize = 2000
    temp_min = 10
    temp_max = 25
    temp_unit = 'K'

    #Temperature og
    Ts_og_title = f'$T$ of {plot_title} with nside 32'
    hp.gnomview(Ts[:,dist], rot=[longitude, latitude], xsize=xsize, ysize=ysize,title=Ts_og_title, nest=True,min=temp_min,max=temp_max, unit='K', 
                notext=True)
    cbar = plt.gcf().axes[-1] #make cbar larger, have to manually do this too
    cbar.tick_params(labelsize=15)
    plt.title(Ts_og_title, fontsize = 16)
    plt.savefig(image_path + "T_og.png")
    plt.close()

    #Temperature new
    Ts_new_title = f'$T$ of {plot_title} with nside 1024'

    plot_map_region(Ts_new, dist, longitude, latitude, xsize, ysize, temp_min, temp_max, Ts_new_title, temp_unit)
    plt.title(Ts_new_title, fontsize = 16)
    plt.savefig(image_path + "T_new.png")
    plt.close()

    #dEBV
    dEBV_title = f'$E(B-V)$ of {plot_title}'
    dEBV_bounds = [0, 1]
    dEBV_unit = 'dEBV'
    plot_map_region(dEBV, dist, longitude, latitude, xsize, ysize, dEBV_bounds[0], dEBV_bounds[1], dEBV_title, dEBV_unit)
    plt.title(dEBV_title, fontsize = 16)
    plt.savefig(image_path + "dEBV.png")
    plt.close()


    for f_index in range(3):
    
        #Temperature Tracer Map at each frequency (i.e Gaussian)
        Ttracer_title = f'$T$ Tracer of {plot_title} at {frequency[f_index]} GHz'
        Ttracer_bounds = [None, None]
        Ttracer_unit = 'None'
        plot_map_region(temptracer[f_index], dist, longitude, latitude, xsize, ysize, Ttracer_bounds[0], Ttracer_bounds[0], 
                        Ttracer_title, Ttracer_unit )
        plt.title(Ttracer_title, fontsize = 16)
        plt.savefig(image_path + f"temptracer_{frequency[f_index]}.png")
        plt.close()

        #Normalized Density x Temperature Map at each frequency
        densxtemp_title = 'Normalized $T$ Tracer X $E(B-V)$ of' + '\n' + f' {plot_title} at {frequency[f_index]} GHz'
        densxtemp_bounds = [0, 0.2]
        densxtemp_unit = 'None'
        plot_map_region(densxtemp[f_index], dist, longitude, latitude, xsize, ysize, densxtemp_bounds[0], densxtemp_bounds[1], 
                        densxtemp_title, densxtemp_unit)
        plt.title(densxtemp_title, fontsize = 16)
        plt.savefig(image_path + f"dxT_norm_{frequency[f_index]}.png", bbox_inches='tight', pad_inches=0.1)
        plt.close()

    ##Now create the panel
    

    #Creating panel of images of Cepheus LMC
    panel = plt.figure(figsize=(22, 34))
    gs = gridspec.GridSpec(5, 3, panel)

    #Get image paths and names to not repeat code

    img_names = ['T_og', 'T_new', 'dEBV', f'temptracer_{frequency[0]}', f'temptracer_{frequency[1]}', f'temptracer_{frequency[2]}', 
                 f'dxT_norm_{frequency[0]}', f'dxT_norm_{frequency[1]}', f'dxT_norm_{frequency[2]}', 'R', 'G', 'B', 'rgb_hist_scaled', 
                 'rgb_hist_ratios', 'rgb_scaled'] #array of names

    for i in range(5):  # Assuming 5 rows
        for j in range(3):  # Assuming 3 columns
            panel.add_subplot(gs[i, j]) #add subplots
            img_index = i * 3 + j  #Can avoid iterating over the image index too
            img = Image.open(image_path + img_names[img_index] + '.png') #Have to string them all together

            if img_index < len(img_names): #Making a check 
                plt.imshow(img)
            plt.axis('off')  

            if img_names[img_index] == 'rgb_scaled': #This image has no title so have to manually add title
                plt.title(f'RGB Image of {plot_title}', fontsize = 16)

    plt.subplots_adjust(wspace=0.05, hspace=0.05, top=0.97)
    plt.suptitle(f'{plot_title} Panel', fontsize = 18)

    plt.savefig(image_path + 'panel.png', facecolor='white', edgecolor='none')

    plt.show()

def visualize_unknown_features(features_list,map, n_distslices, distslices, path, title, save_name, overplot=True, 
                               overplot_list = None, one_dist = False):
    #Loop to create panels
    for ds_index in range(n_distslices):
        dist_path = sm.join_path(path, f'Distance_{ds_index}')
        
        # Get the number of features at this distance slice
        num_features = len(features_list[ds_index])

        if num_features == 0:
            print(f'No unknown features at distance slice {ds_index}')
            continue

        if num_features > 5:
            nrows = 2

        else:
            nrows = 1

        ncols = int((num_features + nrows- 1) // nrows)
        
        # Create a figure with multiple subplots for the current distance slice
        fig, axes = plt.subplots(nrows, ncols, figsize=(10*ncols, 10*nrows))
        fig.patch.set_facecolor('white')

        if nrows > 1:
            axes = axes.flatten()
        
        for i in range(num_features):
            features_atdist = features_list[ds_index]

            # Get positions
            long_full = features_atdist[:, 0]
            lat_full = features_atdist[:, 1]

            # Format for title
            long = '{:.2f}'.format(long_full[i])
            lat = '{:.2f}'.format(lat_full[i])

            unknown_feature_title = 'Unknown Feature at ' + long + ', ' + lat + '\n' + 'with Hot and Cold Regions'

            # Select the appropriate axis for the current feature
            ax = axes[i] if num_features > 1 else axes
            
            # Activate the axis
            plt.sca(ax)

            # Plot using gnomonic projection directly into the provided axis
            hp.gnomview(map[ds_index], rot=[long_full[i], lat_full[i]], xsize=600, ysize=600, title=unknown_feature_title, nest=True, unit='dEBV', 
                        hold=True)
            if overplot == True:
                for plot in overplot_list:
                    if one_dist == False:
                        make_scatter_plot(plot[ds_index])
                    else:
                        make_scatter_plot(plot)
            plt.title(unknown_feature_title, fontsize=16)
            cbar = plt.gcf().axes[-1]
            cbar.tick_params(labelsize=15)
        
        # Save the entire figure for the current distance slice
        fig.suptitle(title + 'Distance ' + '{:.2f}'.format(distslices[ds_index]) + ' kpc', fontsize=20)
        fig.tight_layout(rect=[0, 0.03, 1, 0.95])
        fig.savefig(dist_path + save_name + f'Distance{ds_index}.png', bbox_inches='tight', pad_inches=0.1)
        plt.close(fig)



