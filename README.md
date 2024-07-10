# CITA_StarFormationRegions
Repository containing all work done in the SURF program at CITA. 

**Goal**: Creating an algorithm that can create and identify star formation regions in a 3D Dust Temperature Map of the Milky Way

Notes:

- DataFiles contains all the csv files required in PlottingSFtracers.ipynb. Store these files in csv_files_path from settings so that it can be accessed. 

- dustmap_cropped.png is required in comparison done in 2D_RGBmap.ipynb, but this is not necessary to actually make the images. 

Links for data files can be obtain upon request from Ioana Zelko. Ensure that they are stored in data_files_path from settings. 

### Scripts
Scripts containing functions and variables used in the notebooks in this repository

#### settings_template.py:
Initializes paths to read/store images, data and csv files. Use can change this and save it as "settings.py" for their own personal use.

#### constants.py:
Contains important constants. These include c, h (Planck constant) and k (Stefan-Boltzmann constant) for Wien's Law and finding the Planck function, as well as the longitude (l), latitude (b) and distance (d) of some local molecular clouds LMCs. 

No package/script requirements

#### loading_data.py:
Loads dust map data.

Packages Required:
- h5py

#### calculations.py:
Contains functions used to perform calculations on data. This includes manipulating the arrays (increasing resolution, multiplying, normalizing, extracting parts), finding RGB color channels, and performing calculations like Wien's Law. 

Packages/Scripts Required:
- numpy
- PyEphem
- healpy
- matplotlib.pyplot
- **constants.py**

#### visualizing_funcs.py:
Contains functions used to visualize data (either through matplotlib or images) or manipulate visualizations (i.e change brighteness of an image)

Packages/Scripts Required:
- numpy
- PyEphem
- healpy 
- matplotlib (pyplot and gridspec)
- os
- Image and ImageOps from PIL

### PlottingSFtracer.ipynb:
Jupyter notebook used to plot star formation tracers (YSOs and SFCs) over the temperature maps. It uses the csv files found in the folder "DataFiles". 

Packages/Scripts Required:
- numpy
- healpy
- h5py
- matplotlib (pyplot)
- pandas
- **visualizing_func.py**, **calculations.py**, **loading_data.py**, and **settings.py**


It produces a map of all the sf tracers. It then assigns a distance slice to any with distance and plots them over the temperature map at that distance. 

### 2D_RGBmap.ipynb
Jupyter notebook that tries to create and RGB image of the 2D dust and temperature maps. 

Packages/Scripts Required:
- numpy
- healpy
- h5py
- matplotlib (pyplot and gridspec)
- os
- plotly.express (optional, only for interactive plot)
- PIL (import Image and ImageOps)
- IPython.display (import display) *Only necessary because of how my computer outputs image so may not be required
- **visualizing_func.py**, **calculations.py**, **loading_data.py**, and **settings.py**

The notebook creates a new temperature only map with a higher resolution, then creates a planck function at 3 different frequencies (that will correspond to R, G and B). The planck function requires the script "constants.py". Then the new temperature emissivity (planck function) maps are multiplied by the dust density dEBV maps. From this, the R, G and B channels are extracted and used to create images for the whole sky, the Cepheus molecular cloud and the Orion A and B molecular clouds. 

### flagging_properties.ipynb

Notebook that manually flags properties such as: Regions with high density, cold temperature (high in R channel), hot temperature (high in B channel), a mix of both hot and cold regions, and temperature gradient. It produces panels to provide a visualization of how well the flagging is working. 

Packages/Scripts Required:
- numpy
- os
- healpy
- matplotlib.pyplot
- from PIL, Image
- **constants.py**, **loading_data.py**, **calculations.py**, **visualizing_functions.py**, and **settings.py**

