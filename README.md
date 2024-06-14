# CITA_StarFormationRegions
Repository containing all work done in the SURF program at CITA. 



### PlottingSFtracer.ipynb:
Jupyter notebook used to plot star formation tracers (YSOs and SFCs) over the temperature maps. It uses the csv files founnd in the folder "DataFiles". 

Packages Required:
    - numpy
    - healpy
    - h5py
    - matplotlib (pyplot)
    - pandas

It produces a map of all the sf tracers. It then assigns a distance slice to any with distance and plots them over the temperature map at that distance. 

### 2D_RGBmap.ipynb
Jupyter notebook that tries to create and RGB image of the 2D dust and temperature maps. 

Packages required:
    - numpy
    - healpy
    - h5py
    - matplotlib (pyplot and gridspec)
    - os
    - plotly.express (only for interactive plot)
    - PIL (import Image and ImageOps)
    - IPython.display (import display) *Only necessary because of how my computer outputs image so may not be required

The notebook creates a new temperature only map with a higher resolution, then creates a planck function at 3 different frequencies (that will correspond to R, G and B). The planck function requires the script "constant.py". Then the new temperature emissivity (planck function) maps are multiplied by the dust density dEBV maps. From this, the R, G and B channels are extracted and used to create images for the whole sky, the Cepheus molecular cloud and the Orion A and B molecular clouds. 

