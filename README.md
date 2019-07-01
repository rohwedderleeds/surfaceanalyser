# surfaceanalyser
ImageJ plugin for analysing object surfaces
Surface and Direction Detector
Copy the files into the analysis folder within the plugin the pluginsfolder of ImageJ and restart imageJ activate the plugin.

The ImageJ plugin takes the current image/image stack and generates an internal copy for the analysis. In case of an image stack each slice will be evaluated individually. In a first step a black and white image is 

generated using auto-threshold. To remove small single pixel ("noise") a rank filter is applied to the image using median and bright_outliers filter. The function ParticleAnalyser then is used for isolating the 

cells/cell cluster. This is achieved with the options EXCLUDE_EDGE_PARTICLES, requiring objects not touching the frame, INCLUDE_HOLES, resulting in the detection of the outer objects border only, CENTER_OF_MASS for 

analysing of lateral shifts. Size limits are calculated in percentage from the minimum coverage value from the menu and a reduced width and height value for maximum.
The particle analyzer will try to find an object based on those values first, when it fails it starts reducing the minimum size value until an object is identified. When an object is found, the resulting ROI is 

stored for analysis as a linear array. To calculate the excentricity of the cell is calculated using a circular outer perimeter vs. the maximal inner perimeter using averaging of the ROI values. Differences are 

printed out in a log window.
The inner centre is used to calculate the radius length (yn) of each point at the ROI to the centre.
