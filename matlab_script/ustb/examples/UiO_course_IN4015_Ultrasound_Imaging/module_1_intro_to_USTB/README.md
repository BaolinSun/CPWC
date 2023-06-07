# Module 1 Introduction to the USTB: Getting familiar with the USTB

The lab assignments for module 1 consists of using existing USTB examples to set ut the tools necessary for the rest of the exercises.

+ Get to know the UltraSound ToolBox by running and getting familiar with multiple examples 
	+ Pure USTB examples. Run at least five of these
		+ examples/uff/CPWC_UFF_Alpinion.m
		+ examples/uff/CPWC_UFF_Verasonics.m
		+ examples/uff/FI_UFF_phased_array.m
		+ examples/picmus/experimental_contrast_speckle.m
		+ examples/picmus/experimental_resolution_distortion.m
		+ examples/picmus/carotid_cross.m
		+ examples/picmus/carotid_long.m
		+ examples/acoustical_radiation_force_imaging/ARFI_UFF_Verasonics.m
        + examples/UiO_course_IN4015_Ultrasound_Imaging/module_1_intro_to_USTB/minimal_example.m
        + examples/UiO_course_IN4015_Ultrasound_Imaging/module_1_intro_to_USTB/maximal_example.m
	+ USTB + K-wave examples. You need to install and add the ultrasound simulator k-wave (http://www.k-wave.org/) to your MATLAB path. See simple description below.
		+ examples/kWave/CPWC_linear_array_cyst.m
	+ USTB + Field II examples. You need to install and add ultrasound simulator Field II (https://field-ii.dk/) to your MATLAB path. See simple description below.
		+ examples/field_II/STAI_L11_resolution_phantom.m
		
## How to download k-wave and Field II

### K-wave

Register a new profile (http://www.k-wave.org/forum/register.php) with the necessary information, go to download, and download the k-wave zip-file.
Next, extract the folder on your preferred spot, and add to path (see below).

### Field II

Go to the download site at (https://field-ii.dk/), and download the file for your operating system. Extract the folder (twice) at your preferred location.
Next, add to path (see below).

## How to add to path

There are multiple ways to add to path:

1. Create a new file called "startup.m" in your working directory, and add the path in this file by writing
    + addpath *filepath, i.e. Documents/MATLAB/k-wave*
2. Pathtool: open pathtool in Matlab by writing pathtool in the Command Window in Matlab. Click "Add Folder" to add k-wave and field-ii to path, and press save. To use this method, you have to have Admin rights on the computer you are using.
3. Use the command "addpath()" and possibly "genpath()" in the MATLAB command window.