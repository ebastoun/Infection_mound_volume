# Infection_mound_volume

This folder contains: (1) a representative confocal stack of host MDCK nuclei (stained with DAPI) centered around a *Listeria monocytogenes*-infection mound; (2) the code to convert the raw data into a stack of binary images; and (3) the resulting binarized mound obtained as a result of running the code. It also contains a background image that is subtracted from the nuclei image during the process of binarization and a flatfield image that is used to correct for uneven illumination. The second section of the matlab file contains the code to convert the binarized image into alpha shapes and to obtain a mound volume by adding the area of the alpha shape for each slice and multiplying by the increment between the slices.


The file titled **1control_1.tif** contains the image stack obtained from the microscope, which can be read into the matlab file titled **mound_volume_calculation.m** to output the file titled **tile1.TIFF**, which represents the binarization of the infection mound. Because of its size we were not able to upload it here but one can download it from: https://gitlab.com/theriot_lab/theriot-toolbox/infection-mound-volume. The background image titled **MED_60X_1X_background_3_MMStack_Pos0_1.ome.tif** and **flatfield image titled MED_Concatenated Stacks.tif** are also read into the matlab file titled mound_volume_calculation to produce the binary image stack. The correction for the brightness of the flatfield image must be determined empirically by the user to match the specifications of their own data. Also, we account for the decrease in brightness of the nuclei at higher values of z by periodically altering the extent of background subtracted (analogous to applying a manual threshold). The binary image stack is then read into the second section of the matlab file titled mound_volume_calculation to output a volume calculation. 

________________________________________________________________________________________________
In case you need help on performing background and flatfield correction on the images of the host cell nuclei we have also added below detailed instructions on how to perform that.

**Background Subtraction:** 

Background subtraction removes uneven illumination of the background, which aids in segmentation. 

1.	Prepare a dummy coverslip to be imaged for background subtraction. This is a coverslip that has been treated similar to an experimental coverslip, but which has no foreground objects (in this case, cell nuclei). 
2.	The background image is specific to the objective, optovar, fluorescence intensity, exposure length, camera gain, and excitation wavelength. If any of these parameters are altered while acquiring images of the experimental sample, a new background image must be obtained. 
3.	Find the focal plane of the dummy coverslip using a confocal microscope and acquire an image stack spanning the average height of an infection mound, by specifying the same number of z slices as would be used to image such a mound. Maintain the same interval between z slices as is used to image an infection mound. 
4.	Acquire 10-100 such stacks at multiple positions. 
5.	Create a maximum intensity projection of the stack at each position and generate a median projection across all positions. 
6.	Remove extremely bright pixels by passing a median filter through the image. 

**Flatfield Correction:**

Flatfield correction corrects for uneven illumination of objects in the foreground, such as nuclei at the center of the field of view appearing brighter than nuclei at the corners, to improve accuracy of segmentation. 

1.	To correct for dapi-stained nuclei, first prepare a concentrated dye slide by adding 5μL of 50mg/mL coumarin dye stock solution to a Micro90-cleaned slide. Cover with a clean coverslip and seal using nail polish. Appropriate dyes for excitation wavelengths other than 405 nm include sodium fluorescein for 488 nm, rose bengal for 561 nm, and acid blue for 642 nm light. 
2.	Find the focal plane of the concentrated dye slide using a confocal microscope with the Perfect Focus feature or with the aid of bubbles in the concentrated dye slide. Acquire an image stack by imaging both above and below this focal plane and image sufficient z slices such that the fluorescence light emitted by dye is no longer visible at the most extreme z slices. 
3.	Obtain similar image stacks at various positions on the concentrated dye slide.
4.	Acquire a dark image to account for noise in the camera as well as dust along the optical path. To do so, divert light to the eyepiece or close the camera shutter. Then, acquire 100 images in a time series using the same intensity and exposure settings that was used to obtain the flatfield images. 
5.	Perform a median projection of the dark images. 
6.	Pass each slice of the flatfield images through a median filter to remove very bright pixels. 
7.	Subtract the dark image from each slice of the flatfield images. 
8.	Perform a maximum intensity projection of each dark image-subtracted flatfield image stack. 
9.	Generate a median projection across all positions. 

