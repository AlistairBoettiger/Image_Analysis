# MatlabPhotoshop

Alistair Boettiger

boettiger.alistair@gmail.com

## Introduction
The goal of this project is to duplicate some of the common functions I use in Photoshop into matlab, in a way that allows the data to also become readily available for further processing in Matlab.

In addition this interface and its child functions are intended to allow for easy batch processing for some numerous common image tasks I have.

### Notes
* As written, some functions in here may require `CheckParameter.m` and other functions from my `matlab-storm` directory. 
* This is currently a work in progress directory.  Some updates may cause new bugs.  Not all functions/menu options may yet be programmed and active.  

## GUI
* `ImageShop.fig`
* `ImageShop.m`

## Functions

###   Read/view/exctract LSM files 
(from Zeiss Zen Black 2012 software format)  

*  `LSMextract.m` - Takes a .lsm file, extracts all z-positions and all stage-positions into multicolor images, and also saves a max-projection image for the Z-stack.  Optionally allows user to specify a sub-range of Z for each channel and an aribtrary subset of the positions.  
* `LSMread.m` - creates a Datas structure from an LSM file, needed for other .lsm processing functions.  
*  `readLSMmatfile.m` - Extract a matlab image file from a .lsm file using info from the Datas structure.  This replaces `read_planeT.m`, much, much more cleaned up, removed a ton of unnecessary code.    Not backwards compatible as we dramatically simplified the input and output options.  
* `LSM2image.m` - Takes a Datas file or a filename.lsm, returns an image for a specified z-position and stage-position (defaults to chose the first of each).  Called by LSMextract in saving the images as RGB tifs and in building the max projections.  



