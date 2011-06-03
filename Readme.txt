README
Alistair's Open Image Analysis Repository
http://alistairboettiger.info 


This repository contains my open source Matlab codes for processing and analyzing molecular imaging data.  I work primarily with whole mount blastoderm stage Drosophila embryos, though many of the scripts here should work equally well on different types of data.

My image processing toolbox is primarily organized through a series of GUIs that perform specific functions.  These rely on a common set of custom codes.  

The primary analysis scripts used to count single cytoplasmic constructs as presented at the Annual Drosophila Meeting on March 2nd 2011 are still command-line based (see below).  

Graphical User Interfaces:

Descriptions

im_quant_sna_v3.fig
Description: Begun on 09/25/08, this is the latest version of code to quantify the width of the snail expression pattern of embryos in cc14.  
Algorithm:  A multicolor image of a Drosophila embryo stained for snail is thresholded, the intensity of staining is measured in transects across the pattern, from which width and sharpness measurements can be made.  The nuclei are then segmented, neighbors identified, and the shortest path from one side of the expression domain to the other (through neighboring nuclei) is computed.  This data is exported in a target folder as: 
save([folder,'/',savename,'_data.mat'],'sna_grad','nuc_map','tops',...
      'bots','norm_grads','nuc_tot','cnt_nuc','flp'); 
This is respectively: the raw gradients spanning the expression region, the nuclear map, the indices of nuclei along the top of the pattern, indices along the bottom, length normalized gradients, total nuclei, the counts of nuclei spanning the image, and any horizontal or vertical flips.  


im_singlemolecule.fig
Description: Begun 01/21/11.  This script counts the number of mRNA molecules in every cell from 3D confocal stacks.  Default is Zeiss .lsm formatted inputs.  
Updates: 04/04/11 This script is currently out of date.  The new CheckDotUpDown.m provides a watershed implementation for separating vertically stacked dots that is more precise.  Due to memory issues and processing speed issues the GUI version has been replaced for the moment with UnsupervisedDotFinding.m.   
Algorithm: Step zero loads the data, step 1 assigns channels and max-projects the nuclear channel. Step 2 segments the nuclei, step 3 segments the cells (assigns all space to nearest nuclei), step 4 counts all dots using difference-of-Gaussians, watershedding, and size and min intensity filters.  step 5 compares neighboring levels, removes any dot in current level already found in previous.  step 6 repeats step 4 for chn 2, step 7 repeats step 5 for channel 2, step 8 exports data. 
save([fout,fname],...
        'mRNA_cnt1','mRNA_den1','mRNA_ind1','mRNA_sadj1','DotData1',...
        'mRNA_cnt2','mRNA_den2','mRNA_ind2','mRNA_sadj2','DotData2',...
        'NucLabeled','nuc_cents','nuc_area','In','conn_map','Cell_bnd'); 
This is the raw counts, the densities, the cell size adjusted counts, and the raw dot data centroids (before redundancy check).  Then the labeled matrix for the nuclei, the centroids, area per cell, the projected nuclear map, the connectivity map, and the cell boundary map.  

imviewer_lsm.fig
Description: Begun 02/12/11.  This file reads lsm 5 format (.lsm) files from Zeiss software and saves the data in the Z-stacks as independent .tif files appended as _z#.tif in a target folder.  A max-project image is also saved.  An autocycle function will loop through multiple positions from a multi-position acquisition and automatically save all.  
Algorithm: Two simple steps, one to read in the data, one to max project and write it out.  Max-project is done sequentially to save memory and not need to load all layers into memory simultaneously.  


mRNA Counting 
(as presented at Annual Drosophila Conference, March 2nd 2011).
Data analysis is currently done in 3 parts:
imviewer_lsm.fig: Convert .lsm files to separate tifs and a max projection
imnuc_seg.fig: ID and segment nuclei (this is the only supervised part of this code).  
Unsupervised_DotFinding.m - Find 3d centroids of all puncta (mRNA molecules/diffraction limited spheres) present in image stack.  
Primary methods are in two separate functions.  (1) dotfinder.m, which implements the 2D Gaussian filtering, max-object segmentation algorithm, and watershed separation for fused dots.  (2) CheckDotsUpDown.m which groups the centroids of dots based on an overlap criteria and uses a watershed algorithm to separate stacked dots.  A variety of tricks are necessary to efficiently handle the association of order a million positions with each other.  Data is saved in a target directory as _slidedata.mat.  

 

Installation:
In addition to simply copying the files to a common directory and firing them up from Matlab, there are a few quick setup steps you may need to follow:

If you are running on a Windows, the file separator will need to be changed to backslash '\' instead of forwardslash '/'.  A find all change all should suffice.  
If you haven't copied the .mat files  called _pars#.mat for the respective GUI, find the files called _pars#.mat (where # is the step number of the analysis script).  Execute the commented out pars = {...} save pars command.  *Why it's written this way. 



*All these GUIs write and read the current set of parameters from a .mat file on the hard-drive.  This allows you to find an optimal parameter set for processing your data and have that easily become the program default.


05/31/11
Added AlistairMobile to repository.

06/02/11
Jacques co-editor of repository


