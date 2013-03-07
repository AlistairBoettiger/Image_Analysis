## Introduction

This folder contains code used in my publication:

> Boettiger, A. N. & Levine, M. Rapid Transcription Fosters Coordinate snail Expression in the Drosophila Embryo. Cell reports 3, 8–15 (2013).

These files were used to process the mRNA counting results produced by the code described in the `ImageAnalysis` folder.  They depend upon raw data files produced by that code to run successfully.  These `.mat` files are too large to include in the repository, but you may write to me to request files.

alistair.boettiger@gmail.com

## Files
* MeasureSnaGradient.m - orients snail gradients along a common axis and computes spatial averages for later plotting. 
* MakeFigure1.m - plots the images in Figure 1
* MakeFigure2.m - plots the data in Figure 2
* MakeFigure3.m  - plots the data in Figure 3
* MakeFigure4.m - plots the data in Figure 4

## Functions
* plotstackdots.m - function required by MakeFigure1.m
* Functions in the `ImageAnalysis` directory are also called by scripts and functiosn in this directory, so you will want both in you filepath.




*Note: This readme is written in markdown.  It should be text-readable as well, but will look better if you view it with a markdown compiler.  Github automatically compiles and renders markdown.  
