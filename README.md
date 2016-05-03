# MATLAB
A repository for all MATLAB code

This repository will hold the bare minimum matlab code to run the data from calcium imaging experiments

I am tentatively planning on splitting it into 

Workflow - This will consist of the minimal files necessary to go from movies to rasters
Pattern Recognition - The code for identifying synfire chains
Correlation Analysis - Correlation functions and surrogate generation
Network Analysis - Things such as cc and surrogate generation
Paper figures - Scripts used to generate certain figures

This code is written with the assumption of mapping the sohal server as a network drive on windows with letter Z

Functions are largely adapted from...

> <cite> Mukamel et al. 2009. Automated analysis of cellular signals from large-scale calcium imaging data. <cite>

Only major changes are adapting function to deal with single tiffs as tiff stacks are limited to 4GB files due to 32bit indexing. Future releases might make use of either HDF5 or BigTIFF formats.

* * Code written by Francisco Luongo between 2012 and 2015 in the lab of Vikaas Sohal at UCSF
