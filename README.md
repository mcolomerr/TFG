# Analysis of NMII spatial distribution

In this site you can find all the files I have coded for my Bachelor Thesis titled "Mechanosensitive single molecule dynamics of myosin II" that I have carried out at Lab Wieser (ICFO). 
This repository contains exemplary configuration files and custom-written scripts for the analysis of data. It consists on the following:
- intensity_macro.ijm: batch-processing macro for the program Fiji to track the intensity of a ROI over all the frames.
- myosin_custom.ipynb: code for analysing the tracks created using the custom-written Fiji macro.
- myosin_trackmate.ipynb: code for analysing the tracks created using the Trackmate plugin from Fiji.
- myosin_DBSCAN.ipynb: code for analysing the clustering of the molecules.
- myosin_cytosim.ipynb: code for analysing the text files created using the program Cytosim. It analyses the scatter distribution of myosin in the cell over time and extracts the expansion rate of the contractile simulated cells.
- config.cym: example of a configuration file created for the Cytosim. Combined with the use of batch-processing scripts, several simulations can be performed at the same time.

