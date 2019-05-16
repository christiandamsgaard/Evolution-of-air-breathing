# Evolution-of-air-breathing
Data set and R code for running stochastic character mapping to reconstruct the evolution of air-breathing and amphibious behaviour in fishes. 

## Before running the scripts, you should do the following

1. Make an new folder that will be the working directory and specify its path here
setwd("/Users/christiandamsgaard/Dropbox/Projects/Air-breathing review_full/")                    

2. Download, unzip and paste "Complete “all-taxon assembled” chronograms" and "PFC_short_classification.csv" from "https://fishtreeoflife.org/downloads/" into the working directory. 

3. Download and paste "Data Set 1.xlsx" into the working directory. 

4. Within the working directory folder make a folder called "Figures".

5. Install the packages specified under ## LOAD PACKAGES ## using the "Tools" --> "Install Packages..." option. 

To prepare running the model, first run the lines from the sections #### LOAD PACKAGES, COLOR SCHEMES, FUNCTIONS, IMPORT TREES, IMPORT DATA #### 
From there, you can run modelling within section 1-2 on air-breathing evolution or section 3-4 on amphibious evolution.
Sections 1-2 and sections 3-4 are followed by lines that code the figures. 
Sections 1-2 and sections 3-4 each take up ~ 10GB ram, so if your computer has < 16GB memory,
consider running Sections 1-2, plot the associated data, and shut down R to clear the memory. 
Then re-open R, and run sections 3-4, plot the associated data. 


