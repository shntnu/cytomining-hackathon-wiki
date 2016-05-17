#Input files:
DMSO-normalized variance-feature-extracted single-cell data for 103 samples. Data for each sample is stored in tab-delimit txt format. The first row is the names of 29 selected features. Starting from the second row, every row is the data for one individual cell. In this example dataset, file names contains sample information: the mechanism of action (MOA) and compound/concentration for each sample. 


#Matlab functions needed:
1. SPADE v3 package developed by Peng Qiu (http://pengqiu.gatech.edu/software/SPADE/) 
1. Fast bh_tSNE implementation by Laurens van der Maaten (https://lvdmaaten.github.io/tsne/)

#Matlab script for visualization pipeline:
pipeline_single_cell_vis_to_sample_level.m

1. Read the single-cell data for each sample, and pool them together. (optional: if the total number of cells is too large, we can downsample the data to reduce the number of cells, but we also need to upsample in all subsequent visualizations. How to upsample will be different for each visualization algorithm. Examples are not provided here.)
2. PCA visualization of single cells. 
<img src="https://cloud.githubusercontent.com/assets/18299367/15181524/b785f0be-1755-11e6-9703-a59a4c51e47b.png" width="300">

Visualization of single cells (PCA, tSNE, SPADE). Input of this step is the pooled data matrix. The output is a visualization (map) that visualizes either the data points/cells (PCA and tSNE), or a tree representation of the structure/skeleton of the data (SPADE). (2.1) We can color the visualization by a particular marker. The resulting colored visualization tell us which part of the visualization is high for the marker, which part is low for the marker. Doing this for all the markers one-by-one will show us which part of the visualization is positive for what features, and help us to understand the morphologies corresponding to various parts of the visualization. (2.2) We can color the visualization by distribution of cells in a particular sample. The resulting colored visualization tell us which part of the visualization is occupied by cells in the sample, or in other words, what morphologies are present in the sample with what probability.

Sample distance matrix based on SPADE. The SPADE tree and visualization (2.2) gives us a cell distribution on the visualization for each sample. After obtain the cell distributions of the samples, we can compute pairwise distance of distributions of each pair of samples. The distance metric can be Euclidean, Correlation Distance, or Earth Mover's Distance which takes the SPADE tree structure into consideration.

Visualization of sample similarities. The sample distance matrix can be visualized by heatmap and hierarchical clustering, to reveal clustering pattern of the samples. Alternatively, the distance matrix can be visualized by Minimum Spanning Tree, leading to a SPADE-like visualizatioin where we can observe both clustering pattern of the samples and similarity among the clusters.

