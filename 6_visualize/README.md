# Visualization

@Harmanjit
@holgerhennig
@jnasser3
@prees
@pqiu

Visualization, clustering and classification analysis is typically performed after normalization and feature selection of single-cell level data. 

Data input: cells vs. features matrix for each sample, stored in CSV files or tab-delimit txt files. Standardized file names would be preferred.

Typical data size for a mid-sized primary compound screen: 1000 cells with 20-30 selected features, 10,000 compounds. Makes 10^7 x 30 data points

Preprocessing for data visualization:

1. Pooling / concentration (optional: downsample the data if the total number of cells is too large).  Output of this step is a pooled data matrix. It is a tall thin data matrix, with many many rows corresponding to cells across all samples, and relatively small number of columns corresponding to the features. 

2. Visualization of single cells (PCA, tSNE, SPADE). Input of this step is the pooled data matrix. The output is a visualization (map) that visualizes either the data points/cells (PCA and tSNE), or a tree representation of the structure/skeleton of the data (SPADE). 
(2.1) We can color the visualization by a particular marker. The resulting colored visualization tell us which part of the visualization is high for the marker, which part is low for the marker.  Doing this for all the markers one-by-one will show us which part of the visualization is positive for what features, and help us to understand the morphologies corresponding to various parts of the visualization. 
(2.2) We can color the visualization by distribution of cells in a particular sample. The resulting colored visualization tell us which part of the visualization is occupied by cells in the sample, or in other words, what morphologies are present in the sample with what probability. 

3. Sample distance matrix based on SPADE. The SPADE tree and visualization (2.2) gives us a cell distribution on the visualization for each sample. After obtain the cell distributions of the samples, we can compute pairwise distance of distributions of each pair of samples. The distance metric can be Euclidean, Correlation Distance, or Earth Mover's Distance which takes the SPADE tree structure into consideration. 

4. Visualization of sample similarities. The sample distance matrix can be visualized by heatmap and hierarchical clustering, to reveal clustering pattern of the samples. Alternatively, the distance matrix can be visualized by Minimum Spanning Tree, leading to a SPADE-like visualizatioin where we can observe both clustering pattern of the samples and similarity among the clusters. 
