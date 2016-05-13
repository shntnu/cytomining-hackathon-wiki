# Visualization

@Harmanjit
@holgerhennig
@jnasser3
@prees
@pqiu



Visualization, clustering and classification analysis is typically performed after normalization and feature selection of single-cell level data. 

Data input: cells vs. features matrix for each sample

Typical data size for a mid-sized primary compound screen: 1000 cells with 20-30 selected features, 10,000 compounds. Makes 10^7 x 30 data points

Preprocessing for data visualization:

1. Pooling / concentration (optional: downsample the data if the total number of cells is too large).  Output of this step is a pooled data matrix. It is a tall thin data matrix, with many many rows corresponding to cells across all samples, and relatively small number of columns corresponding to the features. 

2. Visualization of single cells (PCA, tSNE, SPADE). Input of this step is the pooled data matrix. The output is a visualization (map) that visualizes either the data points/cells (PCA and tSNE), or a tree representation of the structure/skeleton of the data (SPADE). 

3. Distance metric based on SPADE distribution.  (L2, Corr, EMD)

4. Visualization of samples heatmap/hierarchical clustering, MST
