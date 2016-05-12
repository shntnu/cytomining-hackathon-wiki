# Visualization

@Harmanjit
@holgerhennig
@jnasser3
@prees
@pqiu


Visualization, clustering and classification analysis is typically performed after normalization and feature selection of single-cell level data. 

Data input: cells vs. features matrix for each sample

Typical data size for a mid-sized primary compound screen: 1000 cells with 20-30 selected features, 10,000 compounds. Makes 10^7 x 30 data points

Preprocessing for data visualization 
1. Pooling / concentration 
2. Visualization of single cells (PCA, tSNE, SPADE)
3. Distance metric based on SPADE distribution (L2, Corr, EMD)
4. Visualization of samples heatmap/hierarchical clustering, MST
