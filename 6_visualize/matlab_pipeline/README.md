#Input files:
DMSO-normalized variance-feature-extracted single-cell data for 103 samples. Data for each sample is stored in tab-delimit txt format. The first row is the names of 29 selected features. Starting from the second row, every row is the data for one individual cell. In this example dataset, file names contains sample information: the mechanism of action (MOA) and compound/concentration for each sample. 


#Matlab functions needed:
1. SPADE v3 package developed by Peng Qiu (http://pengqiu.gatech.edu/software/SPADE/) 
1. Fast bh_tSNE implementation by Laurens van der Maaten (https://lvdmaaten.github.io/tsne/)

#Matlab script for visualization pipeline:
(pipeline_single_cell_vis_to_sample_level.m)

#####Read and pool single-cell data for all samples
Read data for each sample. Pool them together to form a "union" sample that contain all cells from all samples. (Optional: if the total number of cells is too large, we can downsample the data to reduce the number of cells, but we also need to upsample in all subsequent visualizations. How to upsample will be different for each visualization algorithm. Examples are not provided here.)

#####PCA visualization of single cells 
Take the "union" sample and perform PCA to reduce this 29-dimensional dataset to 2D. The resulting 2D data can be visualized by the scatter plot and contour plot below. The two axes correspond to the two principle components. In the scatter plot, each dot represents one cell. The contour plot shows the density of cells in this 2D visualization.
<img src="https://cloud.githubusercontent.com/assets/18299367/15306752/5e8fbaae-1b9b-11e6-8e1a-f8e4a1aa8283.png" width="700" height="300">

The PCA visualization of the "union" sample can be used as a background landscape to visualize each individual sample. Below shows 6 selected samples out of the total of 103 samples. In each sub-figure, the green dots represent all cells from all samples, whereas the blue dots correspond to celles in one particular sample. These 6 samples belong two MOAs. The first three are Aurora Kinase inhibitors (Aur), and the last three are labeled as general Kinase inhibitors (Ki). The difference between the two MOAs is not clearly observed. 
<img src="https://cloud.githubusercontent.com/assets/18299367/15306756/61f952ea-1b9b-11e6-85d2-6aff1a9c0b73.png" width="1000"> 


#####tSNE visualization of single cells 
tSNE can be applied to reduce the 29-dimensional "union" sample into 2D. The resulting 2D data is visualized by scatter plot and contour plot. Similar to PCA, the tSNE visualization appears to be one smear. No clear subpopulations are observed. 
<img src="https://cloud.githubusercontent.com/assets/18299367/15306757/65d47a34-1b9b-11e6-91f2-18e78c571e12.png" width="700" height="300">

However, the tSNE visualization is helpful in separating sampels with different MOAs. Below shows the same set of 6 selected samples, using the tSNE visualization as the background landscape. We can clearly observe that the first 3 samples (Aurora Kinase inhibitors) are similar to each other, the last 3 samples (Kinase inhibitors) are similar to each other, and the two groups are different. This result indicates that tSNE is more informative than PCA in this particular dataset. 
<img src="https://cloud.githubusercontent.com/assets/18299367/15306760/68c8c8e4-1b9b-11e6-84c9-74eabd22f59c.png" width="1000"> 


#####SPADE visualization of single cells 
SPADE can also be applied to visualize the single-cell data. SPADE views the 29-dimensional "union" sample as a 29-dimensional point cloud, and aims to derive a tree representation to approximate the skeleton of the point cloud. In the SPADE tree, each node is a cluster of cells, and the edges form a minimun spanning tree to connect the clusters. If we color the tree according to the 29 features one at a time, we will observe that different parts/branches of tree correspond to different morphological phenotype of cells (not shown here). The structure of the tree is shown below. 
<img src="https://cloud.githubusercontent.com/assets/18299367/15309058/430dca1c-1bb1-11e6-825e-a7d7d094a7e5.png" width="600"> 

Similar as before, the SPADE visualizaiton can be used as the backgroud landscape to visualize individual samples. In each of the 6 sub-figures below, one sample is visualized, and nodes are colored by the precentage/frequence of cells in the sample belong to each node. Basically, each sub-figure shows which part of the SPADE tree is occupied by cells in a particular sample. Again, we can clearly observed that the 6 samples form two groups (MOAs). 
<img src="https://cloud.githubusercontent.com/assets/18299367/15306765/703bbc76-1b9b-11e6-8d06-b1164693855d.png" width="1000"> 


#####Heatmap and MST based on SPADE distributions and Correlation distance
The SPADE visualization of individual samples can be viewed as probability distributions. For each sample, the node colors form a probability distribution of cells on the SPADE tree nodes. Given 103 distributions for the 103 samples in this dataset, we can derive a pairwise sample-to-sample distance metric to visualize the similarity and differences among the samples. 
Pairwise correlation distance is computed and visualizaed by hierarchial clustering heatmap. The colorbar below the heatmap 
<img src="https://cloud.githubusercontent.com/assets/18299367/15306766/73ad5aae-1b9b-11e6-883f-658dea425f56.png" width="1200"> 


visualization. (2.2) We can color the visualization by distribution of cells in a particular sample. The resulting colored visualization tell us which part of the visualization is occupied by cells in the sample, or in other words, what morphologies are present in the sample with what probability.

Sample distance matrix based on SPADE. The SPADE tree and visualization (2.2) gives us a cell distribution on the visualization for each sample. After obtain the cell distributions of the samples, we can compute pairwise distance of distributions of each pair of samples. The distance metric can be Euclidean, Correlation Distance, or Earth Mover's Distance which takes the SPADE tree structure into consideration.

Visualization of sample similarities. The sample distance matrix can be visualized by heatmap and hierarchical clustering, to reveal clustering pattern of the samples. Alternatively, the distance matrix can be visualized by Minimum Spanning Tree, leading to a SPADE-like visualizatioin where we can observe both clustering pattern of the samples and similarity among the clusters.

