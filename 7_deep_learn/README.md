# Deep learning 

@blakeweb
@jccaicedo
@HackerMD
@mvictor212
@okraus

Here we provide a very preliminary example Deep Learning workflow for cell profiling. To download the data, run

`wget https://www.dropbox.com/s/lhf6f9v22v69xdw/az_dataset.h5?dl=1`.

The dataset is about 10GB, and contains the 16-bit integer image arrays, indexed by compound, concentration, replicate, site, and channel.

The `load.py` module facilitates loading up minibatches of the data, based on a passed metadata frame. The full metadata
is available here as `metadata.csv`.

This is by no means definitive, but rather an early prototype to test various DL methods on microscopy data.

There is also an AMI publicly available that already has this dataset, the CUDA installation complete, as well as the
necessary python dependencies: `ami-cafd01aa`

Launching a GPU EC2 instance from that image should make it nearly trivial to get started with deep learning on this
morphological profiling task.
