# Overview
[] 

[] If you do 


`haze-video-dataset` contains code used to generate the haze video dataset. In order to use the dehaze evaluator, instructions in this folder must be followed to download the haze video dataset and raw KITTI data.

`dehaze-systems` contains the code for multiple dehazing methods to be evaluated.

`dehaze-evaluation` contains code to evaluate the dehazing systems.


For more information see the README files inside each folder 

[] this will automatically perform the following steps -
1) Download the required subset of the raw KITTI dataset into the `KITTI` folder
2) Clone the necessary repositories for generating sky masks and depth maps ([SkyAR](https://github.com/jiupinjia/SkyAR) and [ManyDepth](https://github.com/nianticlabs/manydepth) respectively)
3) Generate sky masks, depth maps and hazy sequences.
4) Run the evaluation code on the present dehazing systems to compare their dehazing performance on the newly generated hazy sequences.

