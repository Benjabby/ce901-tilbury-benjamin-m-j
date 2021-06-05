# Haze Video Dataset

This folder contains code []

The "Haze Video Dataset" is contained at [??? Need somewhere to upload this to. Is very large. ???]. It comprises of x sequences of hazy road scenes and corresponding transmission maps for a total of y colour hazy images and x grayscale transmission maps.
Each sequence is derived from a sequence of the KITTI dataset. Only a subset of sequences from the full KITTI dataset is used. For simplicity the dataset follows the same folder structure as the KITTI dataset.
Each sequence has a different, randomly chosen fog density and colour, this information is stored in the `haze/props.json` file.

The dataset can be downloaded from [??? somewhere ???] and should be unziped into `dataset` within this folder.
```
dataset
--- [drive 1]
------ [drive 1 sequence 1]
--------- haze
------------ props.json
------------ image
--------------- 000000000.png
--------------- 000000001.png
--------------- .............
------------ transmission
--------------- 000000000.png
--------------- 000000001.png
--------------- .............
------ [drive 1 sequence 2]
------ [drive 1 sequence 3]
------ ....................
--- [drive 2]
--- [drive 3]
--- ........
```

In addition to the haze video dataset, the raw KITTI data is needed to run the evaluation code, or to generate further hazy sequences.

The raw KITTI data can be downloaded from their [website](http://www.cvlibs.net/datasets/kitti/raw_data.php). As the KITTI folder structure is used you simply need to extract all downloaded files directly into `dataset`. For the following sequences, the 'synced+rectified data' and 'calibration' must be downloaded and merged with the dataset folder.


Alternatively if on a UNIX machine you can run the [] script to download the required data automatically.


## KITTIterator
`kittiterator.py` is a helper script for iterating through sequences of the raw KITTI data and performing actions. This script can be used to generate sky masks (requires [SkyAR](https://github.com/jiupinjia/SkyAR)), depth maps (requires [ManyDepth](https://github.com/nianticlabs/manydepth), [MonoDepth2](https://github.com/nianticlabs/monodepth2), or [LEAStereo](https://github.com/XuelianCheng/LEAStereo)) and combine these with the original data to generate hazy image sequences.

`simple_wrappers` contains wrappers to execute the aforementioned components, provided their repositories are present in `prep_tools`

To replicate the generation of the haze video dataset as it exists at [??? wherever I can upload dataset ???], [SkyAR](https://github.com/jiupinjia/SkyAR) and [ManyDepth](https://github.com/nianticlabs/manydepth) must be cloned into `prep_tools`, then the following command should be run:

`python kittierator.py -haze`

This will generate and save sky masks (if not already present), depth maps (if not already present), hazy images and transmission for all raw KITTI sequences present in the `dataset` folder.
Different seed values can be used to randomize the fog density and colour. 
