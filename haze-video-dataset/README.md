# Haze Video Dataset

This folder contains code []

[The "Haze Video Dataset" can be downloaded here](https://essexuniversity.box.com/s/7hzugcfiv14b6kevpb02ppp1z4d8jfrs) It comprises of 18 sequences of hazy road scenes and corresponding transmission maps for a total of 15,719 colour hazy images, 15,719 grayscale transmission maps, 15,719 sky masks, and 18 haze property files.
Each sequence is derived from a sequence of the KITTI dataset. Only a subset of sequences from the full KITTI dataset is used. For simplicity the dataset follows the same folder structure as the KITTI dataset.
Each sequence has a different, randomly chosen fog density and colour, this information is stored in the `haze/props.json` file for each sequence.

The downloaded dataset should be unziped into `dataset` within this folder.

In addition to the haze video dataset, the raw KITTI data is needed to run the evaluation code, or to generate further hazy sequences.

The raw KITTI data can be downloaded from their [website](http://www.cvlibs.net/datasets/kitti/raw_data.php). As the KITTI folder structure is used you simply need to extract all downloaded files directly into `dataset`. For the following sequences, the 'synced+rectified data' and 'calibration' must be downloaded and merged with the dataset folder.
Alternatively if on a UNIX machine you can run the `dataset/download_kitti.sh` script. This will download all the required sequences and place them in the `dataset` folder.

Note that not all of the data included in KITTI sequences is needed for the generation of the dataset or for the evaluation.

The KITTIterator (`kittiterator.py`) has many function, one of which is to remove unncessary files of the dataset if needed:

Use `python kittierator.py --clean 0` to remove everything except the files needed to generate arbitrary hazy versions of the dataset. This keeps the following folders for each sequence: `image_02`, `image_03`, `haze`, `sky_mask`, `depth`, as well as the required calibration files.

Use `python kittierator.py --clean 1` to remove everything except the files needed only to run the evaluation. This keeps the following folders for each sequence `image_02`, `haze`, `sky_mask`, as well as the required calibration files.

If the dataset is not located in `./dataset/` use the `--path` argument to provide the KITTIterator the path.

The KITTIterator can also can be used to generate sky masks (requires [SkyAR](https://github.com/jiupinjia/SkyAR)), depth maps (requires [ManyDepth](https://github.com/nianticlabs/manydepth), [MonoDepth2](https://github.com/nianticlabs/monodepth2), or [LEAStereo](https://github.com/XuelianCheng/LEAStereo)) and combine these with the original data to generate hazy image sequences.

`simple_wrappers` contains wrappers to execute the aforementioned components, provided their repositories are present in `prep_tools`

To replicate the generation of the haze video dataset [SkyAR](https://github.com/jiupinjia/SkyAR) and [ManyDepth](https://github.com/nianticlabs/manydepth) must be cloned into `prep_tools`, then the following command should be run:

`python kittierator.py -haze`

This will generate and save sky masks (if not already present), ManyDepth depth maps (if not already present), hazy images and transmission for all raw KITTI sequences present in the `dataset` folder.
Different seed values can be used to randomize the fog density and colour. 
