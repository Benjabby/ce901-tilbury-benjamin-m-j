# Dehaze Evalution

This folder contains:

dehaze-evaluation						
├─ 3rd-party	
│   ├─ SpEED-QA							Folder containing 3rd-Party code for the SpEED-QA evaluator
│   ├─ FADE								Folder containing 3rd-Party code for the FADE evaluator
│   └─ loadCalibrationCamToCam.m		3rd-Party code for reading KITTI calibration data
├─ old									Folder containing various old, currently unused code. Kept for attribution of effort
├─ evaluateDehazers.m					Code to evaluate dehazing systems.
├─ groupedResults.mat					Results on the 11 evaluated dehazing systems described in the dissertaion, averaged and grouped by haze category.
├─ groupResults.m						Function to group results by category.
├─ imageCollage.m						Function to create a collage of dehazed images for given dehazing systems
├─ labelScatter.m						Function to plot labeled scatter graphs.
├─ normalisedScore.m					Function to normalise results.
├─ plotScoring.m						Function to normalise and plot results.
├─ results.mat							Results on the 11 evaluated dehazing systems described in the dissertaion, grouped by sequence
└─ saveVideos.m							Function to save the video of a given dehazing system's results on the entire dataset.


In order to run the evaluation the complete dataset is required. The hazy images, transmission, and generated sky masks can be downloaded from https://essexuniversity.box.com/s/7hzugcfiv14b6kevpb02ppp1z4d8jfrs 
The KITTI raw data for these sequences must also be downloaded and merged with the dataset folder. The script at `../haze-video-dataset/dataset/download_kitti.sh` can automatically download the required raw KITTI data.
Alternatively 