Repository for my masters dissertation project.

# Introduction
This repository is comprised of three main sections;

`haze-video-dataset` contains code used to generate the haze video dataset. In order to use the dehaze evaluator, instructions in this folder must be followed to download the haze video dataset and raw KITTI data.

`dehaze-systems` contains the code for the dehazing systems.

`dehaze-evaluation` contains code to evaluate the dehazing systems.

For more information see the README files inside each folder. 

# General Requirements

*  MATLAB 2021a (Older versions may work) 
	* The following MATLAB Toolboxes:
		* Control System Toolbox
		* Signal Processing Toolbox
		* Mapping Toolbox
		* Image Processing Toolbox
		* Statistics and Machine Learning Toolbox
	* The following MATLAB Add-ons:
		* [Text Wait/Progress Bar by Xiaoxuan He](https://www.mathworks.com/matlabcentral/fileexchange/71638-text-wait-progress-bar)
* Python 3
* Windows OS is recommended, as `ChenDehazer` and `LiDehazer` require the Windows only .NET Framework. Everything else can be done on a Unix OS
# Code Overview and References
````
Ref		Structure							Description
-----------------------------------------------------------------------------------
		dehaze-evaluation					Folder containing  code for the evaluation of dehazing systems
		├─ 3rd-party						Folder containing  3rd party tools for evaluation.
[1, 2]		│   ├─ SpEED-QA						3rd-Party code for the SpEED-QA evaluator
[3, 4]		│   ├─ FADE						3rd-Party code for the FADE evaluator
[5]		│   └─ loadCalibrationCamToCam.m			3rd-Party code for reading KITTI calibration data
		├─ old							Folder containing  various old, currently unused code. Kept for attribution of effort
		├─ example-results					Folder containing  small sample of image results for each dehazing
		├─ empiricalComplexity.m				Plots the runtime of dehazers for various size inputs.
		├─ evaluateDehazers.m					Code to evaluate dehazing systems.
		├─ groupedResults.mat					Results on the evaluated dehazing systems, averaged and grouped by haze category.
		├─ groupResults.m					Function to group results by category.
		├─ imageCollage.m					Function to create a collage of dehazed images for given dehazing systems
		├─ labelScatter.m					Function to plot labelled scatter graphs.
		├─ normalisedScore.m					Function to normalise results.
		├─ plotScoring.m					Function to normalise and plot results.
		├─ results.mat						Full results on the evaluated dehazing systems, grouped by sequence
		└─ saveVideos.m						Function to save the video of a given dehazing system's results on the entire dataset.
		dehaze-systems						Folder containing  code for the implemented dehazing systems
		├─ abstracts						Folder containing  abstract base classes for dehazing systems
[6, 7]		│   ├─ BaseDehazer.m					Abstract base class for all dehazing system. Includes guided filter  and window sum filter code.
		│   └─ ExternalPythonDehazer.m				Abstract base class for external python dehazing systems.  (Windows Only)
		├─ neural-nets						Folder containing  code for the neural network methods
[8, 9]		│   ├─ AOD						Folder containing model definitions for AOD neural network
[10, 11]	│   ├─ PMHLD						Folder containing model definitions for PMHLD neural network
		│   └─ pyinterfacer.py					Python program for interfacing with MATLAB to run the two neural network methods
		├─ old							Folder containing  various old, currently unused code. Kept for attribution of effort
		├─ optimise-params					Folder containing  code to optimise the TilburyDehazer
		│   ├─ optimise-tilbury.m				Code for running Bayesian optimisation on the TilburyDehazer
		│   └─ optims.mat					Results of the Bayesian optimisations
[12-14]		├─ BermanDehazer.m					Dehazing system implementing Berman et al.'s method
[15, 16]	├─ CaiDehazer.m						Dehazing system implementing Cai et al.'s method
		├─ ChenDehazer.m					Dehazing system wrapper for Chen et al.'s PMHLD method (Requires TensorFlow) (Windows Only)
[17, 18]	├─ HeDehazer.m						Dehazing system implementing He et al.'s DCP method using the guided filter.
		├─ LiDehazer.m						Dehazing system wrapper for Li et al.'s AOD neural network (Requires PyTorch) (Windows Only)
[19, 20]	├─ TarelDehazer.m					Dehazing system implementing Tarel et al.'s method
		├─ TilburyDehazer.m					Proposed dehazing system containing different implementations of the smooth-extremum filter
		├─ TsaiDehazer.m					Dehazing system implementing Tsai et al.'s method
[21, 22]	├─ ZhuDehazer.m						Dehazing system implementing Zhu et al.'s CAP method
		└─ UnrefinedDehazer.m					Dehazing system implementing He et al.'s DCP method without refinement.
		haze-video-dataset					Folder containing  code for generating the haze video dataset.
		├─ dataset						Default location for the haze dataset and KITTI data to be downloaded to. 
[23]		│   └─ download_kitti.sh				3rd-Party script to download the required raw KITTI sequences automatically (Requires Bash shell)
		├─ prep-tools						Default location to put repositories of components needed in the generation of the dataset.
		├─ simplewrappers					Folder containing  wrappers to interface with the 3rd-Party components.
		│   ├─ simpleleastereo.py				Wrapper for running LEAStereo to generate depth maps
		│   ├─ simplemanydepth.py				Wrapper for running ManyDepth to generate depth maps
		│   ├─ simplemonodepth2.py				Wrapper for running MonoDepth2 to generate depth maps
		│   └─ simpleskyar.py					Wrapper for running SkyAR to generate sky masks
		├─ hazeutils.py						Code to generate hazy images
		├─ hazegenerator.mlap					MATLAB app to preview hazy sequences
		└─ kittiterator.py					A multifunctional program to iterate through the KITTI data. Generates the haze dataset.
		figures 						Folder containing  figures used in the dissertation paper.
		haze-subset.txt						List of KITTI sequences in the haze dataset. By default, any other sequences encountered are ignored.
````

[1]	C. G. Bampis, P. Gupta, R. Soundararajan, and A. C. Bovik. (2017, Accessed: 26/08/2021). SpEED-QA Software Release. Available: https://github.com/christosbampis/SpEED-QA_release

[2]	C. G. Bampis, P. Gupta, R. Soundararajan, and A. C. Bovik, "SpEED-QA: Spatial efficient entropic differencing for image and video quality," IEEE signal processing letters, vol. 24, no. 9, pp. 1333-1337, 2017.

[3]	L. K. Choi, J. You, and A. C. Bovik. (2015, Accessed: 26/08/2021). FADE Software Release. Available: http://live.ece.utexas.edu/research/fog/FADE_release.zip

[4]	L. K. Choi, J. You, and A. C. Bovik, "Referenceless prediction of perceptual fog density and perceptual image defogging," IEEE Transactions on Image Processing, vol. 24, no. 11, pp. 3888-3901, 2015.

[5]	A. Geiger, P. Lenz, C. Stiller, and R. Urtasun, "Vision meets robotics: The kitti dataset," The International Journal of Robotics Research, vol. 32, no. 11, pp. 1231-1237, 2013.

[6]	K. He. (2015, Accessed: 26/08/2021). Fast Guided Filter MATLAB Code. Available: http://kaiminghe.com/eccv10/fast-guided-filter-code-v1.rar

[7]	K. He and J. Sun, "Fast guided filter," arXiv preprint arXiv:1505.00996, 2015.

[8]	M. Singal. (2018, Accessed: 26/08/2021). PyTorch-Image-Dehazing. Available: https://github.com/MayankSingal/PyTorch-Image-Dehazing

[9]	B. Li, X. Peng, Z. Wang, J. Xu, and D. Feng, "Aod-net: All-in-one dehazing network," in Proceedings of the IEEE international conference on computer vision, 2017, pp. 4770-4778.

[10]	W.-T. Chen, H.-Y. Fang, J.-J. Ding, and S.-Y. Kuo, "PMHLD: Patch Map-Based Hybrid Learning DehazeNet for Single Image Haze Removal," IEEE Transactions on Image Processing, vol. 29, pp. 6773-6788, 2020.

[11]	W.-T. Chen, H.-Y. Feng, J.-J. Ding, and S.-Y. Kuo. (2020, Accessed: 26/08/2021). PMHLD-Patch-Map-Based-Hybrid-Learning-DehazeNet-for-Single-Image-Haze-Removal. Available: https://github.com/weitingchen83/Dehazing-PMHLD-Patch-Map-Based-Hybrid-Learning-DehazeNet-for-Single-Image-Haze-Removal-TIP-2020

[12]	D. Berman. (2017, Accessed: 26/08/2021). non-local-dehazing. Available: https://github.com/danaberman/non-local-dehazing

[13]	D. Berman, S. Avidan, and T. Treibitz, "Non-local image dehazing," in Proceedings of the IEEE conference on computer vision and pattern recognition, 2016, pp. 1674-1682.

[14]	D. Berman, S. Avidan, and T. Treibitz, "Air-light estimation using haze-lines," 2017 IEEE International Conference on Computational Photography (ICCP), pp. 1-9, 2017.

[15]	B. Cai, X. Xu, and D. Tao, "Real-time video dehazing based on spatio-temporal mrf," in Pacific Rim conference on multimedia, 2016, pp. 315-325: Springer.

[16]	B. Cai, X. Xu, and D. Tao. (2016, Accessed: 26/08/2021). Real-time Video Dehazing based on Spatio-temporal MRF Software Release. Available: https://github.com/caibolun/ST-MRF/

[17]	K. He, J. Sun, and X. Tang, "Single image haze removal using dark channel prior," IEEE transactions on pattern analysis and machine intelligence, vol. 33, no. 12, pp. 2341-2353, 2010.

[18]	S. Tierney. (2014, Accessed: 26/08/2021). MATLAB implementation of "Single Image Haze Removal Using Dark Channel Prior". Available: https://github.com/sjtrny/Dark-Channel-Haze-Removal

[19]	J.-P. Tarel. (2013, Accessed: 26/08/2021). Vision Enhancement in Homogeneous and Heterogeneous Fog MATLAB Source Code. Available: http://perso.lcpc.fr/tarel.jean-philippe/visibility/visibresto2.zip

[20]	J.-P. Tarel, N. Hautiere, L. Caraffa, A. Cord, H. Halmaoui, and D. Gruyer, "Vision enhancement in homogeneous and heterogeneous fog," IEEE Intelligent Transportation Systems Magazine, vol. 4, no. 2, pp. 6-20, 2012.

[21]	Q. Zhu, J. Mai, and L. Shao. (2015, Accessed: 26/08/2021). Color-Attenuation-Prior-Dehazing. Available: https://github.com/JiamingMai/Color-Attenuation-Prior-Dehazing

[22]	Q. Zhu, J. Mai, and L. Shao, "A fast single image haze removal algorithm using color attenuation prior," IEEE transactions on image processing, vol. 24, no. 11, pp. 3522-3533, 2015.

[23]	pean1128. (2018, Accessed: 26/08/2021). KITTI-download-scipt  [sic]. Available: https://github.com/pean1128/KITTI-download-scipt

