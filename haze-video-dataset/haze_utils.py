import os
import sys
import cv2
import numpy as np
from scipy.io import loadmat

def create_haze(img_path, depth_path, out_path, a, vis, mask_path=None, double_mask=True, depth_scale=1):
    img = cv2.imread(img_path, cv2.COLOR_BGR2RGB).astype(np.float32)/255.0
    depth = loadmat(depth_path)['depth']
    depth = depth[...,None]
    depth*=depth_scale
    depth/=1000

    disp = 1.0/depth
    
    sky_mask = None
    if mask_path is not None:
        sky_mask = 1.0-cv2.imread(mask_path, cv2.IMREAD_ANYDEPTH).astype(np.float32)/65535
        sky_mask = sky_mask[...,None]
        disp = disp*sky_mask
        
    
    beta = 3.912/vis;
    beta = np.array([[[beta, beta, beta]]])
    transmission = np.exp(-beta*depth)
    if double_mask and sky_mask is not None: transmission = transmission*sky_mask

    if a.ndim < 3: a = a.reshape([1,1,3])

    hazy = img*transmission+a*(1.0-transmission)
    hazy = (hazy*255).astype(np.uint8)
    cv2.imwrite(out_path, hazy)
