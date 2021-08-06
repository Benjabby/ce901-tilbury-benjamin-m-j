import os
import sys
import cv2
import numpy as np
from scipy.io import loadmat

def create_haze(img_path, depth_path, out_dir, name, a, vis, mask_path=None, double_mask=True, depth_scale=1):
    img = cv2.imread(img_path, cv2.COLOR_BGR2RGB).astype(np.float32)/255.0
    depth = loadmat(depth_path)['depth']
    depth = depth[...,None]

    if depth.shape[0:2] != img.shape[0:2]:
        h, w = img.shape[0:2]
        x = int(img.shape[1]/2 - w/2)
        y = int(img.shape[0]/2 - h/2)

        depth = depth[y:y+h, x:x+w]
    
    depth*=depth_scale
    depth/=1000

    disp = 1.0/depth
    
    sky_mask = None
    if mask_path is not None:
        sky_mask = 1.0-cv2.imread(mask_path, cv2.IMREAD_ANYDEPTH).astype(np.float32)/65535
        sky_mask = sky_mask[...,None]
        disp = disp*sky_mask
        
    #depth = 1.0/disp
    beta = 3.912/vis;
    beta = np.array([[[beta, beta, beta]]])
    transmission = np.exp(-beta*depth)
    if double_mask and sky_mask is not None: transmission = transmission*sky_mask

    if a.ndim < 3: a = a.reshape([1,1,3])

    image_path = os.path.join(out_dir,"image",name)

    hazy = img*transmission+a*(1.0-transmission)
    hazy = (hazy*255).astype(np.uint8)
    cv2.imwrite(image_path, hazy)

    transmission = (transmission*65536).astype(np.uint16);
    trans_path = os.path.join(out_dir,"transmission",name);
    cv2.imwrite(trans_path, transmission)
