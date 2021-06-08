import sys

import numpy as np
import matplotlib.pyplot as plt
import cv2
import os
import glob
import torch
import builtins
from cv2.ximgproc import guidedFilter
import imageio
import importlib
from scipy.ndimage import uniform_filter

def print(*args, **kwargs):
    builtins.print("SkyAR         > ", *args, **kwargs)

def special_guided(guide, target, radius, eps):
    m, n = guide.shape

    avgDenom = window_sum(np.ones_like(guide), radius);
    mean_g = window_sum(guide, radius) / avgDenom;
    corr_gg = window_sum(guide * guide, radius) / avgDenom;
    var_g = corr_gg - mean_g * mean_g;
    mask = mean_g-var_g;

    mean_t = window_sum(target, radius) / avgDenom;
    corr_gt = window_sum(guide * target, radius) / avgDenom;
    cov_gt = corr_gt - mean_g * mean_t;
    a = cov_gt / (var_g + eps);
    b = mean_t - a * mean_g;
    mean_a = window_sum(a, radius) / avgDenom;
    mean_b = window_sum(b, radius) / avgDenom;

    q = mean_a * guide + mean_b;

    qd = np.copy(q);
    qd[q>mask]=1;

    return q, qd


    
def window_sum(image, r):
    h, w = image.shape
    sumImg = np.zeros_like(image)

    im_cum = np.cumsum(image, axis=0);
    sumImg[:r+1, :] = im_cum[r:2*r+1, :]
    sumImg[r+1:h-r, :] = im_cum[2*r+1:, :] - im_cum[:h-2*r-1, :];
    sumImg[h-r:, :] = np.tile(im_cum[h-1, :], [r, 1]) - im_cum[h-2*r-1:h-r-1, :];

    im_cum = np.cumsum(sumImg,axis=1);

    sumImg[:, :r+1] = im_cum[:, r:2*r+1];
    sumImg[:, r+1:w-r] = im_cum[:, 2*r+1:] - im_cum[:, :w-2*r-1];
    sumImg[:, w-r:] = np.tile(im_cum[:, w-1], [r, 1]).T - im_cum[:, w-2*r-1:w-r-1];

    return sumImg


class SimpleSkyFilter():

    def __init__(self, args):

        # Janky fix for importing avoiding name clashes with other components
        if os.path.abspath(args.skyar_path) not in sys.path:
            sys.path.insert(1, os.path.abspath(args.skyar_path))
        import networks as skyarnetworks
        sys.modules['skyarnetworks'] = sys.modules['networks']
        del sys.modules['networks']
        skyarnetworks.print = print
        sys.path.remove(os.path.abspath(args.skyar_path))

        self.verbose = args.sky_masks_verbose
            
        self.device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
        ckptdir = os.path.join(args.skyar_path,"checkpoints_G_coord_resnet50")
        self.net_G = skyarnetworks.define_G(input_nc=3, output_nc=1, ngf=64, netG="coord_resnet50").to(self.device)
        
        print('Loading SkyAR Model...')
        checkpoint = torch.load(os.path.join(ckptdir, 'best_ckpt.pt'),map_location=None if torch.cuda.is_available() else self.device)
        self.net_G.load_state_dict(checkpoint['model_G_state_dict'])
        self.net_G.to(self.device)
        self.net_G.eval()

        self.in_size_w, self.in_size_h = 384,384

    def skymask_refinement(self, G_pred, img):

        r, eps = 20, 0.01
        refined_skymask = guidedFilter(img[:,:,2], G_pred, r, eps)
        
        return np.clip(refined_skymask, a_min=0, a_max=1)

    def skymask_improved_refinement(self, G_pred, img):

        r, eps = 20, 0.001
        gray = np.mean(img, axis=-1)
        
        basic, helper = special_guided(gray, G_pred, r, eps)
        
        refined_skymask, _ = special_guided(basic, helper, r, eps)
        
        return np.clip(refined_skymask, a_min=0, a_max=1)


    def synthesize(self, img_HD):

        h, w, c = img_HD.shape

        img = cv2.resize(img_HD, (self.in_size_w, self.in_size_h))

        img = np.array(img, dtype=np.float32)
        img = torch.tensor(img).permute([2, 0, 1]).unsqueeze(0)

        with torch.no_grad():
            G_pred = self.net_G(img.to(self.device))
            G_pred = torch.nn.functional.interpolate(G_pred, (h, w), mode='bicubic', align_corners=False)
            G_pred = G_pred[0, :].permute([1, 2, 0])
            G_pred = torch.cat([G_pred, G_pred, G_pred], dim=-1)
            G_pred = np.array(G_pred.detach().cpu())
            G_pred = np.clip(G_pred, a_max=1.0, a_min=0.0)

        G_pred = G_pred[...,0]
        skymask = self.skymask_refinement(G_pred, img_HD)
        improved = self.skymask_improved_refinement(G_pred, img_HD)

        return G_pred, skymask, improved


    def cvtcolor_and_resize(self, img_HD):

        img_HD = cv2.cvtColor(img_HD, cv2.COLOR_BGR2RGB)
        img_HD = np.array(img_HD / 255., dtype=np.float32)
        #img_HD = cv2.resize(img_HD, outsize)

        return img_HD

    def run_partial(self, data_dir, output_dir, files):
        print("run_partial currently unimplemented. Defaulting to run_full")
        self.run_full(data_dir, output_dir)
        
    def run_full(self, data_dir, output_dir):

        img_names = os.listdir(data_dir)

        for idx in range(len(img_names)):

            this_dir = os.path.join(data_dir, img_names[idx])
            img_HD = cv2.imread(this_dir, cv2.IMREAD_COLOR)
            img_HD = self.cvtcolor_and_resize(img_HD)

            G_pred, skymask, improved  = self.synthesize(img_HD)

            raw_path = os.path.join(output_dir, "raw","img", img_names[idx])
            im = (G_pred*65535)
            imageio.imwrite(raw_path, im.astype(np.uint16))
            
            refined_path = os.path.join(output_dir, "default_refined","img", img_names[idx])
            im = (skymask*65535)
            imageio.imwrite(refined_path, im.astype(np.uint16))

            improved_path = os.path.join(output_dir, "improved","img", img_names[idx])
            im = (improved*65535)
            imageio.imwrite(improved_path, im.astype(np.uint16))

            if self.verbose: print("    > ",'Processed: %d / %d ...' % (idx+1, len(img_names)))

