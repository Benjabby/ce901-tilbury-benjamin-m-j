import os
import sys
import json
import argparse
import numpy as np
from PIL import Image
import matplotlib as mpl
import matplotlib.cm as cm
import builtins
from scipy.io import savemat 
from functools import partial
import cv2

from argparse import Namespace
import torch
from torch.autograd import Variable
from torchvision import transforms

from inspect import getmembers, isfunction

#Why must you be so awkward with your printing.
def prefix_print(old_print, *args, **kwargs):     # >
    if '\n' in str(args[0]):
        sc = args[0].split('\n')
        for s in sc[:-1]:
            old_print("LEAStereo     > ",s, *args[1:], **kwargs)
    else:
        old_print("LEAStereo     > ", *args, **kwargs)

def prefix(f):
    def closure(*args, **kwargs):
        old = builtins.print
        builtins.print = partial(prefix_print,old)
        v = f(*args, **kwargs)
        builtins.print = old
        if v is not None:
            return v
    return closure

@prefix
def test_transform(temp_data, crop_height, crop_width):
    _, h, w=np.shape(temp_data)

    if h <= crop_height and w <= crop_width: 
        # padding zero 
        temp = temp_data
        temp_data = np.zeros([6, crop_height, crop_width], 'float32')
        temp_data[:, crop_height - h: crop_height, crop_width - w: crop_width] = temp    
    else:
        start_x = int((w - crop_width) / 2)
        start_y = int((h - crop_height) / 2)
        temp_data = temp_data[:, start_y: start_y + crop_height, start_x: start_x + crop_width]
    left = np.ones([1, 3,crop_height,crop_width],'float32')
    left[0, :, :, :] = temp_data[0: 3, :, :]
    right = np.ones([1, 3, crop_height, crop_width], 'float32')
    right[0, :, :, :] = temp_data[3: 6, :, :]
    return torch.from_numpy(left).float(), torch.from_numpy(right).float(), h, w

@prefix
def load_data(leftname, rightname):
    left = Image.open(leftname)
    right = Image.open(rightname)
    size = np.shape(left)
    height = size[0]
    width = size[1]
    temp_data = np.zeros([6, height, width], 'float32')
    left = np.asarray(left)
    right = np.asarray(right)
    r = left[:, :, 0]
    g = left[:, :, 1]
    b = left[:, :, 2]
    temp_data[0, :, :] = (r - np.mean(r[:])) / np.std(r[:])
    temp_data[1, :, :] = (g - np.mean(g[:])) / np.std(g[:])
    temp_data[2, :, :] = (b - np.mean(b[:])) / np.std(b[:])
    r = right[:, :, 0]
    g = right[:, :, 1]
    b = right[:, :, 2]	
    #r,g,b,_ = right.split()
    temp_data[3, :, :] = (r - np.mean(r[:])) / np.std(r[:])
    temp_data[4, :, :] = (g - np.mean(g[:])) / np.std(g[:])
    temp_data[5, :, :] = (b - np.mean(b[:])) / np.std(b[:])
    return temp_data


class SimpleLEAStereo():
    @prefix
    def __init__(self, args):

        if os.path.abspath(args.lea_path) not in sys.path:
            sys.path.insert(1, os.path.abspath(args.lea_path))
        
        from retrain.LEAStereo import LEAStereo
        sys.path.remove(os.path.abspath(args.lea_path))

        stargs = Namespace()
        stargs.kitti2015 = 1
        stargs.cuda = False
        stargs.maxdisp = 192
        stargs.data_path = os.path.join(args.lea_path,'dataset','kitti2015','testing')
        stargs.test_list = os.path.join(args.lea_path,'dataloaders','lists','kitti2015_test.list')
        stargs.save_path = os.path.join(args.lea_path,'predict','kitti2015','images')
        stargs.fea_num_layers = 6
        stargs.mat_num_layers = 12
        stargs.fea_filter_multiplier = 8
        stargs.fea_block_multiplier = 4
        stargs.fea_step = 3 
        stargs.mat_filter_multiplier = 8
        stargs.mat_block_multiplier = 4
        stargs.mat_step = 3
        stargs.net_arch_fea = os.path.join(args.lea_path,'run','sceneflow','best','architecture','feature_network_path.npy') 
        stargs.cell_arch_fea = os.path.join(args.lea_path,'run','sceneflow','best','architecture','feature_genotype.npy')
        stargs.net_arch_mat = os.path.join(args.lea_path,'run','sceneflow','best','architecture','matching_network_path.npy')
        stargs.cell_arch_mat = os.path.join(args.lea_path,'run','sceneflow','best','architecture','matching_genotype.npy')

        self.model = LEAStereo(stargs)

        self.model = torch.nn.DataParallel(self.model)

        self.cuda = torch.cuda.is_available()

        if self.cuda:
            self.model.cuda()
            checkpoint = torch.load(os.path.join(args.lea_path,'run','Kitti15','best','best.pth'))
        else:
##            print("NOTE: Unfortunately the native LEAStereo code requires a CUDA complition of torch and cannot be run on CPU")
##            print("See the README file for how to edit LEAStereo to get this code running on CPU")
##            print("If you have already made these changes, please comment out this code")
##            exit(0)

            
            checkpoint = torch.load(os.path.join(args.lea_path,'run','Kitti15','best','best.pth'),map_location=torch.device('cpu'))

        self.model.load_state_dict(checkpoint['state_dict'], strict=True)

        self.model.eval()
        
        self.save_heatmaps = not args.skip_heatmaps
        self.verbose = args.lea_verbose

    @prefix
    def run_partial(self, data_dir, output_dir, drive_dir, missing):
        print("run_partial not implemented, running full")
        self.run_full(data_dir, output_dir, drive_dir)

    @prefix
    def run_full(self, data_dir, output_dir, drive_dir):

        img_names = os.listdir(data_dir)
        img_names.sort()
        with torch.no_grad():
            prev_img = None
            for idx, img_name in enumerate(img_names):
                
                left_dir = os.path.join(data_dir, img_name)
                right_dir = left_dir.replace("image_02","image_03")
                
                left, right, height, width = test_transform(load_data(left_dir, right_dir), 384, 1248)

                left = Variable(left, requires_grad = False)
                right = Variable(right, requires_grad = False)

                if self.cuda:
                    input1 = input1.cuda()
                    input2 = input2.cuda()

                pred = self.model(left,right)
                pred = pred.cpu().numpy().squeeze()

                name_dest_mat = os.path.join(output_dir, "mat", img_name[:-4]+'.mat')
                depth = 1.0/pred
                #np.save(name_dest_npy, depth)
                savemat(name_dest_mat,{'depth':depth})

                if self.save_heatmaps:
                    vs = 1.0-pred
                    normalizer = mpl.colors.Normalize(vmin=vs.min(), vmax=vs.max())
                    mapper = cm.ScalarMappable(norm=normalizer, cmap='magma')
                    colormapped_im = (mapper.to_rgba(vs)[:, :, :3] * 255).astype(np.uint8)
                    im = Image.fromarray(colormapped_im)
                    name_dest_im = os.path.join(output_dir, "img", img_name)
                    im.save(name_dest_im)

                if self.verbose: print("Processed {:d} of {:d} images".format(idx + 1, len(img_names)))


