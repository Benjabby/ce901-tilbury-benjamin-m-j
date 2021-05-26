import os
import sys
import json
import argparse
import numpy as np
import PIL.Image as pil
import matplotlib as mpl
import matplotlib.cm as cm
import builtins
from scipy.io import savemat 

import torch
from torchvision import transforms

from inspect import getmembers, isfunction

def print(*args, **kwargs):     # >
    builtins.print("ManyDepth     > ", *args, **kwargs)

def read_calib_file(filepath):
    """Read in a calibration file and parse into a dictionary."""
    data = {}

    with open(filepath, 'r') as f:
        for line in f.readlines():
            key, value = line.split(':', 1)
            # The only non-float values in these files are dates, which
            # we don't care about anyway
            try:
                data[key] = np.array([float(x) for x in value.split()])
            except ValueError:
                pass

    return data

def load_calib_cam_to_cam(cam_to_cam_file):
    # We'll return the camera calibration as a dictionary
    data = {}

    # Load and parse the cam-to-cam calibration data
    filedata = read_calib_file(cam_to_cam_file)

    # Create 3x4 projection matrices
    P_rect_00 = np.reshape(filedata['P_rect_00'], (3, 4))
    P_rect_10 = np.reshape(filedata['P_rect_01'], (3, 4))
    P_rect_20 = np.reshape(filedata['P_rect_02'], (3, 4))
    P_rect_30 = np.reshape(filedata['P_rect_03'], (3, 4))

    # Compute the camera intrinsics
    data['K_cam0'] = P_rect_00[0:3, 0:3]
    data['K_cam1'] = P_rect_10[0:3, 0:3]
    data['K_cam2'] = P_rect_20[0:3, 0:3]
    data['K_cam3'] = P_rect_30[0:3, 0:3]

    data['b00'] = P_rect_00[0, 3] / P_rect_00[0, 0]
    data['b10'] = P_rect_10[0, 3] / P_rect_10[0, 0]
    data['b20'] = P_rect_20[0, 3] / P_rect_20[0, 0]
    data['b30'] = P_rect_30[0, 3] / P_rect_30[0, 0]

    return data


def load_and_preprocess_image(image_path, resize_width, resize_height):
    image = pil.open(image_path).convert('RGB')
    original_width, original_height = image.size
    image = image.resize((resize_width, resize_height), pil.LANCZOS)
    image = transforms.ToTensor()(image).unsqueeze(0)
    if torch.cuda.is_available():
        return image.cuda(), (original_height, original_width)
    return image, (original_height, original_width)


def load_and_preprocess_intrinsics(intrinsics_path, original_width, original_height, resize_width, resize_height):
    K = np.eye(4)
    K[:3, :3] = load_calib_cam_to_cam(intrinsics_path)['K_cam2']
    K[0,:] /= original_width
    K[1,:] /= original_height

    # Convert normalised intrinsics to 1/4 size unnormalised intrinsics.
    # (The cost volume construction expects the intrinsics corresponding to 1/4 size images)
    K[0, :] *= resize_width // 4
    K[1, :] *= resize_height // 4

    invK = torch.Tensor(np.linalg.pinv(K)).unsqueeze(0)
    K = torch.Tensor(K).unsqueeze(0)

    if torch.cuda.is_available():
        return K.cuda(), invK.cuda()
    return K, invK


class SimpleManyDepth():
    def __init__(self, args):

        if os.path.abspath(args.manyd_path) not in sys.path:
            sys.path.insert(1, os.path.abspath(args.manyd_path))
        from manydepth import networks as manydepthnetworks
        from manydepth.layers import transformation_from_parameters
        self.tfp = transformation_from_parameters
        sys.path.remove(os.path.abspath(args.manyd_path))

        
        model_path = os.path.join(args.manyd_path,"models","KITTI_HR")
        device = torch.device("cuda") if torch.cuda.is_available() else torch.device("cpu")
        print("Loading model from ", model_path)

        # Loading pretrained model
        print("Loading pretrained encoder")
        self.encoder_dict = torch.load(os.path.join(model_path, "encoder.pth"), map_location=device)
        self.encoder = manydepthnetworks.ResnetEncoderMatching(18, False,
                                                 input_width=self.encoder_dict['width'],
                                                 input_height=self.encoder_dict['height'],
                                                 adaptive_bins=True,
                                                 min_depth_bin=self.encoder_dict['min_depth_bin'],
                                                 max_depth_bin=self.encoder_dict['max_depth_bin'],
                                                 depth_binning='linear',
                                                 num_depth_bins=96)

        filtered_dict_enc = {k: v for k, v in self.encoder_dict.items() if k in self.encoder.state_dict()}
        self.encoder.load_state_dict(filtered_dict_enc)

        print("Loading pretrained decoder")
        self.depth_decoder = manydepthnetworks.DepthDecoder(num_ch_enc=self.encoder.num_ch_enc, scales=range(4))

        loaded_dict = torch.load(os.path.join(model_path, "depth.pth"), map_location=device)
        self.depth_decoder.load_state_dict(loaded_dict)

        print("Loading pose network")
        pose_enc_dict = torch.load(os.path.join(model_path, "pose_encoder.pth"), map_location=device)
        pose_dec_dict = torch.load(os.path.join(model_path, "pose.pth"), map_location=device)

        self.pose_enc = manydepthnetworks.ResnetEncoder(18, False, num_input_images=2)
        self.pose_dec = manydepthnetworks.PoseDecoder(self.pose_enc.num_ch_enc, num_input_features=1, num_frames_to_predict_for=2)

        self.pose_enc.load_state_dict(pose_enc_dict, strict=True)
        self.pose_dec.load_state_dict(pose_dec_dict, strict=True)

        # Setting states of networks
        self.encoder.eval()
        self.depth_decoder.eval()
        self.pose_enc.eval()
        self.pose_dec.eval()
        if torch.cuda.is_available():
            self.encoder.cuda()
            self.depth_decoder.cuda()
            self.pose_enc.cuda()
            self.pose_dec.cuda()


        self.save_heatmaps = not args.skip_heatmaps
        self.verbose = args.manyd_verbose

    def run_partial(self, data_dir, output_dir, drive_dir, missing):
        print("run_partial not implemented, running full")
        self.run_full(data_dir, output_dir, drive_dir)

    def run_full(self, data_dir, output_dir, drive_dir):

        img_names = os.listdir(data_dir)
        img_names.sort()

        K, invK = None, None
        
        with torch.no_grad():
            prev_img = None
            for idx, img_name in enumerate(img_names):
                
                this_dir = os.path.join(data_dir, img_name)
                current_img, original_size = load_and_preprocess_image(this_dir, resize_width=self.encoder_dict['width'], resize_height=self.encoder_dict['height'])
                if prev_img is None:
                    prev_img = current_img
                    K, invK = load_and_preprocess_intrinsics(os.path.join(drive_dir,'calib_cam_to_cam.txt'), original_width = current_img.shape[1], original_height=current_img.shape[0], resize_width=self.encoder_dict['width'], resize_height=self.encoder_dict['height'])
        

                pose_inputs = [prev_img, current_img]
                pose_inputs = [self.pose_enc(torch.cat(pose_inputs, 1))]
                axisangle, translation = self.pose_dec(pose_inputs)
                pose = self.tfp(axisangle[:, 0], translation[:, 0], invert=True)

                output, lowest_cost, _ = self.encoder(current_image=current_img, lookup_images=prev_img.unsqueeze(1), poses=pose.unsqueeze(1), K=K, invK=invK, min_depth_bin=self.encoder_dict['min_depth_bin'], max_depth_bin=self.encoder_dict['max_depth_bin'])

                output = self.depth_decoder(output)

                disp = output[("disp", 0)]
                disp_resized = torch.nn.functional.interpolate(disp, original_size, mode="bilinear", align_corners=False)
                disp_resized_np = disp_resized.cpu().numpy().squeeze()

                name_dest_mat = os.path.join(output_dir, "mat", img_name[:-4]+".mat")
                depth = 1.0/disp_resized_np
                #np.save(name_dest_npy, depth)
                savemat(name_dest_mat,{'depth':depth})

                if self.save_heatmaps:
                    vs = 1.0-disp_resized_np
                    normalizer = mpl.colors.Normalize(vmin=vs.min(), vmax=vs.max())
                    mapper = cm.ScalarMappable(norm=normalizer, cmap='magma')
                    colormapped_im = (mapper.to_rgba(vs)[:, :, :3] * 255).astype(np.uint8)
                    im = pil.fromarray(colormapped_im)
                    name_dest_im = os.path.join(output_dir, "img", img_name)
                    im.save(name_dest_im)

                if self.verbose: print("Processed {:d} of {:d} images".format(idx + 1, len(img_names)))

                prev_img = current_img


