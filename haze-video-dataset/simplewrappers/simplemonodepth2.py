from __future__ import absolute_import, division, print_function

import os
import sys
import glob
import argparse
import numpy as np
import PIL.Image as pil
import matplotlib as mpl
import matplotlib.cm as cm
import builtins
import importlib.util
from scipy.io import savemat 

import torch
from torchvision import transforms, datasets

def print(*args, **kwargs):
    builtins.print("MonoDepth2    > ", *args, **kwargs)

class SimpleMonoDepth2():
    def __init__(self, args):

        model_name = "mono+stereo_1024x320"

        # Janky fix for importing avoiding name clashes with other components
        if os.path.abspath(args.mono2_path) not in sys.path:
            sys.path.insert(1, os.path.abspath(args.mono2_path))   
        import networks as monodepth2networks
        sys.modules['monodepth2networks'] = sys.modules['networks']
        del sys.modules['networks']
        monodepth2networks.print = print
        sys.path.remove(os.path.abspath(args.mono2_path))

        if torch.cuda.is_available():
            self.device = torch.device("cuda")
        else:
            self.device = torch.device("cpu")

        model_path = os.path.join(args.mono2_path,"models", model_name)
        print("Loading model from ", model_path)
        encoder_path = os.path.join(model_path, "encoder.pth")
        depth_decoder_path = os.path.join(model_path, "depth.pth")

        # LOADING PRETRAINED MODEL
        print("Loading pretrained encoder")
        self.encoder = monodepth2networks.ResnetEncoder(18, False)
        self.loaded_dict_enc = torch.load(encoder_path, map_location=self.device)

        # extract the height and width of image that this model was trained with
        self.feed_height = self.loaded_dict_enc['height']
        self.feed_width = self.loaded_dict_enc['width']
        self.filtered_dict_enc = {k: v for k, v in self.loaded_dict_enc.items() if k in self.encoder.state_dict()}
        self.encoder.load_state_dict(self.filtered_dict_enc)
        self.encoder.to(self.device)
        self.encoder.eval()

        print("Loading pretrained decoder")
        self.depth_decoder = monodepth2networks.DepthDecoder(num_ch_enc=self.encoder.num_ch_enc, scales=range(4))

        self.loaded_dict = torch.load(depth_decoder_path, map_location=self.device)
        self.depth_decoder.load_state_dict(self.loaded_dict)

        self.depth_decoder.to(self.device)
        self.depth_decoder.eval()

        self.save_heatmaps = not args.skip_heatmaps
        self.verbose = args.mono2_verbose

    def run_partial(self, data_dir, output_dir, files):
           print("run_partial currently unimplemented. Defaulting to run_full")
           self.run_full(data_dir, output_dir)

    def run_full(self, data_dir, output_dir):

        img_names = os.listdir(data_dir)

        # PREDICTING ON EACH IMAGE IN TURN
        with torch.no_grad():
            for idx, img_name in enumerate(img_names):

                this_dir = os.path.join(data_dir, img_name)
                # Load image and preprocess
                input_image = pil.open(this_dir).convert('RGB')
                original_width, original_height = input_image.size
                input_image = input_image.resize((self.feed_width, self.feed_height), pil.LANCZOS)
                input_image = transforms.ToTensor()(input_image).unsqueeze(0)

                # PREDICTION
                input_image = input_image.to(self.device)
                features = self.encoder(input_image)
                outputs = self.depth_decoder(features)

                disp = outputs[("disp", 0)]
                disp_resized = torch.nn.functional.interpolate(disp, (original_height, original_width), mode="bilinear", align_corners=False)
                disp_resized_np = disp_resized.squeeze().cpu().numpy()

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
