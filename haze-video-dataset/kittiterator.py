import sys
import hashlib
import numpy as np
import matplotlib.pyplot as plt
import cv2
import os
import glob
import argparse
from argparse import RawDescriptionHelpFormatter
import shutil
import json
from functools import reduce
from datetime import datetime
from PIL import Image
from haze_utils import create_haze

# Note, I use a slightly different convention when it comes to arguments that I find far more informative:
# -x is for flags (i.e those requiring the form '[...] -x [...]')
# --x is for option arguments (i.e. those requiring the form '[...] --x [operand] [...]')

# General arguments
parser = argparse.ArgumentParser(description='Iterates through KITTI dataset and does various things\n\
Note; I use a slightly different convention when it comes to arguments for simple scripts, that I find far more informative:\n\
Single dash { -x} is used for flags (i.e. "| -xxx | -yyy | -zzz |")\n\
Double dash {--x} is used for parameter arguments (i.e. "| --xxx 20 | --yyy something | --zzz something/else |")', formatter_class=RawDescriptionHelpFormatter)
parser.add_argument('--path', type=str, default='./KITTI', metavar='str', help='KITTI data path')
#parser.add_argument('--ftype',type=str, default='both',choices=['img','mat','both'], help='Whether to save using images, matlab files or both')
#parser.add_argument("-mat", help='If a required map is not found as an .mat file a new one will be generated instead of loading from an image', action='store_true')
parser.add_argument("-clean", help='Deletes unnecessary KITTI files', action='store_true')
parser.add_argument("-override", help='Replace any existing masks or depth maps', action='store_true')
parser.add_argument("-skip_heatmaps", help="Don't save KITTI style disparity images for any depth maps", action='store_true')
parser.add_argument("-verbose", help='Run all components as verbose', action='store_true')

# Sky Mask Generation arguments
parser.add_argument("-sky_masks", help='Generate and save sky masks', action='store_true')
parser.add_argument('--skyar_path', type=str, default='./prep_tools/SkyAR', metavar='str', help='SkyAR repository path')
parser.add_argument("-sky_masks_verbose", help='Verbose output of SkyAR', action='store_true')

# MonoDepth2 Generation arguments
parser.add_argument("-mono2", help='Generate and save MonoDepth2 depth maps', action='store_true')
parser.add_argument('--mono2_path', type=str, default='./prep_tools/monodepth2', metavar='str', help='MonoDepth2 repository path')
parser.add_argument("-mono2_verbose", help='Use verbose output for MonoDepth2', action='store_true')

# ManyDepth Generation arguments
parser.add_argument("-manyd", help='Generate and save ManyDepth depth maps', action='store_true')
parser.add_argument('--manyd_path', type=str, default='./prep_tools/manydepth', metavar='str', help='ManyDepth repository path')
parser.add_argument("-manyd_verbose", help='Use verbose output for ManyDepth', action='store_true')

# LEAStereo Generation arguments
parser.add_argument("-lea", help='Generate and save LEAStereo depth maps', action='store_true')
parser.add_argument('--lea_path', type=str, default='./prep_tools/LEAStereo', metavar='str', help='LEAStereo repository path')
parser.add_argument("-lea_verbose", help='Use verbose output for LEAStereo', action='store_true')

# Haze Image Generation         [DEPTH MAP + SKY MASK]
parser.add_argument("-haze", help='Generate hazy frames from KITTI sequences and depth maps', action='store_true')
parser.add_argument('--haze_depth_type',type=str, default='manyd',choices=['mono2','manyd','lea'], help='Which depth map to use for haze generation')
parser.add_argument("--haze_mask_type",type=str, default='double',choices=['normal','double','none'], help='Method of sky masking for haze generation')
parser.add_argument("--haze_seed",type=int, default=0, help='Random seed for haze generation')
parser.add_argument("-haze_verbose", help='Use verbose output for haze frame generation', action='store_true')
#parser.add_argument('--masked_depth_mask',type=str, default='default',choices=['default','volume'], help='What sky mask should be used to mask a depth map')

############## [UNUSED]
# Sky Mask Refinement arguments                 [SKY MASK]

# Sky Augmented LiDAR Map Generation arguments  [SKY MASK]

# Full Augmented LiDAR Map Generation arguments [UNGUIDED DEPTH MAP + SKY MASK]

# Augmented Depth Generation arguments          [AUGMENTED LIDAR MAP [UNGUIDED DEPTH MAP + SKY MASK]]

############## KITTI Data
# |OK| implicit  | Colour Image Data            | [sequence path]/image_02/data
#               [LIDAR MAPS]
# |OK| implicit  | LiDAR Maps                   | [sequence path]/proj_depth/groundtruth/image_02
# |  | generated | Sky Augmented LiDAR Maps     | [sequence path]/proj_depth/sky_aug
# |  | generated | Full Augmented LiDAR Maps    | [sequence path]/proj_depth/full_aug
#               [SKY MASKS]
# |  | generated | Raw Sky Masks                | [sequence path]/sky_mask/raw/img
# |  | generated | GF Refined Sky Masks         | [sequence path]/sky_mask/default_refined/img
# |  | generated | Volume GF Refined Sky Masks  | [sequence path]/sky_mask/volume_refined/img
#               [DEPTH MAPS]
# |  | generated | MonoDepth2 Depth Map         | [sequence path]/depth/monodepth/[img|mat]
# |  | generated | ManyDepth Depth Map          | [sequence path]/depth/manydepth/[img|mat]
# |  | generated | LEAStereo Depth Map          | [sequence path]/depth/leastereo/[img|mat]
# |  | generated | Augmented Depth Map          | [sequence path]/depth/augmented/[img|mat]


## Component initializers
def init_skyar(args):
    print("SkyAR         > ","Initializing...")
    if not os.path.exists(args.skyar_path):
        print("SkyAR         > ","FATAL ERROR: SkyAR repository not present. Please ensure SkyAR is cloned into the skyar_path ({})".format(args.skyar_path))
        exit(0);
    from simplewrappers import SimpleSkyFilter
    skyar = SimpleSkyFilter(args)
    
    return skyar

##def init_sky_mask_volume_refine(args):
##    
##    pass

def init_monodepth(args):
    print("MonoDepth2    > ","Initializing...")
    if not os.path.exists(args.mono2_path):
        print("MonoDepth2    > ","FATAL ERROR: MonoDepth2 repository not present. Please ensure MonoDepth2 is cloned into the mono2_path ({})".format(args.mono2_path))
        exit(0);
    from simplewrappers import SimpleMonoDepth2
    mono = SimpleMonoDepth2(args)
    
    return mono

def init_manydepth(args):
    print("ManyDepth     > ","Initializing...")
    if not os.path.exists(args.manyd_path):
        print("ManyDepth     > ","FATAL ERROR: ManyDepth repository not present. Please ensure ManyDepth is cloned into the manyd_path ({})".format(args.manyd_path))
        exit(0);
    from simplewrappers import SimpleManyDepth
    manyd = SimpleManyDepth(args)
    
    return manyd

def init_leastereo(args):
    print("LEAStereo     > ","Initializing...")
    if not os.path.exists(args.lea_path):
        print("ManyDepth     > ","FATAL ERROR: LEAStereo repository not present. Please ensure LEAStereo is cloned into the lea_path ({})".format(args.lea_path))
        exit(0);
    from simplewrappers import SimpleLEAStereo
    lea = SimpleLEAStereo(args)
    
    return lea

def init_haze(args):
    return None

# Component execution
# should be executed in this order
def run_skyar(instance, args):
    if args.override:
        args.current_missing['sky masks default'] = 'all'
        args.current_missing['sky masks raw'] = 'all'


    if args.current_missing['sky masks default'] == 'none' and args.current_missing['sky masks raw'] == 'none':
        print("SkyAR         > ", "Skipping {}: Already complete".format(args.current_folder))
        return

    image_dir = os.path.join(args.current_path,"image_02","data")
    output_dir = os.path.join(args.current_path,"sky_mask")

    os.makedirs(os.path.join(args.current_path,"sky_mask","raw","img"),exist_ok=True)
    os.makedirs(os.path.join(args.current_path,"sky_mask","default_refined","img"),exist_ok=True)
    
    if (args.current_missing['sky masks default'] == 'none' and not args.current_missing['sky masks raw'] == 'none') or (not args.current_missing['sky masks default'] == 'all' and args.current_missing['sky masks raw'] == 'all'):
        miss = args.current_missing['sky masks raw']
    elif (not args.current_missing['sky masks default'] == 'none' and args.current_missing['sky masks default'] == 'none') or (args.current_missing['sky masks default'] == 'all' and not args.current_missing['sky masks raw'] == 'all'):
        miss = args.current_missing['sky masks default']
    elif args.current_missing['sky masks default'] == 'all' and args.current_missing['sky masks raw'] == 'all':
        miss = 'all'
    else:
        miss = args.current_missing['sky masks raw'].union(args.current_missing['sky masks default'])

    if miss=='all':
        print("SkyAR         > ","Predicting on {:d} images in {}".format(len(os.listdir(image_dir)),args.current_folder))
        instance.run_full(image_dir, output_dir)
    else:
        print("SkyAR         > ","Predicting on {:d} of {:d} images in {}".format(len(miss),len(os.listdir(image_dir)),args.current_folder))
        instance.run_partial(image_dir, output_dir, miss)

##def run_sky_mask_volume_refine(instance, args):
##    pass

def run_monodepth(instance, args):
    miss = args.current_missing['depth maps mono']
    
    if args.override:
        miss = 'all'

    if miss == 'none':
        print("MonoDepth2    > ", "Skipping {}: Already complete".format(args.current_folder))
        return

    image_dir = os.path.join(args.current_path,"image_02","data")
    output_dir = os.path.join(args.current_path,"depth","monodepth")

    os.makedirs(os.path.join(args.current_path,"depth","monodepth","mat"),exist_ok=True)    
    if not args.skip_heatmaps: os.makedirs(os.path.join(args.current_path,"depth","monodepth","img"), exist_ok=True)


    if miss=='all':
        print("MonoDepth2    > ","Predicting on {:d} images in {}".format(len(os.listdir(image_dir)),args.current_folder))
        instance.run_full(image_dir, output_dir)
    else:
        print("MonoDepth2    > ","Predicting on {:d} of {:d} images in {}".format(len(miss),len(os.listdir(image_dir)),args.current_folder))
        instance.run_partial(image_dir, output_dir, miss)


def run_manydepth(instance, args):
    miss = args.current_missing['depth maps many']
    
    if args.override:
        miss = 'all'

    if miss == 'none':
        print("ManyDepth     > ", "Skipping {}: Already complete".format(args.current_folder))
        return

    image_dir = os.path.join(args.current_path,"image_02","data")
    output_dir = os.path.join(args.current_path,"depth","manydepth")

    os.makedirs(os.path.join(args.current_path,"depth","manydepth","mat"),exist_ok=True)    
    if not args.skip_heatmaps: os.makedirs(os.path.join(args.current_path,"depth","manydepth","img"), exist_ok=True)

    if miss=='all':
        print("ManyDepth     > ","Predicting on {:d} images in {}".format(len(os.listdir(image_dir)),args.current_folder))
        instance.run_full(image_dir, output_dir, args.current_drive)
    else:
        print("ManyDepth     > ","Predicting on {:d} of {:d} images in {}".format(len(miss),len(os.listdir(image_dir)),args.current_folder))
        instance.run_partial(image_dir, output_dir, args.current_drive, miss)

def run_leastereo(instance, args):
    miss = args.current_missing['depth maps many']
    
    if args.override:
        miss = 'all'

    if miss == 'none':
        print("LEAStereo     > ", "Skipping {}: Already complete".format(args.current_folder))
        return

    image_dir = os.path.join(args.current_path,"image_02","data")
    output_dir = os.path.join(args.current_path,"depth","leastereo")

    os.makedirs(os.path.join(args.current_path,"depth","leastereo","mat"),exist_ok=True)    
    if not args.skip_heatmaps: os.makedirs(os.path.join(args.current_path,"depth","leastereo","img"), exist_ok=True)

    if miss=='all':
        print("LEAStereo     > ","Predicting on {:d} images in {}".format(len(os.listdir(image_dir)),args.current_folder))
        instance.run_full(image_dir, output_dir, args.current_drive)
    else:
        print("LEAStereo     > ","Predicting on {:d} of {:d} images in {}".format(len(miss),len(os.listdir(image_dir)),args.current_folder))
        instance.run_partial(image_dir, output_dir, args.current_drive, miss)    

def run_haze(instance, args):
    seed = args.haze_seed ^ (int(hashlib.sha1(args.current_folder.encode("utf-8")).hexdigest(), 16) % (2**32))
    seed = np.array(seed).astype(np.uint32)
    np.random.seed(seed)
    vis = np.random.uniform(low=0.1,high=1)
    a0 = np.random.uniform(low=0.7,high=1)
    a = np.array([a0,a0,a0])

    out_dir = os.path.join(args.current_path,"haze")
    os.makedirs(out_dir,exist_ok=True)
    os.makedirs(os.path.join(out_dir,"image"),exist_ok=True)
    os.makedirs(os.path.join(out_dir,"transmission"),exist_ok=True)

    props = {}
    props['visibility'] = vis
    props['A'] = a.tolist()
    props['depth_source'] = args.haze_depth_type
    props['mask_type'] = args.haze_mask_type
    props['gen_timestamp'] = datetime.now().isoformat()
    
    with open(os.path.join(out_dir,"props.json"),'w') as f:
              json.dump(props,f)

    if args.haze_depth_type=="mono2":
        depth_dir = os.path.join(args.current_path, paths['depth maps mono'])
        scale=1
    elif args.haze_depth_type=="manyd":
        depth_dir = os.path.join(args.current_path, paths['depth maps many'])
        scale=6
    else:
        depth_dir = os.path.join(args.current_path, paths['depth maps lea'])
        scale=1000

    if args.haze_mask_type!="none":
        mask_dir = os.path.join(args.current_path, paths['sky masks default'])
    

    image_dir = os.path.join(args.current_path,"image_02","data")
    print("HazeGen     > ","Running on {:d} images in {}".format(len(os.listdir(image_dir)),args.current_folder))

    
    img_names = os.listdir(image_dir)
    img_names.sort()

    for img_name in img_names:
        
        img_path = os.path.join(image_dir,img_name)
        depth_path = os.path.join(depth_dir, img_name[:-4]+".mat")

        #out_path = os.path.join(out_dir, img_name)
        
        if args.haze_mask_type=="none":
            create_haze(img_path,depth_path,out_dir,img_name,a,vis,scale=scale)
        else:
            double = args.haze_mask_type=="double"
            create_haze(img_path,depth_path,out_dir,img_name,a,vis,mask_path=os.path.join(mask_dir,img_name),double_mask=double,depth_scale=scale)
            

    

helper = {
    'init':{'skyar':init_skyar,
##            'volrefine':init_sky_mask_volume_refine,
##            'maskdepth':init_depth_masker,
            'manyd':init_manydepth,
            'mono2':init_monodepth,
            'lea':init_leastereo,
            'haze':init_haze},
    
    'run':{'skyar':run_skyar,
##           'volrefine':run_sky_mask_volume_refine,
##           'maskdepth':run_depth_masker,
           'manyd':run_manydepth,
           'mono2':run_monodepth,
           'lea':run_leastereo,
           'haze':run_haze}
    }

paths = {}
paths['sky masks default']              = os.path.join('sky_mask','default_refined','img')
paths['sky masks raw']                  = os.path.join('sky_mask','raw','img')
paths['sky masks volume']               = os.path.join('sky_mask','volume_refined','img')
paths['depth maps mono']                = os.path.join('depth','monodepth','mat')
paths['depth maps many']                = os.path.join('depth','manydepth','mat')
paths['depth maps lea']                 = os.path.join('depth','lea','mat')
paths['depth maps augmented']           = os.path.join('depth','augmented','mat')
paths['lidar sky aug']                  = os.path.join('proj_depth','sky_aug')
paths['lidar full aug']                 = os.path.join('proj_depth','full_aug')

all_maps = list(paths.keys())

clamped_range = {m:False for m in all_maps}
clamped_range['depth maps augmented']           = True
clamped_range['lidar sky aug']                 = True
clamped_range['lidar full aug']                = True


def build_requirements_list(args):

    requirements = {m:False for m in all_maps}
    requirements['sky masks default']           = args.sky_masks or (args.haze and args.haze_mask_type!="none")
##    requirements['sky masks raw']                   = False
##    requirements['sky masks volume']                = False
    requirements['depth maps mono']    = args.mono2 or (args.haze and args.haze_depth_type=="mono2")
    requirements['depth maps many']    = args.manyd or (args.haze and args.haze_depth_type=="manyd")
    requirements['depth maps lea']    = args.lea or (args.haze and args.haze_depth_type=="lea")
##    requirements['depth maps augmented']  = False
##    requirements['lidar sky aug']                 = False
##    requirements['lidar full aug']                = False
    requirements['haze']                = args.haze
    
    missings = {}

    drives = os.listdir(args.path)
    for drive in drives:
        drive_path = os.path.join(data_path,drive)

        folders = os.listdir(drive_path)
        for folder in folders:
            folder_path = os.path.join(drive_path,folder)
            if os.path.isfile(folder_path): continue
            
            contents = os.listdir(folder_path)
            if 'image_02' not in contents:
                print("Images missing for {}".format(folder))
                continue

            folder_missings = {m:'all' for m in all_maps}

            img_names = os.listdir(os.path.join(folder_path,'image_02','data'))

            ids = set([int(x[:10]) for x in img_names])
            t = len(ids)
            clamped_ids = ids - {0,1,2,3,4, t-1, t-2, t-3, t-4, t-5}

            for m in all_maps:
                t_path = os.path.join(folder_path,paths[m])
                #print(t_path)
                if(os.path.exists(t_path)):
                    #print("Apple")
                    t_ids = set([int(x[:10]) for x in os.listdir(t_path)])
                    #print(t_ids)
                    idset = clamped_ids if clamped_range[m] else ids
                    miss = idset - t_ids
                    if not miss:    folder_missings[m] = 'none'
                    elif not t_ids: folder_missings[m] = 'all'
                    else:           folder_missings[m] = miss

            missings[folder] = folder_missings

        #print(folder_missings)


    #print(requirements)
    #print()
    fully_present = {m:False if args.override else reduce(lambda a, b: a and b, [v[m]=='none' for v in missings.values()], True) for m in all_maps}
    #print(fully_present)
    
    components = []
    # The order in which these are added matters.
    if (requirements['sky masks default'] and not fully_present['sky masks default']) or (requirements['sky masks raw'] and not fully_present['sky masks raw']):
        components.append('skyar')
##    elif args.sky_masks:
##        print("SkyAR > ", "Skipping as all sky masks are present")
    
    if (requirements['depth maps mono'] and not fully_present['depth maps mono']):
        components.append('mono2')
##    elif args.mono2:
##        print("MonoDepth2 > ", "Skipping as all MonoDepth2 depth maps present")
    
    if (requirements['depth maps many'] and not fully_present['depth maps many']):
        components.append('manyd')

    if (requirements['depth maps lea'] and not fully_present['depth maps lea']):
        components.append('lea')
##    elif args.manyd:
##        print("ManyDepth  > ", "Skipping as all ManyDepth depth maps present")
        

    if (args.haze):
        components.append('haze')

    return components, missings

if __name__ == '__main__':

    
    args = parser.parse_args()
    data_path = args.path

    if args.verbose:
        args.sky_masks_verbose = True
        args.mono2_verbose = True
        args.manyd_verbose = True
        args.lea_verbose = True
        args.haze_verbose = True
    
    components, missings = build_requirements_list(args)

    if not components:
        print("Nothing to generate")
        if args.clean:
            print("Cleaning only")
        else:
            exit(0)

    instances = {}
    for component in components:
        instances[component] = helper['init'][component](args)

    if args.clean:
        print("Cleaning will delete files, are you sure? [y]")
        if input().lower() != 'y':
            exit(0)

    #print(missings)
    
    drives = os.listdir(data_path)
    for drive in drives:
        drive_path = os.path.join(data_path,drive)
        args.current_drive = drive_path
        
        folders = os.listdir(drive_path)
        for folder in folders:
            folder_path = os.path.join(drive_path,folder)

            if os.path.isfile(folder_path): continue

            args.current_path = folder_path
            args.current_folder = folder
            args.current_missing = missings[folder]
            
            contents = os.listdir(folder_path)

            if args.clean:
                for content in contents:
                    if content in ['image_02','image_03','proj_depth','sky_mask','depth']: continue
                    #print("{} would be removed".format(os.path.join(folder_path,content)))
                    shutil.rmtree(os.path.join(folder_path,content))

            if not missings[folder]: continue
            
            for component in components:
                helper['run'][component](instances[component], args)



