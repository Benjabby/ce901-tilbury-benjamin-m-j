import cv2
import numpy as np
import os
import sys
import matlab
import tensorflow as tf
import time

import traceback


try:

    matlab_out = sys.stdout
    sys.stdout = sys.stderr

    model_type = sys.argv[1]

    print("######### SHELL WINDOW FOR MATLAB MANAGED PYTHON NEURAL NETWORK - TYPE " + model_type + " #########")
    print("DO NOT CLOSE\n\n")


    this_path = os.path.dirname(os.path.abspath(__file__))

    if model_type == "AOD":
        import torch
        import torch.nn.parallel
        import torchvision.transforms as transforms

        from AOD import net

        
        device = torch.device('cuda') if torch.has_cuda else torch.device('cpu')
        model = net.dehaze_net()
        model.load_state_dict(torch.load(os.path.join(this_path,"AOD","model.pth"),map_location=device))

        if torch.has_cuda:
            model = model.cuda()

    else:
        from PMHLD.model.model import build_combine_model
        model = build_combine_model(this_path)
            

    while(True):
        # Get next command
        command = int.from_bytes(sys.stdin.buffer.read(1),sys.byteorder)

        # 0 - Exit
        # 1 - New Image Incoming.

        if command==0:
            exit(0)
        elif command==1:
            # Next two bytes are height of upcoming image.
            h = sys.stdin.buffer.read(2)
            h = int.from_bytes(h,sys.byteorder)
            # Next two bytes are width of upcoming image.
            w = sys.stdin.buffer.read(2)
            w = int.from_bytes(w,sys.byteorder)

            byte_size_in = h*w*3

            # Next 'byte_size_in' bytes are the image
            byte = sys.stdin.buffer.read(byte_size_in)
            img = np.frombuffer(byte,dtype=matlab.uint8);
            
            img = img.reshape(h,w,3)
            img = img/255
            
            t = time.time()
            if model_type == "AOD":
                img = torch.from_numpy(img).float()
                img = img.permute(2,0,1)
                img = img.unsqueeze(0)
                
                if torch.has_cuda: img = img.cuda()

                pred = model(img).data.cpu().numpy().squeeze().transpose((1,2,0))
                
            else:
                img = img.astype(matlab.double)
                
                img = cv2.resize(img, (640, 480), interpolation=cv2.INTER_CUBIC)
                pred = model.predict(img[None,...])
                pred = cv2.resize(pred[0], (w, h), interpolation=cv2.INTER_CUBIC)
                

            elapsed = time.time()-t
            
            pred = pred.astype(matlab.double) # just in case the model doesn't output np.float64 for whatever reason.

            elapsed = np.array(elapsed,dtype=matlab.single)
            matlab_out.buffer.write(elapsed.tobytes()) # Write the time taken
            matlab_out.buffer.flush()
        
            matlab_out.buffer.write(pred.tobytes(order='F')) # Write the output image
            matlab_out.buffer.flush()
    
except Exception as e:
    print("\nACTUAL ERROR\n")
    print(e)
    traceback.print_exc();
    time.sleep(100)
