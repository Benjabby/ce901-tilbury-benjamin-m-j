path = "E:\dataset\2011_09_26\2011_09_26_drive_0009_sync\haze\image";

IPath = path;
ISize = length(dir(IPath))-2;

L = lines.S2011_09_26_drive_0009_sync;

Is = zeros(375,1242,3,ISize);

for frameID = 0:ISize-1
    framePath = num2str(frameID,'%010.f')+".png";
    I = imread(fullfile(IPath,framePath));
    I = im2double(I);

    Is(:,:,:,frameID+1) = I;
    
    Is(L(frameID+1),:,:,frameID+1) = repmat(reshape([1,0,1], 1,1,3,1),1,1242,1,1);
    
end