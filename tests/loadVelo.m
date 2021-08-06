function points = loadVelo(path, frame)
fid = fopen(sprintf('%s/velodyne_points/data/%010d.bin',path,frame),'rb');
points = fread(fid,[4 inf],'single')';
fclose(fid);
end