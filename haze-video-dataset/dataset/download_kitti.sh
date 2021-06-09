#!/bin/bash


#KITTI download script from https://github.com/pean1128/KITTI-download-scipt

files=(2011_09_26_calib.zip
2011_09_26_drive_0001
2011_09_26_drive_0009
2011_09_26_drive_0015
2011_09_26_drive_0022
2011_09_26_drive_0028
2011_09_26_drive_0036
2011_09_26_drive_0051
2011_09_26_drive_0061
2011_09_26_drive_0064
2011_09_26_drive_0093
2011_09_26_drive_0096
2011_09_26_drive_0101
2011_09_26_drive_0107
2011_09_28_calib.zip
2011_09_28_drive_0002
2011_09_29_calib.zip
2011_09_29_drive_0071
2011_09_30_calib.zip
2011_09_30_drive_0028
2011_10_03_calib.zip
2011_10_03_drive_0042
2011_10_03_drive_0047)
 

########filter out the downloaded files#######
#rootpath=/dataset/
dir=$(ls $1)
for dirname in $dir; do  
	if [ -d $dirname ]
	then
   		 #filepath=$rootpath$dirname
		 filepath=$dirname
	fi
	subfiles=$(ls $filepath)
	for f in $subfiles; do
		if [ ${f:(-3)} != "txt" ]
		then
			#array=$array" "${f:0:21}
			array+=${f:0:21}
		else
			continue
		fi
	done
done	
echo $array

time=0
for arr in ${array[@]}; do
	time=$[time+1]
	#echo $arr
done
echo $time
	
for f in ${files[@]}; do
	echo $array | grep -q $f 
	if [ $? -eq 0 ]
	then
		echo $"file "$f$" has been downloaded"
		continue
	#fi
	#if [ $? -ne 0 ]
	#then
	else
		echo $f
		if [ ${f:(-3)} != "zip" ]
       		then
                	shortname=$f'_sync.zip'
                	fullname=$f'/'$f'_sync.zip'
        	else
                	shortname=$f
                	fullname=$f
        	fi
		echo "Downloading: "$shortname
        	wget -c 'https://s3.eu-central-1.amazonaws.com/avg-kitti/raw_data/'$fullname
        	unzip -o $shortname
        	rm $shortname
	fi
done