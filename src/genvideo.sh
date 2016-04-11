#!/bin/bash
# script called if -v flag is present. Stitches together all frames after they have been rendered
cd ../frames;
echo -n "Generating video. This may take a while...  "
rename 's/\d+/sprintf("%05d",$&)/e' *.bmp
ffmpeg -y -framerate 30 -i %05d.bmp -c:v libx264 ../out.mp4 &> /dev/null
echo "done."