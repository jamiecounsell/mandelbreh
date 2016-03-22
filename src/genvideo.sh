#!/bin/bash
cd ../frames;
echo -n "Generating video. This may take a while...  "
rename 's/\d+/sprintf("%05d",$&)/e' *.bmp
ffmpeg -y -framerate 10 -i %05d.bmp -c:v libx264 ../out.mp4 &> /dev/null
echo "done."