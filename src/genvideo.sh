cd ../frames;
rename 's/\d+/sprintf("%05d",$&)/e' *.bmp
ffmpeg -y -framerate 10 -i %05d.bmp -c:v libx264 ../out.mp4
