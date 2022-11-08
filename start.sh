PROJECT_NAME="./project"
BASE_DIR=$(pwd)
THREADS=16

# Clean old files and logs
rm $PROJECT_NAME
rm output/*.xml
rm output/*.svg
rm output/*.png
rm output/*.mat
rm cells.csv

# Build new program files
cd build
make
ret=$?
cd ..

if (($ret == 0));
then
	# Run PhysiCell
	$PROJECT_NAME
	ret=$?
	echo "Ok"
else
	echo "Error in Makefile: Do not continue"
fi

if (($ret == 0));
then
	# convert from svg to png in parallel
	parallel -j $THREADS convert -size 842x901 {} {.}.png ::: $(ls output/snapshot*.svg)

	# Create movie
	# ffmpeg -y -threads 40 -r 15 -f image2 -pattern_type glob -i 'output/*.png' -vcodec libx264 -pix_fmt yuv444p -strict -2 -tune animation movie.mp4
	# ffmpeg -y -pattern_type glob -i 'output/*.png' -c:v libx264 -b:v 15000k -minrate 5000k -maxrate 8000k -pass 1 -pix_fmt yuv444p -an -f mp4 /dev/null && \
	# ffmpeg -pattern_type glob -i 'output/*.png' -c:v libx264 -b:v 15000k -minrate 5000k -maxrate 8000k -pass 2 -pix_fmt yuv444p -c:a aac -b:a 192k -movflags faststart output.mp4
	ffmpeg -y -pattern_type glob -i 'output/*.png' -c:v libx264 -b:v 15000k -minrate 5000k -maxrate 15000k -pass 1 -pix_fmt yuv422p -an -f mp4 /dev/null && \
	ffmpeg -pattern_type glob -i 'output/*.png' -c:v libx264 -b:v 15000k -minrate 20000k -maxrate 30000k -pass 2 -pix_fmt yuv422p -c:a aac -b:a 192k -movflags faststart output.mp4
else
	echo "Error in Execution of $PROJECT_NAME: Do not continue"
fi
