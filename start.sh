PROJECT_NAME="./project"
BASE_DIR=$(pwd)

# Clean old files and logs
# rm $PROJECT_NAME
rm output/*.xml
rm output/*.svg
rm output/*.png
rm output/*.mat

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
else
	echo "Error in Makefile: Do not continue"
fi

if (($ret == 0));
then
	# convert from svg to png in parallel
	parallel -j 40 convert -size 842x901 {} {.}.png ::: $(ls output/snapshot*.svg)

	# Create movie
	ffmpeg -y -threads 40 -r 15 -f image2 -pattern_type glob -i 'output/*.png' -vcodec libx264 -pix_fmt yuv444p -strict -2 -tune animation movie.mp4
else
	echo "Error in Execution of $PROJECT_NAME: Do not continue"
fi
