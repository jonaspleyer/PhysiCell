binary_name=secretion_project

# File names and settings for individual video generation
framerate=20
output_file_svg=out_left.mp4
output_file_pic=out_right.mp4

# Combined videos/pictures
output_png=PhysiCell-PatternSimulation
output_combined=PhysiCell-PatternSimulation-combined

dir=$(pwd)
cd ..
latex_dir=$(pwd)/src/media/
cd simulation

# First we want to compile and run the simulation
cd build
make
cd ..
rm output/*.xml
rm output/*.svg
rm output/*.png
rm output/*.mat
./$binary_name

# Now generate the pictures
rm output/*.png
python main.py

# Now generate the movie for substrate plots
cd pictures
ffmpeg -y -threads 40 -r $framerate -f image2 -pattern_type glob -i '*.png' -vcodec libx264 -pix_fmt yuv444p -strict -2 -tune animation ../$output_file_pic
cp plot_substr_activator_00000000.png ..

# Now generate movie from PhysiCell files
# First convert svg to png files
cd ..
cd output
parallel -j 40 convert -size 842x901 {} {.}.png ::: $(ls snapshot*.svg)
ffmpeg -y -threads 40 -r $framerate -f image2 -pattern_type glob -i '*.png' -vcodec libx264 -pix_fmt yuv444p -strict -2 -tune animation $dir/$output_file_svg
cp snapshot00000000.png $dir
cd $dir

# Now combine both movies
cd $dir
ffmpeg -y -i out_left.mp4 -i out_right.mp4 -filter_complex hstack $output_combined.mp4

# Also create the montage as cover for the video embedding later in LaTeX
montage snapshot00000000.png plot_substr_activator_00000000.png -tile 2x1 -geometry +0+0 $output_combined.png

# Copy all files to the LaTeX folder
mv snapshot00000000.png $latex_dir/$output_png-left.png
mv plot_substr_activator_00000000.png $latex_dir/$output_png-right.png
cp $output_combined.png $output_combined.mp4 $latex_dir

# Now combine with the generated video in the latex dir
cd $latex_dir
cd generated
./generate_anim.sh
cd ..
./combine_anims.sh
