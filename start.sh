ODIR_VALIDITY="../Validity_Testing/PDE-Single-Point-Source/output"
PROJECT_NAME="./secretion_project"
BASE_DIR=$(pwd)

# Clean old files and logs
rm $PROJECT_NAME
rm output/logs*

# Build new program files
cd build
make
cd ..

# Run PhysiCell
$PROJECT_NAME

#cd build
#make movie
#cd ..
#vlc output/out.mp4
