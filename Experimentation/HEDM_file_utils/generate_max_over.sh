#!/bin/bash
# HEDM utilities for Linux BASH.
#
# Harshad Paranjape (hparanja@mines.edu)
# MechanicsNext
#
# Generate a max-over-frames GE2/GE3 file from all files in a directory.
# Scroll to the bottom of the script to specify which folders to process.
# This file must have 'execute' permissions
# $>chmod 744 generate_max_over.sh
#
# Utility function. These need to be defined *before* they are used.
#
# Function: Max over all frames
function max_over_frames() {
    for var in "$@"
    do
        filename=$(basename "$var")
        extension="${filename##*.}"
        filename="${filename%.*}"

        if [ ${extension^^} == "TIF" ] || [ ${extension^^} == "TIFF" ]
        then
            echo "Processing : max over : $extension : $var"
            convert "$var" -evaluate-sequence max "$filename"-max."$extension"
        elif [ ${extension^^} == "GE2" ] || [ ${extension^^} == "GE3" ]
        then
            echo "Processing : max over : $extension : $var"
            convert -endian LSB -depth 16 -size 2048x2048+8192 gray:"$var" \
                -endian LSB -depth 16 -size 2048x2048 -evaluate-sequence max gray:"$filename"."$extension".max
            dd if="$filename"."$extension".max of="$filename"-max."$extension" obs=8192 seek=1
            rm -f "$filename"."$extension".max
        fi
    done
}
#
# Function: Combine multiple GE2
function combine_ge() {
    infile_counter=1

    for var in "$@"
    do
        filename=$(basename "$1")
        extension="${filename##*.}"
        filename="${filename%.*}"

        if [ $infile_counter -eq 1 ]
        then
            outname="$filename-combined.ge2"
            echo "Combining multiple inputs to $outname"
        fi

        dd bs=8192 skip=1 if="$var" of="$outname.tmp" oflag=append conv=notrunc
        infile_counter=$[ $infile_counter + 1 ]
    done

    dd if="$outname.tmp" of="$outname" obs=8192 seek=1
    rm -f "$outname.tmp"
}

# Calculate max-over-frames GE2/GE3 file from all files in a folder.
# All folders in FFDIRS will be processed. One max-over from each folder
# will be saved in the current directory.
FFDIRS=(10 12 14 16 18 29 31 35 37 39 50 52 54 56 58 61 63 65 67 69 77 79 81 83 85 95 97 99 101 103 106 108 110 112 114)

for II in "${FFDIRS[@]}"
do
    mkdir temp
    cd temp
    max_over_frames /FULL/PATH/TO/CHESS/DATA/$II/ff/*.ge2
    combine_ge *.ge2
    max_over_frames *max-combined.ge2
    mv *max-combined-max.ge2 ../
    cd ..
    rm -rf temp
done

