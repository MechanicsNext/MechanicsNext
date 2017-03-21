# HEDM utilities for Linux BASH.
# Include this file in your ~/.bashrc by adding the following line in ~/.bashrc
# source ~/bash-hedm-utils.sh
#
# Harshad Paranjape (hparanja@mines.edu)
# MechanicsNext
#
# Convert TIFF/PNG etc to GE2 (single or multi frame, works on multiple files)
function tiff2ge() {
    for var in "$@"
    do
	echo "Processing : tiff2ge : $var"
        filename=$(basename "$var")
        extension="${filename##*.}"
        filename="${filename%.*}"
        convert "$var" -endian LSB -depth 16 -size 2048x2048 gray:"$filename".ge2.tmp
        dd if="$filename".ge2.tmp of="$filename".ge2 obs=8192 seek=1
        rm -f "$filename.ge2.tmp"
    done
}

# Convert GE2 image to tiff
function ge2tiff() {
    for var in "$@"
    do
	echo "Processing : ge2tiff : $var"
        filename=$(basename "$var")
        extension="${filename##*.}"
        filename="${filename%.*}"
        convert -endian LSB -depth 16 -size 2048x2048+8192 gray:"$var" "$filename".tiff
    done
}

# Max over all frames
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

# Sum over all frames
function sum_over_frames() {
    for var in "$@"
    do
        filename=$(basename "$var")
        extension="${filename##*.}"
        filename="${filename%.*}"

        if [ ${extension^^} == "TIF" ] || [ ${extension^^} == "TIFF" ]
        then
            echo "Processing : sum over : $extension : $var"
            convert "$var" -evaluate-sequence add "$filename"-sum."$extension"
        elif [ ${extension^^} == "GE2" ] || [ ${extension^^} == "GE3" ]
        then
            echo "Processing : sum over : $extension : $var"
            convert -endian LSB -depth 16 -size 2048x2048+8192 gray:"$var" \
                -endian LSB -depth 16 -size 2048x2048 -evaluate-sequence add gray:"$filename"."$extension".sum
            dd if="$filename"."$extension".sum of="$filename"-sum."$extension" obs=8192 seek=1
            rm -f "$filename"."$extension".sum
        fi
    done
}

# Mean over all frames
function mean_over_frames() {
    for var in "$@"
    do
        filename=$(basename "$var")
        extension="${filename##*.}"
        filename="${filename%.*}"

        if [ ${extension^^} == "TIF" ] || [ ${extension^^} == "TIFF" ]
        then
            echo "Processing : max over : $extension : $var"
            convert "$var" -evaluate-sequence mean "$filename"-mean."$extension"
        elif [ ${extension^^} == "GE2" ]
        then
            echo "Processing : max over : $extension : $var"
            convert -endian LSB -depth 16 -size 2048x2048+8192 gray:"$var" \
                -endian LSB -depth 16 -size 2048x2048 -evaluate-sequence mean gray:"$filename"."$extension".mean
            dd if="$filename"."$extension".mean of="$filename"-mean."$extension" obs=8192 seek=1
            rm -f "$filename"."$extension".mean
        fi
    done
}

# Extract a specific frame from multi-frame image
function extract_frame() {

    filename=$(basename "$1")
    extension="${filename##*.}"
    filename="${filename%.*}"

    if [ ${extension^^} == "TIF" ] || [ ${extension^^} == "TIFF" ]
    then
        echo "Processing : extract frame $2 : $extension : $1"
        convert "$1"[$2] "$filename"_"$2"."$extension"
    elif [ ${extension^^} == "GE2" ]
    then
        echo "Processing : extract frame $2 : $extension : $1"
        convert -endian LSB -depth 16 -size 2048x2048+8192 gray:"$1"["$2"] \
            -endian LSB -depth 16 -size 2048x2048 -evaluate-sequence max gray:"$filename"_"$2"."$extension".tmp
        dd if="$filename"_"$2"."$extension".tmp of="$filename"_"$2"."$extension" obs=8192 seek=1
        rm -f "$filename"_"$2"."$extension".tmp
    fi 
}

# Split a multi-frame ge2/3 into chuncks of 240 frames max
function split_ge_240() {

    filename=$(basename "$1")
    extension="${filename##*.}"
    filename="${filename%.*}"

    TEMP=( $( ls -ln "$1" ) )
    filebytes=${TEMP[4]}
    framesfile=$((($filebytes - 8192) / (2048*2048*2)))

    counter=0
    filecounter="$2"
    while [ "$counter" -lt $framesfile ]; do
        echo "Writing frames $counter - `expr $counter + 239` to ff_`printf %05d $filecounter`.$extension"
        counterend=`expr $counter + 239`
        convert -endian LSB -depth 16 -size 2048x2048+8192 gray:"$1"["$counter"-"$counterend"] \
            -endian LSB -depth 16 -size 2048x2048 gray:"ff_`printf %05d $filecounter`"."$extension".tmp
        dd if="ff_`printf %05d $filecounter`"."$extension".tmp of="ff_`printf %05d $filecounter`"."$extension" obs=8192 seek=1
        rm -f "ff_`printf %05d $filecounter`"."$extension".tmp
        counter=`expr $counter + 240`
        filecounter=`expr $filecounter + 1`
    done
}

# Subtract a dark frame from an image (not tested for multiframe)
function subtract_dark() {
    filename=$(basename "$1")
    extension="${filename##*.}"
    filename="${filename%.*}"

    if [ ${extension^^} == "TIF" ] || [ ${extension^^} == "TIFF" ]
    then
        echo "Processing : subtract dark : $extension : $2 from $1"
        convert $1 null: $2  -compose difference -layers composite  "$filename"-dark-subtracted."$extension"
        #composite -compose difference $1 $2 "$filename"-dark-subtracted."$extension"
    elif [ ${extension^^} == "GE2" ]
    then
        echo "Processing : subtract dark : $extension : $2 from $1"
        convert -endian LSB -depth 16 -size 2048x2048+8192 gray:$1 \
            -endian LSB -depth 16 -size 2048x2048+8192 null: gray:$2 \
            -compose difference -layers composite  -endian LSB -depth 16 -size 2048x2048+8192 "$filename"-dark-subtracted."$extension".tmp
        #composite -compose difference -endian LSB -depth 16 -size 2048x2048+8192 gray:"$1" \
        #    -endian LSB -depth 16 -size 2048x2048+8192 gray:"$2" \
        #    -endian LSB -depth 16 -size 2048x2048+8192  "$filename"-dark-subtracted."$extension".tmp
        dd if="$filename"-dark-subtracted."$extension".tmp of="$filename"-dark-subtracted."$extension" obs=8192 seek=1
        rm -f "$filename"-dark-subtracted."$extension".tmp
    fi
}

# Cmbine GE2
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
# Display file numbers from an ff directory
function show_ff_data() {
    # Find all directories with name ff and loop over them
    find "$@" -not -empty -name "ff" -print | while read f; do
	# Show directory name
	echo "$f"
	# Get size of all contents in the directory (in bytes)
	tot_byte_size=$(du -b "$f")
	# By default du prints the folder name too. Strip that.
	tot_byte_size=($tot_byte_size)
	tot_byte_size=${tot_byte_size[0]}
	# Estimate number of FF frames assuming one frame = 2048*2048*2 bytes
	num_frames_estimate=$( echo $tot_byte_size / 8388608 | bc )
	# Get size of all contents in a human readable unit (M, G, T)
	tot_file_size=$(du -h "$f")
	tot_file_size=($tot_file_size)
	# Print total size and frame number estimate
	echo Total size = ${tot_file_size[0]}, Estimated frames = $num_frames_estimate
	# List the contents of the folders
	filenames1=$(eval ls \'"$f"\' | tr '\n' ' ')
	filelist=''
	# Clean up the filenames assuming ff_%5d.ge2 format
	for filename in $filenames1
	do
	    filenum1=${filename:3:5}
	    filenum=$(echo $filenum1 | sed 's/^0*//')
	    filelist+=", $filenum"
	done
	filelist1=${filelist:2}
	# Print comma separated list of file numbers
	echo $filelist1
    done
}
