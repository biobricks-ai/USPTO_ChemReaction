#!/usr/bin/env bash

# Script to process different types of files and build output files

# Get local path
localpath=$(pwd)
echo "Local path: $localpath"

# Set raw path
rawpath="$localpath/raw"
echo "Raw path: $rawpath"

# Create output directory
outputpath="$localpath/output"
mkdir -p $outputpath
echo "Output path: $outputpath"

process_xml() {
    local infile=$1
    local filename=$(basename "${infile%.*}")
    local outfile="$outputpath/$filename.parquet"

    echo "Processing XML file: $infile"
    echo "Output file: $outfile"
    python3 stages/xml2parquet.py "$infile" "$outfile"
}

export -f process_xml

# Process XML files
find "$rawpath" -type f -name '*.xml' -print0 | xargs -0 -P14 -n1 -I{} bash -c 'process_xml "$@"' _ {}
