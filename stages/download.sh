#!/usr/bin/env bash

# Script to download data from Figshare

# Define the Figshare URL for the dataset
figshare_url="https://ndownloader.figshare.com/articles/5104873/versions/1"

# Get local path
localpath=$(pwd)
echo "Local path: $localpath"

# Create the download directory
downloadpath="$localpath/download"
echo "Download path: $downloadpath"
mkdir -p "$downloadpath"

# Change to the download directory
cd $downloadpath;

# Download the file from Figshare
echo "Downloading data from Figshare..."
wget -nc "$figshare_url" -O chemical_reactions.zip

echo "Download complete."
