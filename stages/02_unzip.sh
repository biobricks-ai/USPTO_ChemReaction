#!/usr/bin/env bash
extract() {
    local sourcefile="$1"
    local destdir="$2"

    # Determine file type
    local filetype=$(file --brief --mime-type "$sourcefile")

    # Extract based on filetype
    case "$filetype" in
        application/zip|application/x-zip-compressed)
            unzip -o "$sourcefile" -d "$destdir"
            ;;
        application/x-7z-compressed)
            7z x "$sourcefile" -o"$destdir"
            ;;
        *)
            echo "Unsupported file type: $sourcefile"
            return
            ;;
    esac
}

# Function to find and extract nested archives
extract_nested() {
    local searchdir="$1"
    find "$searchdir" -mindepth 1 -type f \( -name '*.zip' -o -name '*.7z' \) -print0 | while IFS= read -r -d '' file; do
        extract "$file" "$(dirname "$file")"
    done
}

# Get local path
localpath=$(pwd)
echo "Local path: $localpath"

# Set download path
downloadpath="$localpath/download"
echo "Download path: $downloadpath"

# Create raw path
rawpath="$localpath/raw"
mkdir -p $rawpath
echo "Raw path: $rawpath"

# Main execution
primaryfile="chemical_reactions.zip"
primarypath="$downloadpath/$primaryfile"

# Extract primary file
if [[ -f "$primarypath" ]]; then
    extract "$primarypath" "$rawpath"
    extract_nested "$rawpath"
else
    echo "Primary file not found: $primarypath"
fi
