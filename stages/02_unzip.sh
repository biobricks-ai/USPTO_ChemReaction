#!/usr/bin/env bash
extract() {
    local sourcefile="$1" destdir="$2"
    local filetype=$(file --brief --mime-type "$sourcefile")

    case "$filetype" in
        application/zip|application/x-zip-compressed) unzip -o "$sourcefile" -d "$destdir" ;;
        application/x-7z-compressed) python3 extract_7z.py "$sourcefile" "$destdir" ;;
        *) echo "Unsupported file type: $sourcefile" ; return ;;
    esac
}

extract_nested() {
    local searchdir="$1"
    find "$searchdir" -type f \( -name '*.zip' -o -name '*.7z' \) -print0 |
        while IFS= read -r -d '' file; do extract "$file" "$(dirname "$file")"; done
}

mkdir -p raw
extract download/chemical_reactions.zip raw
extract_nested raw
