import sys
from py7zr import unpack_7zarchive
import shutil

def main():
    source_file = sys.argv[1]
    destination_dir = sys.argv[2]
    shutil.register_unpack_format('7zip', ['.7z'], unpack_7zarchive)
    shutil.unpack_archive(source_file, destination_dir)

if __name__ == '__main__':
    main()
