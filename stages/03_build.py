# TODO: output to brick with two df parquets, change yaml too.
# TODO: get the reaction product thing working.
# TODO: blog post
import os
import pandas as pd
from tqdm import tqdm

def process_files(start_path, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for dirpath, dirnames, filenames in os.walk(start_path):
        for filename in tqdm(filenames, desc="Processing files", unit="file"):
            if filename.endswith('.rsmi'):
                file_path = os.path.join(dirpath, filename)
                try:
                    df = pd.read_csv(file_path, sep='\t', dtype={'ParagraphNum': str})
                    output_filename = os.path.join(output_dir, filename.replace('.rsmi', '.parquet'))
                    df.to_parquet(output_filename)
                    print(f"Data saved to {output_filename}")
                except pd.errors.ParserError as e:
                    print(f"Error reading file {file_path}: {e}")

process_files('raw', 'brick')

    