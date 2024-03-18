import pandas as pd
import xml.etree.ElementTree as ET
import sys
import pickle

def parse_xml(xml_file):
    tree = ET.parse(xml_file)
    root = tree.getroot()
    
    # Example of extracting data - you'll need to adapt this to your specific XML structure
    data = []
    for reaction in root.findall('.//reactionSmiles'):  # Adapt XPath according to your needs
        data.append({
            'reactionSmiles': reaction.text,
            # Extract other elements as needed
        })

    return pd.DataFrame(data)

def convert_to_parquet(input_xml, output_parquet):
    try:
        print(f"Converting file {input_xml}")
        df = parse_xml(input_xml)
        df.to_parquet(output_parquet)
        print(f"Conversion successful for {input_xml}")
    except Exception as e:
        print(f"Error processing file {input_xml}: {e}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python xml2parquet.py [input XML file] [output Parquet file]")
    else:
        in_file_name = sys.argv[1]
        out_file_name = sys.argv[2]
        convert_to_parquet(in_file_name, out_file_name)
