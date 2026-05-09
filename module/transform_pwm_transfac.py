import sys
import argparse

def convert_homer_to_transfac(input_path, output_path):
    try:
        with open(input_path, 'r') as f_in, open(output_path, 'w') as f_out:
            f_out.write("ID  RNA_Motif\n")
            f_out.write("BF  Species\n")
            f_out.write("P0      A      C      G      U\n")
            line_count = 1
            for line in f_in:
                if line.startswith('>') or not line.strip():
                    continue
                
                parts = line.split()
                if len(parts) != 4:
                    continue
                
                vals = [float(x) * 100 for x in parts]
                f_out.write(f"{line_count:02d}  {vals[0]:>5.1f}  {vals[1]:>5.1f}  "
                            f"{vals[2]:>5.1f}  {vals[3]:>5.1f}\n")
                line_count += 1
            
            f_out.write("XX\n//\n")
        print(f"Successfully converted! Output path: {output_path}")
        
    except Exception as e:
        print(f"Failed to convert! Reason: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="convert homer motif matrix to the WebLogo TRANSFAC format")
    parser.add_argument("-i", "--input", required=True, help="input homer motif matrix(txt/motif)")
    parser.add_argument("-o", "--output", required=True, help="output TRANSFAC file path")
    
    args = parser.parse_args()
    convert_homer_to_transfac(args.input, args.output)
