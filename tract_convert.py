import argparse
import nibabel as nib

def convert_trk_to_tck(trk_file, tck_file):
    # Load the .trk file
    trk_tractogram = nib.streamlines.load(trk_file)

    # Save the Tractogram object as a .tck file
    nib.streamlines.save(trk_tractogram.tractogram, tck_file)
    print(f"Converted {trk_file} to {tck_file}")

def main():
    parser = argparse.ArgumentParser(description="Convert .trk tractography file to .tck format.")
    parser.add_argument('trk_file', type=str, help='Path to the .trk file')
    parser.add_argument('tck_file', type=str, help='Path to save the .tck file')
    args = parser.parse_args()

    convert_trk_to_tck(args.trk_file, args.tck_file)

if __name__ == "__main__":
    main()