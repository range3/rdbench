import argparse
from pathlib import Path
import struct
import math
import numpy as np
from PIL import Image, ImageDraw
from tqdm import tqdm
# import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap


def read_dat(file: Path):
    ar = np.array(struct.unpack(f'{file.stat().st_size // 8}d', file.read_bytes()))
    return ar.reshape((int(math.sqrt(ar.shape[0])), -1))

def create_frame_image(data, cmap = "Blues", frame_number = None):
    """フレームデータを画像に変換し、カラーマップを適用"""
    # Normalize data for visualization
    normalized_data = (data - data.min()) / (data.max() - data.min())
    
    # Apply colormap using matplotlib
    if cmap == 'PastelReds':  # Custom pastel-like Reds colormap
        colors = [(1, 1, 1), (1, 0.8, 0.8), (1, 0, 0)]  # white -> pink -> red
        colormap = LinearSegmentedColormap.from_list("PastelReds", colors, N=256)
    else:
        colormap = cm.get_cmap(cmap)

    colored_data = colormap(normalized_data)  # This returns an (M, N, 4) array with RGBA values
    
    # Convert the RGBA colormap array to RGB (dropping the alpha channel)
    rgb_data = (colored_data[:, :, :3] * 255).astype(np.uint8)
    
    # Convert numpy array to RGB image
    img = Image.fromarray(rgb_data, 'RGB')
    
    # Draw the frame number on the image
    # draw = ImageDraw.Draw(img)
    # draw.text((10, 10), f'Frame {frame_number}', fill=(255, 255, 255))  # Use white color for text in RGB mode
    return img

def generate_gif(input_dir, output_file, initial_duration, frame_duration, cmap, verbose):
    input_files = sorted(list(Path(input_dir).glob('*.bin')))
    
    if not input_files:
        print(f"No binary files found in directory: {input_dir}")
        return
    
    # Initialize the list to store the images (frames)
    images = []
    
    if verbose:
        print(f"Processing {len(input_files)} files...")
    
    # Process each file
    for i, file in enumerate(tqdm(input_files, desc="Processing files", disable=not verbose)):
        data = read_dat(file)
        img = create_frame_image(data, cmap, i)
        images.append(img)
    
    # Save the images as a GIF
    if verbose:
        print(f"Saving GIF to {output_file}...")
    
    images[0].save(
        output_file, save_all=True, append_images=images[1:], 
        duration=[initial_duration] + [frame_duration] * (len(images) - 1),
        loop=0
    )

    if verbose:
        print(f"GIF saved successfully: {output_file}")

def main():
    parser = argparse.ArgumentParser(description="Create an animated GIF from binary data files.")
    parser.add_argument("input", type=str, help="Input directory containing .bin files")
    parser.add_argument("output", type=str, help="Output GIF file")
    parser.add_argument("--initial_duration", type=int, default=1000, help="Duration of the initial frame in milliseconds")
    parser.add_argument("--frame_duration", type=int, default=100, help="Duration of each frame in milliseconds")
    parser.add_argument("--cmap", type=str, default="Blues", help="Colormap to use for the frames (default: Blues)")
    parser.add_argument("--verbose", "-v", action="store_true", help="Enable verbose output with progress")

    args = parser.parse_args()

    generate_gif(args.input, args.output, args.initial_duration, args.frame_duration, args.cmap, args.verbose)

if __name__ == "__main__":
    main()
