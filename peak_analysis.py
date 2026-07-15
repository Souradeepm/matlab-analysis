#!/usr/bin/env python3
"""
Extract peak frequencies from batch result files and generate RGB visualizations
"""

import os
import re
from pathlib import Path
from collections import defaultdict
import numpy as np
from PIL import Image

# Define datasets
datasets = [
    ('S2022Sap', 's2022sap'),
    ('S2022Al', 's2022al'),
    ('S2222Sap', 's2222sap'),
    ('S2222Al', 's2222al'),
    ('S2302Sap', 's2302sap'),
    ('S2302Al', 's2302al'),
    ('S2322Sap', 's2322sap'),
    ('S2322Al', 's2322al'),
    ('S2332Sap', 's2332sap'),
    ('S2422Sap', 's2422sap'),
    ('S2422Al', 's2422al'),
]

methods = ['Bayes-DRT', 'Paper Method', 'Residual Method']
method_tags = {
    'Bayes-DRT': 'bayes_drt_matlab2011',
    'Paper Method': 'paper_method',
    'Residual Method': 'residual_method'
}

# Color palette for peak groups (RGB, normalized 0-1)
color_palette = np.array([
    [1.0, 0.0, 0.0],    # Red - 1 peak
    [0.0, 1.0, 0.0],    # Green - 2 peaks
    [0.0, 0.0, 1.0],    # Blue - 3 peaks
    [1.0, 1.0, 0.0],    # Yellow - 4 peaks
    [1.0, 0.0, 1.0],    # Magenta - 5 peaks
    [0.0, 1.0, 1.0],    # Cyan - 6 peaks
    [1.0, 0.5, 0.0],    # Orange - 7+ peaks
], dtype=np.float32)

# Create output directory
os.makedirs('peak_visualizations', exist_ok=True)

print("Peak Frequency Analysis and RGB Visualization")
print("=" * 50)
print()

# Collect peak statistics
peak_stats = defaultdict(list)  # (dataset, method) -> [(temp, peak_count), ...]
all_peak_counts = []

print("Scanning batch files for peak data...")
print()

for ds_name, ds_tag in datasets:
    # Scan for batch files for each method
    for method, method_tag in method_tags.items():
        pattern = f"{ds_tag}_{method_tag}_b*.txt"
        batch_files = sorted(Path('.').glob(pattern))
        
        if batch_files:
            print(f"  {ds_name} - {method}: ", end='')
            temp_count = 0
            
            for batch_file in batch_files:
                try:
                    with open(batch_file, 'r') as f:
                        content = f.read()
                        
                    # Find all temperature and peak count pairs
                    # Pattern: lines with format like "5.100000,1.00e-04,...,2,..."
                    lines = content.split('\n')
                    for line in lines:
                        if ',' in line and not line.startswith('Temperature'):
                            parts = line.split(',')
                            try:
                                # Column indices: 0=Temp, 6=PeakCount
                                if len(parts) >= 7:
                                    temp = float(parts[0].strip())
                                    peak_count = int(float(parts[6].strip()))
                                    
                                    if temp > 0 and peak_count > 0:
                                        peak_stats[(ds_name, method)].append((temp, peak_count))
                                        all_peak_counts.append(peak_count)
                                        temp_count += 1
                            except (ValueError, IndexError):
                                pass
                except Exception as e:
                    print(f"Error reading {batch_file}: {e}")
            
            if temp_count > 0:
                print(f"{temp_count} temperatures")
            else:
                print("0 temperatures")

print()
print("Peak Groups Identified:")
print("-" * 50)

# Identify unique peak groups
unique_peaks = sorted(set(all_peak_counts))
print(f"Peak counts found: {unique_peaks}")
print(f"Number of groups: {len(unique_peaks)}")
print()

# Generate visualizations
print("Generating visualizations...")
print()

for (ds_name, method), data in sorted(peak_stats.items()):
    if not data:
        continue
    
    temps, peaks = zip(*sorted(data, key=lambda x: x[0]))
    n_temps = len(temps)
    
    print(f"  {ds_name} - {method}: {n_temps} temperatures")
    
    # Create RGB image: each row is a temperature, colored by peak group
    img_height = n_temps
    img_width = max(10, len(unique_peaks) + 2)  # Minimum width
    
    # Create RGB image (PIL uses RGB format)
    rgb_array = np.ones((img_height, img_width, 3), dtype=np.uint8) * 240  # Light gray background
    
    for t in range(n_temps):
        peak_count = int(peaks[t])
        
        # Find color for this peak group
        group_idx = unique_peaks.index(peak_count) if peak_count in unique_peaks else -1
        
        if group_idx >= 0 and group_idx < len(color_palette):
            color = (color_palette[group_idx] * 255).astype(np.uint8)
        else:
            color = np.array([128, 128, 128], dtype=np.uint8)  # Gray for unknown
        
        # Fill the entire row with this color
        rgb_array[t, :] = color
    
    # Create image and save
    img = Image.fromarray(rgb_array, 'RGB')
    
    # Scale up for better visibility
    scale_factor = 2
    img = img.resize((img_width * scale_factor, img_height * scale_factor), Image.NEAREST)
    
    filename = f'peak_visualizations/{ds_name.replace(" ", "_")}_{method.replace(" ", "_")}_peak_dist.png'
    img.save(filename)

print()
print("Creating summary report...")
print()

# Write summary report
with open('peak_analysis_summary.txt', 'w') as f:
    f.write("PEAK FREQUENCY ANALYSIS AND RGB VISUALIZATION\n")
    f.write("=" * 50 + "\n\n")
    
    f.write("Peak Groups and RGB Mapping:\n")
    f.write("-" * 50 + "\n")
    for i, pc in enumerate(unique_peaks):
        if i < len(color_palette):
            rgb = color_palette[i]
            count = all_peak_counts.count(pc)
            f.write(f"Group {i+1}: {pc}-Peak Process\n")
            f.write(f"  RGB Color: ({rgb[0]:.1%}, {rgb[1]:.1%}, {rgb[2]:.1%})\n")
            f.write(f"  Occurrences: {count}\n\n")
    
    f.write("\nDataset-Method Peak Profiles:\n")
    f.write("-" * 50 + "\n\n")
    
    for ds_name, _ in datasets:
        f.write(f"\n{ds_name}:\n")
        for method in methods:
            key = (ds_name, method)
            if key in peak_stats and peak_stats[key]:
                _, peak_list = zip(*peak_stats[key])
                unique_in_ds = sorted(set(peak_list))
                f.write(f"  {method}: Peak counts = {unique_in_ds}\n")

print("Summary saved to: peak_analysis_summary.txt")
print(f"Visualizations saved to: peak_visualizations/")
print(f"Total dataset-method combinations: {len(peak_stats)}")
print(f"Total temperature records: {sum(len(v) for v in peak_stats.values())}")
print(f"Total unique peak groups: {len(unique_peaks)}")
print()
print("Analysis complete!")
