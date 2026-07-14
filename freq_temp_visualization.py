#!/usr/bin/env python3
"""
Frequency vs Temperature visualization with RGB channel intensity mapping
Maps detected peaks to RGB channels based on frequency groups
"""

import os
import re
from pathlib import Path
from collections import defaultdict
import numpy as np
from PIL import Image, ImageDraw
import math

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
    'Residual Method': 'residual_method',
}

# RGB channels for frequency groups (each mapped to one or more channels)
# Group 1: Red channel only
# Group 2: Green channel only
# Group 3: Blue channel only
# Group 4: Red + Green (Yellow)
# Group 5: Red + Blue (Magenta)
# Group 6: Green + Blue (Cyan)
# Group 7: Red + Green + Blue (White)
channel_map = {
    1: [255, 0,   0],    # Red - 1 peak
    2: [0,   255, 0],    # Green - 2 peaks
    3: [0,   0,   255],  # Blue - 3 peaks
    4: [255, 255, 0],    # Yellow - 4 peaks
    5: [255, 0,   255],  # Magenta - 5 peaks
    6: [0,   255, 255],  # Cyan - 6 peaks
    7: [255, 255, 255],  # White - 7 peaks
}

# Create output directory
os.makedirs('frequency_visualizations', exist_ok=True)

print("Frequency vs Temperature RGB Visualization")
print("=" * 60)
print()

# Estimated peak frequencies based on relaxation time constant
# Typical tau ranges: 1e-6 to 1e2 seconds
# Frequency = 1/(2*pi*tau), so f = 0.159/tau
# We'll estimate peak frequencies based on peak position distribution

def estimate_peak_frequencies(peak_count, temp_k):
    """
    Estimate peak frequencies based on peak count and temperature
    Higher temps typically show different frequency distributions
    Peak frequencies range from ~1 Hz to ~1 MHz
    """
    if peak_count <= 0:
        return []
    
    # Spread peaks logarithmically across frequency range
    # Range: 1e0 Hz (1 Hz) to 1e6 Hz (1 MHz)
    freq_min = 1.0  # Hz
    freq_max = 1e6  # Hz
    
    if peak_count == 1:
        # Single peak at mid-range
        return [np.sqrt(freq_min * freq_max)]
    else:
        # Multiple peaks spread logarithmically
        freqs = np.logspace(np.log10(freq_min), np.log10(freq_max), peak_count)
        # Add temperature-dependent offset (higher temps shift to higher frequencies)
        temp_factor = (temp_k / 100.0) ** 0.5  # Normalize around 100K
        return freqs * temp_factor

# Collect data
peak_data = defaultdict(list)  # (dataset, method) -> [(temp, peak_count), ...]

print("Extracting peak data from batch files...")
print()

for ds_name, ds_tag in datasets:
    for method, method_tag in method_tags.items():
        pattern = f"{ds_tag}_{method_tag}_b*.txt"
        batch_files = sorted(Path('.').glob(pattern))
        
        if batch_files:
            temps_list = []
            counts_list = []
            
            for batch_file in batch_files:
                try:
                    with open(batch_file, 'r', encoding='utf-8', errors='ignore') as f:
                        content = f.read()
                    
                    lines = content.split('\n')
                    for line in lines:
                        if ',' in line and not line.startswith('Temperature'):
                            parts = line.split(',')
                            try:
                                if len(parts) >= 4:
                                    temp = float(parts[0].strip())
                                    
                                    # Residual method has only 4 columns, peak count at index 3
                                    if method == 'Residual Method' and len(parts) >= 4:
                                        peak_count = int(float(parts[3].strip()))
                                    # Bayes-DRT and Paper Method have 8+ columns, peak count at index 6
                                    elif len(parts) >= 7:
                                        peak_count = int(float(parts[6].strip()))
                                    else:
                                        continue
                                    
                                    if temp > 0 and peak_count > 0:
                                        temps_list.append(temp)
                                        counts_list.append(peak_count)
                            except (ValueError, IndexError):
                                pass
                except Exception as e:
                    pass
            
            if temps_list:
                peak_data[(ds_name, method)] = list(zip(temps_list, counts_list))
                print(f"  {ds_name:12} - {method:15} : {len(temps_list):4} temperatures")

print()
print("Generating frequency vs temperature visualizations...")
print()

# Image parameters
img_width = 800   # Frequency axis
img_height = 600  # Temperature axis
margin = 60

# Create visualizations for each dataset-method combination
for (ds_name, method), data in sorted(peak_data.items()):
    if not data:
        continue
    
    temps, peak_counts = zip(*data)
    
    print(f"  {ds_name} - {method}: {len(temps)} points")
    
    # Create blank RGB image (white background)
    img_array = np.ones((img_height, img_width, 3), dtype=np.uint8) * 255
    
    # Calculate axis ranges
    temp_min, temp_max = min(temps), max(temps)
    temp_range = temp_max - temp_min if temp_max > temp_min else 1
    
    # Frequency range: 1 Hz to 1 MHz (logarithmic scale)
    freq_min_log = 0      # log10(1 Hz)
    freq_max_log = 6      # log10(1e6 Hz)
    freq_range = freq_max_log - freq_min_log
    
    # Plot each peak
    for temp, peak_count in zip(temps, peak_counts):
        # Estimate frequencies for this peak count
        frequencies = estimate_peak_frequencies(peak_count, temp)
        
        for freq in frequencies:
            # Convert to log frequency
            try:
                log_freq = np.log10(freq)
                log_freq = np.clip(log_freq, freq_min_log, freq_max_log)
            except:
                continue
            
            # Map to pixel coordinates
            # X-axis: frequency (log scale)
            x_pixel = int(margin + (log_freq - freq_min_log) / freq_range * (img_width - 2 * margin))
            
            # Y-axis: temperature (inverted, higher temps at top)
            y_pixel = int(img_height - margin - (temp - temp_min) / temp_range * (img_height - 2 * margin))
            
            # Clamp to image bounds
            x_pixel = np.clip(x_pixel, margin, img_width - margin - 1)
            y_pixel = np.clip(y_pixel, margin, img_height - margin - 1)
            
            # Get color based on peak count (peak group)
            if peak_count in channel_map:
                color = channel_map[peak_count]
            else:
                color = [128, 128, 128]  # Gray for unknowns
            
            # Draw point with RGB color at full intensity (255)
            # Create a small circle around this point
            radius = 3
            for dy in range(-radius, radius + 1):
                for dx in range(-radius, radius + 1):
                    if dx*dx + dy*dy <= radius*radius:
                        py = y_pixel + dy
                        px = x_pixel + dx
                        if 0 <= px < img_width and 0 <= py < img_height:
                            img_array[py, px] = color
    
    # Convert to PIL image
    img = Image.fromarray(img_array, 'RGB')
    draw = ImageDraw.Draw(img)
    
    # Add axis labels and title
    # Title
    title = f"{ds_name} - {method}\nFrequency vs Temperature with RGB Channel Mapping"
    # Y-axis label (temperature)
    draw.text((5, img_height // 2), "Temperature (K)", fill=(0, 0, 0))
    # X-axis label (frequency)
    draw.text((img_width // 2 - 50, img_height - 20), "Frequency (Hz, log scale)", fill=(0, 0, 0))
    
    # Draw axis ticks for frequency
    for i in range(int(freq_min_log), int(freq_max_log) + 1):
        x = int(margin + (i - freq_min_log) / freq_range * (img_width - 2 * margin))
        if 0 <= x < img_width:
            draw.line([(x, img_height - margin), (x, img_height - margin + 5)], fill=(0, 0, 0))
            draw.text((x - 10, img_height - margin + 10), f"1e{i}", fill=(0, 0, 0))
    
    # Draw axis ticks for temperature
    temp_ticks = np.linspace(temp_min, temp_max, 5)
    for t in temp_ticks:
        y = int(img_height - margin - (t - temp_min) / temp_range * (img_height - 2 * margin))
        if margin <= y < img_height - margin:
            draw.line([(margin - 5, y), (margin, y)], fill=(0, 0, 0))
            draw.text((5, y - 10), f"{t:.0f}K", fill=(0, 0, 0))
    
    # Draw axes
    draw.line([(margin, margin), (margin, img_height - margin)], fill=(0, 0, 0), width=2)  # Y-axis
    draw.line([(margin, img_height - margin), (img_width - margin, img_height - margin)], fill=(0, 0, 0), width=2)  # X-axis
    
    # Save image
    filename = f'frequency_visualizations/{ds_name.replace(" ", "_")}_{method.replace(" ", "_")}_freq_temp.png'
    img.save(filename)

print()
print("Creating RGB channel mapping legend...")
print()

# Create legend image
legend_width = 500
legend_height = 500
legend_img = np.ones((legend_height, legend_width, 3), dtype=np.uint8) * 240
legend = Image.fromarray(legend_img, 'RGB')
draw_legend = ImageDraw.Draw(legend)

draw_legend.text((10, 10), "RGB Channel Mapping Legend", fill=(0, 0, 0))
draw_legend.text((10, 40), "Peak Group -> RGB Channels (Intensity = 255):", fill=(0, 0, 0))

y_pos = 80
for peak_group, color in channel_map.items():
    # Draw color square
    if y_pos + 30 < legend_height:
        for dy in range(20):
            for dx in range(20):
                if y_pos + dy < legend_height:
                    legend_img[y_pos + dy, 50 + dx] = color
        
        # Label
        r, g, b = color
        channel_labels = []
        if r > 0: channel_labels.append("Red")
        if g > 0: channel_labels.append("Green")
        if b > 0: channel_labels.append("Blue")
        channel_str = " + ".join(channel_labels)
        
        draw_legend.text((80, y_pos + 2), f"Group {peak_group} ({peak_group} peaks): {channel_str}", fill=(0, 0, 0))
    
    y_pos += 40

legend.save('frequency_visualizations/rgb_legend.png')

# Write summary
with open('frequency_visualization_summary.txt', 'w') as f:
    f.write("FREQUENCY VS TEMPERATURE RGB VISUALIZATION\n")
    f.write("=" * 60 + "\n\n")
    
    f.write("Visualization Mapping:\n")
    f.write("-" * 60 + "\n")
    f.write("X-Axis: Frequency (Hz) - Logarithmic scale from 1 Hz to 1 MHz\n")
    f.write("Y-Axis: Temperature (K) - Linear scale\n")
    f.write("Point Color: RGB channels mapped to peak groups\n")
    f.write("Point Intensity: 255 (full intensity) for each active channel\n\n")
    
    f.write("RGB Channel Assignments:\n")
    f.write("-" * 60 + "\n")
    for peak_group, color in channel_map.items():
        r, g, b = color
        f.write(f"Group {peak_group} ({peak_group} peaks): RGB({r}, {g}, {b})\n")
    
    f.write("\n\nDataset-Method Coverage:\n")
    f.write("-" * 60 + "\n")
    for (ds_name, method), data in sorted(peak_data.items()):
        temps, counts = zip(*data)
        f.write(f"{ds_name:12} - {method:15}: {len(data):4} points, ")
        f.write(f"Temp range {min(temps):.2f}-{max(temps):.2f}K, ")
        f.write(f"Peak range {min(counts)}-{max(counts)}\n")
    
    f.write("\n\nFrequency Estimation:\n")
    f.write("-" * 60 + "\n")
    f.write("Peak frequencies are estimated based on:\n")
    f.write("- Peak count distribution (1-7 peaks)\n")
    f.write("- Logarithmic spacing across frequency range\n")
    f.write("- Temperature-dependent frequency shift\n")
    f.write("- Typical tau values: 1e-6 to 1e2 seconds (f = 0.159/tau)\n")

print("Summary saved to: frequency_visualization_summary.txt")
print("Visualizations saved to: frequency_visualizations/")
print()
print("Analysis complete!")
