##########################################
# CALCULATING VAN VLECK SECOND MOMENT M2 #
#  J. H. Van Vleck Phys. Rev. 74, 1168   #
# -------------------------------------- #
# Ettore Bartalucci, 06.11.2023, Aachen  #
##########################################

# import modules
import numpy as np
import matplotlib.pyplot as plt
import os
import math

# 19F-19F homonuclear distances have been calculated from the crystal structures 
# of the TLa racemic and optical pure 

# Cutoff 0 to 10 Angstrom
# All distances have been included for the F1, F2 and F3 for both inter and intramolecular contacts

# ----------- Section 1 ----------- #
# Function to extract 19F-19F distances from file
# molecules are flagged with '$'
def extract_distances(file_path):
    optical_pure_distances = []
    racemic_distances = []
    current_case = None

    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip() # strip file
            
            # Check if the line denotes a new case
            if line.startswith('$'): # flag for molecule
                current_case = line[1:]
                continue

            # Skip comments and empty lines
            if line.startswith('#') or not line: # flag for comment
                continue

            # Split the line into components (assuming space-separated values)
            pairs = line.split(':') # each distance pair features the ':'
            
            # Extracting data (assuming the format is Nucleus1-Nucleus2: Distance)
            nuclei, distance = pairs[0], float(pairs[1]) # pairs[0] is nuclei tuple and pairs[1] are distances

            # Determine the case and store the distance accordingly
            if current_case == 'optical_pure':
                optical_pure_distances.append((nuclei, distance))
            elif current_case == 'racemic':
                racemic_distances.append((nuclei, distance))

    return optical_pure_distances, racemic_distances


# ----------- Section 2 ----------- #
# Plot distance histograms for comparison between distance distributions
def plot_distance_histogram(optical_pure, racemic):
    # Extracting nuclei and distances for optical pure and racemic cases
    opt_nuclei, opt_distances = zip(*optical_pure)
    rac_nuclei, rac_distances = zip(*racemic)

    # Creating a dictionary to store distances for each nucleus pair
    distances_dict = {nuclei: {'optical_pure': [], 'racemic': []} for nuclei in set(opt_nuclei + rac_nuclei)}

    # Populating the dictionary with distances
    for nuclei, distances, case in zip([opt_nuclei, rac_nuclei], [opt_distances, rac_distances], ['optical_pure', 'racemic']):
        for nucleus, distance in zip(nuclei, distances):
            distances_dict[nucleus][case].append(distance)

    # Define colors for optical pure and racemic
    color_mapping = {'optical_pure': 'blue', 'racemic': 'red'}

    # Plotting horizontal bar charts in the same figure with pairwise colors
    for i, (nucleus, distances) in enumerate(distances_dict.items()):

        plt.barh(nucleus, distances['optical_pure'], color=color_mapping['optical_pure'], alpha=0.5, label='Optical Pure')
        plt.barh(nucleus, distances['racemic'], color=color_mapping['racemic'], alpha=0.5, label='Racemic')

    plt.title(f'19F-19F Intramolecular Distances')
    plt.xlabel('Distance')
    plt.ylabel('Nucleus Pair')
    plt.legend()
    plt.show()



# ----------- RUN CODE ----------- #
# Get current working directory
working_dir = os.getcwd()

# Get text file with distances
tla_distances = os.path.join(working_dir, '19F_distances_TLa.txt')

# Run distance extraction
optical_pure, racemic = extract_distances(tla_distances)

# Print values: Optical pure
print("Optical Pure Distances:")
for nuclei, distance in optical_pure:
    print(f"{nuclei}: {distance}")

# Print values: Racemic
print("\nRacemic Distances:")
for nuclei, distance in racemic:
    print(f"{nuclei}: {distance}")

# Plot distance histograms for nuclei tuples for comparison
plot_distance_histogram(optical_pure, racemic)

# Calculate van Vleck second moment (M2)

# Plot M2 histograms for comparison
