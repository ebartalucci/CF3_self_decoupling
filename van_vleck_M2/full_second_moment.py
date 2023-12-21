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
import pandas as pd
import seaborn as sns


# 19F-19F homonuclear distances have been calculated from the crystal structures 
# of the TLa racemic and optical pure 

# Cutoff 0 to 10 Angstrom
# All distances have been included for the F1, F2 and F3 for both inter and intramolecular contacts

# ------ Section 1 ------ #
def extract_correlated_distances(file_path):
    
    # Read the CSV file into a DataFrame
    df = pd.read_csv(file_path, delimiter=':').rename(columns=lambda x: x.strip())
    
    # Check if everything is alright
    print(df)

    # Create a dictionary to store correlated distances for each nucleus in 'Atom 1'
    correlated_distances = {}

    # Iterate through the DataFrame and populate the dictionary
    for index, row in df.iterrows():
        atom_1 = row['Atom 1']
        atom_2 = row['Atom 2']
        distance = row['d 1,2 [A]']

        # Check if the atom_1 is already in the dictionary, if not, add it
        if atom_1 not in correlated_distances:
            correlated_distances[atom_1] = []

        # Add the correlated distance for atom_1 and atom_2
        correlated_distances[atom_1].append((atom_2, distance))

    return correlated_distances


# ------ Section 2 ------ #
def plot_distance_distribution(df):
    # Get unique nuclei in 'Atom 1' column
    nuclei_list = df['Atom 1'].unique()

    # Set up subplots
    num_plots = len(nuclei_list)
    num_cols = 2  # You can adjust the number of columns as per your preference
    num_rows = (num_plots + num_cols - 1) // num_cols

    fig, axes = plt.subplots(num_rows, num_cols, figsize=(15, 5 * num_rows))

    # Flatten the axes array for easy iteration
    axes = axes.flatten()

    # Plot for each nucleus
    for i, nucleus in enumerate(nuclei_list):
        # Filter rows where 'Atom 1' is equal to the specified nucleus
        nucleus_df = df[df['Atom 1'] == nucleus]

        # Extract distances from the 'd 1,2 [A]' column
        distances = nucleus_df['d 1,2 [A]'].str.replace(',', '.').astype(float)

        # Plot the frequency distribution
        sns.histplot(distances, bins=20, kde=True, color='skyblue', ax=axes[i])
        axes[i].set_title(f'Distance Distribution for {nucleus}')
        axes[i].set_xlabel('Internuclear Distance (Angstroms)')
        axes[i].set_ylabel('Frequency')

    # Adjust layout
    plt.tight_layout()
    plt.show()



# ------ Section 3 ------ #
# compute homonuclear van Vleck

# ------ Section 4 ------ #
# plot comparison of homonuclear van vleck as histogram



# ------ Run code ------ #

# Get current working directory
working_dir = os.getcwd()

# Get text file with distances
tla_distances = os.path.join(working_dir, 'distances_racemic_TLa.csv')

result = extract_correlated_distances(tla_distances)

# Print the result
for nucleus, distances in result.items():
    print(f"Nucleus {nucleus}: {distances}")


df = pd.read_csv(tla_distances, delimiter=':').rename(columns=lambda x: x.strip())
plot_distance_distribution(df, 'F1')