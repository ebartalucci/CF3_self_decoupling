##########################################
# CALCULATING VAN VLECK SECOND MOMENT M2 #
#  J. H. Van Vleck Phys. Rev. 74, 1168   #
# -------------------------------------- #
# Ettore Bartalucci, 21.12.2023, Aachen  #
##########################################

# Compute van Vleck equation: 
# For a powder made of crystallites of random orientations, we can average
# the (1 - 3cos^2theta_jk)^2 over all directions (see page 112, abragam)
# which leads to the following equation:
# M2 = 3/5 gamma^4 h_bar^2 I(I+1) SUM_k 1/r_jk^6

# Extract distances:
# Homonuclear 19F-19F distances are extracted from the respective crystal structures using
# Mercury and setting a cutoff to 10 Angstroms.
# All distances have been included for the F1, F2 and F3 for both inter and intramolecular contacts
# Files are:
# - distances_racemic_TLa.csv is the racemic file
# - distances_optical_pure_TLa.csv is the optical pure file

# Working with this script:
# STEP 1: convert .csv to .txt by prompting the path of the file to convert. Text file is written in the same folder
# STEP 2:

# Import modules
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import math


# Set constants

# STEP 1: Convert CSV to TXT file
def convert_csv_to_txt(csv_file, txt_file):
    """
    Convert a CSV file to a text file.
    :param csv_file: Input CSV file to convert
    :param txt_file: Output textfile
    """
    with open(csv_file, 'r') as csv_input, open(txt_file, 'w') as txt_output:
        for line in csv_input:
            # Replace semicolon with spaces
            line = line.replace(';',' ')
            # write to file
            txt_output.write(line)


# Compute M2: M2 = 3/5 gamma^4 h_bar^2 I(I+1) SUM_k 1/r_jk^6
def compute_M2(gamma,h_bar,I,r_jk):
    """
    Compute van Vleck second moment based on the given parameters and the homonuclear 19F-19F pair distances.
    
    :param gamma: gyromagnetic ratio in rad s^-1 T^-1
    :param h_bar: The reduced Planck constant 
    :param I: spin of resonant nuclei
    :param r_jk: An array of 19F-19F homonuclear distances (r_jk). In this case: cutoff to 10 Angstrom. 
    :return: The computed M2 value.
    """

    # handle eventual division by zero
    try:
        M2 = (3/5) * gamma**4 * h_bar**2 * I*(I+1) * sum(1 / r_jk**6 for r_jk in r_jk) #summation over pair distances

        return M2
    
    except ZeroDivisionError:
        raise ValueError("Division by zero detected! Please adjust input data")


def main():
    # define constants


    # STEP 1
    # input file
    input_csv = input("Enter the path to the input CSV file: ")
    
    # Use the same name as the input CSV file for the output text file
    base_name = os.path.splitext(os.path.basename(input_csv))[0]
    output_text = os.path.join(os.getcwd(), f'{base_name}.txt')
    
    convert_csv_to_txt(input_csv, output_text)

    print(f"Conversion complete. CSV file '{input_csv}' converted to text file '{output_text}'.")


    # Van Vleck second moment calculations
    gamma_19F = 251.815  # rad s^-1 T^-1
    h_bar = 1.05457266 * 10e-34
    I_19F = 1/2  # spin of resonant nuclei
    r_jk_homo = [1, 2, 3, 4, 1, 2, 3, 4] # to change with the actual values
    try: 
        second_moment = compute_M2(gamma_19F, h_bar, I_19F, r_jk_homo)
        print(f"The computed van Vleck Homonuclear 19F-19F second moment is: {second_moment}")
        
    except ValueError as error:
        print(f"Error: {error}")

if __name__ == "__main__":
    main()