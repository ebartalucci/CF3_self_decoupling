

# Ettore Bartalucci
# Compute van Vleck equation 
# For a powder made of crystallites of random orientations, we can average
# the (1 - 3cos^2theta)^2 over all directions (see page 112, abragam)

import numpy as np

import numpy as np

def compute_homo_van_vleck(homonuclear_distances):
    
    # Define variables
    gamma_19F = 251.815  # rad s^-1 T^-1
    h_bar = 1.05457266 * 10e-34
    I_19F = 1/2  # spin of resonant nuclei

    # initialize M2
    M2 = 0

    for r_ij in homonuclear_distances:
        # Compute homonuclear van Vleck second moment
        M2 += 3/5 * gamma_19F**4 * h_bar**2 * I_19F * (I_19F + 1) * 1/r_ij

    # Return the computed value after the loop
    return M2


homonuclear_distances = [3.005, 3.602, 3.587, 3.682, 3.587, 4.690, 3.792, 3.005, 3.682]

result = compute_homo_van_vleck(homonuclear_distances)

print('second moments are:', result)

