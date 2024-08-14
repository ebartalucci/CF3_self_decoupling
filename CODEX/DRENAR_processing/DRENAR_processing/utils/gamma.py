"""
List of gyromagnetic ratios of common NMR nuclei
Unit in rad*s^-1*T^-1

Source: scipy.constants.physical_constants; Harris et al., Pure. Appl. Chem., 2001, 73, 1795-1818 (https://publications.iupac.org/pac/2001/pdf/7311x1795.pdf)
Author: Yufei Wu
Date: 2024-01-25
"""

from scipy.constants import physical_constants

e = physical_constants['electron gyromag. ratio'][0]
n = physical_constants['neutron gyromag. ratio'][0]
p = physical_constants['proton gyromag. ratio'][0]

H = 267522128
D = 41066279.1
C = 67282840
Li7 = 103977013
B11 = 85847044
N14 = 19337792
N = -27126180.4
O = -36280800
F = 251814800
Na = 70808493
Al = 69762715
Si = -53190000
P = 108394000
Fe = 8680624
Cu63 = 71117890
Zn = 16766880
Pt = 58385000
Pb = 55804600