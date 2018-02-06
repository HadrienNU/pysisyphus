from pysisyphus.constants import ANG2BOHR


# Taken from periodictable-1.5.0
MASS_DICT = {
    'n': 14.0067, 'h': 1.00794, 'he': 4.002602, 'li': 6.941, 'be': 9.012182,
    'b': 10.811, 'c': 12.0107, 'o': 15.9994, 'f': 18.9984032, 'ne': 20.1797,
    'na': 22.98977, 'mg': 24.305, 'al': 26.981538, 'si': 28.0855, 'p': 30.973761,
    's': 32.065, 'cl': 35.453, 'ar': 39.948, 'k': 39.0983, 'ca': 40.078,
    'sc': 44.95591, 'ti': 47.867, 'v': 50.9415, 'cr': 51.9961, 'mn': 54.938049,
    'fe': 55.845, 'co': 58.9332, 'ni': 58.6934, 'cu':63.546, 'zn': 65.409,
    'ga': 69.723, 'ge': 72.64, 'as': 74.9216, 'se': 78.96, 'br': 79.904,
    'kr': 83.798,'rb': 85.4678, 'sr': 87.62, 'y': 88.90585, 'zr': 91.224,
    'nb': 92.90638, 'mo': 95.94, 'tc': 98, 'ru': 101.07, 'rh': 102.9055,
    'pd': 106.42, 'ag': 107.8682, 'cd': 112.411, 'in': 114.818, 'sn': 118.71,
    'sb': 121.76, 'te': 127.6, 'i': 126.90447, 'xe': 131.293, 'cs': 132.90545,
    'ba': 137.327, 'la': 138.9055, 'ce': 140.116, 'pr': 140.90765, 'nd': 144.24,
    'pm': 145, 'sm': 150.36, 'eu': 151.964, 'gd': 157.25, 'tb': 158.92534,
    'dy': 162.5, 'ho': 164.93032, 'er': 167.259, 'tm': 168.93421, 'yb': 173.04,
    'lu': 174.967, 'hf': 178.49, 'ta': 180.9479, 'w': 183.84, 're': 186.207,
    'os': 190.23, 'ir': 192.217, 'pt': 195.078, 'au': 196.96655, 'hg': 200.59,
    'tl': 204.3833, 'pb': 207.2, 'bi': 208.98038, 'po': 209, 'at': 210, 'rn': 222,
    'fr': 223, 'ra': 226, 'ac': 227, 'th': 232.0381, 'pa': 231.03588, 'u': 238.02891,
    'np': 237, 'pu': 244, 'am': 243, 'cm': 247, 'bk': 247, 'cf': 251, 'es': 252,
    'fm': 257, 'md': 258, 'no': 259, 'lr': 262, 'rf': 261, 'db': 262, 'sg': 266,
    'bh': 264, 'hs': 277, 'mt': 268, 'ds': 281, 'rg': 272, 'cn': 285, 'nh': 286,
    'fl': 289,'mc': 289, 'lv': 293, 'ts': 294, 'og': 294
}

# {el.symbol: el.covalent_radius for el  in periodictable.elements
#  if el.covalent_radius}
# In Angstroem
COVALENT_RADII_ANG = {
    'n': 0.71, 'h': 0.31, 'he': 0.28, 'li': 1.28, 'be': 0.96, 'b': 0.84,
    'c': 0.76, 'o': 0.66, 'f': 0.57, 'ne': 0.58, 'na': 1.66, 'mg': 1.41,
    'al': 1.21, 'si': 1.11, 'p': 1.07, 's': 1.05, 'cl': 1.02, 'ar': 1.06,
    'k': 2.03, 'ca': 1.76, 'sc': 1.7, 'ti': 1.6, 'v': 1.53, 'cr': 1.39,
    'mn': 1.39, 'fe': 1.32, 'co': 1.26, 'ni': 1.24, 'cu': 1.32, 'zn': 1.22,
    'ga': 1.22, 'ge': 1.2, 'as': 1.19, 'se': 1.2, 'br': 1.2, 'kr': 1.16,
    'rb': 2.2, 'sr': 1.95, 'y': 1.9, 'zr': 1.75, 'nb': 1.64, 'mo': 1.54,
    'tc': 1.47, 'ru': 1.46, 'rh': 1.42, 'pd': 1.39, 'ag': 1.45, 'cd': 1.44,
    'in': 1.42, 'sn': 1.39, 'sb': 1.39, 'te': 1.38, 'i': 1.39, 'xe': 1.4,
    'cs': 2.44, 'ba': 2.15, 'la': 2.07, 'ce': 2.04, 'pr': 2.03, 'nd': 2.01,
    'pm': 1.99, 'sm': 1.98, 'eu': 1.98, 'gd': 1.96, 'tb': 1.94, 'dy': 1.92,
    'ho': 1.92, 'er': 1.89, 'tm': 1.9, 'yb': 1.87, 'lu': 1.87, 'hf': 1.75,
    'ta': 1.7, 'w': 1.62, 're': 1.51, 'os': 1.44, 'ir': 1.41, 'pt': 1.36,
    'au': 1.36, 'hg': 1.32, 'tl': 1.45, 'pb': 1.46, 'bi': 1.48, 'po': 1.4,
    'at': 1.5, 'rn': 1.5, 'fr': 2.6, 'ra': 2.21, 'ac': 2.15, 'th': 2.06,
    'pa': 2.0, 'u': 1.96, 'np': 1.9, 'pu': 1.87, 'am': 1.8, 'cm': 1.69
}
# In Bohr
COVALENT_RADII = {
    'n': 1.3417, 'h': 0.5858, 'he': 0.5291, 'li': 2.4188, 'be': 1.8141,
    'b': 1.5874, 'c': 1.4362, 'o': 1.2472, 'f': 1.0771, 'ne': 1.096,
    'na': 3.1369, 'mg': 2.6645, 'al': 2.2866, 'si': 2.0976, 'p': 2.022,
    's': 1.9842, 'cl': 1.9275, 'ar': 2.0031, 'k': 3.8361, 'ca': 3.3259,
    'sc': 3.2125, 'ti': 3.0236, 'v': 2.8913, 'cr': 2.6267, 'mn': 2.6267,
    'fe': 2.4944, 'co': 2.3811, 'ni': 2.3433, 'cu': 2.4944, 'zn': 2.3055,
    'ga': 2.3055, 'ge': 2.2677, 'as': 2.2488, 'se': 2.2677, 'br': 2.2677,
    'kr': 2.1921, 'rb': 4.1574, 'sr': 3.685, 'y': 3.5905, 'zr': 3.307,
    'nb': 3.0992, 'mo': 2.9102, 'tc': 2.7779, 'ru': 2.759, 'rh': 2.6834,
    'pd': 2.6267, 'ag': 2.7401, 'cd': 2.7212, 'in': 2.6834, 'sn': 2.6267,
    'sb': 2.6267, 'te': 2.6078, 'i': 2.6267, 'xe': 2.6456, 'cs': 4.6109,
    'ba': 4.0629, 'la': 3.9117, 'ce': 3.855, 'pr': 3.8361, 'nd': 3.7983,
    'pm': 3.7606, 'sm': 3.7417, 'eu': 3.7417, 'gd': 3.7039, 'tb': 3.6661,
    'dy': 3.6283, 'ho': 3.6283, 'er': 3.5716, 'tm': 3.5905, 'yb': 3.5338,
    'lu': 3.5338, 'hf': 3.307, 'ta': 3.2125, 'w': 3.0614, 're': 2.8535,
    'os': 2.7212, 'ir': 2.6645, 'pt': 2.57, 'au': 2.57, 'hg': 2.4944,
    'tl': 2.7401, 'pb': 2.759, 'bi': 2.7968, 'po': 2.6456, 'at': 2.8346,
    'rn': 2.8346, 'fr': 4.9133, 'ra': 4.1763, 'ac': 4.0629, 'th': 3.8928,
    'pa': 3.7795, 'u': 3.7039, 'np': 3.5905, 'pu': 3.5338, 'am': 3.4015,
    'cm': 3.1936
}

# Taken from mendeleev 0.41
# CRC Handbook of Chem. and Phys 95th ed.
# for e in element(["H", "N", "O", "F", "P", "S", "Cl"]):
#   print(e.symbol, e.vdw_radius)
# pm -> Bohr
# {key: f"{value/100*ANG2BOHR:.4f}" for key, value in VDW_RADII.items()}
# VDW_RADII = { # In pm
#    "h": 110.0000,
#    "n": 155.0,
#    "o": 152.0,
#    "f": 147.0,
#    "p": 180.0,
#    "s": 180.0,
#    "cl": 175.0,
# }
# In Bohr, for hydrogen bonds
VDW_RADII = {
    'h': '2.0787', 'n': '2.9291', 'o': '2.8724', 'f': '2.7779',
    'p': '3.4015', 's': '3.4015', 'cl': '3.3070'
}
