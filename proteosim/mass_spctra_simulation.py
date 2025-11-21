amino_acid_mass_dalton = {
    'A': 71.08, 'R': 156.19, 'N': 114.10, 'D': 115.09,
    'C': 103.15, 'E': 129.12, 'Q': 128.13, 'G': 57.05,
    'H': 137.14, 'I': 113.16, 'L': 113.16, 'K': 128.17,
    'M': 131.19, 'F': 147.18, 'P': 97.12, 'S': 87.08,
    'T': 101.11, 'W': 186.21, 'Y': 163.18, 'V': 99.13,
}

def calculate_mol_mass(peptide_seq, amino_acid_mass_dict):
    """
    Calculating the molecular mass of single peptide sequences by summing the amino-acid masses.

    Parameters
    ----------
    peptide_seq : str
        Peptide sequence consisting of one-letter amino acid codes.

    amino_acid_mass_dict : dict
        Dictionary mapping amino-acid one-letter codes (str) to their monoisotopic or average masses (float).

    Returns
    -------
    mol_mass : dict
        Dictionary containing the molecular mass for the peptide sequence.
    """
    mol_mass = {}
    mass = 0
    for aa in peptide_seq:
        mass += amino_acid_mass_dict[aa]

    mol_mass[peptide_seq] = mass
    return mol_mass

def calculate_mol_mass_collection(peptides, amino_acid_mass_dict):
    """
    Calculating the molecular mass of single peptide sequences by summing the amino-acid masses.

    Parameters
    ----------
    peptides : list of str
        Peptide sequence consisting of one-letter amino acid codes.

    amino_acid_mass_dict : dict
        Dictionary mapping amino-acid one-letter codes (str) to their monoisotopic or average masses (float).

    Returns
    -------
    mol_masses : dict
        Dictionary containing the molecular mass for the peptide sequence.
    """
    mol_masses = {}
    
    for peptide_seq in peptides:
        mass = 0
        for aa in peptide_seq:
            mass += amino_acid_mass_dict[aa]        

        mol_masses[peptide_seq] = mass
    return mol_masses

def calculate_mz_collection(peptide_mass_map, charge=2, proton_mass=1.007):
    """
    Calculate the mass-to-charge ratio (m/z) for a collection of peptides.

    Parameters
    ----------
    peptide_mass_map : dict
        Dictionary mapping peptide sequences (str) to their neutral molecular masses (float).

    charge : int, optional
        Peptide charge state (default = 2).

    proton_mass : float, optional
        Mass of a proton to be added per charge (default = 1.007 Da).

    Returns
    -------
    mass_charge_map : dict
        Dictionary mapping each peptide sequence (str) to its m/z value (float).

    """
    mass_charge_map = {}
    
    for peptide_seq, mass in peptide_mass_map.items():
        mass_charge = (mass + charge * proton_mass) / charge

        mass_charge_map[peptide_seq] = mass_charge
    return mass_charge_map

import numpy as np
import matplotlib.pyplot as plt

def plot_spectrum(mz_values, random_count_range=(0, 30000), seed=42, title = "MS spectrum"):
    """
    Plot a MS spectrum with m/z values and random intensities.

    Parameters
    ----------
    mz_values : list of floats
        The m/z values (x-axis) to plot the spectrum for.

    count_range : tuple of int, optional
        Range (min, max) of random intensity counts for the spectrum.
        Default is (0, 30000).

    seed : int, optional
        Random seed for reproducibility. Default is 42.

    title : str, optional
        Name of the Plot. Default is MS spectrum.

    Returns
    -------
    None
        Displays the spectrum plot.
    """

    np.random.seed(seed)
    intensities = np.random.randint(random_count_range[0], random_count_range[1], size=len(mz_values))

    plt.bar(mz_values, intensities, width=2, color='skyblue', edgecolor='black')
    plt.xlabel("m/z")
    plt.ylabel("Intensity")
    plt.title(title)
    plt.show()

def fragment_peptide(peptide):
    """
    Generate a combined list of b-ion and y-ion fragments for a given peptide.

    Parameters
    ----------
    peptide : str
        Peptide sequence (one-letter amino acid codes).
    Returns
    -------
    fragments : list of str
        List containing all b- and y-ion fragments.
        b-ions come first, followed by y-ions.

    """
    fragmented_peptides = []
    for l in range(1, len(peptide)):
        fragmented_peptides.append(peptide[:l])

    for l in range(1, len(peptide)):
        fragmented_peptides.append(peptide[-l:])

    return fragmented_peptides