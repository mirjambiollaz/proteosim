from proteosim.mass_spctra_simulation import (
    calculate_mol_mass,
    calculate_mol_mass_collection,
    calculate_mz_collection,
    fragment_peptide,
)

def test_calculate_mol_mass():
    aa_mass_dict = {'A': 1, 'B': 2, 'C': 3, 'D': 4, 'E': 5, 'F': 6, 'G': 7,}
    test_peptide = "ABCDEFG"
    expected_mass = {"ABCDEFG": 28}
    calculated_mass = calculate_mol_mass(test_peptide, aa_mass_dict)
    assert calculated_mass == expected_mass

def test_calculate_mol_mass_collection():
    aa_mass_dict = {'A': 1, 'B': 2, 'C': 3, 'D': 4, 'E': 5, 'F': 6, 'G': 7,}
    peptides = {"AAA", "ABC",}
    expected = {"AAA": 3, "ABC": 6}
    actual = calculate_mol_mass_collection(peptides, aa_mass_dict)

    assert actual == expected

def test_calculate_mz_collection():
    peptide_mass_map = {"B": 10, "C": 50}
    actual = calculate_mz_collection(peptide_mass_map, charge=2)
    expected = {"B": 6.007, "C": 26.007}

    assert actual == expected

def test_fragment_peptide():
    peptide = "PEPT"
    expected = fragment_peptide("PEPT")
    actual = ["P", "PE", "PEP", "T", "PT", "EPT"]

    assert actual == expected  