enzyme_cleavage_patterns = {
    'LysC': r'(?<=K)',
    'LysN': r'(?=K)',
    'ArgC': r'(?<=R)',
    'Trypsin': r'(?<=[KR])(?!P)',
}

def digest_protein_sequence(protein_seq, cleave_pattern):
    """
    Function to digest protein sequnece with a specific cleavage pattern

    Parameters
    ----------
    protein_seq : str
        containing sequence of the protein to digest
    
    cleave_pattern : str
        containing specific cleavage pattern
    
    Returns
    -------
    peptides : list
        containing all the fragments after cleaving
    
    """
    peptides = re.split(cleave_pattern, protein_seq)
    return peptides

def digest_protein_collection(protein_map, cleave_pattern, min_pep_len=5, max_pep_len=30):
    """
    Add a short description here.
    protein_seq : str
        containing sequence of the protein to digest
    
    cleave_pattern : dict
        containing specific cleavage pattern
    
    min_pep_len : int, Default = 5, optional
        decides how long the peptide has to be
    
    max_pep_len : int, Default = 30, optional
        decides how long the peptide is allowed to be

    Returns
    -------
    digest_collection : dict
        all Proteins with the belonging peptides
    """
    digest_collection = {}
    for protein_id, protein_seq in protein_map.items():
        if protein_seq is None:
            raise ValueError("Protein sequence is missing (None).")
        if max_pep_len < min_pep_len:
            raise ValueError("max_pep_len must be >= min_pep_len.")
        
        all_peptides = re.split(cleave_pattern, protein_seq)
        peptides_filterd = [p for p in all_peptides if min_pep_len <= len(p) <= max_pep_len]
        digest_collection[protein_id] = peptides_filterd
    return digest_collection

def compute_sequence_coverage(protein_seq, peptides):
    """
    Compute the sequence coverage of a protein given a list of peptides.

    Sequence coverage is the fraction of residues in the protein
    that are represented by at least one peptide. The result
    is returned as a percentage (0–100).

    Parameters
    ----------
    protein_seq : str
        Full protein sequence.
    peptides : list of str
        List of peptides obtained from digestion.

    Returns
    -------
    coverage_percent : float
        Protein sequence coverage as a percentage.
    """
    if not protein_seq:
        return 0.0  # empty protein

    protein_len = len(protein_seq)
    coverage_mask = [False] * protein_len

    for pep in peptides:
        if not pep:
            continue
        start = 0
        while True:
            idx = protein_seq.find(pep, start)
            if idx == -1:
                break
            # mark covered residues
            for i in range(idx, idx + len(pep)):
                coverage_mask[i] = True
            start = idx + 1  # handle overlapping occurrences

    covered_residues = sum(coverage_mask)
    coverage_percent = (covered_residues / protein_len) * 100
    return coverage_percent