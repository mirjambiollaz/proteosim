def read_fasta(filepath):
    """
    To import Fasta files to save the IDs and sequence in a dictonary

    Parameters
    -------
    filepath: str
        Input the path of the fasta file to be read

    Returns
    -------
    protein_map: dict
        Dictonary with IDs + Protein sequence
    """

    protein_map = {}
    current_id = None
    current_sequence = []

    with open(filepath, 'r', encoding='utf-8') as fasta_handle:
        for line in fasta_handle:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith('>'):
                if current_id is not None:
                    protein_map[current_id] = ''.join(current_sequence)
                    current_sequence = []
                current_id = stripped.split('|')[1]             # muss man rausnehmen wenn man keine 
            else:
                current_sequence.append(stripped)

    if current_id is not None:
        protein_map[current_id] = ''.join(current_sequence)
    
    return protein_map