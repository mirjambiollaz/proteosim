from proteosim.protein_digestion import (
    enzyme_cleavage_patterns,
    digest_protein_collection,
    )

def test_digest_protein_collection():
    dummy_proteins = {
        "PROTEIN_ID_1": "ABCDEFGHIJKLMNOPQRSTUVWXYZ",
        "PROTEIN_ID_2": "ETZCGKPHUUVGTRXUHVGKJKUZRPINIVUC"
    }
    cleave_pattern = r'(?<=[KR])(?!P)'                      # keine globalen Variablen, mehrere Schnittstellen
    test_digested_peptides_map = digest_protein_collection(
    dummy_proteins,
    cleave_pattern=cleave_pattern,
    min_pep_len=5,
    max_pep_len=30,
)
    assert test_digested_peptides_map["PROTEIN_ID_1"] == ["ABCDEFGHIJK", "LMNOPQR", "STUVWXYZ"]
    assert test_digested_peptides_map["PROTEIN_ID_2"] == ["ETZCGKPHUUVGTR","XUHVGK","UZRPINIVUC"]

