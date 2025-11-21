from proteosim.liquid_chromatography import (
    predict_lc_retention_times,
    select_retention_time_window,
)

def test_predict_lc_retention_times():
    peptides = ["PEPTIDE"]
    expected = {"PEPTIDE": 7.80}

    actual = predict_lc_retention_times(peptides)
    assert actual == expected

def test_select_retention_time_window():
    peptide_rt_map = {"PEPTIDEA" : 4.5, "PEPTIDEB" : 6, "PEPTIDEC" : -1 }
    selected = select_retention_time_window(peptide_rt_map, lower_ret_time=0, upper_ret_time=5)

    assert selected == {"PEPTIDEA" : 4.5}
