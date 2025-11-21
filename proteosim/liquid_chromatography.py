from pyteomics import achrom

def predict_lc_retention_times(peptides):
    """
    Predict LC retention times for a list of peptide sequences.
    
    Parameters
    ----------
    peptides : list of str
        A list of peptide sequences for which the retention times should be predicted.

    Returns
    -------
    predict_rt : dict
        A dictionary mapping each peptide sequence (str) to its predicted relative LC retention time (float).
    """
    predict_rt = {}

    for peptide in peptides:
        relative_RT = achrom.calculate_RT(peptide, achrom.RCs_guo_ph7_0)
        predict_rt[peptide] = round(float(relative_RT), 2)

    return predict_rt

import matplotlib.pyplot as plt

def plot_retention_time(retention_times, resolution=30):
    """
    Plot a histogram of predicted peptide retention times.

    Parameters
    ----------
    retention_times : list of float
        A sequence of predicted retention times to be visualized.
    
    resolution : int, optional
        Number of histogram bins. Default is 30.
    
    Returns
    -------
    None - Displays the histogram plot.
    """
    
    plt.hist(retention_times, resolution)

    plt.title("Retention Time Distribution ")
    plt.xlabel("Retention Time")
    plt.ylabel("Intensity")

    return plt.show()

def select_retention_time_window(peptide_rt_map, lower_ret_time, upper_ret_time):
    """
    Predicted LC retention times for a list of peptide sequences are filter in a time-window
    
    Parameters
    ----------
    predict_rt : dict
    A dictionary mapping each peptide sequence (str) to its predicted relative LC retention time (float).

    lower_ret_time : int
    lower bound of time window (in min)

    upper_ret_time : int
    upper bound of time window (in min)
    
    Returns
    -------
    rt_time_window : dict
       
    """
    rt_time_window = {}
    
    for seq, rt in peptide_rt_map.items():
        if lower_ret_time <= rt <= upper_ret_time:
            rt_time_window[seq] = rt

    return rt_time_window