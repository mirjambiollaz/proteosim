# Proteosim Course Package

## Overview

The **Proteosim Course Package** provides **Python-based educational tools** and exercises for teaching core concepts in **computational proteomics**. It focuses on **simulating realistic mass spectrometry (MS) experiments**—from protein digestion to chromatographic separation and MS/MS fragmentation—enabling students to explore the full protein identification workflow in a controlled, reproducible environment.
This repository is under **active development**, and interfaces or behavior may change.

## What This Project Does

The package enables you to:

* **Load and process protein sequences** (typically from FASTA files).
* **Simulate enzymatic digestion**, producing theoretical peptides.
* **Model chromatographic retention** and **MS1** / **MS2** acquisition behavior.
* **Generate synthetic but realistic mass spectra**, including:

  * *m/z tables*
  * *fragment ion spectra*
  * *chromatographic peak shapes*
* **Build end-to-end MS experiments** that can be used for teaching or prototyping algorithms.

These components provide a fully simulated pipeline that mirrors real proteomics data flow—from proteins to spectra.

## Installation

You can install the package and its dependencies using:

```bash
pip install -r requirements.txt
pip install .
```

This will install all required libraries as well as the local package itself.

## Expected Inputs

Typical inputs include:

* **FASTA files** containing protein sequences
* Optional configuration parameters for:

  * enzyme specificity
  * digestion rules
  * chromatographic simulation settings
  * MS scanning parameters

## Intermediate Processing Steps

The internal workflow usually proceeds as:

1. **Protein digestion**
   Applies enzyme rules to generate peptide sequences.

2. **Peptide filtering and modification assignment**
   (e.g., missed cleavages, optional modifications).

3. **Chromatography simulation**
   Assigns retention times and optionally generates peak shapes.

4. **Mass spectrometry simulation**

   * MS1 isotope envelopes
   * MS2 fragmentation (b/y ions, charge states, intensities)

## Typical Outputs

The package can generate:

* **m/z tables** for MS1 and MS2 spectra
* **Fragment spectra** (theoretical or simulated)
* **Chromatographic profiles**
* **Intermediate peptide tables** (sequence, mass, retention time, charge, etc.)
* **Complete synthetic MS experiments** suitable for teaching or benchmarking

## Module Overview

A quick guide to the core modules:

| Module                       | Description                                         |
| ---------------------------- | --------------------------------------------------- |
| `file_handling.py`           | FASTA loading & basic I/O utilities                 |
| `protein_digestion.py`       | Enzymatic digestion & peptide generation            |
| `liquid_chromatography.py`   | Retention time prediction & chromatography modeling |
| `mass_spectra_simulation.py` | MS1/MS2 simulation & fragmentation models           |

---

## Quickest Way to Try the Package 

The fastest way to run a full experiment is through the included end-to-end notebook:

📍 **`ms_experiment_final.ipynb`**

This notebook walks you through:

* Loading proteins
* Digesting sequences
* Simulating chromatography
* Generating MS spectra
* Inspecting and analyzing outputs

It demonstrates the **complete proteomics workflow** from start to finish.