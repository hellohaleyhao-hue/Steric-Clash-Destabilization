# SSCREE
Steric Strain Calculation and Repulsive Energy Estimation (SSCREE)

## Introduction
Estimating the impact of specific missense gene mutations on protein structure, thus function, comes from various sources of evidence such as changes to bonding, positional shifts, and evolutionary conservation. Steric clashes occur when the atoms of the mutated amino acid are in close proximity to surrounding atoms.

Generally, the more clashes there are, the less stable the variant will make of the protein. Since natural selection optimizes native protein chains to pack tightly without severe destabilizing overlaps in the electron cloud, this tool, which calculates the total steric overlap from all atoms in the sidechain and surrounding atoms as a strain value, can be used to predict the pathogenicity of that variant, especially variants of uncertain significance (VUS).

## Prerequisites
### PyMOL
The script is used and tested on **PyMOL v3.1.8** (Schrödinger, LLC)

* **Tested version**: PyMOL 3.1.8
* **Compatibility**: Compatible with PyMOL 2.x and 3.x
* **Dependencies**: Standard python libraries, no additional dependencies needed

### Protein model
A protein model is needed in the script, represented by the `pathway` variable, for analysis.

* Requires a `.pdb` or `.pse` file containing the target protein model for clash analysis.
* Should include all relevant objects (e.g. protein chains, ligands, cofactors) intended for strain value evaluation.
* Remove extraneous or non-target objects or else they'll be included in calculation.

## Implementation & Usage
The script takes a variant as input and outputs a table displaying statistical frequency and strain values (amount of overlap) for various groups (e.g. nucleic acids, inorganic, protein) for each rotamer. It calculates strain score through the summation of the squares of Van der Waals Overlaps as a repulsive steric potential estimation: $\text{Total Strain} = \sum (\text{Overlap})^2$. All distance units are in Angstroms (Å).

Run the following command in the PyMOL command line:
```pymol
run /path/to/script.py
```
Further instructions are provided on the console.
