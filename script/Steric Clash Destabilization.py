from pymol import cmd, stored
import math
import time
import re
# ==========================================================
# Enter "run <pathway to this file>" into the PyMOL command line to run this script
# All distance values are in Angstroms (Å)

# ==========================================================
#              User Configuration Variables
# ==========================================================
# Pathway to the PyMOL session file or PDB file that contains the protein structure to be analyzed.
pathway = "/Users/apple/Desktop/CCIR Academy/Transthyretin Amyloidosis Research/TTR monomer (pdb).pse"

# If set to True, hydrogen atoms will be included in the steric clash calculations.
# If set to False, hydrogen atoms will be excluded from the calculations.
Include_hydrogen_atoms = False
# ==========================================================
#                 Calculation Parameters
# ==========================================================
# Maximum search radius for detecting possible steric interactions.
# This is a computational cutoff rather than a biological parameter.
search_distance = 5.0 

# Distance tolerance subtracted from the VDW-overlap criterion
buffer = 0.4
# ==========================================================

print("NOTE: Starting this calculate_strain_value script reinitializes everything.\n")
print("Type \"calculate_strain <mutation (e.g. Val30Met)>, <chain (e.g. chain C, default: chain A)>\"\n"
      "     to calculate strain value for that amino acid.")

def calculate_strain_value(mutation: str, chain: str = "chain A") -> list:
    # Tracks execution start time
    start_time = time.perf_counter()

    def reset():
        cmd.reinitialize()
        cmd.load(pathway)

        if Include_hydrogen_atoms:
            cmd.h_add()
    reset()

    amino_acids = ('ALA','ARG','ASP','ASN','CYS','GLU','GLN','GLY','HIS','ILE',
                'LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL')
    strain_value_data = []

    h_filter = "" if Include_hydrogen_atoms else " and not elem H"

    # Checks if the mutation string matches the expected format
    match = re.fullmatch(r"([A-Za-z]{3})(\d+)([A-Za-z]{3})", mutation)
    if not match:
        print(f"[!] Error: Invalid format for mutation '{mutation}'. Please use the format <Wild-Type Residue><Residue Position><Mutant Residue> "
              "(e.g., Val30Met) with residue names being the three letter code, case insensitive")
        return []

    # Extracts information from the mutation string
    position = match.group(2)
    wild_type_residue = match.group(1).upper()
    mutant_residue = match.group(3).upper()
    residue = f"{chain} and resi {position}"

    # Validates the wild-type and mutant residues against the list of standard amino acids
    if wild_type_residue not in amino_acids or mutant_residue not in amino_acids:
        print(f"[!] Error: Invalid amino acid in mutation '{mutation}'. Please use three-letter codes for standard amino acids.")
        return []

    # Get the residue name at the specified position
    stored.residue_name = []
    cmd.iterate(f"{residue} and name CA", "stored.residue_name.append(resn)")
    residue_name = stored.residue_name[0] if stored.residue_name else None

    # Checks for other potential errors
    if len(stored.residue_name) > 1:
        print(f"[!] Error: Something went wrong. Program detected {len(stored.residue_name)} CA atoms at position {position} in {chain}.")
        return []
    if not residue_name:
        print(f"[!] Error: Resi {position} not found in {chain}.")
        return []
    if residue_name != wild_type_residue:
        print(f"[!] Error: Wild-type residue mismatch at position {position} in {chain}. Expected {wild_type_residue}, but found {residue_name}.")
        return []

    # Opens the PyMOL Mutagenesis Wizard and gets the total number of rotamers
    cmd.wizard("mutagenesis")
    cmd.get_wizard().do_select(residue)
    cmd.get_wizard().set_mode(mutant_residue)
    num_rotamers = cmd.count_states("mutation")
    cmd.set_wizard()

    if num_rotamers < 1:
        print(f"[!] Error: No rotamers found for mutation '{mutation}'.")
        return []

    # Loops each rotamer in all the possible rotamer configurations
    for rotamer_number in range(1, num_rotamers + 1):

        # Opens PyMOL wizard mutagenesis and applies the mutation
        cmd.wizard("mutagenesis")
        cmd.get_wizard().do_select(residue)
        cmd.get_wizard().set_mode(mutant_residue)
        cmd.frame(rotamer_number)
        Statistical_frequency = cmd.get_title("mutation", rotamer_number)
        cmd.get_wizard().apply()
        cmd.set_wizard()

        Total_protein_strain = 0
        Total_nucleic_strain = 0
        Total_organic_strain = 0   # (Ligands/Cofactors)
        Total_inorganic_strain = 0 # (Ions/Metals)
        Total_solvent_strain = 0   # (Water)
        Total_other_strain = 0     # (Other compounds)

        # Get atoms in the mutated residue's side chain (R group)
        model_resi = cmd.get_model(f"{residue} and sidechain and not name CA {h_filter}")

        # Loops each atom in the R group
        for atom_resi in model_resi.atom:
            x_i, y_i, z_i = atom_resi.coord

            # Get surrounding atoms (excluding the residue's own atoms)
            cur_atom = f"model {atom_resi.model} and index {atom_resi.index}"
            model_surr = cmd.get_model(
                f"({cur_atom} around {search_distance}) "
                f"and not (byres ({cur_atom})) {h_filter}"
            )
            # Loops through all the surrounding atoms for each R group atom
            for atom_sur in model_surr.atom:
                x_j, y_j, z_j = atom_sur.coord

                Distance = math.sqrt((x_i - x_j) ** 2 + (y_i - y_j) ** 2 + (z_i - z_j) ** 2)
                # Calculates the overlap distance
                Overlap = max(0, atom_resi.vdw + atom_sur.vdw - Distance - buffer)
                # Approximate potential destabilization energy through quadratic function of overlap distance,
                # which is a common approach in molecular modeling to estimate steric strain.
                Strain = Overlap ** 2

                # Uniquely scope the selection to both the specific model/object and atom index
                target_atom = f"model {atom_sur.model} and index {atom_sur.index}"

                # Check which PyMOL category the surrounding atom belongs to:
                if cmd.count_atoms(f"{target_atom} and polymer.protein"):
                    Total_protein_strain += Strain
                elif cmd.count_atoms(f"{target_atom} and polymer.nucleic"):
                    Total_nucleic_strain += Strain
                elif cmd.count_atoms(f"{target_atom} and organic"):
                    Total_organic_strain += Strain
                elif cmd.count_atoms(f"{target_atom} and inorganic"):
                    Total_inorganic_strain += Strain
                elif cmd.count_atoms(f"{target_atom} and solvent"):
                    Total_solvent_strain += Strain
                else:
                    Total_other_strain += Strain

        # Append all strain values to the data array:
        Total_strain = (
            Total_protein_strain
            + Total_nucleic_strain
            + Total_organic_strain
            + Total_inorganic_strain
            + Total_solvent_strain
            + Total_other_strain
        )
        strain_value_data.append([
            rotamer_number,
            Statistical_frequency,
            Total_strain,
            Total_protein_strain, 
            Total_nucleic_strain, 
            Total_organic_strain, 
            Total_inorganic_strain, 
            Total_solvent_strain,
            Total_other_strain
        ])
        reset()

    # Forces python to stop until all processes are completed
    cmd.sync()
    cmd.refresh()
    # Tracks execution end time
    end_time = time.perf_counter()

    # Print the results in a table
    print(f"\nStrain results for {mutation} (run_time: {end_time - start_time:.2f}s)")
    print("=" * 111)
    print(
        f"{'Rotamer':<8} | {'Frequency':<10} | {'Total':<10} | "
        f"{'Protein':<10} | {'Nucleic':<10} | {'Organic':<10} | "
        f"{'Inorganic':<10} | {'Solvent':<10} | {'Other':<10}"
    )
    print("-" * 111)
    for row in strain_value_data:
        print(
            f"{row[0]:<8} | {row[1]:<10} | {row[2]:<10.4f} | "
            f"{row[3]:<10.4f} | {row[4]:<10.4f} | {row[5]:<10.4f} | "
            f"{row[6]:<10.4f} | {row[7]:<10.4f} | {row[8]:<10.4f}"
        )
    print("=" * 111)

    return strain_value_data

# Connects to PyMOL command line
cmd.extend("calculate_strain", calculate_strain_value)
