import re
from collections import deque

#Extract perturbation data
def extract_perturbation(log_file):
    with open(log_file, "r", encoding="utf-8") as file:
        lines = file.readlines()
    perturbation_data = []
    record = False
    for line in lines:
        if "SECOND ORDER PERTURBATION THEORY ANALYSIS OF FOCK MATRIX IN NBO BASIS" in line:
            record = True
        if record:
            if "NATURAL BOND ORBITALS (Summary):" in line:
                break
            if "RY" in line:
                continue
            perturbation_data.append(line.strip())
    output_file = log_file.replace(".log", "_Perturbation.txt")
    with open(output_file, "w", encoding="utf-8") as out_perturb:
        out_perturb.writelines("\n".join(perturbation_data))
    print(f"Perturbation data saved to {output_file}")
    return perturbation_data

#Extract contribution data
def extract_contribution(filename):
    contribution_data = []
    record = False
    with open(filename, "r", encoding="utf-8") as file:
        lines = file.readlines()
    for line in lines:
        if "(Occupancy)" in line:
            record = True 
            continue
        if record:
            if "81." in line: 
                break
            if any(c.isalpha() for c in line):
                contribution_data.append(line.strip())
    output_file = filename.replace(".log", "_Contribution.txt")
    with open(output_file, "w", encoding="utf-8") as out_con:
        out_con.writelines("\n".join(contribution_data))
    print(f"Contribution data saved to {output_file}")
    return contribution_data

bond_types = {
    "BD(1)": "σ",  
    "BD*(1)": "σ*", 
    "BD(2)": "π",  
    "BD*(2)": "π*", 
    "LP": "n"      
}

#Classify interactions
def classify_interactions(perturbation_data, log_file):
    classified_data = {
        "nσ*": [], "nπ*": [], "σσ*": [], "σπ*": [], "ππ*": []
    }
    energy_totals = {}

    for line in perturbation_data:
        if "RY" in line or not line.strip():
            continue

        if "LP" in line:
            donor = "n"
        elif "BD ( 1)" in line:
            donor = "σ"
        elif "BD ( 2)" in line:
            donor = "π"
        else:
            continue

        if "BD*( 1)" in line:
            acceptor = "σ*"
        elif "BD*( 2)" in line:
            acceptor = "π*"
        else:
            continue

        interaction = donor + acceptor
        if interaction not in classified_data:
            continue

        try:
            energy = float(line.strip().split()[-3])
        except Exception:
            energy = 0.0

        classified_data[interaction].append((line.strip(), energy))
        energy_totals[interaction] = energy_totals.get(interaction, 0.0) + energy

    classification_file = log_file.replace(".log", "_Classification.txt")
    with open(classification_file, "w", encoding="utf-8") as out_class:
        for category, items in classified_data.items():
            if items:
                out_class.write(f"{category} Interactions\n")
                for entry, energy in items:
                    out_class.write(f"{entry}\n")
                out_class.write(f"Total {category} energy: {energy_totals[category]:.2f} kcal\n\n")

    print(f"Classified interactions saved to {classification_file}")
    return classified_data


#Read PDB connectivity
def read_pdb_connectivity(pdb_file):
    connectivity = {}
    with open(pdb_file, 'r') as file:
        for line in file:
            if line.startswith("CONECT"):
                parts = line.split()
                first_atom = parts[1]
                if first_atom not in connectivity:
                    connectivity[first_atom] = []
                for bonded_atom in parts[2:]:
                    connectivity.setdefault(bonded_atom, []).append(first_atom)
                    connectivity[first_atom].append(bonded_atom)
    return connectivity

#Find shortest bond path
def find_bond_path_with_distance(connectivity, start_atom, target_atom):
    queue = deque([(start_atom, [start_atom])])
    visited = set()
    while queue:
        current_atom, path = queue.popleft()
        if current_atom == target_atom:
            return len(path) - 1, path
        visited.add(current_atom)
        for neighbor in connectivity.get(current_atom, []):
            if neighbor not in visited:
                queue.append((neighbor, path + [neighbor]))


def extract_atoms_from_line(line):
    atoms = re.findall(r'([A-Z]{1,2})\s+(\d+)', line)
    if len(atoms) >= 2:
        return atoms[0][1], atoms[-1][1]
    digits = re.findall(r'\d+', line)
    if len(digits) >= 4:
        return digits[2], digits[4]
    return None, None

#Classify interactions with bond range
def classify_interactions_with_bond_range(classified_data, connectivity, log_file):
    new_file = log_file.replace(".log", "_BondRangeClassification.txt")
    with open(new_file, "w", encoding="utf-8") as out:
        for category, interactions in classified_data.items():
            total_energy_short = 0.0
            total_energy_long = 0.0
            short_range_entries = []
            long_range_entries = []
            
            for line, energy in interactions:
                atom1, atom2 = extract_atoms_from_line(line)
                if atom1 and atom2:
                    bond_count, path = find_bond_path_with_distance(connectivity, atom1, atom2)
                    if bond_count == -1:
                        continue
                                        
                    range_type = "short range" if bond_count < 3 else "long range"
                    clean_line = "  ".join(line.split()[:-2]) 
                    
                    if range_type == "short range":
                        short_range_entries.append(f"{clean_line}  {bond_count} bonds  {range_type}")
                        total_energy_short += energy
                    elif range_type == "long range":
                        long_range_entries.append(f"{clean_line}  {bond_count} bonds  {range_type}")
                        total_energy_long += energy

            if short_range_entries:
                out.write(f"{category} Short Range Interactions\n")
                out.writelines("\n".join(short_range_entries) + "\n")
                out.write(f"Total {category} short range energy: {total_energy_short:.2f} kcal\n\n")
                
            if long_range_entries:
                out.write(f"{category} Long Range Interactions\n")
                out.writelines("\n".join(long_range_entries) + "\n")
                out.write(f"Total {category} long range energy: {total_energy_long:.2f} kcal\n\n")
    print(f"Bond range classification saved to {new_file}")

#Classify long-range interactions only
def classify_long_range_interactions_only(classified_data, connectivity, log_file):
    new_file = log_file.replace(".log", "_LongRangeOnly.txt")
    with open(new_file, "w", encoding="utf-8") as out:
        for category, interactions in classified_data.items():
            total_energy_long = 0.0
            long_range_entries = []
            
            for line, energy in interactions:
                atom1, atom2 = extract_atoms_from_line(line)
                if atom1 and atom2:
                    bond_count, path = find_bond_path_with_distance(connectivity, atom1, atom2)
                    if bond_count == -1:
                        continue
                    
                    range_type = "long range" if bond_count >= 3 else None
                    if range_type == "long range":
                        clean_line = "  ".join(line.split()[:-2]) 
                        long_range_entries.append(f"{clean_line}  {bond_count} bonds  {range_type}")
                        total_energy_long += energy
            
            if long_range_entries:
                out.write(f"{category} Long Range Interactions\n")
                out.writelines("\n".join(long_range_entries) + "\n")
                out.write(f"Total {category} long range energy: {total_energy_long:.2f} kcal\n\n")

    print(f"Long-range only classification saved to {new_file}")

#input
log_file = input("Enter the path to the log file: ").strip()
pdb_file = input("Enter the path to the PDB file: ").strip()

connectivity = read_pdb_connectivity(pdb_file)

perturbation_data = extract_perturbation(log_file)
contribution_data = extract_contribution(log_file)
classified_data = classify_interactions(perturbation_data, log_file)
classify_interactions_with_bond_range(classified_data, connectivity, log_file)
classify_long_range_interactions_only(classified_data, connectivity, log_file)
