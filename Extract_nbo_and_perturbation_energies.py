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

def extract_occupancy(filename):
    occupancy_data = []
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
                occupancy_data.append(line.strip())

    output_file = filename.replace(".log", "_Occupancy.txt")
    with open(output_file, "w", encoding="utf-8") as out_occ:
        out_occ.writelines("\n".join(occupancy_data))

    print(f"Occupancy data saved to {output_file}")
    return occupancy_data 

bond_types = {
    "BD(1)": "σ",
    "BD*(1)": "σ*",
    "BD(2)": "π",
    "BD*(2)": "π*",
    "LP": "n"
}

def classify_interactions(perturbation_data, log_file):
    classified_data = {
        "nσ": [], "nσ*": [], "nπ": [], "nπ*": [],
        "σσ*": [], "σπ": [], "σπ*": [], "ππ*": []
    }

    for line in perturbation_data:
        words = line.split()
        if len(words) < 5:
            continue  
        donors = []
        for key in bond_types.keys():
            if key.replace(" ", "") in line.replace(" ", ""):  
                donors.append(bond_types[key])
        for donor in donors:
            for acceptor in donors:
                if donor != acceptor: 
                    interaction = f"{donor}{acceptor}"
                    if interaction in classified_data:
                        classified_data[interaction].append(line)

    classification_file = log_file.replace(".log", "_Classification.txt")
    with open(classification_file, "w", encoding="utf-8") as out_class:
        for category, interactions in classified_data.items():
            if interactions:  
                out_class.write(f" {category} Interactions \n")
                out_class.writelines("\n".join(interactions))
                out_class.write("\n\n")  

    print(f"Classified interactions saved to {classification_file}")
    return classified_data

def filter_interactions_by_indices(classified_data):
    atom_index = input("Enter the atom no. to find all interactions: ").strip()
    print(f"\nInteraction involving atom {atom_index}:")
    found = False

    for category, interactions in classified_data.items():
        for interaction in interactions:
            if atom_index in interaction:
                parts = interaction.split()
                energy = parts[-3]
                output_interaction = " ".join(parts[:-3])
                print(f"{category}: {output_interaction} {energy} kcal")
                found = True
    if not found:
        print("No interactions found for the given atom index")

log_file = input("Enter the Gaussian log file: ").strip()
perturbation_data = extract_perturbation(log_file)
occupancy_data = extract_occupancy(log_file)
classified_data = classify_interactions(perturbation_data, log_file)
filter_interactions_by_indices(classified_data)











#Based on orbital indices
#def filter_interactions_by_indices(classified_data):
    #index1 = input("Enter the first orbital no.: ").strip()
    #index2 = input("Enter the second orbital no.: ").strip()
    #print(f"\nInteraction between orbitals {index1} and {index2}:")
    #found = False

    #for category, interactions in classified_data.items():
        #for interaction in interactions:
            #if index1 in interaction and index2 in interaction:
                #parts = interaction.split()
                #energy = parts[-3]
                #output_interaction = " ".join(parts[:-3])
                #print(f"{category}: {output_interaction} {energy} kcal")
                #found = True
                #break
        #if found:
            #break

    #if not found:
        #print("No interactions found between these orbitals.")
