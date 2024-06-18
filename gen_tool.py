import sys
from generator_io import GenToolIO

def generate_condensation_reactions(data):

    
    species = data["species"][1:]
    n_s, n_d, v = map(float, data["conds"])
    condensation_reactions = {'reactions': [], 'v': v}
    for i in range(len(species)):
        for j in range(len(species)):
            reagent_1 = species[i][0]  # Reimpostiamo il prefisso R-
            reagent_2 = species[j][0]
            if len(reagent_1) == n_s and len(reagent_2) == n_d:
                condensation_reactions['reactions'].append([reagent_1, reagent_2])  # "V" rappresenta la velocità associata alla reazione

    return condensation_reactions


def generate_cleavage_reactions(data):

    species = data["species"][1:]
    Nts = list(map(int, [nt for nt in data["clls"][0]]))
    v = float(data["conds"][1])
    cleavage_reactions = {'reactions': [], 'v': v}
    for nt in Nts:
        for specie in species:
            specie_name = specie[0]
            subspecies = [specie_name[i:i+nt] for i in range(0, len(specie_name)) if len(specie_name[i:i+nt]) == nt]
            cleavage_reactions['reactions'] += subspecies  # Concatena i risultati alla lista cleavage_reactions

    #flat
    cleavage_reactions_set = set(cleavage_reactions['reactions'])
    cleavage_reactions['reactions'] = list(cleavage_reactions_set)

    return cleavage_reactions



if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Error!")
        print("Correct usage: python3 basic_gen.py <file-name.txt>")
        sys.exit(1)

    file_path = sys.argv[1]
    genTool = GenToolIO(file_path)
    try:
        parsed_data = genTool.parse_data()
        parsed_data["gen-conds"] = generate_condensation_reactions(parsed_data)
        parsed_data["gen-clls"] = generate_cleavage_reactions(parsed_data)

        genTool.write_data(parsed_data)
        
    except Exception as e:
        print("An error occurred:", str(e))
        sys.exit(1)