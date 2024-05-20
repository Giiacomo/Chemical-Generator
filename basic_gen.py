import sys
import random
from generator_io import GeneratorIO

#generate_catalyzers() returns a dictionary that contains the condensation and cleavage catalyzers randomly chosen among the species 
#PRIMA RIGA = <lunghezza minima di una specie chimica per diventare catalizzatore>
#SECONDA RIGA = <numero di specie chimiche catalizzatrici di una classe di condensazioni>
#TERZA RIGA = <numero di specie chimiche catalizzatrici di una classe di cleavage>
#QUARTA RIGA = <catalizzatore di condensazione e cleavage: ON/OFF>
def generate_catalyzers(data):
    catalyzer_params = data["catalyzer_params"]
    species = data["species"]

    min_length, max_length = catalyzer_params[0]
    num_cond_catalyzers = catalyzer_params[1]
    num_cll_catalyzers = catalyzer_params[2]
    both_on = catalyzer_params[3] == 'ON'

    # Filtra le specie che soddisfano la lunghezza minima
    eligible_species = [specie[0] for specie in species if len(specie[0]) >= min_length and len(specie[0]) <= max_length]


    if len(eligible_species) < (num_cond_catalyzers + num_cll_catalyzers):
        raise ValueError("Error! Not enough eligible species to satisfy the catalyzer requirements or catalysis is OFF.")

    random.shuffle(eligible_species)

    if both_on:
        catalyzer_candidates = eligible_species[:num_cond_catalyzers + num_cll_catalyzers] #3 + 2 = 5
        cond_catalyzers = catalyzer_candidates[:num_cond_catalyzers] # i primi 2
        cll_catalyzers = catalyzer_candidates[num_cond_catalyzers:]
    else:
        cond_catalyzers = eligible_species[:num_cond_catalyzers]
        cll_catalyzers = eligible_species[num_cond_catalyzers:num_cond_catalyzers + num_cll_catalyzers]

    catalyzer = {
        "cond": cond_catalyzers,
        "cll": cll_catalyzers
    }

    return catalyzer


def generate_condensation_reactions(data):
    if len(data["catalyzers"]["cond"]) == 0:
        return []
    

    species = data["species"]
    reactions = data["reactions"]["conds"]

    condensation_reactions = []

    for i in range(len(species)):
        for j in range(len(species)):
            reagent_1 = species[i][0]
            reagent_2 = species[j][0]

            for reaction in reactions:
                if reagent_1.endswith(reaction[0]) and reagent_2.startswith(reaction[1]):
                    new_reaction = [reagent_1+reagent_2 ,reagent_1, reagent_2, reaction[2]] #es: AB, BA, 0.1 
                    condensation_reactions.append(new_reaction)

    return condensation_reactions


def generate_cleavage_reactions(cll_catalyzers, species, reactions):
    if len(cll_catalyzers) == 0:
        return []

    cleavage_reactions = []

    for specie in species:

        specie_name = specie[0]
        for i in range(1, len(specie_name)):
            cleavage_1 = specie_name[:i]
            cleavage_2 = specie_name[i:]
            for reaction in reactions:
                reactant = reaction[0]
                reactant_length = len(reactant)
                reactant_half_length = reactant_length // 2 #Es: R-AB-R -> we divide this reactant at the half of its "known" part. In this case AB -> R-A + B-R
                #This condition checks if the first part of R-AB-R (R-A) ends in A AND the last part (B-R) starts with B 
                if cleavage_1.endswith(reactant[:reactant_half_length]) and cleavage_2.startswith(reactant[reactant_half_length:]): 
                    new_reaction = [specie_name, cleavage_1, cleavage_2, reaction[1]]
                    cleavage_reactions.append(new_reaction)

    return cleavage_reactions


def generate_new_species(data):
    new_species = data["cond_reactions"]
    cleavage_products = []

    
    # Itero fino a quando ci sono nuovi prodotti generati dai cleavage, che non erano presenti nelle specie
    while True:
        new_cleavage_products = generate_cleavage_reactions(data["catalyzers"]["cll"], 
                                                            new_species, 
                                                            data["reactions"]["clls"])

        cleavage_products.extend(new_cleavage_products)

        new_species_copy = data["species"]

        cleavage_products = [product for product in cleavage_products if product[0] not in [species[0] for species in data["species"]]]
        
        new_species_to_add = list(set((product[0], '1.00E-15', '0.') for product in cleavage_products))
        new_species_to_add = [list(item) for item in new_species_to_add]
        
        data["species"].extend(new_species_to_add)
        
        if new_species_copy == data["species"]:
            break

    print(data["species"])

    return new_species




#I dati vengono inseriti in un dizionario con questa forma:
# {
# 'species':  [[<nomespecie>, <concentrazione>, <contributo>]],
# 'catalyzer_params': [<prob_catalyzer>, <prob_cond>],  
# 'reactions': {'conds': [<specie>], 'clls': [<specie>] },
# 'catalyzers': {'cond': [<specie>], 'cll': [<specie>],
#                'cond_reactions': [<reactant_1>, <reactant_2>, <v>]
#                'cll_reactions': [<specie>, <cleavage_1>, <cleavage_2>, <v>]}
# }

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Error!")
        print("Correct usage: python3 basic_gen.py <file-name.txt>")
        sys.exit(1)

    file_path = sys.argv[1]
    generatorIO = GeneratorIO(file_path)
    try:
        parsed_data = generatorIO.parse_data()
        parsed_data["catalyzers"] = generate_catalyzers(parsed_data)
        print(parsed_data["catalyzers"])
        # parsed_data["cond_reactions"] = generate_condensation_reactions(parsed_data)
        # parsed_data["cll_reactions"] = generate_cleavage_reactions(parsed_data["catalyzers"]["cll"],
        #                                                            parsed_data["species"],
        #                                                            parsed_data["reactions"]["clls"]
        #                                                           )
        # generate_new_species(parsed_data)
        # generatorIO.write_data(parsed_data)

    except Exception as e:
        print("An error occurred:", str(e))
        sys.exit(1)