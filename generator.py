import sys
import random
import argparse
from generator_io import GeneratorIO

def generate_catalyzers(data):
    catalyzer_params = data["catalyzer_params"]
    species = data["species"]

    min_length, max_length = catalyzer_params[0]
    num_cond_catalyzers = catalyzer_params[1]
    num_cll_catalyzers = catalyzer_params[2]
    both_on = catalyzer_params[3] == 'ON'

    eligible_species = [specie[0] for specie in species[1:] if min_length <= len(specie[0]) <= max_length]

    if len(eligible_species) < (num_cond_catalyzers + num_cll_catalyzers):
        raise ValueError("Error! Not enough eligible species to satisfy the catalyzer requirements.")

    random.shuffle(eligible_species)
    cond_catalyzers = eligible_species[:num_cond_catalyzers]
    cll_catalyzers = eligible_species[num_cond_catalyzers:num_cond_catalyzers + num_cll_catalyzers]
    remaining_species = eligible_species[num_cond_catalyzers + num_cll_catalyzers:]

    for specie in remaining_species:
        if random.choice([True, False]):
            cond_catalyzers.append(specie)
        else:
            cll_catalyzers.append(specie)

    catalyzers = {'eligible_cond_species': cond_catalyzers, 'eligible_cll_species': cll_catalyzers, 
                  'both_on': both_on, 'num_cond_catalyzers': num_cond_catalyzers, 'num_cll_catalyzers': num_cll_catalyzers}

    return catalyzers

def get_reaction_catalyzer(catalyzers, reaction_type):
    eligible_cond_species = catalyzers['eligible_cond_species']
    eligible_cll_species = catalyzers['eligible_cll_species']
    both_on = catalyzers['both_on']
    num_cond_catalyzers = catalyzers['num_cond_catalyzers']
    num_cll_catalyzers = catalyzers['num_cll_catalyzers']
    
    if reaction_type == "condensation":
        if both_on:
            catalyzers = random.sample(eligible_cond_species + eligible_cll_species, num_cond_catalyzers)
        else:
            catalyzers = random.sample(eligible_cond_species, num_cond_catalyzers)
    elif reaction_type == "cleavage":
        if both_on:
            catalyzers = random.sample(eligible_cond_species + eligible_cll_species, num_cll_catalyzers)
        else:
            catalyzers = random.sample(eligible_cll_species, num_cll_catalyzers)
    else:
        raise ValueError("Invalid reaction type. It should be 'condensation' or 'cleavage'.")
    
    return catalyzers

def generate_condensation_reactions(data):
    species = data["species"][1:]
    reactions = data["reactions"]["conds"]
    condensation_reactions = []

    for i in range(len(species)):
        for j in range(len(species)):
            reagent_1 = species[i][0]
            reagent_2 = species[j][0]

            for reaction in reactions:
                if reagent_1.endswith(reaction[0]) and reagent_2.startswith(reaction[1]):
                    new_reaction = [reagent_1 + reagent_2, reagent_1, reagent_2, reaction[2]]
                    catalyzer = get_reaction_catalyzer(data["catalyzers"], "condensation") 
                    new_reaction.append(catalyzer)
                    condensation_reactions.append(new_reaction)

    return condensation_reactions

def generate_cleavage_reactions(catalyzers, species, reactions):
    cleavage_reactions = []

    for specie in species[1:]:
        specie_name = specie[0]
        for i in range(1, len(specie_name)):
            cleavage_1 = specie_name[:i]
            cleavage_2 = specie_name[i:]
            for reaction in reactions:
                reactant = reaction[0]
                reactant_length = len(reactant)
                reactant_half_length = reactant_length // 2
                if cleavage_1.endswith(reactant[:reactant_half_length]) and cleavage_2.startswith(reactant[reactant_half_length:]): 
                    new_reaction = [specie_name, cleavage_1, cleavage_2, reaction[1]]
                    catalyzer = get_reaction_catalyzer(catalyzers, "cleavage")
                    new_reaction.append(catalyzer)
                    cleavage_reactions.append(new_reaction)

    return cleavage_reactions

def generate_new_species(data):
    new_species = data["cond_reactions"]
    cleavage_products = []

    while True:
        new_cleavage_products = generate_cleavage_reactions(data["catalyzers"], 
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

    return new_species

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate species and reactions.")
    parser.add_argument("file_path", help="The path to the input file.")
    parser.add_argument("-o", "--output", help="The name of the output file.")
    parser.add_argument("-debug", action="store_true", help="Enable debug mode.", default=False)
    args = parser.parse_args()

    debug = args.debug
    file_path = args.file_path
    output_file = args.output

    generatorIO = GeneratorIO(file_path, debug, output_file)
    try:
        parsed_data = generatorIO.parse_data()
        parsed_data["catalyzers"] = generate_catalyzers(parsed_data)
        parsed_data["cond_reactions"] = generate_condensation_reactions(parsed_data)
        parsed_data["cll_reactions"] = generate_cleavage_reactions(parsed_data["catalyzers"],
                                                                   parsed_data["species"],
                                                                   parsed_data["reactions"]["clls"])
        generate_new_species(parsed_data)
        generatorIO.write_data(parsed_data)
    except Exception as e:
        print("An error occurred:", str(e))
        sys.exit(1)