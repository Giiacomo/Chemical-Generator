import sys
import random
import argparse
from generator_io import GeneratorIO

def assign_catalyzers(num_catalyzers, eligible_species, reactions):
    catalyzer_list = []
    species_pool = eligible_species[:]

    if num_catalyzers > len(reactions):
        for reaction in reactions:
            if not species_pool:
                raise ValueError("Error! Not enough unique species to cover all required catalyzers.")
            chosen = random.choice(species_pool)
            catalyzer_list.append({'catalyzer_specie': chosen, 'reaction': reaction})
            species_pool.remove(chosen)
    else:
        assigned_reactions = random.sample(reactions, num_catalyzers)
        for reaction in assigned_reactions:
            if not species_pool:
                species_pool = eligible_species[:]
            chosen = random.choice(species_pool)
            catalyzer_list.append({'catalyzer_specie': chosen, 'reaction': reaction})
            species_pool.remove(chosen)

    while len(catalyzer_list) < num_catalyzers:
        if not species_pool:
            species_pool = eligible_species[:]
            
        chosen = random.choice(species_pool)
        available_reactions = [reaction for reaction in reactions if any(c['reaction'] == reaction for c in catalyzer_list)]
        
        if not available_reactions:
            continue

        chosen_reaction = random.choice(available_reactions)
        catalyzer_list.append({'catalyzer_specie': chosen, 'reaction': chosen_reaction})

    return catalyzer_list

def generate_catalyzers(data):
    catalyzer_params = data["catalyzer_params"]
    species = data["species"]
    cond_reactions = [tuple(r) for r in data["reactions"]["conds"]]
    cll_reactions = [tuple(r) for r in data["reactions"]["clls"]]

    min_length, max_length = catalyzer_params[0]
    num_cond_catalyzers = catalyzer_params[1]
    num_cll_catalyzers = catalyzer_params[2]
    both_on = catalyzer_params[3] == 'ON'

    # Scelgo le specie che rispettano il criterio
    eligible_species = [specie[0] for specie in species[1:] if min_length <= len(specie[0]) <= max_length]
    # Se i potenziali catalizzatori sono meno di quelli richiesti, da errore
    if len(eligible_species) < (num_cond_catalyzers + num_cll_catalyzers):
        raise ValueError("Error! Not enough eligible species to satisfy the catalyzer requirements.")
    
    catalyzers = {
        'eligible_cond_species': [],
        'eligible_cll_species': [],
        'both_on': both_on,
        'num_cond_catalyzers': num_cond_catalyzers,
        'num_cll_catalyzers': num_cll_catalyzers
    }

    catalyzers['eligible_cond_species'] = assign_catalyzers(num_cond_catalyzers, eligible_species, cond_reactions)

    if both_on:
        cll_candidates = eligible_species
    else:
        cll_candidates = [s for s in eligible_species if s not in [c['catalyzer_specie'] for c in catalyzers['eligible_cond_species']]]

    catalyzers['eligible_cll_species'] = assign_catalyzers(num_cll_catalyzers, cll_candidates, cll_reactions)

    return catalyzers

def get_reaction_catalyzer(catalyzers, reaction):
    reaction_tuple = tuple(reaction)
    return [c for c in catalyzers['eligible_cond_species'] + catalyzers['eligible_cll_species'] if c['reaction'] == reaction_tuple]

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
                    catalyzer = get_reaction_catalyzer(data["catalyzers"], reaction) 
                    new_reaction.append(catalyzer)
                    condensation_reactions.append(new_reaction)

    return condensation_reactions
    
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
                    catalyzer = get_reaction_catalyzer(data["catalyzers"], reaction) 
                    new_reaction.append(catalyzer)
                    condensation_reactions.append(new_reaction)

    return condensation_reactions

def generate_cleavage_reactions(catalyzers, species, reactions):
    cleavage_reactions = []

    species = [specie for specie in species if specie != "Cont"]

    for specie_name in species:
        for i in range(1, len(specie_name)):
            cleavage_1 = specie_name[:i]
            cleavage_2 = specie_name[i:]
            for reaction in reactions:
                reactant = reaction[0]
                reactant_length = len(reactant)
                reactant_half_length = reactant_length // 2
                if cleavage_1.endswith(reactant[:reactant_half_length]) and cleavage_2.startswith(reactant[reactant_half_length:]): 
                    new_reaction = [specie_name, cleavage_1, cleavage_2, reaction[1]]
                    catalyzer = get_reaction_catalyzer(catalyzers, reaction)
                    new_reaction.append(catalyzer)
                    cleavage_reactions.append(new_reaction)

    return cleavage_reactions


def generate_new_species(data):
    cond_species = {reaction[0] for reaction in data["cond_reactions"]}
    cll_species = {product for reaction in data["cll_reactions"] for product in reaction[1:3]}
    
    new_species = list(cond_species | cll_species)
    
    cleavage_products = []

    while True:
        new_cleavage_products = generate_cleavage_reactions(
            data["catalyzers"], 
            new_species, 
            data["reactions"]["clls"]
        )
        
        cleavage_products.extend(new_cleavage_products)
        
        new_species_copy = data["species"]
        
        cleavage_products = [
            product for product in cleavage_products 
            if product[0] not in {species[0] for species in data["species"]}
        ]
        
        new_species_to_add = list(set((product[0], '1.00E-15', '0.') for product in cleavage_products))
        new_species_to_add = [list(item) for item in new_species_to_add]
        
        data["species"].extend(new_species_to_add)
        
        if new_species_copy == data["species"]:
            break

    return new_species

def eliminate_duplicate_reactions(reactions):
    reaction_dict = {}
    
    for reaction in reactions:
        reagents_products = (reaction[1], reaction[2], reaction[0])
        reagents_products = tuple(sorted(reagents_products[:2])) + (reagents_products[2],)
        
        if reagents_products not in reaction_dict:
            reaction_dict[reagents_products] = [reaction]
        else:
            reaction_dict[reagents_products].append(reaction)
    
    unique_reactions = []
    for key in reaction_dict:
        unique_reactions.append(random.choice(reaction_dict[key]))
    
    return unique_reactions


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
        species_names = [specie[0] for specie in parsed_data["species"]]
        parsed_data["cll_reactions"] = generate_cleavage_reactions(parsed_data["catalyzers"],
                                                           species_names,
                                                           parsed_data["reactions"]["clls"])
        
        parsed_data["cond_reactions"] = eliminate_duplicate_reactions(parsed_data["cond_reactions"])
        parsed_data["cll_reactions"] = eliminate_duplicate_reactions(parsed_data["cll_reactions"])


        generate_new_species(parsed_data)
        generatorIO.write_data(parsed_data)
    except Exception as e:
        print("An error occurred:", str(e))
        sys.exit(1)