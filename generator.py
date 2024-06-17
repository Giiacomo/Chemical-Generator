import sys
import random
import argparse
from generator_io import GeneratorIO

class ReactionGenerator:
    def __init__(self, data):
        self.species = data["species"]
        self.reactions = data["reactions"]
        self.catalyzer_params = data["catalyzer_params"]
        self.catalyzers = None
        self.cond_reactions = None
        self.cll_reactions = None

    def assign_catalyzers(self, num_catalyzers, eligible_species, reactions):
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

    def generate_catalyzers(self):
        catalyzer_params = self.catalyzer_params
        species = self.species
        cond_reactions = [tuple(r) for r in self.reactions["conds"]]
        cll_reactions = [tuple(r) for r in self.reactions["clls"]]

        min_length, max_length = catalyzer_params[0]
        num_cond_catalyzers = catalyzer_params[1]
        num_cll_catalyzers = catalyzer_params[2]
        both_on = catalyzer_params[3] == 'ON'

        eligible_species = [specie[0] for specie in species[1:] if min_length <= len(specie[0]) <= max_length]
        
        if len(eligible_species) < (num_cond_catalyzers + num_cll_catalyzers):
            raise ValueError("Error! Not enough eligible species to satisfy the catalyzer requirements.")
        
        self.catalyzers = {
            'eligible_cond_species': [],
            'eligible_cll_species': [],
            'both_on': both_on,
            'num_cond_catalyzers': num_cond_catalyzers,
            'num_cll_catalyzers': num_cll_catalyzers
        }

        self.catalyzers['eligible_cond_species'] = self.assign_catalyzers(num_cond_catalyzers, eligible_species, cond_reactions)

        if both_on:
            cll_candidates = eligible_species
        else:
            cll_candidates = [s for s in eligible_species if s not in [c['catalyzer_specie'] for c in self.catalyzers['eligible_cond_species']]]

        self.catalyzers['eligible_cll_species'] = self.assign_catalyzers(num_cll_catalyzers, cll_candidates, cll_reactions)

    def get_reaction_catalyzer(self, reaction):
        reaction_tuple = tuple(reaction)
        return [c for c in self.catalyzers['eligible_cond_species'] + self.catalyzers['eligible_cll_species'] if c['reaction'] == reaction_tuple]

    def generate_condensation_reactions(self):
        species = self.species[1:]
        reactions = self.reactions["conds"]
        condensation_reactions = []

        for i in range(len(species)):
            for j in range(len(species)):
                reagent_1 = species[i][0]
                reagent_2 = species[j][0]

                for reaction in reactions:
                    if reagent_1.endswith(reaction[0]) and reagent_2.startswith(reaction[1]):
                        new_reaction = [reagent_1 + reagent_2, reagent_1, reagent_2, reaction[2]]
                        catalyzer = self.get_reaction_catalyzer(reaction)
                        new_reaction.append(catalyzer)
                        condensation_reactions.append(new_reaction)

        return condensation_reactions

    def generate_cleavage_reactions(self):
        species = [specie[0] for specie in self.species[1:]]
        reactions = self.reactions["clls"]
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
                        catalyzer = self.get_reaction_catalyzer(reaction)
                        new_reaction.append(catalyzer)
                        cleavage_reactions.append(new_reaction)

        return cleavage_reactions

    def generate_new_species(self, cond_reactions, cll_reactions):
        cond_species = {reaction[0] for reaction in cond_reactions}
        cll_species = {product for reaction in cll_reactions for product in reaction[1:3]}
        
        new_species = list(cond_species | cll_species)

        return new_species

    def eliminate_duplicate_reactions(self, reactions):
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

    def run_generation(self):
        self.generate_catalyzers()
        self.cond_reactions = self.generate_condensation_reactions()
        self.cll_reactions = self.generate_cleavage_reactions()

        self.cond_reactions = self.eliminate_duplicate_reactions(self.cond_reactions)
        self.cll_reactions = self.eliminate_duplicate_reactions(self.cll_reactions)

        new_species = self.generate_new_species(self.cond_reactions, self.cll_reactions)
        current_species = {specie[0] for specie in self.species}
        for species in new_species:
            if species not in current_species:
                self.species.append([species, '1.00E-15', '0.'])

        generated_data = {
            "catalyzers": self.catalyzers,
            "cond_reactions": self.cond_reactions,
            "cll_reactions": self.cll_reactions,
            "species": self.species
        }

        return generated_data

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
        generator = ReactionGenerator(parsed_data)
        generated_data = generator.run_generation()

        generatorIO.write_data(generated_data)
    except Exception as e:
        print("An error occurred:", str(e))
        sys.exit(1)
