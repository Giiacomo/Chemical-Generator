import sys

class BaseIO:
    def __init__(self, file_path):
        self.input_file = file_path
        self.output_file = file_path.replace(".txt", "") + "_generated.txt"

    def parse_data(self):
        data = {}
        current_section = None
        catalyzer_params_counter = 0
        try:
            with open(self.input_file, 'r') as file:
                for line in file:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    data, current_section, catalyzer_params_counter = self._process_line(line, data, current_section, catalyzer_params_counter)
        except FileNotFoundError:
            print("Error: File not found.")
            sys.exit(1)
        except ValueError as e:
            print(e)
            sys.exit(1)
        except Exception as e:
            print("Error:", str(e))
            sys.exit(1)
        return data

    def _parse_species(self, line):
        parts = line.split()
        if len(parts) != 3:
            raise ValueError("Error!\nSpecies correct form:\n<speciename> <concentration> <contribution>")
        return parts

    def _parse_catalyzer_param(self, line, catalyzer_params_counter):
        try:
            if catalyzer_params_counter == 0:
                value = list(map(int, line.split(',')))
                if int(value[0]) < 1:
                    raise ValueError("Error!\nMinimum length of a chemical species to become a catalyst must be at least 1.")
            elif catalyzer_params_counter in [1, 2]:
                value = int(line)
                if value < 0:
                    raise ValueError("Error!\nThe number of catalyst species must be non-negative.")
            elif catalyzer_params_counter == 3:
                if line not in ['ON', 'OFF']:
                    raise ValueError("Error!\nCondensation and cleavage catalyst must be either 'ON' or 'OFF'.")
                value = line
            else:
                raise ValueError("Error!\nUnexpected number of parameters in the CATALYZERS section.")
        except ValueError:
            raise ValueError("Error!\nCatalyzer values must be integers for the first three lines.")
        return value

class GeneratorIO(BaseIO):

    def __init__(self, file_path, debug=False, output_file=None):
        super().__init__(file_path)  
        self.debug = debug
        if output_file:
            self.output_file = output_file


    def _process_line(self, line, data, current_section, catalyzer_params_counter):
        if line.startswith('SPECIES'):
            current_section = 'species'
            data[current_section] = []
        elif line.startswith('CATALYZER_PARAMS'):
            current_section = 'catalyzer_params'
            data[current_section] = []
        elif line.startswith('REACTIONS'):
            current_section = 'reactions'
            data[current_section] = {"conds": [], "clls": []}
        else:
            data, catalyzer_params_counter = self._read_data(line, data, current_section, catalyzer_params_counter)
        return data, current_section, catalyzer_params_counter

    def _read_data(self, line, data, current_section, catalyzer_params_counter):
        if current_section == 'species':
            data[current_section].append(self._parse_species(line))
        elif current_section == 'catalyzer_params':
            data[current_section].append(self._parse_catalyzer_param(line, catalyzer_params_counter))
            catalyzer_params_counter += 1
        elif current_section == 'reactions':
            self._parse_reactions(line, data)
        return data, catalyzer_params_counter

    def _parse_reactions(self, line, data):
        #info
        self.initial_species_count = len(data["species"])

        parts = [part.replace("R-", "").replace("-R", "") for part in line.split()]
        if len(parts) == 3:
            data["reactions"]["conds"].append(parts)
        elif len(parts) == 2:
            data["reactions"]["clls"].append(parts)
        else:
            raise ValueError("Error!\nCondensation correct form:\n<reactant_1> <reactant_2> <reaction_speed>\nCleavage correct form:\n<reactant> <reaction_speed>")

    def write_data(self, data):
        with open(self.output_file, 'w') as file:
            counter_species = 0
            for specie in data["species"]:
                file.write(" ".join(specie) + "\n")
                counter_species += 1

            file.write("\n")

            self.counter_cond = 0
            for r in data["cond_reactions"]:  
                for catalyzer in r[4]:
                    file.write(r[1] + " + " + r[2] + " + " + catalyzer['catalyzer_specie'] + " > " + r[0] + " + " + catalyzer['catalyzer_specie'] + " ; " + r[3] + "\n")
                    self.counter_cond += 1

            file.write("\n")

            self.counter_cll = 0
            for r in data["cll_reactions"]:
                for catalyzer in r[4]:
                    file.write(r[0] + " + " + catalyzer['catalyzer_specie'] + " > " + r[1] + " + " + r[2] + " + " + catalyzer['catalyzer_specie'] + " ; " + r[3] + "\n")
                    self.counter_cll += 1

            #Debug info
            if self.debug:
                self.print_info(data)
    
    def print_info(self, data):
        print("The chemical file has been generated. Here's some info!")

        new_species_count = len(data["species"]) - self.initial_species_count
        print(f"{new_species_count} new species have been generated:")
        print(", ".join([species[0] for species in data["species"][self.initial_species_count:]]))
        print()
        print(f"{self.counter_cond} condensation reactions have been generated")
        print(f"{self.counter_cll} cleavage reactions have been generated")
        print()
        print("Condensation catalyzers for this chemical are: " )
        out = ""
        for catalyzer in data["catalyzers"]["eligible_cond_species"]:
            out += f"{catalyzer['catalyzer_specie']}\t"
        print(f'{out}\n')

        print("Cleavage catalyzers for this chemical are: " )
        out = ""
        for catalyzer in data["catalyzers"]["eligible_cll_species"]:
            out += f"{catalyzer['catalyzer_specie']}\t"
        print(f'{out}\n')

        print("Assigned reactions for each catalyzer:")
        print("Condensation reactions:")
        for catalyzer in data["catalyzers"]["eligible_cond_species"]:
            reaction = catalyzer['reaction']
            print(f"{catalyzer['catalyzer_specie']} is assigned to reaction:\t R-{reaction[0]} + {reaction[1]}-R")

        print("\nCleavage reactions:")
        for catalyzer in data["catalyzers"]["eligible_cll_species"]:
            reaction = catalyzer['reaction']
            print(f"{catalyzer['catalyzer_specie']} is assigned to reaction:\t R-{reaction[0]}-R")




# Class to handle AutoTool input-output
class AutoToolIO(BaseIO):

    def _process_line(self, line, data, current_section, catalyzer_params_counter):
        
        if line.startswith('SPECIES'):
            current_section = 'species'
            data[current_section] = []
        elif line.startswith('CATALYZER_PARAMS'):
            current_section = 'catalyzer_params'
            data[current_section] = []
        elif line.startswith('CONDS'):
            current_section = 'conds'
            data[current_section] = []
        elif line.startswith('CLLS'):
            current_section = 'clls'
            data[current_section] = []
        else:
            data, catalyzer_params_counter = self._read_data(line, data, current_section, catalyzer_params_counter)
        return data, current_section, catalyzer_params_counter

    def _read_data(self, line, data, current_section, catalyzer_params_counter):
        if current_section == 'species':
            data[current_section].append(self._parse_species(line))
        elif current_section == 'catalyzer_params':
            data[current_section].append(self._parse_catalyzer_param(line, catalyzer_params_counter))
            catalyzer_params_counter += 1
        elif current_section in ['conds', 'clls']:
            self._parse_conditions_or_cleavage(line, data, current_section)
        return data, catalyzer_params_counter

    def _parse_conditions_or_cleavage(self, line, data, current_section):
        parts = line.split()
        if current_section == 'conds' and len(parts) == 3:
            data[current_section] = parts
        elif current_section == 'clls' and len(parts) == 2:
            nt_s = parts[0].split(',')
            data[current_section] = [nt_s, parts[1]]
        else:
            raise ValueError(f"Error!\n{current_section} correct form:\n"
                             f"CONDS: <reactant_1> <reactant_2> <reaction_speed>\n"
                             f"CLLS: <reactant> <reaction_speed>")

    def write_data(self, data):
        with open(self.output_file, 'w') as file:
            file.write("SPECIES\n")
            for specie in data["species"]:
                file.write(" ".join(specie) + "\n")
            print(data["catalyzer_params"])

            file.write("\nCATALYZER_PARAMS\n")
            first_cata_line = ",".join(map(str, data['catalyzer_params'][0]))
            file.write(first_cata_line + "\n")  # Write the first catalyzer parameter line
            for param in data['catalyzer_params'][1:]:
                file.write(str(param) + "\n")       


            file.write("\nREACTIONS\n")

            for r in data["gen-conds"]["reactions"]:
                file.write(f'R-{r[0]}\t{r[1]}-R\t{data["gen-conds"]["v"]}\n')
            file.write("\n")

            for r in data["gen-clls"]["reactions"]:
                file.write(f'R-{r}-R\t{data["gen-clls"]["v"]}\n')
            print(f"The chemical {self.output_file} has been generated and is ready to be used as input to the actual chemical generator!")