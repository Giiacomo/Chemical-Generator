import sys

class BaseIO:
    def __init__(self, file_path):
        self.input_file = file_path
        self.output_file = file_path.replace(".txt", "") + "_output.txt"

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

            self._validate_system_params(data.get('system', {}))
        
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

    def _validate_system_params(self, system_data):
        required_params = {'ML', 'CLL_ML_ACTIVE', 'D_CONCENTRATION', 'D_CONTRIB'}
        if not required_params.issubset(system_data.keys()):
            missing = required_params - system_data.keys()
            raise ValueError(f"Error!\n Missing system parameters: {', '.join(missing)}")

    def _parse_system_param(self, line, system_data):
        parts = line.split()
        possible_params = ['ML', 'CLL_ML_ACTIVE', 'D_CONCENTRATION', 'D_CONTRIB']
        if len(parts) != 2 or parts[0] not in possible_params:
            raise ValueError("Error!\nInvalid parameter. Check the documentation to understand more about system parameters!")

        system_data[parts[0]] = parts[1]

    def _parse_new_species_params(self, line):
        parts = line.split()
        if len(parts) != 4:
            raise ValueError("Error!\nNEW_SPECIES_PARAMS correct form:\n<class lengths> <p_cata_cond> <p_cata_cll> <specificity>")
        classes = parts[0].split(',')
        p_cata_cond = float(parts[1])
        p_cata_cll = float(parts[2])
        specificity = int(parts[3])
        return {cls: {'p_cata_cond': p_cata_cond, 'p_cata_cll': p_cata_cll, 'specificity': specificity} for cls in classes}

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
        elif line.startswith('SYSTEM'):
            current_section = 'system'
            data[current_section] = {}
        elif line.startswith('NEW_SPECIES_PARAMS'):
            current_section = 'new_species_params'
            data[current_section] = {}
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
        elif current_section == 'system':
            self._parse_system_param(line, data['system'])
        elif current_section == 'new_species_params':
            new_params = self._parse_new_species_params(line)
            data[current_section].update(new_params)
        return data, catalyzer_params_counter


    def _parse_reactions(self, line, data):
        self.initial_species_count = len(data["species"])

        parts = [part.replace("R-", "").replace("-R", "") for part in line.split()]

        if len(parts) != 3:
            raise ValueError(
                "Error!\nCondensation correct form:\n<reactant_1> <reactant_2> <reaction_speed>\nCleavage correct form:\n<reactant> <reaction_speed> <n_split>"
            )

        try:
            reaction_speed = float(parts[1])
            
            try:
                n_split = int(parts[2])
                if n_split >= len(parts[0]):
                    raise RuntimeError("Error!\nThe value of the third parameter for cleavage reaction class must be less than the size of the defined string\nFor example R-AABBA-R, n_split must be < 5!")
                data["reactions"]["clls"].append(parts)
            except ValueError:
                raise ValueError(
                    "Error!\nCleavage correct form:\n<reactant> <reaction_speed> <n_split>"
                )
        except ValueError:
            data["reactions"]["conds"].append(parts)
            return



    


    def write_data(self, data):
        with open(self.output_file, 'w') as file:
            counter_species = 0
            max_lengths = [0] * len(data["species"][0])
            for specie in data["species"]:
                for i, item in enumerate(specie):
                    max_lengths[i] = max(max_lengths[i], len(str(item)))

            for specie in data["species"]:
                formatted_specie = [str(item).ljust(max_lengths[i]) for i, item in enumerate(specie)]
                file.write("\t".join(formatted_specie) + "\n")
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

            if self.debug:
                self.print_info(data)
                file.write("\n")
                file.write("# Debug Summary:\n")
                file.write(f"# Total new species generated: {self.new_species_count}\n")
                file.write(f"# Total cond reaction lines written: {self.counter_cond}\n")
                file.write(f"# Total cll reaction lines written: {self.counter_cll}\n")
    
    def print_info(self, data):
        print("The chemical file has been generated. Here's some info!")

        self.new_species_count = len(data["species"]) - self.initial_species_count
        print(f"{self.new_species_count} new species have been generated:")
        #print(", ".join([species[0] for species in data["species"][self.initial_species_count:]]))
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
            print(f"\t- {catalyzer['catalyzer_specie']} is assigned to reaction\t⟶\tR-{reaction[0]} + {reaction[1]}-R")

        print("\nCleavage reactions:")
        for catalyzer in data["catalyzers"]["eligible_cll_species"]:
            reaction = catalyzer['reaction']
            print(f"\t- {catalyzer['catalyzer_specie']} is assigned to reaction\t⟶\tR-{reaction[0]}-R")




class GenToolIO(BaseIO):

    def _process_line(self, line, data, current_section, catalyzer_params_counter):
        if line.startswith('SYSTEM'):
            current_section = 'system'
            data[current_section] = {}
        elif line.startswith('SPECIES'):
            current_section = 'species'
            data[current_section] = []
        elif line.startswith('NEW_SPECIES_PARAMS'):
            current_section = 'new_species_params'
            data[current_section] = {}
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
        if current_section == 'system':
            self._parse_system_param(line, data['system'])
        elif current_section == 'species':
            data[current_section].append(self._parse_species(line))
        elif current_section == 'new_species_params':
            new_params = self._parse_new_species_params(line)
            data[current_section].update(new_params)
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

            file.write("SYSTEM\n")
            for key, value in data['system'].items():
                file.write(f"{key}\t{value}\n")
            file.write("\n")

            file.write("SPECIES\n")
            for specie in data["species"]:
                file.write(" ".join(specie) + "\n")

            file.write("\nCATALYZER_PARAMS\n")
            first_cata_line = ",".join(map(str, data['catalyzer_params'][0]))
            file.write(first_cata_line + "\n")  
            for param in data['catalyzer_params'][1:]:
                file.write(str(param) + "\n")       


            file.write("\nREACTIONS\n")

            for r in data["gen-conds"]["reactions"]:
                file.write(f'R-{r[0]}\t{r[1]}-R\t{data["gen-conds"]["v"]}\n')
            file.write("\n")

            for r in data["gen-clls"]["reactions"]:
                file.write(f'R-{r}-R\t{data["gen-clls"]["v"]}\t 1\n')
            print(f"\nThe chemical {self.output_file} has been generated and is ready to be used as input to the actual chemical generator!\n")