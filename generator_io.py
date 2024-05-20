import sys

class BaseIO:
    def __init__(self, file_path):
        self.input_file = file_path
        self.output_file = file_path.replace(".txt", "") + "_generated.txt"

class GeneratorIO(BaseIO):
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


    def _parse_reactions(self, line, data):
        parts = [part.replace("R-", "").replace("-R", "") for part in line.split()]
        if len(parts) == 3:
            data["reactions"]["conds"].append(parts)
        elif len(parts) == 2:
            data["reactions"]["clls"].append(parts)
        else:
            raise ValueError("Error!\nCondensation correct form:\n<reactant_1> <reactant_2> <reaction_speed>\nCleavage correct form:\n<reactant> <reaction_speed>")


    def _parse_reactions(self, line, data):
        parts = [part.replace("R-", "").replace("-R", "") for part in line.split()]
        if len(parts) == 3:
            data["reactions"]["conds"].append(parts)
        elif len(parts) == 2:
            data["reactions"]["clls"].append(parts)
        else:
            raise ValueError("Error!\nCondensation correct form:\n<reactant_1> <reactant_2> <reaction_speed>\nCleavage correct form:\n<reactant> <reaction_speed>")

    def write_data(self, data):
        with open(self.output_file, 'w') as file:
            for specie in data["species"]:
                file.write(" ".join(specie) + "\n")

            #Non ci sono catalizzatori
            if (len(data["catalyzers"]) < 1):
                return

            file.write("\n")

            for catalyzer in data["catalyzers"]["cond"]:
                for r in data["cond_reactions"]:  
                    file.write(r[1] + " + " + r[2] + " + " + catalyzer + " > " + r[0] + " + " + catalyzer + " ; " + r[3] + "\n")

            file.write("\n")

            for catalyzer in data["catalyzers"]["cll"]:
                for r in data["cll_reactions"]:  
                    file.write(r[0] + " + " + catalyzer + " > " + r[1] + " + " + r[2] + " + " + catalyzer + " ; " + r[3] + "\n")
            print(f"The chemical {self.output_file} has been generated, with {len(data['cond_reactions'])*len(data['catalyzers']['cond'])} condensation reactions and {len(data['cll_reactions'])*len(data['catalyzers']['cll'])} cleavage reactions!")
            
            print("\nCondensation catalyzers:")
            if len(data["catalyzers"]["cond"]) > 0:
                print(", ".join(data["catalyzers"]["cond"]))
            else:
                print("None")

            print("\nCleavage catalyzers:")
            if len(data["catalyzers"]["cll"]) > 0:
                print(", ".join(data["catalyzers"]["cll"]))
            else:
                print("None")


# Class to handle AutoTool input-output

class AutoToolIO(BaseIO):

    def _parse_species_section(self, line, data):
        parts = line.split()
        if len(parts) != 3:
            raise ValueError("Error! Species correct form: <speciename> <concentration> <contribution>")
        data['species'].append(parts)

    def _parse_probs_section(self, line, data):
        prob = float(line)
        if not 0 <= prob <= 1:
            raise ValueError("Error! Works with probability values in [0, 1]")
        data['probs'].append(prob)

    def _parse_conds_section(self, line, data):
        parts = line.split()
        if len(parts) != 3:
            raise ValueError("Error! Condensation correct form: <N_s> <N_d> <reaction_speed>")
        data['conds'] = parts

    def _parse_clls_section(self, line, data):
        parts = line.split()
        n_ts = parts[0].split(',')
        if len(parts) != 2:
            raise ValueError("Error! Cleavage correct form: <N_t1,...N_tn> <reaction_speed>")
        data['clls'] = [n_ts, parts[1]]

    def parse_data(self):
        data = {'species': [], 'probs': [], 'conds': [], 'clls': []}
        current_section = None
        
        section_parsers = {
            'SPECIES': self._parse_species_section,
            'Cont': self._parse_species_section,
            'PROBS': self._parse_probs_section,
            'CONDS': self._parse_conds_section,
            'CLLS': self._parse_clls_section
        }
        
        try:
            with open(self.input_file, 'r') as file:
                for line in file:
                    line = line.strip()
                    
                    if not line or line.startswith('#'):
                        continue
                    
                    if line.startswith(('SPECIES', 'Cont', 'PROBS', 'CONDS', 'CLLS')):
                        current_section = line.split()[0]
                        if current_section not in section_parsers:
                            raise ValueError("Error! Unknown section: {}".format(current_section))
                        continue
                    
                    if current_section:
                        section_parsers[current_section](line, data)
                    
        except FileNotFoundError:
            print("Error: File not found.")
            sys.exit(1)
        except ValueError as e:
            print(e)
            sys.exit(1)
        except Exception as e:
            print("Error:", e)
            sys.exit(1)

        return data

    def write_data(self, data):
        with open(self.output_file, 'w') as file:
            file.write("SPECIES\n")
            for specie in data["species"]:
                file.write(" ".join(specie) + "\n")

            file.write("\nPROBS\n")
            file.write(f"{data['probs'][0]}\n{data['probs'][1]}")

            file.write("\nREACTIONS\n")

            for r in data["gen-conds"]["reactions"]:
                file.write(f'R-{r[0]}\t{r[1]}-R\t{data["gen-conds"]["v"]}\n')

            file.write("\n")

            for r in data["gen-clls"]["reactions"]:
                file.write(f'R-{r}-R\t{data["gen-clls"]["v"]}\n')

            print(f"The chemical {self.output_file} has been generated and is ready to be used as input to the actual chemical generator!")
            
            