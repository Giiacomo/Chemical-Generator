# Chemical Generator


## Usage 🚀

The basic version of this chemical generator is implemented in the script `generator.py`. You can execute it using the following command:

```sh
python3 generator.py <input-filename>.txt [-debug] [-o <output-filename>.txt]
```

where:
- **-debug** enables debug info printing
- **-o <filename>** enables the choice of output file name

The script executes in order the following functions:
```py
    parsed_data = parse_data_from_file(file_path)
    parsed_data["catalyzers"] = generate_catalyzers(parsed_data)
    parsed_data["cond_reactions"] = generate_condensation_reactions(parsed_data)
    parsed_data["cll_reactions"] = generate_cleavage_reactions(parsed_data)
    generate_new_species(parsed_data)
    generatorIO.write_data(parsed_data)
```

At the end of the script execution, `parsed_data` takes the following structure:
```py
{
'species':  [[<nomespecie>, <concentrazione>, <contributo>]],
'catalyzer_params': [[<range>], <n_cond_catalyzers>, <n_cll_catalyzers>, <both_on>],  
'reactions': {'conds': [<specie>], 'clls': [<specie>] },
'catalyzers': {'cond': [<specie>], 'cll': [<specie>] },
'cond_reactions': [<reactant_1>, <reactant_2>, <v>, [<catalyzers>]]
'cll_reactions': [<specie>, <cleavage_1>, <cleavage_2>, <v>, [<catalyzers>]]
}
```


### Usage Example ✏️

```
SPECIES
#FORMA 
#PRIMA RIGA = <CONTENITORE>
#PRIMA COLONNA = <NOME SPECIE>
#SECONDA COLONNA = <CONCENTRAZIONE SPECIE IN PROTOCELLA>
#TERZA COLONNA = <CONTRIBUTO AL CONTENITORE> 
Cont	1.35E-16	0				
A	    1.00E-15	0.
B	    1.00E-15	0.
AA	    1.00E-15	0.
AB	    1.00E-15	0.
BA	    1.00E-15	0.
BB	    1.00E-15	0.
AABB	1.00E-15	1E-3
BABA    1.00E-15	1E-3
AAA	    1.00E-15	1E-3
BAB	    1.00E-15	1E-3
BBBB    1.00E-15	1E-3
AAAB    1.00E-15	1E-3
AAB	    1.00E-15	1E-3

CATALYZER_PARAMS
#FORMA
#PRIMA RIGA = <range della lunghezza di una specie chimica per diventare catalizzatore, se non c'è limite non inserire nulla>
#SECONDA RIGA = <numero di specie chimiche catalizzatrici di una classe di condensazioni>
#TERZA RIGA = <numero di specie chimiche catalizzatrici di una classe di cleavage>
#QUARTA RIGA = <catalizzatore di condensazione e cleavage: ON/OFF (both_on)>
1,2
3
2
ON

REACTIONS
#FORMA CONDS
#REAGENTE_1 = <stringa qualsiasi che termina in AB>
#REAGENTE_2 = <stringa qualsiasi che inizia con AB> 
#V_i = <velocità associata alla condensazione>
R-AB    BA-R    0.1 
R-BA    AB-R    0.2

#FORMA CLLS
#REAGENTE = <stringa generica, pattern, stringa generica>
#V_i = <velocità associata alla condensazione>
R-AB-R  0.3
R-BB-R  0.5
```

## Dev Notes 🛠️

The generator is composed of a class `ReactionGenerator`, that has a constructor that saves the data parsed by the GeneratorIO class. It also has a `run_generation` method that executes all the step in order to obtain the generated reactions and the new species.

Let's see in detail how it works!

First of all, it generates the catalyzers following these steps:
- In order for a specie to be acatalyzers it needs to be in the length range defined in the input file
- It first chooses the catalyzers for condensation:
    1. If there are enough catalyzer species for each class of reaction, every class of reaction will happen. For each reaction a catalyzer specie is extracted randomly to be the chosen catalyzer for that reaction, and is removed from the pool of extractable catalyzer species.
    2. otherwise **n** class of reactions will be randomly chosen to happen (where **n** is the number of species that serve as catalyzers), then the same happens for this case, so for each chosen reaction, a catalyzer specie will be assigned.
    3. For both cases this will happen: if the chosen catalyzers are less than 
    