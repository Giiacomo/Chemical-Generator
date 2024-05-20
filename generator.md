# Chemical Generator

## Basic Version

The basic version of this chemical generator is implemented in the script `basic_gen.py`. You can execute it using the following command:

```sh
python3 basic_gen.py <filename>.txt
```

The script executes in order the following functions:
```py
    parsed_data = parse_data_from_file(file_path)
    parsed_data["catalyzers"] = generate_catalyzers(parsed_data)
    parsed_data["cond_reactions"] = generate_condensation_reactions(parsed_data)
    parsed_data["cll_reactions"] = generate_cleavage_reactions(parsed_data)
    write_data_to_file(file_path=file_to_write, data=parsed_data)
```

At the end of the script execution, `parsed_data` takes the following structure:
```py
{
'species':  [[<nomespecie>, <concentrazione>, <contributo>]],
'probs': [<prob_catalyzer>, <prob_cond>],  
'reactions': {'conds': [<specie>], 'clls': [<specie>] },
'catalyzers': {'cond': [<specie>], 'cll': [<specie>],
               'cond_reactions': [<reactant_1>, <reactant_2>, <v>]
               'cll_reactions': [<specie>, <cleavage_1>, <cleavage_2>, <v>]}
}
```


## Usage Example

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

PROBS
#FORMA
#PRIMA RIGA = <probabilità che una specie chimica diventi catalizzatore>
#SECONDA RIGA = <probabilità che una reazione sia una condensazione>
0.1 
0.5

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