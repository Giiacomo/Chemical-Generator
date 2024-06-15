## GenTool Documentation

### Introduction
GenTool is a tool designed to facilitate the generation of files based on input data of a particular format (example at the end of this file). These files serve as inputs to the `basic_gen` generator.

### Species Format
The input files for GenTool have the following format:
- chemical species, their concentrations, contributions to the membrane
- parameters for condensation and cleavage catalyzers
- conditions for condensation, conds speed
- cleavage windows, clls speed

#### Chemical Species
- Each species is represented by a symbol (e.g.: A, B, AA) along with its concentration and contribution to the membrane.
- Comments can be added using '#' at the beginning of the line.

#### Catalyzers Section
- **first line**: range of the length of a chemical species to become a catalyst, if there is no limit do not include anything
- **second line**:number of chemical species catalyzing a class of condensations
- **third line**: number of chemical species catalyzing a class of cleavage
- **fourth line**: catalyst for condensation and cleavage: ON/OFF

#### Condensation Section
- Allows users to define the number of symbols on the left (N_s) and right (N_d) sides for condensation reactions.
- Allows users to define condensation reactions speed

#### Cleavage Section
- Defines windows of a certain length (N_t, which is a even integer) used for cleavage reactions.
- Allows users to define cleavage reactions speed

### Usage Example
```SPECIES
Cont 1.35E-16 0
A 1.00E-15 0.
B 1.00E-15 0.
AA 1.00E-15 0.
AB 1.00E-15 0.
BA 1.00E-15 0.
BB 1.00E-15 0.
AABB 1.00E-15 1E-3
BABA 1.00E-15 1E-3
AAA 1.00E-15 1E-3
BAB 1.00E-15 1E-3
BBBB 1.00E-15 1E-3
AAAB 1.00E-15 1E-3
AAB 1.00E-15 1E-3

CATALYZER_PARAMS
1,2
3
2
ON

CONDS   
1   2   0.2

#Windows (CN_t) must be even number! (Working on CN_s, CN_t, so that if CN_s = 3, CN_d = 1 and specie ABBA, we'll have ABB + A as cleavage product)
CLLS
4,3   0.1
```