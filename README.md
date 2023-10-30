# Bridging Scales

This is the code base for the publication:
**Bridging scales: A hybrid model to simulate vascular tumor growth and treatment response**
by Duswald, Lima, Oden, and Wohlmuth (2023) available [here](https://doi.org/10.1016/j.cma.2023.116566).

![05-treatment](https://github.com/TobiasDuswald/angiogenesis/assets/44771875/bb8312a7-0b5d-486b-b28c-d57dd7dd22b7)

## Software Stack and Depndencies

*  BioDynaMo. To install, consult the 
  [BDM dev-guide](https://biodynamo.org/docs/devguide/build/). Install the 
  latest master as described in the link. Following the dev-guide allows us to
  easily upgrade and extend BioDynaMo if necessary.

## Running the simulations

Building on the BioDynaMo platform, the commands to running the simulations
assume a functional installation of BioDynaMo and are
```
source <path>/thisbdm.sh  # source biodynamo 
bdm build                 # compile application code
bdm run                   # execute application code
bdm view                  # open results in ParaView
```

## Reproducing the results of the publication

Reproducing the results of the different experiments of the publication
requires setting the correct parameters but also minimal changes to code, i.e.,
possibly adding a few lines of code to adequately set the initial conditions.
This process it automated in a second repository that can be found at
`TobiasDuswald/bdm-angiogenesis-reproducer`. Please consult the
[README.md](https://github.com/TobiasDuswald/bdm-angiogenesis-reproducer/blob/main/README.md)
of that [repository](https://github.com/TobiasDuswald/bdm-angiogenesis-reproducer)
for further instructions.
