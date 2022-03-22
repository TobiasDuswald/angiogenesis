# Angiogenesis

## Project description

* Goal: simulate cancer and the related phenomenon angiogenesis with a hybrid
  (ABM+PDE) model.
* Stages: 
  1) Scenario: few hypoxic cancer cells surrounding a vessel. Cells secrete VEGF
     to trigger sprouting and apical growth. The vessel then grows towards the 
     cells and if the nutrient concentration at the hypoxic cells is beyond a 
     threshold (or they are simply close to a vessel), then the cell become 
     quiescent. Vessels are supposed to follow the gradient.
  2) Possible extensions:
     * Include growth of primary tumor
     * Possibly include complicated initial vessel configurations from image 
       data
     * Use more complicated continuum models
* Data: yet to be defined
* Metrics: yet to be defined

## Software stack

*  BioDynaMo. To install, consult the 
  [BDM dev-guide](https://biodynamo.org/docs/devguide/build/). Install the 
  latest master as described in the link. Following the dev-guide allows us to
  easily upgrade and extend BioDynaMo if necessary.

## Running the simulation

BioDynaMo needs to be installed and sourced, e.g.
```bash
source <path_to_biodynamo>/biodynamo/build/bin/thisbdm.sh
```
The simulation can then be executed by simply running 
```bash
bdm run
```
in your project directory. You can verify component of the model by running a 
very similar command
```bash
bdm test
``` 
We outsourced all parameters steering the simulation to `src/sim_param.h` and 
one can modify the parameters in `bdm.json` without the need to recompile the 
simulation. Further ways to feed parameters to the simulation are the following:

1. ```cpp
    // Define parameters
    auto set_param = [&](Param* param) {
      auto* sparam = param->Get<SimParam>();
      param->statistics = true;
      param->calculate_gradients = false;
      sparam->some_value = some_value
    };

    // Initialize the simulation
    Simulation simulation(argc, argv, set_param);
    ```
2. We can use BioDynaMo's CLI, e.g. 
   ```
   ./build/angiogenesis --inline-config '{ "bdm::SimParam": { "x": 6.28 } }'
   ```

The output of the simulation is, by default, available in 
`output/angiogenesis`.

## Visualization

Assuming that one has successfully sourced and executed the simulation, the user
simply follows up with the paraview command"
```bash
source <path>/thisbdm.sh
bdm run
paraview
```
The place where we call `paraview` is indeed important. This opens a ParaView
window. Follow with
```
File >> Load State... >> "output/angiogenesis/angiogenesis.pvsm" >> \
Use File Name from States
```
and the simulation output is ready to view. Simply click play.
