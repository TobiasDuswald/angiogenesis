# Paper Experiments

This README file gives an overview of the different experiments and how to
replicate the paper results. To reproduce, you will need to clone three
repositories and they need to be structured in a specific way. The three
repositories are:

1. `biodynamo/biodynamo`
2. `LukasBreitwieser/bdm-paper-examples`
3. `TobiasDuswald/angiogenesis`

## Project structure

To set up the project structure, proceed as follows. Go to some `<DIR>` of
choice and executed
```bash
cd <DIR>
mkdir angiogenesis-simulation && cd angiogenesis-simulation
git clone https://github.com/BioDynaMo/biodynamo.git
cd biodynamo
git checkout <fix-some-commit> #ToDo
cd ..
git clone https://github.com/LukasBreitwieser/bdm-paper-examples.git
cd bdm-paper-examples/bdm/
rm -rf *
git clone https://github.com/TobiasDuswald/angiogenesis.git
```

## Docker environment

BioDynaMo is the simulation backend used for all experiments, bdm-paper-examples
provides access to a docker container that allows us to provide the right
environment for the simulations. Having `docker` is, consequently, a mandatory
prerequisite. It is equally well possible to run the simulations bare metal
but you'll have to deal with dependencies yourself.

To launch the docker container:
```bash
cd <DIR>/bdm-paper-examples
docker/build.sh # only once
docker/run.sh
docker/exec.sh
```

## Reproducing the results

After starting your docker container, simply execute the scripts in
`scripts/experiments` to run the simulations:
```bash
pwd # should end with `.../bdm-paper-examples`
cd bdm/angiogenesis
./scripts/experiments/run_<experiment>.sh
```

## Available experiments

1. **Avascular tumor growth:** experiment corresponding to Fig. XY in the paper.
  Run with `./scripts/experiments/run_avascular_spheroid.sh`
2. **Porous tumor core:** experiment corresponding to Fig. XY in the paper.
  Run with `./scripts/experiments/run_porous_tumor_core.sh`
3. **Treatment of spheroid:** experiment corresponding to Fig. XY in the paper.
  Run with `./scripts/experiments/run_spheroid_treatment.sh`
4. **Simplified growth model:** experiment corresponding to Fig. XY in the paper.
  Run with `./scripts/experiments/run_simplified_growth.sh`
5. **Vessels growing to center:** experiment corresponding to Fig. XY in the paper.
  Run with `./scripts/experiments/run_vessel_to_center.sh`
6. **Vessel coupling:** experiment corresponding to Fig. XY in the paper.
  Run with `./scripts/experiments/run_vessel_coupling.sh`
7. **Full scale model:** experiment corresponding to Fig. XY in the paper.
  Run with `./scripts/experiments/run_full_scale_model.sh`
