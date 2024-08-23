# Roadmap of ToddGPT

## TeraChem
- [ ] Set up terachem input file
- [x] Run terachem input file
- [ ] Get data from terachem output using qcio
- [ ] May have to make own parsers since qcio may not support all of terachem output styles

## Optimization
- [ ] Optimize from Frank-Condon point
- [ ] Find MECI
- [ ] Perform NEB from FC minima to MECI Geom to get barrier height
- [ ] Calculate lifetime using boltzmann/RRKM?

## Suggest Methods for Excited State Dynamics
- [ ] Analyze literature references to find methods used for similar systems
- [ ] Run multiple optimizations on FC minima and MECI geom
- [ ] Run NEB on new FC minima and MECI geom using various methods
- [ ] Compare results and suggest methods

## FMS24
- [ ] Once FMS24 is up and running, create necessary input files for FMS24