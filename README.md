# andy_spie_journal

### `build` contains Dockerfile for running meep simulations

### `src` contains python scripts for generating meep data
    _3x3Pillars.py contains the class which serves as a wrapper for meep, specifically for running our simulations.
    `utils` contains the parameter manager, which manages all of our meta atom parameters, which are initialized in config.yaml.
    
### To generate data, `src/3x3_pillar_sims/gen_data.py` line 9 contains the command for generating a single sample. 
