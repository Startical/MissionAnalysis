# MissionAnalysis

Mission Analysis Tools for Spacecraft and Constellations


## Getting started with the development environment

Check [README](https://github.com/Startical/.github-private/blob/main/profile/README.md) in Startical's organization landing page

## Repo structure

The `MissionAnalysis` source folder contains: 
- OrbitTools - Utilities for orbit definition, frame transformation and analysis of orbits
- Model - Main classes used for defining the constellation and Spacecraft
- Scripts - source files of the analysis supported by the tool
- MissionAnalysis.py - the main script to run the analysis.


## Run the sample script

Run `MissionAnalysis.py` in the source folder

Currently two kind of analysis are supported: 
- Parametric analysis, to analyse key drivers for the constellation geometry and analyse geometric coverage
- Coverage analysis, starting from the definition of constellation candidates, analysing the coverage at the target flight level

Both analysis generate an automatic report that is stored in the `results` folder (will be created automatically by the script if not existing).


