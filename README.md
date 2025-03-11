# MissionAnalysis

Mission Analysis Tools for Spacecraft and Constellations


## Getting started with the development environment

<details>
<summary>Configuration for Visual Studio</summary>
<br>

- STEP 1 - Add your GitHub account to Visual Studio; go to File -> Account Settings -> Add
![Add GitHub account to VS](/doc/resources/VS_config_1.png)

- STEP 2 - Configure the git settings, go to Account options / Source Control

| | |
|:-- |:--|
|![Configure git in VS 1](/doc/resources/VS_config_2.png)|![Configure git in VS 2](/doc/resources/VS_config_3.png)|

(make sure that your user name in the second image is the same as your GitHub user name).

- STEP 3 - Clone the repository, go to File -> Clone repository ...
![Clone repository (1)](/doc/resources/VS_config_4.png)

- Enter the repository URL and the local path where you want to clone the repository
![Clone repository (2)](/doc/resources/VS_config_5.png)

```
https://github.com/Startical/MissionAnalysis.git
```


</details>

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


