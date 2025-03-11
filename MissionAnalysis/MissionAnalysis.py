from Model.ConstellationDesign import Constellation
import numpy as np
from Scripts.ParametricAnalysis import *
from Scripts.ConstellationCoverageAnalysis import *

if __name__ == "__main__":

    ## Parametric analysis
    runParametricAnalysis();

    ## Evaluate constellations
    constellations = []
    constellations.append(Constellation("case0", 650, 80*np.pi/180, 10, 28, "star", 65*np.pi/180, 18))
    constellations.append(Constellation("case1a", 650, 80*np.pi/180, 9, 20, "star", 65*np.pi/180, 18))
    constellations.append(Constellation("case1b", 650, 80*np.pi/180, 6, 20, "star", 65*np.pi/180, 18))
    constellations.append(Constellation("case2a", 700, 80*np.pi/180, 6, 15, "star", 65*np.pi/180, 18))
    constellations.append(Constellation("case2b", 700, 80*np.pi/180, 6, 12, "star", 65*np.pi/180, 18))

    evaluate_constellations(constellations)





    
    



