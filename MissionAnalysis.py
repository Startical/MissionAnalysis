from Model.ConstellationDesign import Constellation
import numpy as np
from Scripts.ParametricAnalysis import *
from Scripts.ConstellationCoverageAnalysis import *

if __name__ == "__main__":

    folder = 'results'
    save_report = True

    ## Parametric analysis
    runParametricAnalysis(save_report, folder)
    
    ## Evaluate constellations
    constellations = []
    #constellations.append(Constellation("case0",  650, 80*np.pi/180, 10, 28, "star", 65*np.pi/180, 18))
    #constellations.append(Constellation("case1a", 650, 80*np.pi/180, 9, 20, "star", 65*np.pi/180, 18))
    #constellations.append(Constellation("case1b", 650, 80*np.pi/180, 6, 20, "star", 65*np.pi/180, 18))
    #constellations.append(Constellation("case2a", 700, 80*np.pi/180, 6, 15, "star", 65*np.pi/180, 18))
    #constellations.append(Constellation("case2b", 700, 80*np.pi/180, 6, 12, "star", 65*np.pi/180, 18))

    constellations.append(Constellation("case0",  650, 45*np.pi/180, 6, 20, "delta", 65*np.pi/180, 18))
    constellations.append(Constellation("case1a", 650, 45*np.pi/180, 6, 15, "star" , 65*np.pi/180, 18))
    #constellations.append(Constellation("case1b", 650, 45*np.pi/180, 6, 15, "delta", 65*np.pi/180, 18))
    #constellations.append(Constellation("case2a", 700, 45*np.pi/180, 6, 12, "star" , 65*np.pi/180, 18))
    #constellations.append(Constellation("case2b", 700, 45*np.pi/180, 6, 12, "delta", 65*np.pi/180, 18))
    #constellations.append(Constellation("case3a", 700, 45*np.pi/180, 4, 12, "star" , 65*np.pi/180, 18))
    #constellations.append(Constellation("case3b", 700, 45*np.pi/180, 4, 12, "delta", 65*np.pi/180, 18))

    evaluate_constellations(constellations, save_report, folder)





    
    



