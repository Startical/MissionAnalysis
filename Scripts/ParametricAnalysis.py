import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import webbrowser

from CommonTools.Constants import EARTH_RADIUS
from OrbitTools.CoverageAnalysis import antenna_swath

def runParametricAnalysis(save_report, folder):
    # Parametric inputs
    H = 18; # flight level in km
    h = np.arange(500, 1250, 50) #altitude in kms
    beta = np.arange(70,55,-2.5)*np.pi/180 # antenna aperture in radians
    
    ## STEP 1 - Geometrical coverage
    creq = 1500; # coverage requirement in km
    fig1 = geometricalCoverage(h, beta, H, creq)

    ## STEP 2 - Global coverage
    fig2 = inclination_vs_altitude(h, H)
    table2 = number_of_planes(h, beta, H)
    
    ## STEP 3
    N1=2
    N2=3
    table3 = multiCoverage(h, beta, H, N1, N2)

    ## Print report
    if save_report:
        printReport(folder, 'Constellation_sizing-Parametric_analysis', creq, H, fig1, fig2, table2, table3, N1, N2)

def geometricalCoverage(h, beta, H = 0, creq = 1500):
    """Parametric analysis of the geometrical coverage at a given flight level (H)"""

    # Initialize parameters
    clim = np.zeros(len(h))
    betamax = np.zeros(len(h))
    cmax = np.zeros(len(h))
    betamin = np.zeros(len(h))

    # loop
    for i in range(len(h)):

        p = (EARTH_RADIUS+h[i])/(EARTH_RADIUS+H)

        # lim - antenna tangent to flight level
        betamax[i] = np.arcsin(1/p)
        clim[i] = (EARTH_RADIUS+H)*(np.pi/2-betamax[i])

        # maximum coverage
        alpha, cmax[i] = antenna_swath(h[i],H, np.pi/2)


    # Plot 
    fig = plt.figure()
    
    for j in range(len(beta)):
        c = np.zeros(len(h))
        
        for i in range(len(h)):
            alpha, c[i] = antenna_swath(h[i],H, beta[j])
            
        plt.plot(c,h,label=f'beta = {beta[j]*180/np.pi:.1f} deg')
        
    plt.legend(loc='upper left')
    plt.plot(clim,h,'k--', label='lim')
    plt.plot(cmax,h,'k--', label='cmax')
    plt.plot([creq, creq], [h[0], h[-1]], 'r--', label='req')
    plt.text(clim[-1]-100, h[-1], 'lim', fontsize=12, verticalalignment='bottom')
    plt.text(cmax[-1]-150, h[-1], 'cmax', fontsize=12, verticalalignment='bottom')
    plt.xlabel('Coverage (km)')
    plt.ylabel('Altitude (km)')
    plt.title(f'Geometrical Coverage Analysis (FL {H/3*100:.0f})')
    plt.grid()
    
    return fig;

def inclination_vs_altitude(h, H=0):

    ## Inclination
    fig = plt.figure()

    inc = np.zeros(len(h))

    for i in range(len(h)):

        inc1 = np.pi/2
        res = 1000

        hh = (EARTH_RADIUS+h[i])/(EARTH_RADIUS+H);

        inc[i] = np.arcsin(1/hh)

    plt.plot(h,inc*180/np.pi)

    plt.ylabel('Inclination (deg)')
    plt.xlabel('Altitude (km)')
    plt.title(f'Minimum inclination to ensure global coverage (FL {H/3*100:.0f})')
    plt.grid()

    return fig

def number_of_planes(h, beta, H=0):

    ## Number of planes    
    beta_grid, h_grid = np.meshgrid(beta,h);
    m1 = np.zeros_like(h_grid);
    m2 = np.zeros_like(h_grid);
    m = np.empty_like(m1, dtype=object)

    for j in range(len(beta)):
        for i in range(len(h)):
        
            alpha, c = antenna_swath(h[i],H, beta[j])
            if c>1500:
                m1[i,j] = np.ceil(np.pi*(EARTH_RADIUS+H)/c)
                m2[i,j] = np.ceil(np.pi*(EARTH_RADIUS+H)/2/c)
            else:
                m1[i,j] = 0
                m2[i,j] = 0
            
            m[i,j] = f'{m2[i,j]}/{m1[i,j]}'

            if m1[i,j]==0 or m2[i,j]==0:
                m[i,j] = 'N/A'

    # Create DataFrame
    line_names = [f'{h} km' for h in h]  # String array for line names'
    column_names = [f'{b*180/np.pi:.1f} deg' for b in beta]  # String array for column names
    
    table = pd.DataFrame(m, index=line_names, columns=column_names)
    
    return table




def multiCoverage(h, beta, H=0, N1 = 3, N2 = 0):

    # Initialize mesh
    beta_grid, H_grid = np.meshgrid(beta, h)
    n1 = np.zeros_like(H_grid)
    n2 = np.zeros_like(H_grid)
    n = np.empty_like(n1, dtype=object)
    
    for j in range(len(beta)):
        for i in range(len(h)):
            alpha, c = antenna_swath(h[i],H, beta[j])

            if c>1500:
                n1[i,j] = np.ceil(np.pi/alpha*N1)
                n2[i,j] = np.ceil(np.pi/alpha*N2)
                n[i,j] = f'{n1[i,j]}'
                if N2>0:
                    n[i,j] = f'{n1[i,j]}/{n2[i,j]}'

            else:
                n1[i,j] = 0
                n2[i,j] = 0
                n[i,j] = f'N/A'

    # Create DataFrame
    line_names = [f'{h} km' for h in h]  # String array for line names'
    column_names = [f'{b*180/np.pi:.1f} deg' for b in beta]  # String array for column names
    
    table = pd.DataFrame(n, index=line_names, columns=column_names)
    
    return table

def printReport(folder, file_name, creq, H, fig1, fig2, table2, table3, N1, N2):

    # Create folders
    if not os.path.exists(folder):
        os.makedirs(folder)
    if not os.path.exists(f'{folder}/figures'):
        os.makedirs(f'{folder}/figures')

    # save figures
    fig1.name = 'STEP1_altitude_x_antenna_aperture.png'
    fig1.savefig(f'{folder}/figures/{fig1.name}')

    fig2.name = 'STEP2_altitude_x_inclination.png'
    fig2.savefig(f'{folder}/figures/{fig2.name}')

    # write report
    style = f'''
            <style>
                #toc ul {{
                    list-style-type: none;
                    padding-left: 0;
                }}
                #toc li {{
                    margin-left: 0;
                }}
                #toc .h2 {{
                    margin-left: 20px;
                }}
                #toc .h3 {{
                    margin-left: 40px;
                }}
            </style>
    '''
    script = f'''
            <script>
            // JavaScript to generate the Table of Contents
            document.addEventListener("DOMContentLoaded", function() {{
                const toc = document.getElementById("toc");
                const headers = document.querySelectorAll("h1, h2, h3");
                const tocList = document.createElement("ul");

                headers.forEach(header => {{
                    const listItem = document.createElement("li");
                    listItem.className = header.tagName.toLowerCase();
                    const link = document.createElement("a");
                    link.href = `#${{header.id}}`;
                    link.textContent = header.textContent;
                    listItem.appendChild(link);
                    tocList.appendChild(listItem);
                }});

                toc.appendChild(tocList);
            }});
            </script>
    '''

    html = f'''
    <html>
        <head>
            <title>Constellation Sizing - Parametric Analysis</title>
            {style}
        </head>
        <body>
            <div id="toc"></div> <!-- Table of Contents will be inserted here -->
            <h1 id="default">Constellation Sizing - Parametric Analysis</h1>
            <p>This document presents a parametrical analysis for the preliminary design of Startical constellation.</p>
            <p>The analysis covers only the geometrical parameters of the constellation (altitude, inclination, number of planes, number of satellites per plane).</p>
            <p></p>
            <h2 id="intro">Constellation description</h2>
            <table border="1">
                <thead>
                    <tr>
                        <th>STEP</th>
                        <th>Design driver considered</th>
                        <th>Description</th>
                        <th>Outcome</th>
                    </tr>
                </thead>
                <tbody>
                    <tr>
                        <td>STEP 1</td>
                        <td>[D-2] Minimum footprint > {creq}km @ FL{H/3*100:.0f}</td>
                        <td>The footprint is mainly affected by the altitude and antenna aperture</td>
                        <td>h, beta</td>
                    </tr>
                    <tr>
                        <td>STEP 2</td>
                        <td>[D-2] Minimum footprint > {creq}km @ FL{H/3*100:.0f}<br>[D-1] Global coverage</td>
                        <td>
                            STEP 2a – Minimum inclination to reach any latitude (poles), considering also the altitude and antenna aperture (h, beta)<br>
                            STEP 2b – Considering the footprint at equator, given by (h, beta), it is possible to define the maximum separation between planes, giving m
                        </td>
                        <td>inc, m</td>
                    </tr>
                    <tr>
                        <td>STEP 3</td>
                        <td>[D-3] {N2}x coverage</td>
                        <td>Considering (h, beta) and the overlap between adjacent planes (inferred from m), the number of satellites per plane is calculated (n)</td>
                        <td>n</td>
                    </tr>
                </tbody>
            </table>
            <p></p>
            <h2 id="section1">STEP 1: ALTITUDE VS. ANTENNA APERTURE TO ENSURE MINIMUM FOOTPRINT</h2>
            <p>The footprint (c) at a given flight level (H) is primarily affected by the altitude of the constellation (h) and the antenna aperture (beta).</p>
            <img src=figures/{fig1.name} width="700">
            <h2 id="section2">STEP2: INCLINATION AND NUMBER OF PLANES TO ENSURE GLOBAL COVERAGE</h2>
            <h3 id="section2a">STEP2a: Inclination vs. altitude</h3>
            <p>The inclination required to provide global coverage of the Earth could be calculated with a trigonometric construction considering the altitude (h), the target flight level (H).</p>
            <p>The effect of antenna aperture is not considered in the following image (higher antenna aperture allows to ensure global coverage with lower inclination)</p>
            <img src=figures/{fig2.name} width="700">
            <h3 id="section2b">STEP2b: Number of planes</h3>
            <p>The footprint requirement together with the global coverage need allows to define the maximum allowable separation between adjacent planes.</p>
            <p>The number of planes needed varies depending on the coverage overlap. In the following table the extreme cases of no overlap / complete overlap are presented, sizing the min and max cases.</p>
            <p>Number of planes (min/max) for different altitudes (h) and antenna apertures (beta):</p>
            {table2.to_html()}
            <h2 id="section2">STEP3: NUMBER OF SATELLITES IN PLANE</h2>
            <p>To ensure the {N2}x coverage requirement two alternative solutions will be explored: </p>
            <ul>
                <li>Case 2: {N2}x coverage with {N1} satellites in the same plane and one in an adjacent plane</li>
                <li>Case 1: {N2}x coverage with satellites in the same plane</li>
            </ul>
            <p>The following table presents the number of satellites per plane considering the two cases described above (min / max need).</p>
            <p>Number of satellites in plane (min/max) for different altitudes (h) and antenna apertures (beta):</p>
            {table3.to_html()}
            {script}
        </body>
    </html>
    '''


    # save report
    file_path = f'{folder}/{file_name}.html'
    print(f'Saving report to {file_path}')
    with open(f'{file_path}', 'w') as f:
        f.write(html)

    webbrowser.open('file://' + os.path.realpath(file_path))







