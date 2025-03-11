from Model.ConstellationDesign import Constellation
import numpy as np
import os
import webbrowser


from matplotlib import pyplot as plt

def evaluate_constellations(constellations):
    
    # Create the results folder if it does not exist
    folder = 'results'
    if not os.path.exists(folder):
        os.makedirs(folder)
    figs_folder = f'{folder}/constellations'
    if not os.path.exists(figs_folder):
        os.makedirs(figs_folder)

    constellation_reports = ''
    constellation_definitions = ''

    for constellation in constellations:
        print(f'Running analysis: {constellation.idFullContext()}')
        constellation.initializeSpacecraft()

        fig1 = constellation.plot_constellation_groundTrack()
        fig2 = constellation.plot_constellation_coverage()

        fig1.savefig(f'{figs_folder}/{fig1.name}.png')
        fig2.savefig(f'{figs_folder}/{fig2.name}.png')

        # report
        constellation_params = f'''
        <tr>
            <td><a href="#{constellation.idFullContext()}">{constellation.id}</a></td>
            <td>{constellation.walker_option}</td>
            <td>{constellation.h}</td>
            <td>{constellation.inc*180/np.pi:.0f}</td>
            <td>{constellation.m}x{constellation.n} = {constellation.m*constellation.n}</td>
            <td>{constellation.antenna_aperture*180/np.pi:.0f}</td>
        </tr>
        '''
        constellation_definitions+= f'''{constellation_params}'''
        constellation_reports += f'''
                <h2 id={constellation.idFullContext()}>{constellation.id}</h2>
                <table border="1">
                <thead>
                    <tr>
                        <th>ID</th>
                        <th>Option</th>
                        <th>h (km)</th>
                        <th>inc (deg)</th>
                        <th>Number of satellites</th>
                        <th>Antenna aperture (deg)</th>
                    </tr>
                </thead>
                <tbody>
                    {constellation_params}
                </tbody>
                </table>
                <p>Geometrical aspect of the constellation</p>
                <img src=constellations/{fig1.name}.png width="700">
                <p>Coverage analysis at FL{constellation.H/3*100:.0f}</p>
                <img src=constellations/{fig2.name}.png width="700">
        '''
    

    # create the report
    style, script = reportStyle()

    html = f'''
    <html>
        <head>
            <title>Constellation Sizing - Coverage Analysis</title>
            {style}
        </head>
        <body>
            <div id="toc"></div> <!-- Table of Contents will be inserted here -->
            <h1 id="default">Constellation Sizing - Coverage Analysis</h1>
            <p>This document presents the coverage analysis of constellation candidates for Startical Space Based Communication, Navigation and Surveillance (SB-CNS) System.</p>
            <p>The analysis covers only the geometrical parameters of the constellation (altitude, inclination, number of planes, number of satellites per plane).</p>
            <p></p>
            <h2 id="intro">Constellation options</h2>
            <table border="1">
            <thead>
                <tr>
                    <th>ID</th>
                    <th>Option</th>
                    <th>h (km)</th>
                    <th>inc (deg)</th>
                    <th>Number of satellites</th>
                    <th>Antenna aperture (deg)</th>
                </tr>
            </thead>
            <tbody>
                {constellation_definitions}
            </tbody>
            </table>
            {constellation_reports}
            {script}
        </body>
    </html>
    '''


    # save report
    file_name = 'Constellation_sizing-Coverage_analysis'
    file_path = f'{folder}/{file_name}.html'
    print(f'Saving report to {file_path}')
    with open(f'{file_path}', 'w') as f:
        f.write(html)

    webbrowser.open('file://' + os.path.realpath(file_path))

def reportStyle():
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
    return style, script
