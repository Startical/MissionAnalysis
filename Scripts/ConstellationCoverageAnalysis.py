from OrbitTools.ConstellationPlotUtils import *
import numpy as np
import os
import webbrowser
import pandas as pd

def evaluate_constellations(constellations, save_report = True, folder = 'tmp'):
    
    figs_folder = f'{folder}/constellations'
    
    # Create the results folder if it does not exist
    if not os.path.exists(folder):
        os.makedirs(folder)
        
    if not os.path.exists(figs_folder):
        os.makedirs(figs_folder)

    constellation_reports = ''
    constellation_definitions = ''
    all_tables = []

    for constellation in constellations:
        print(f'Running analysis: {constellation.idFullContext()}')

        if len(constellation.spacecraft)==0:
            constellation.initializeSpacecraft()

        fig0_name = plot_constellation3D(constellation,save_fig=True, save_folder=figs_folder)
        fig1_name, *_ = plot_constellation_groundTrack(constellation,save_fig=True, save_folder=figs_folder)
        fig2_name, full_table, *_ = plot_constellation_coverage(constellation,save_fig=True, save_folder=figs_folder)

        table = full_table.iloc[0:1] # keep only the first row of the table

        # Set the index of the table to the constellation identifier
        table.index = [f'{constellation.name}']

        # Collect all tables
        all_tables.append(table)

        # report
        constellation_params = f'''
        <tr>
            <td><a href="#{constellation.idFullContext()}">{constellation.name}</a></td>
            <td>{constellation.walker_option}</td>
            <td>{constellation.h}</td>
            <td>{constellation.inc*180/np.pi:.0f}</td>
            <td>{constellation.m}x{constellation.n} = {constellation.m*constellation.n}</td>
            <td>{constellation.antenna_aperture*180/np.pi:.0f}</td>
        </tr>
        '''
        constellation_definitions+= f'''{constellation_params}'''
        constellation_reports += f'''
                <h2 id={constellation.idFullContext()}>{constellation.name}</h2>
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
                <div style="display: inline-block; text-align: center;">
                    <img src=constellations/{fig0_name}.png alt="Big Figure" width="700">
                </div>
                <div style="display: inline-block; text-align: center;">
                    <div>
                        <img src=constellations/{fig0_name}_XY.png alt="Small Figure 1" width="350">
                    </div>
                    <div>
                        <img src=constellations/{fig0_name}_YZ.png alt="Small Figure 2" width="350">
                    </div>
                </div>
                <p></p>
                <img src=constellations/{fig1_name}.png width="700">
                <p>Coverage analysis at FL{constellation.H/3*100:.0f}</p>
                <img src=constellations/{fig2_name}.png width="700">
                <p>Satellites visible at FL{constellation.H/3*100:.0f}</p>
                {table.to_html()}
        '''

    # Concatenate all tables
    concatenated_table = pd.concat(all_tables)

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
            <h2>Final summary</h2>
            {concatenated_table.to_html()}
            {script}
        </body>
    </html>
    '''


    # save report
    if save_report:
        file_name = 'Constellation_sizing-Coverage_analysis'
        file_path = f'{folder}/{file_name}.html'
        print(f'Saving report to {file_path}')
        with open(f'{file_path}', 'w') as f:
            f.write(html)

        webbrowser.open('file://' + os.path.realpath(file_path))

    else:
        plt.show()

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
