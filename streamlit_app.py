import streamlit as st
import os

# Path to the flag file
flag_file = 'setup_completed.flag'

# Check if the setup has already been completed
if not os.path.exists(flag_file):
    print("Clearing cache data...")
    st.cache_data.clear()
    print("Installing system dependencies...")
    os.system('apt-get update && apt-get install -y libxrender1')

    # Create the flag file to indicate setup completion
    with open(flag_file, 'w') as f:
        f.write('Setup completed')

    print("Done! starting app...")

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from Model.ConstellationDesign import Constellation, FL_TO_ALT, ALT_TO_FL




def section_model_description():

    st.markdown('''
    The model assumes circular orbits and a symmetrical disposition of the orbital planes following a [Walker](https://www.researchgate.net/publication/359647076_NGSO_Constellation_Design_for_Global_Connectivity) *Start* or *Delta* pattern.
    ''')

    st.image("doc/figures/walker_constellations.png")
    st.markdown('''
    Image from [NGSO Constellation Design for Global Connectivity](https://www.researchgate.net/publication/359647076_NGSO_Constellation_Design_for_Global_Connectivity).
    ''')

    st.markdown('''
    Under these assumptions, the constellation is described by the following parameters:
    
    | Parameter | Description |
    | :-- | :-- |
    | h (km) | Altitude | 
    | inc (deg) | Inclination | 
    | m | Number of planes | 
    | n | Number of satellites per plane | 
    | walker option | *Star* or *Delta* | 

    Additionally, for the study of constellation coverage, it is assumed that the satellites are nominally pointed towards Nadir, and the field of view is characterized by the aperture angle of the instrument (half cone angle):
    
    | Parameter | Description |
    | :-- | :-- |
    | beta (deg) | Aperture angle of the instrument (half cone angle) | 

    The user can opt to calculate the coverage at sea level (ground coverage) or at a given altitude (specified in km above sea level or alterntively using the Flight Level).

    ''')

def remove_folder(folder):
    for filename in os.listdir(folder):
            file_path = os.path.join(folder, filename)
            try:
                # Check if it's a file and remove it
                if os.path.isfile(file_path) or os.path.islink(file_path):
                    os.unlink(file_path)
                # Check if it's a directory and remove its contents
                elif os.path.isdir(file_path):
                    remove_folder(file_path) 
                    os.rmdir(file_path)
            except Exception as e:
                print(f'Failed to delete {file_path}. Reason: {e}')
    

if __name__ == "__main__":

    folder = 'tmp'
    if not os.path.exists(folder):
        os.makedirs(folder)
    else:
        remove_folder(folder) # remove the folder contents (keep the folder)
    
    

    # Streamlit app
    st.title("Constellation sizing analysis tool")

    # Explanatory text
    st.markdown('''
    This app implements a geometric model for analyzing coverage provided by spacecraft constellations.
    ''')

    st.markdown("## Model description")
    
    with st.expander("Model parameters"):
        section_model_description()

    ################# User inputs
    st.sidebar.header("User Inputs")

    st.sidebar.markdown('### Constellation label')
    constellation_id = st.sidebar.text_input("Id", value="New Constellation")
    constellation_id = constellation_id.replace(" ","_")

    st.sidebar.markdown('### Constellation geometry')

    h = st.sidebar.number_input("Altitude (km)", value=650)
    inc = st.sidebar.number_input("Inclination (deg)", value=80)
    m = st.sidebar.number_input("Number of planes", value=10)
    n = st.sidebar.number_input("Number of satellites per plane", value=20)
    
    st.sidebar.markdown(f'''Total number of satellites N = {m*n}''')

    walker_option = st.sidebar.selectbox("Walker pattern:", ("star", "delta"))

    st.sidebar.markdown('### Instrument features')
    beta = st.sidebar.number_input("Aperture half-cone angle (deg)", value=65)

    #################  Plot and display results
    st.markdown("## Analysis section")

    # Create a placeholder for the info message
    info_placeholder = st.empty()
    info_placeholder.info("Check user inputs (bar on the left) and run analysis!")

    # Run analysis button
    st.sidebar.markdown('### Run analysis')

    input_type = st.sidebar.radio("Choose input type for the air/ground terminal altitudes (to evaluate coverage):", ("Flight Level (FL)", "Altitude (km)"))
    
    if input_type == "Flight Level (FL)":
        FL = st.sidebar.number_input("Flight level", value=450)
        H = FL*FL_TO_ALT
    else:
        H = st.sidebar.number_input("Altitude (km)", value=450*ALT_TO_FL)

    if st.sidebar.button("Run analysis"):

        st.markdown(f'''
        #### User inputs

        Constellation:
        | Altitude | Inclination | Number of satellites | Walker option |
        |:-- | :-- | :-- | :-- |
        |{h} km | {inc} deg | {m}x{n}={m*n} | {walker_option} |        
        
        Instrument:
        | Aperture angle |
        |:-- |
        | {beta} deg |

        Targets:
        | Altitude / Flight Level |
        |:--  |
        | {np.round(10*H)/10} km / {np.round(10*H*ALT_TO_FL)/10}|
        ''')
        
        info_placeholder.empty()
        info_placeholder.info("Creating the constellation (this may take a couple of minutes)...")

        # Perform calculations
        constellation = Constellation(constellation_id,  h, inc*np.pi/180, m, n, walker_option, beta*np.pi/180, np.round(10*H)/10)
        constellation.initializeSpacecraft()

        # Display results
        fig0_name = constellation.plot_constellation3D(True,folder)

        st.markdown('''
        #### 3D Plot 
        Geometrical aspect of the constellation
        ''')

        st.image(f'{folder}/{fig0_name}.png')
        
        # Create columns for side-by-side display
        col1, col2 = st.columns(2)

        # Display images in the columns
        with col1:
            st.image(f'{folder}/{fig0_name}_XY.png', caption='View from Z (above North Pole)', use_container_width =True)
        with col2:
            st.image(f'{folder}/{fig0_name}_YZ.png', caption='View from X (above Equator)', use_container_width =True)

        fig1_name, fig1 = constellation.plot_constellation_groundTrack(False,folder)
        
        st.markdown('''
        #### Ground track
        Ground projection of the constellation, including instrument field of view 
        ''')

        st.pyplot(fig1)

        fig2_name, table, fig2 = constellation.plot_constellation_coverage(False,folder)

        st.markdown(f'''
        #### Constellation coverage
        Coverage analysis, number of satellites visible from any location at altitude {np.round(10*H)/10} km
        ''')
        st.pyplot(fig2)

        st.markdown('''
        Percentage of Earth surface covered by *N* satellites simultaneously
        ''')
        st.dataframe(table)

        # Clear the info message
        info_placeholder.empty()