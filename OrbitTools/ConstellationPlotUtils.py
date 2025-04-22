
from Model.Spacecraft import Spacecraft
from OrbitTools.Constants import EARTH_RADIUS
from OrbitTools.PlotUtils import initialize_Earth_2D_plot, plot_circle_projection
from OrbitTools.CoverageAnalysis import antenna_swath, check_visibility_lower_upper_terminal

import numpy as np

import matplotlib
#matplotlib.use('Qt5Agg')  # or 'Qt5Agg' if you prefer

import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import ListedColormap, BoundaryNorm, LinearSegmentedColormap

import cartopy.crs as ccrs
import cartopy.feature as cfeature

import pandas as pd

from OrbitTools.FrameTransformations import FrameTransformations as frames

import vtk

# FLIGHT LEVEL TO ALTITUDE
ALT_TO_FL = 0.03284*1000  # ALT in km
FL_TO_ALT = 1/ALT_TO_FL   # ALT in km


def plot_constellation3D(constellation, save_fig = False, save_folder = '', enable_interaction = False):
    """Plot the Constellation 3D geometry"""

    # Create a sphere source to represent the Earth
    earthSource = vtk.vtkSphereSource()
    earthSource.SetRadius(EARTH_RADIUS)
    earthSource.SetThetaResolution(100)
    earthSource.SetPhiResolution(100)

    # Create a mapper and actor for the Earth
    earthMapper = vtk.vtkPolyDataMapper()
    earthMapper.SetInputConnection(earthSource.GetOutputPort())
    earthActor = vtk.vtkActor()
    earthActor.GetProperty().SetColor(0.5, 0.5, 0.5)  # Grey color
    earthActor.SetMapper(earthMapper)

    # Create a renderer
    renderer = vtk.vtkRenderer()
        
    # Add the Earth actor to the scene
    renderer.AddActor(earthActor)

    # Create axes
    axes = vtk.vtkAxesActor()
    axes.SetTotalLength(2 * EARTH_RADIUS, 2 * EARTH_RADIUS, 2 * EARTH_RADIUS)  # Set the length of the axes
    axes.SetShaftTypeToCylinder()  # Optional: Set the shaft type
    axes.SetCylinderRadius(0.02)  # Optional: Set the radius of the cylinder
    axes.SetConeRadius(0.1)  # Optional: Set the radius of the cone
    axes.SetSphereRadius(0.1)  # Optional: Set the radius of the sphere

    # Set the color of the axes to black
    axes.GetXAxisShaftProperty().SetColor(0, 0, 0)
    axes.GetYAxisShaftProperty().SetColor(0, 0, 0)
    axes.GetZAxisShaftProperty().SetColor(0, 0, 0)
    axes.GetXAxisTipProperty().SetColor(0, 0, 0)
    axes.GetYAxisTipProperty().SetColor(0, 0, 0)
    axes.GetZAxisTipProperty().SetColor(0, 0, 0)
        
    # Set the label text properties to smaller size
    axes.GetXAxisCaptionActor2D().GetCaptionTextProperty().SetFontSize(1)
    axes.GetYAxisCaptionActor2D().GetCaptionTextProperty().SetFontSize(1)
    axes.GetZAxisCaptionActor2D().GetCaptionTextProperty().SetFontSize(1)

    # Add the axes to the renderer
    renderer.AddActor(axes)

    # Define a colormap
    colormap = cm.get_cmap('viridis', constellation.m)

    # Plot satellites
    x_plane = []
    y_plane = []
    z_plane = []

    # Create satellite line and points
    satellite_lines = []
    satellite_points = []

    for sc in constellation.spacecraft:
        xyz = sc.xyz.get_data_at_time(0)
        dt, xyz1 = sc.xyz.get_data_at_index(1)

        dxyz = (xyz1-xyz)/np.linalg.norm((xyz1-xyz))

        i_index = int(sc.id.split('_(')[1].split(',')[0])  # Extract the i index from the spacecraft ID
        j_index = int(sc.id.split('_(')[1].split(',')[1].split(')')[0])  # Extract the j index from the spacecraft ID

        color = colormap(i_index / constellation.m)  # Get the color from the colormap
        color_r, color_g, color_b, _ = color


        # Create a sphere source to represent the satellite
        satelliteSource = vtk.vtkSphereSource()
        satelliteSource.SetRadius(100)
        satelliteSource.SetCenter(xyz[0], xyz[1], xyz[2])

        # Create a mapper and actor for the satellite
        satelliteMapper = vtk.vtkPolyDataMapper()
        satelliteMapper.SetInputConnection(satelliteSource.GetOutputPort())
        satelliteActor = vtk.vtkActor()
        satelliteActor.GetProperty().SetColor(color_r, color_g, color_b)
        satelliteActor.SetMapper(satelliteMapper)

        # Add the satellite actor to the scene
        renderer.AddActor(satelliteActor)

        # Create an arrow source to represent the direction
        arrowSource = vtk.vtkArrowSource()
        arrowSource.SetTipLength(0.1)
        arrowSource.SetTipRadius(0.05)
        arrowSource.SetShaftRadius(0.02)

        # Create a transform to position the arrow
        transform = vtk.vtkTransform()
        transform.Translate(xyz[0], xyz[1], xyz[2])
            
        velocity_vector = dxyz / np.linalg.norm(dxyz)
        default_arrow = np.array([1, 0, 0])
        rotation_axis = np.cross(default_arrow, velocity_vector)

        if np.linalg.norm(rotation_axis) > 1e-6:
            rotation_axis = rotation_axis / np.linalg.norm(rotation_axis)
            rotation_angle = np.degrees(np.arccos(np.dot(default_arrow, velocity_vector)))
        else:
            rotation_axis = np.array([0, 0, 1])
            rotation_angle = 0
            
        transform.RotateWXYZ(rotation_angle, rotation_axis[0], rotation_axis[1], rotation_axis[2])
        transform.Scale(1000, 1000, 1000)  # Set the arrow length

        # Create a transform filter
        transformFilter = vtk.vtkTransformPolyDataFilter()
        transformFilter.SetTransform(transform)
        transformFilter.SetInputConnection(arrowSource.GetOutputPort())
        transformFilter.Update()

        # Create a mapper and actor for the arrow
        arrowMapper = vtk.vtkPolyDataMapper()
        arrowMapper.SetInputConnection(transformFilter.GetOutputPort())
        arrowActor = vtk.vtkActor()
        arrowActor.GetProperty().SetColor(color_r, color_g, color_b)
        arrowActor.SetMapper(arrowMapper)

        # Add the arrow actor to the scene
        renderer.AddActor(arrowActor)

        # Create a circle to represent the antenna coverage
        r,long,lat = frames.cartesian2spherical(xyz)
        alpha, c = antenna_swath(constellation.h, constellation.H, constellation.antenna_aperture)
        R = (r-(EARTH_RADIUS+constellation.H)*np.cos(alpha))*np.tan(constellation.antenna_aperture)
        XYZ = frames.spherical2cartesian(EARTH_RADIUS+constellation.H, long, lat)

        circleSource = vtk.vtkRegularPolygonSource()
        circleSource.SetNumberOfSides(100)  # Higher number for smoother circle
        circleSource.SetRadius(R)
        circleSource.SetCenter(0, 0, 0)

        # Create a transform to position the circle
        transform2 = vtk.vtkTransform()
        transform2.Translate(XYZ[0], XYZ[1], XYZ[2])
            
        normal_vector = xyz / np.linalg.norm(xyz)
        default_normal = np.array([0, 0, 1])
        rotation_axis = np.cross(default_normal, normal_vector)
        rotation_axis = rotation_axis / np.linalg.norm(rotation_axis)
        rotation_angle = np.degrees(np.arccos(np.dot(default_normal, normal_vector)))

        transform2.RotateWXYZ(rotation_angle, rotation_axis[0], rotation_axis[1], rotation_axis[2])

        # Create a transform filter
        transformFilter2 = vtk.vtkTransformPolyDataFilter()
        transformFilter2.SetTransform(transform2)
        transformFilter2.SetInputConnection(circleSource.GetOutputPort())
        transformFilter2.Update()

        # Create a mapper
        circleMapper = vtk.vtkPolyDataMapper()
        circleMapper.SetInputConnection(transformFilter2.GetOutputPort())
            
        # Create an actor
        circleActor = vtk.vtkActor()
        circleActor.SetMapper(circleMapper)
            
        # Set color and transparency
        circleActor.GetProperty().SetColor(color_r, color_g, color_b)
        circleActor.GetProperty().SetOpacity(0.2)  # 50% transparency

        # Add the arrow actor to the scene
        renderer.AddActor(circleActor)


    # Add a title using vtkTextActor
    title = vtk.vtkTextActor()
    title.SetInput(f'{constellation.name}\n{constellation.h} km, {constellation.inc*180/np.pi:.0f} deg\n{constellation.m} x {constellation.n} , walker {constellation.walker_option}')
    title.GetTextProperty().SetFontSize(24)
    title.GetTextProperty().SetColor(0.0, 0.0, 0.0)  # Black color
    title.SetPosition(10, 450)  # Position the title (x, y)

    renderer.AddActor(title)

    # Set the background color and camera position
    renderer.SetBackground(1, 1, 1)  # White background
    camera = renderer.GetActiveCamera()
    camera.SetPosition(4 * EARTH_RADIUS, 4 * EARTH_RADIUS, 4 * EARTH_RADIUS)
    camera.SetFocalPoint(0, 0, 0)
    camera.SetViewUp(0,0,1)  # Optional: Set the view up vector
    camera.Zoom(1)  # Optional: Zoom in or out

    # render window, and interactor
    renderWindow = vtk.vtkRenderWindow()
    renderWindow.OffScreenRenderingOn()
    renderWindow.AddRenderer(renderer)
        
    # Set the size of the render window
    renderWindow.SetSize(800, 800)  # Width, Height

    renderWindowInteractor = vtk.vtkRenderWindowInteractor()
    renderWindowInteractor.SetRenderWindow(renderWindow)
    renderWindowInteractor.Initialize()

    # Render the scene
    renderWindow.Render()
    renderer.ResetCamera()
    renderWindow.Render()
        
    fig_name = f'{constellation.idFullContext()}_3D'

    if save_fig:

        # Export the scene as an image
        windowToImageFilter = vtk.vtkWindowToImageFilter()
        windowToImageFilter.SetInput(renderWindow)
        windowToImageFilter.Update()

        writer = vtk.vtkPNGWriter()
        writer.SetFileName(f'{save_folder}/{fig_name}.png')
        writer.SetInputConnection(windowToImageFilter.GetOutputPort())
        writer.Write()

        ## change camera position - View from Z
        camera.SetPosition(0,0, 4 * EARTH_RADIUS)
        camera.SetViewUp(0,1,0)  # Optional: Set the view up vector
            
        renderWindow.Render()
        renderer.ResetCamera()
        renderWindow.Render()
        windowToImageFilter2 = vtk.vtkWindowToImageFilter()
        windowToImageFilter2.SetInput(renderWindow)
        windowToImageFilter2.Update()
        writer.SetFileName(f'{save_folder}/{fig_name}_XY.png')
        writer.SetInputConnection(windowToImageFilter2.GetOutputPort())
        writer.Write()

        ## change camera position - View from X
        camera.SetPosition(4 * EARTH_RADIUS, 0, 0)
        camera.SetViewUp(0,0,1)  # Optional: Set the view up vector
            
        renderWindow.Render()
        renderer.ResetCamera()
        renderWindow.Render()
        windowToImageFilter3 = vtk.vtkWindowToImageFilter()
        windowToImageFilter3.SetInput(renderWindow)
        windowToImageFilter3.Update()
        writer.SetFileName(f'{save_folder}/{fig_name}_YZ.png')
        writer.SetInputConnection(windowToImageFilter3.GetOutputPort())
        writer.Write()


    if enable_interaction:
        # this will block the execution of the script until the user closes the window
        renderWindowInteractor.Start()

    return fig_name


def plot_constellation_groundTrack(constellation, save_fig = False, save_folder = '', enable_interaction = False):
    """Plot the Earth in longitude/latitude coordinates"""
        
    fig, ax = initialize_Earth_2D_plot()
    
    # Define a colormap
    colormap = cm.get_cmap('viridis', constellation.m)
    
    # Plot satellites
       
    for sc in constellation.spacecraft:
        xyz = sc.xyz.get_data_at_time(0)
        dt, xyz1 = sc.xyz.get_data_at_index(1)
    
        R,long,lat = frames.cartesian2spherical(xyz)
        R1, long1, lat1 = frames.cartesian2spherical(xyz1)
    
        i_index = int(sc.id.split('_(')[1].split(',')[0])  # Extract the i index from the spacecraft ID
        j_index = int(sc.id.split('_(')[1].split(',')[1].split(')')[0])  # Extract the j index from the spacecraft ID
            
        color = colormap(i_index / constellation.m)  # Get the color from the colormap
            
        # Plot each satellite
        ax.scatter(long*180/np.pi, lat*180/np.pi, color=color, s=10, transform=ccrs.PlateCarree())
        ax.text(long*180/np.pi, lat*180/np.pi, sc.id, fontsize=4, verticalalignment='bottom')
        arrow_length = 5
        dlong = long1-long
        dlat = lat1-lat
        ax.arrow(long * 180 / np.pi, lat * 180 / np.pi, dlong/np.linalg.norm([dlong,dlat]) * arrow_length, dlat/np.linalg.norm([dlong,dlat])  * arrow_length, color=color, head_width=1, head_length=2, transform=ccrs.PlateCarree())
    
        # Add a disk centered on each point
        alpha, c = antenna_swath(constellation.h, constellation.H, constellation.antenna_aperture)
            
        plot_circle_projection(ax, long*180/np.pi, lat*180/np.pi, 180/np.pi*alpha, color=color, alpha = 0.2)

    
    plt.title(f'Constellation coverage: {constellation.h} km, {constellation.inc*180/np.pi:.0f} deg, {constellation.m} x {constellation.n} , walker {constellation.walker_option}')
        
    fig.name = f'{constellation.idFullContext()}_fig1'

    if save_fig:
        fig.savefig(f'{save_folder}/{fig.name}.png')

    if enable_interaction:
        plt.show()

    return fig.name, fig
    
def plot_constellation_coverage(constellation, save_fig = False, save_folder = '', enable_interaction = False):
    H = constellation.H
    TOTAL_AREA = 4*np.pi*(EARTH_RADIUS+H)**2
    AREA_BELOW_75 = TOTAL_AREA - 2*2*np.pi*(EARTH_RADIUS+H)**2*(1-np.cos(np.pi/2-75*np.pi/180))
    AREA_BELOW_60 = TOTAL_AREA - 2*2*np.pi*(EARTH_RADIUS+H)**2*(1-np.cos(np.pi/2-60*np.pi/180))
    
    deltaLong = np.pi/50
    long = np.arange(0,2*np.pi+deltaLong, deltaLong)
    deltaLat = np.pi/50
    lat = np.arange(-np.pi/2, np.pi/2+deltaLat,deltaLat)
    
    long_grid, lat_grid = np.meshgrid(long, lat)
    n = np.zeros_like(long_grid)

    # Define contour levels
    levels = np.arange(0,6,1)

    # area placeholder
    area = np.zeros(len(levels))

    area_75 = np.zeros(len(levels))
    area_60 = np.zeros(len(levels))
    
    for i in range(len(long)): 
        for j in range(len(lat)):

            # area element
            da = ((EARTH_RADIUS+H)**2)*np.cos(lat[j])*deltaLong*deltaLat
    
            # point position
            target_pos = [H, long[i], lat[j]]    

            for sc in constellation.spacecraft:

                xyz = sc.xyz.get_data_at_time(0)
                r, long_sc, lat_sc = frames.cartesian2spherical(xyz)
                spacecraft_pos = [r-EARTH_RADIUS, long_sc, lat_sc]

                is_visible = check_visibility_lower_upper_terminal(target_pos, spacecraft_pos, lower_terminal_fov=np.pi, upper_terminal_fov=constellation.antenna_aperture)
    
                if is_visible:
                    n[j,i] = n[j,i]+1
        
            idx = int(min(n[j,i],len(area)-1))
            area[idx] = area[idx]+da/TOTAL_AREA

            if abs(lat[j])<=75*np.pi/180:
                area_75[idx] = area_75[idx]+da/AREA_BELOW_75
            if abs(lat[j])<=60*np.pi/180:
                area_60[idx] = area_60[idx]+da/AREA_BELOW_60

    # 
    fig, ax = initialize_Earth_2D_plot()
    
    levels1 = np.append(levels-0.5, np.max(n))
        
    # Create a custom colormap
    coolwarm = plt.get_cmap('coolwarm')
    colors = [(1, 1, 1)] + [coolwarm(i) for i in range(coolwarm.N)]
    cmap = ListedColormap(colors)
    norm = BoundaryNorm(levels1, ncolors=cmap.N, clip=True)
    
    # Plot coverage as contour
    contour = ax.contourf(long_grid*180/np.pi, lat_grid*180/np.pi, n, levels=levels1, cmap=cmap, norm=norm)
    
    # Add colorbar with title
    cbar = plt.colorbar(contour, ax=ax, orientation='horizontal', pad=0.05)
        
    labels = [f'{levels[i]}' for i in levels]
    labels[-1] = labels[-1]+'+'
    cbar.set_ticks(levels)
    cbar.set_ticklabels(labels)
    cbar.set_label('Number of visible satellites')
    
    plt.title(f'Constellation coverage: {constellation.h} km, {constellation.inc*180/np.pi:.0f} deg, {constellation.m} x {constellation.n} , walker {constellation.walker_option}')
    fig.name = f'{constellation.idFullContext()}_fig2'

    if save_fig:
        fig.savefig(f'{save_folder}/{fig.name}.png')

    if enable_interaction:
        plt.show()

    ## Create results table

    index_names = ['% total area covered', '% lat<75deg', '% lat<60deg']
    column_names = [f'{idx}' for idx in levels]  # String array for column names
    column_names[-1] = '>'+column_names[-1]

    table = pd.DataFrame([[f'{a*100:.1f} %' for a in area],[f'{a*100:.1f} %' for a in area_75],[f'{a*100:.1f} %' for a in area_60]], index=index_names, columns=column_names)

    return fig.name, table, fig

     