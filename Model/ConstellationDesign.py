from Model.Spacecraft import Spacecraft
from OrbitTools.Constants import EARTH_RADIUS
from Scripts.ParametricAnalysis import antenna_swath

import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')  # or 'Qt5Agg' if you prefer

import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.patches import Circle, Ellipse
from matplotlib.colors import ListedColormap, BoundaryNorm, LinearSegmentedColormap
from mpl_toolkits.mplot3d import Axes3D

import cartopy.crs as ccrs
import cartopy.feature as cfeature

import pandas as pd

from OrbitTools.FrameTransformations import FrameTransformations as frames

import vtk

# FLIGHT LEVEL TO ALTITUDE
ALT_TO_FL = 0.03284*1000  # ALT in km
FL_TO_ALT = 1/ALT_TO_FL   # ALT in km


class Constellation(object):
    """Constellation class"""

    id = ""
    spacecraft = []

    H = 0; # flight level
    h = 650;
    inc = 90*np.pi/180;
    walker_option = "star";
    n = 16; # number of satellites per plane
    m = 9; # number of orbital planes
    N = n*m;

    antenna_aperture = 60*np.pi/180;

    def __init__(self, id):
        self.id = id;

    def __init__(self, id, h, inc, m, n, walker_option, antenna_aperture, H):
        self.id = id
        self.h = h
        self.inc = inc
        self.m = m
        self.n = n
        self.walker_option = walker_option
        self.antenna_aperture = antenna_aperture
        self.H = H
        self.N = m*n

    def idFullContext(self):
        return f'{self.id}_{self.h:.0f}_{self.inc*180/np.pi:.0f}_{self.m}_{self.n}_{self.walker_option}_{self.antenna_aperture*180/np.pi:.0f}_FL{self.H/3*100:.0f}'

    def antenna_swath(self):

        hh = (EARTH_RADIUS+self.h)/(EARTH_RADIUS+self.H);
        gamma = np.arccos(1/hh);

        while np.sin(gamma)-(hh-np.cos(gamma))*np.tan(self.antenna_aperture) > 1e-6:
            gamma = np.arcsin((hh-np.cos(gamma))*np.tan(self.antenna_aperture));

        return gamma


    def planesRAANdistribution(self):
        raanDistribution = 2*np.pi
        if self.walker_option == "star":
            raanDistribution = np.pi

        return raanDistribution

    def initializeSpacecraft(self):
        
        # Initialize spacecraft parameters
        defaultKeplerParams = np.array([EARTH_RADIUS+self.h, 0, self.inc, 0, 0, 0])
        self.spacecraft = []

        for i in range(self.m):  # loop over planes
            for j in range(self.n):  # loop over satellites in each plane
                sc_id = f"SC_{i,j}"
                kepler_parameters = defaultKeplerParams

                kepler_parameters[3] = self.planesRAANdistribution()/self.m*i
                kepler_parameters[5] = 2*np.pi/self.n*(j+i/2) % (2 * np.pi)

                newSC = Spacecraft(sc_id, kepler_parameters)
                newSC.propagate(newSC.get_orbitalPeriod(),newSC.get_orbitalPeriod()/20)
                self.spacecraft.append(newSC)




    def plot_constellation3D(self, save_fig = False, save_folder = '', enable_interaction = False):
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
        colormap = cm.get_cmap('viridis', self.m)

        # Plot satellites
        x_plane = []
        y_plane = []
        z_plane = []

        # Create satellite line and points
        satellite_lines = []
        satellite_points = []

        for sc in self.spacecraft:
            xyz = sc.xyz.get_data_at_time(0)
            dt, xyz1 = sc.xyz.get_data_at_index(1)

            dxyz = (xyz1-xyz)/np.linalg.norm((xyz1-xyz))

            i_index = int(sc.id.split('_(')[1].split(',')[0])  # Extract the i index from the spacecraft ID
            j_index = int(sc.id.split('_(')[1].split(',')[1].split(')')[0])  # Extract the j index from the spacecraft ID

            color = colormap(i_index / self.m)  # Get the color from the colormap
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
            rotation_axis = rotation_axis / np.linalg.norm(rotation_axis)
            rotation_angle = np.degrees(np.arccos(np.dot(default_arrow, velocity_vector)))
            
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
            alpha, c = antenna_swath(self.h, self.H, self.antenna_aperture)
            R = (r-(EARTH_RADIUS+self.H)*np.cos(alpha))*np.tan(self.antenna_aperture)
            XYZ = frames.spherical2cartesian(EARTH_RADIUS+self.H, long, lat)

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
        title.SetInput(f'{self.id}\n{self.h} km, {self.inc*180/np.pi:.0f} deg\n{self.m} x {self.n} , walker {self.walker_option}')
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
        
        fig_name = f'{self.idFullContext()}_3D'

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


    def plot_constellation_groundTrack(self, save_fig = False, save_folder = '', enable_interaction = False):
        """Plot the Earth in longitude/latitude coordinates"""
        
        fig, ax = plt.subplots(figsize=(10, 5), subplot_kw={'projection': ccrs.PlateCarree()})
    
        # Add coastlines and borders
        ax.coastlines()
        ax.add_feature(cfeature.BORDERS, linestyle=':')
    
        # Add gridlines
        gl = ax.gridlines(draw_labels=True)
        gl.top_labels = False
        gl.right_labels = False
    
        # Define a colormap
        colormap = cm.get_cmap('viridis', self.m)
    
        # Plot satellites
       
        for sc in self.spacecraft:
            xyz = sc.xyz.get_data_at_time(0)
            dt, xyz1 = sc.xyz.get_data_at_index(1)
    
            R,long,lat = frames.cartesian2spherical(xyz)
            R1, long1, lat1 = frames.cartesian2spherical(xyz1)
    
            i_index = int(sc.id.split('_(')[1].split(',')[0])  # Extract the i index from the spacecraft ID
            j_index = int(sc.id.split('_(')[1].split(',')[1].split(')')[0])  # Extract the j index from the spacecraft ID
            
            color = colormap(i_index / self.m)  # Get the color from the colormap
            
            # Plot each satellite
            ax.scatter(long*180/np.pi, lat*180/np.pi, color=color, s=10, transform=ccrs.PlateCarree())
            ax.text(long*180/np.pi, lat*180/np.pi, sc.id, fontsize=4, verticalalignment='bottom')
            arrow_length = 5
            dlong = long1-long
            dlat = lat1-lat
            ax.arrow(long * 180 / np.pi, lat * 180 / np.pi, dlong/np.linalg.norm([dlong,dlat]) * arrow_length, dlat/np.linalg.norm([dlong,dlat])  * arrow_length, color=color, head_width=1, head_length=2, transform=ccrs.PlateCarree())
    
            # Add a disk (ellipse to consider the 2D projection) centered on each point
            alpha, c = antenna_swath(self.h, self.H, self.antenna_aperture)
            r2 = 180/np.pi*alpha
            
            if np.cos(lat)>1e-2:
                r1 = r2/np.cos(lat);
            else:
                r1 = 180;

            ellipse = Ellipse((long*180/np.pi, lat*180/np.pi), width=2*r1, height=2*r2, color=color, alpha=0.2, transform=ccrs.PlateCarree())
            ax.add_patch(ellipse)

    
        plt.title(f'Constellation coverage: {self.h} km, {self.inc*180/np.pi:.0f} deg, {self.m} x {self.n} , walker {self.walker_option}')
        
        fig.name = f'{self.idFullContext()}_fig1'

        if save_fig:
            fig.savefig(f'{save_folder}/{fig.name}.png')

        if enable_interaction:
            plt.show()

        return fig.name, fig
    
    def plot_constellation_coverage(self, save_fig = False, save_folder = '', enable_interaction = False):
        H = self.H
        TOTAL_AREA = 4*np.pi*((EARTH_RADIUS+H)**2)
    
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
    
        for i in range(len(long)): 
            for j in range(len(lat)):

                # area element
                da = ((EARTH_RADIUS+H)**2)*np.cos(lat[j])*deltaLong*deltaLat
    
                # point position
                xyz_P = frames.spherical2cartesian(EARTH_RADIUS + H, long[i], lat[j]);
                R = np.linalg.norm(xyz_P)
                u_P = xyz_P/R
                gamma_lim = np.arcsin(EARTH_RADIUS/R)
    
                for sc in self.spacecraft:
    
                    xyz_SAT = sc.xyz.get_data_at_time(0)
                    u_SAT = xyz_SAT/np.linalg.norm(xyz_SAT)
                    xyz_P_SAT = xyz_SAT - xyz_P
                    u_P_SAT = xyz_P_SAT/np.linalg.norm(xyz_P_SAT)
    
                    # check if satellite is above the horizon
                    if (angle_u_v(-u_P, u_P_SAT) > gamma_lim):
                        # satellite is above the horizon, check if visible within antenna aperture
                        if (angle_u_v(-u_SAT, -u_P_SAT)< self.antenna_aperture):
                            n[j,i] = n[j,i]+1
        
                idx = int(min(n[j,i],len(area)-1))
                area[idx] = area[idx]+da/TOTAL_AREA

        fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
    
        levels1 = np.append(levels-0.5, np.max(n))
        
        # Create a custom colormap
        coolwarm = plt.get_cmap('coolwarm')
        colors = [(1, 1, 1)] + [coolwarm(i) for i in range(coolwarm.N)]
        cmap = ListedColormap(colors)
        norm = BoundaryNorm(levels1, ncolors=cmap.N, clip=True)
    
        # Plot coverage as contour
        contour = ax.contourf(long_grid*180/np.pi, lat_grid*180/np.pi, n, levels=levels1, cmap=cmap, norm=norm)
    
        # Add coastlines and borders
        ax.coastlines()
        ax.add_feature(cfeature.BORDERS, linestyle=':')
    
        # Add gridlines
        gl = ax.gridlines(draw_labels=True)
        gl.top_labels = False
        gl.right_labels = False
    
        # Add colorbar with title
        cbar = plt.colorbar(contour, ax=ax, orientation='horizontal', pad=0.05)
        
        labels = [f'{levels[i]}' for i in levels]
        labels[-1] = labels[-1]+'+'
        cbar.set_ticks(levels)
        cbar.set_ticklabels(labels)
        cbar.set_label('Number of visible satellites')
    
        plt.title(f'Constellation coverage: {self.h} km, {self.inc*180/np.pi:.0f} deg, {self.m} x {self.n} , walker {self.walker_option}')
        fig.name = f'{self.idFullContext()}_fig2'

        if save_fig:
            fig.savefig(f'{save_folder}/{fig.name}.png')

        if enable_interaction:
            plt.show()

        ## Create results table

        index_names = ['% area covered']
        column_names = [f'{idx}' for idx in levels]  # String array for column names
        column_names[-1] = '>'+column_names[-1]

        table = pd.DataFrame([[f'{a*100:.1f} %' for a in area]], index=index_names, columns=column_names)

        return fig.name, table, fig

def angle_u_v(u,v):
    u_unit = u/np.linalg.norm(u)
    v_unit = v/np.linalg.norm(v)

    f = np.dot(u_unit,v_unit)

    if f<-1.0:
        f=-1.0
    if f>1.0:
        f=1.0

    return np.arccos(f)

    
        

