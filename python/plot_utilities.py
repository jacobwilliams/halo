#
# Plot a trajectory using PyVista
#
# https://docs.pyvista.org/version/stable/
#
# Example usage: python ./python/plot_utilities.py solution_20240101120000_L2_S_NREVS=100.txt
#

import sys
import os
import math
import csv
import json
import numpy as np
import pyvista as pv
from pyvista import examples
from pyvista.core.utilities import lines_from_points
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter

# pv.set_plot_theme('paraview')

matplotlib.rcParams["font.family"] = "Serif"
# matplotlib.rcParams["font.size"] = 20
# matplotlib.rcParams["axes.labelsize"] = 25
# matplotlib.rcParams["xtick.labelsize"] = 20
# matplotlib.rcParams["ytick.labelsize"] = 20
# matplotlib.rcParams["legend.fontsize"] = 20

J2000_DATENUM = -946728000.000000  # et of python date epoch (1970-01-01)
SEC2DAY = 3600.0*24.0
PYTHON_DIR = os.path.dirname(__file__)  # directory containing this script and other files

# things to test:
save_as_html = True
plot_starmap = False
subplot_test = True
link_views = False # for subplot_test
center_camera_on_body = False

if plot_starmap:
    font_color = 'white'
    background_color = 'black'
    cylinder_color='cyan'
    cylinder_opacity=0.1
else:
    font_color = 'black'
    background_color = 'white'
    cylinder_color='grey'
    cylinder_opacity=0.05

###################################################################
def generate_eclipse_plot(filename : str, save : bool = True,
                          show : bool = False, ylabel : str = None,
                          showtitle : bool = True, dpi : int = 200):

    """
    Generate a plot of the eclipses for a run. [from the CSV file]
    """

    et = []
    phi = []
    with open(filename, 'r') as f:
        reader = csv.reader(f)
        for i,row in enumerate(reader):
            if i==0: continue # skip header
            et.append(float(row[0]))
            phi.append(float(row[1]))

    fig = plt.figure(figsize=(20,10),facecolor="white",tight_layout=True)
    ax = fig.add_subplot(1, 1, 1)
    ax.grid()
    ax.set_axisbelow(True)
    ax.set_xlabel("Ephemeris time (sec)")
    if not ylabel:
        ylabel = 'Phi (<0 means in Earth shadow)'
    ax.set_ylabel(ylabel)
    if showtitle:
        ax.set_title(filename)
    ax.plot(et,phi,"b-",linewidth=2,label="phi")

    if save:
        plot_filename = os.path.splitext(filename)[0] + '.png'
        plt.savefig(plot_filename, dpi=dpi)

    if show:
        plt.show()

    return fig

###################################################################
def et2date(et):
    """et to python date for x axis plot"""
    return (et - J2000_DATENUM )/SEC2DAY

###################################################################
def generate_eclipse_plot_from_json(filename : str, et0 : bool, etf : bool,
                                    save : bool = True,
                                    show : bool = False, ylabel : str = None,
                                    showtitle : bool = True, dpi : int = 200):

    """
    Generate a plot of the eclipses for a run. [from the JSON file]
    """

    with open(filename) as f:
        data = json.load(f)

    et = data['et']
    phi = data['phi']

    fig = plt.figure(figsize=(20,10),facecolor="white",tight_layout=True)
    ax = fig.add_subplot(1, 1, 1)
    ax.grid()
    ax.set_axisbelow(True)
    ax.set_xlabel("Time")
    if not ylabel:
        ylabel = 'Phi (<0 means in Earth shadow)'
    ax.set_ylabel(ylabel)
    if showtitle:
        ax.set_title(filename)

    ax.plot([et2date(et0), et2date(etf)], [1.0, 1.0],"b-",linewidth=2,label="zero_line")
    for e,p in zip(et, phi):
        ax.plot([et2date(e), et2date(e)], [1.0, p+1.0],"b-",linewidth=2,label="phi")
    date_form = DateFormatter("%m/%d/%y")   # x axis is a date
    ax.xaxis.set_major_formatter(date_form)

    ax.set_ylim([0.0, 1.0])
    if save:
        plot_filename = os.path.splitext(filename)[0] + '.png'
        plt.savefig(plot_filename, dpi=dpi)

    if show:
        plt.show()

    return fig

###################################################################
def generate_defect_plot(filename : str, save : bool = True,
                         show : bool = False, ylabel : str = None,
                         showtitle : bool = True, dpi : int = 200):

    """
    Generate a plot of the constraint violations for a run.
    """

    idx = []
    r = []
    v = []
    with open(filename, 'r') as f:
        reader = csv.reader(f)
        for i,row in enumerate(reader):
            if i==0: continue # skip header
            idx.append(float(i))
            r.append(float(row[0]))
            v.append(float(row[1]))

    fig = plt.figure(figsize=(20,10),facecolor="white",tight_layout=True)
    ax = fig.add_subplot(1, 1, 1)
    ax.grid()
    ax.set_axisbelow(True)
    ax.set_xlabel("Segment interface number")
    if not ylabel:
        ylabel = 'Constraint violation magnitude'
    ax.set_ylabel(ylabel)
    if showtitle:
        ax.set_title(filename.strip('./'))
    ax.plot(idx,r,"b-",linewidth=2,label="r (km)")
    ax.plot(idx,v,"r-",linewidth=2,label="v (km/s)")
    ax.legend()

    if save:
        plot_filename = os.path.splitext(filename)[0] + '.png'
        plt.savefig(plot_filename, dpi=dpi)

    if show:
        plt.show()

    return fig

##########################################################################################################
def add_spline_trajectory_plot(p,x,y,z,
                               name: str = None,
                               color:str = None,
                               n_points:int = 500,
                               radius:float = 200.0,
                               show_scalar_bar:bool = False,
                               scalars:str=None):
    points = np.column_stack((x, y, z))
    spline = pv.Spline(points, n_points).tube(radius=radius)
    # spline.plot(scalars='arc_length', show_scalar_bar=False)
    # p.add_mesh(spline, name='trajectory', scalars='arc_length', show_scalar_bar=False)
    p.add_mesh(spline, name=name, stitle=name, color=color,
               show_scalar_bar=show_scalar_bar,
               scalars=scalars)

def add_line_trajectory_plot(p,x,y,z, scalars = None,
                             show_scalar_bar : bool = False,
                             line_width : int = 2,
                             name = None,
                             cmap = None):
    pts = np.column_stack((x, y, z)).astype(float)
    lines = lines_from_points(pts)  # lines is a pyvista.core.pointset.PolyData
    _ = p.add_mesh(lines, show_edges=True, show_scalar_bar=False,
                      scalars = scalars, flip_scalars = True,
                      line_width=line_width, name=name, cmap=cmap)
    # tube = lines.tube(50.0).elevation()
    # p.add_mesh(tube, show_edges=False, show_scalar_bar=False, scalars = rmag) # doesn't work. wrong number of points for scalars ?

    if show_scalar_bar:
        _ = p.add_scalar_bar(
            title='Radius (km)',
            color=font_color,
            interactive=False,
            vertical=True,
            title_font_size=35,
            label_font_size=30,
            font_family='times',
            outline=False,
            position_x = 0.8,
            fmt='%10.0f' )

    return lines

# this is based on a function from: https://github.com/pyvista/pyvista/blob/main/pyvista/examples/planets.py
def _sphere_with_texture_map(radius=1.0, lat_resolution=50, lon_resolution=100, lon_rotation:float = None,
                             x : float = 0.0, y : float = 0.0, z : float = 0.0):
    """Sphere with texture coordinates.

    Parameters
    ----------
    radius : float, default: 1.0
        Sphere radius.

    lat_resolution : int, default: 100
        Set the number of points in the latitude direction.

    lon_resolution : int, default: 100
        Set the number of points in the longitude direction.

    Returns
    -------
    pyvista.PolyData
        Sphere mesh with texture coordinates.
    """
    # https://github.com/pyvista/pyvista/pull/2994#issuecomment-1200520035
    theta, phi = np.mgrid[0 : np.pi : lat_resolution * 1j, 0 : 2 * np.pi : lon_resolution * 1j]
    x = x + radius * np.sin(theta) * np.cos(phi + lon_rotation)
    y = y + radius * np.sin(theta) * np.sin(phi + lon_rotation)
    z = z + radius * np.cos(theta)
    sphere = pv.StructuredGrid(x, y, z)
    texture_coords = np.empty((sphere.n_points, 2))
    texture_coords[:, 0] = phi.ravel('F') / phi.max()
    texture_coords[:, 1] = theta[::-1, :].ravel('F') / theta.max()
    sphere.active_texture_coordinates = texture_coords
    return sphere.extract_surface(pass_pointid=False, pass_cellid=False)

def plot_body(p, name:str = 'earth', texturemap:str = None,
              x:float = 0.0, y:float = 0.0, z:float = 0.0,
              r:float = 6378.0, draw_axes:float = True,
              axis_color:str = font_color,
              lon_rotation : float = 0.0):

    sphere = _sphere_with_texture_map(radius=r, lat_resolution=50,
                                      lon_resolution=100, lon_rotation = lon_rotation,
                                      x=x, y=y, z=z)
    texture = pv.Texture(texturemap)
    p.add_mesh(sphere, texture=texture, smooth_shading=False)

    if draw_axes:
        x_arrow = pv.Arrow(start=np.array([x,y,z]),
                           direction=np.array([1.0, 0.0, 0.0]), scale=r*2, tip_length=0.1,
                           tip_radius=0.05, tip_resolution=50, shaft_radius=0.02, shaft_resolution=50)
        y_arrow = pv.Arrow(start=np.array([x,y,z]),
                           direction=np.array([0.0, 1.0, 0.0]), scale=r*2, tip_length=0.1,
                           tip_radius=0.05, tip_resolution=50, shaft_radius=0.02, shaft_resolution=50)
        z_arrow = pv.Arrow(start=np.array([x,y,z]),
                           direction=np.array([0.0, 0.0, 1.0]), scale=r*2, tip_length=0.1,
                           tip_radius=0.05, tip_resolution=50, shaft_radius=0.02, shaft_resolution=50)
        p.add_mesh(x_arrow, show_edges=False, color=axis_color)
        p.add_mesh(y_arrow, show_edges=False, color=axis_color)
        p.add_mesh(z_arrow, show_edges=False, color=axis_color)

def read_trajectory_solution_file(filename : str, traj : int = 1):
    """
    Read the trajectory file (.txt or .json)
    """

    ext = os.path.splitext(filename)[1]

    if ext=='.txt':
        # inertial data only in the .txt file (this it mkspk input file)
        if not inertial:
            raise Exception(r'The .txt trajectory file only contains the inertial data. Use .json instead.')
        et = []
        x = []
        y = []
        z = []
        rmag = []
        with open(filename, 'r') as f:
            for line in f:
                data = line.strip().split(';')
                et.append(float(data[0]))
                x.append(float(data[1]))
                y.append(float(data[2]))
                z.append(float(data[3]))
                rmag.append(math.sqrt(float(data[1])**2 + float(data[2])**2 + float(data[3])**2)) # will color line based on norm(r)

    elif ext=='.json':
        # the json file has both inertial and rotating

        with open(filename) as f:
            data = json.load(f)

        if traj==1:
            et_str = 'et'
            x_str = 'x_inertial'
            y_str = 'y_inertial'
            z_str = 'z_inertial'
        elif traj==2:
            et_str = 'et'
            x_str = 'x_rotating'
            y_str = 'y_rotating'
            z_str = 'z_rotating'
        elif traj==3:
            et_str = 'et_se'
            x_str = 'x_se_rotating'
            y_str = 'y_se_rotating'
            z_str = 'z_se_rotating'

        if traj==3:
            # this is being read from the eclipse file
            et = np.array(data[et_str])
            x  = np.array(data[x_str])
            y  = np.array(data[y_str])
            z  = np.array(data[z_str])
            rmag = np.array(data['phi_se'])  # this is really phi (<0 means eclipsed)
            rmag = np.where(rmag >= 0.0, 0.0, rmag)

            # eclipsed_epochs = np.extract(rmag<0.0, rmag)
            # print(f'eclipsed_points = {eclipsed_epochs}')

        else:
            # these are in the file by segment ...
            et   = []
            x    = []
            y    = []
            z    = []
            for seg in data['segs']:
                if seg[et_str][0] > seg[et_str][-1]:
                    # reverse for backward propagated segments
                    seg[et_str].reverse()
                    seg[x_str].reverse()
                    seg[y_str].reverse()
                    seg[z_str].reverse()
                et.extend(seg[et_str])
                x.extend(seg[x_str])
                y.extend(seg[y_str])
                z.extend(seg[z_str])

            # convert to numpy arrays:
            et = np.array(et)
            x  = np.array(x)
            y  = np.array(y)
            z  = np.array(z)
            rmag = np.linalg.norm([x,y,z], axis=0)

    else:
        raise Exception(f'invalid file extension in read_trajectory_solution_file: {ext}')

    return et,x,y,z,rmag

##########################################################################################################

if __name__ == "__main__":

    # the first argument in the file to plot (can be the .txt or .json output file)
    #
    # examples:
    # file_to_plot = '../100-rev-solution/solution_20220704131413_L2_S_NREVS=100.txt'
    # file_to_plot = '../solution_20220704131413_L2_S_NREVS=10.txt'
    # file_to_plot = '../solution_20220704131413_L2_S_NREVS=5.txt'
    # file_to_plot = '../solution_20220704131415_L2_S_NREVS=100.txt'
    # file_to_plot = '../solution_20220627172948_L2_S_NREVS=100.txt'  # 100 rev with no eclipses
    # file_to_plot = '../solution_20220704131415_L2_S_NREVS=20.txt'
    # file_to_plot = '../solution_20220704131415_L2_S_NREVS=100.txt'
    # file_to_plot = '../solution_20240101120000_L2_S_NREVS=100.txt'

    file_to_plot = sys.argv[1]

    if subplot_test:
        n_plots=3
    else:
        n_plots=1

    for iplot in range(n_plots):

        if iplot==0:
            # rotating plot
            p = pv.Plotter(notebook=False, shape=(1, n_plots))
            p.subplot(0, 0)
            # p.add_text('Moon Texturemap', name='banner1', position='upper_left', color=font_color)
            # plot_body( p, name = 'moon', r = 1737.4, texturemap = os.path.join(PYTHON_DIR,'lroc_color_poles_1k.jpg') )
            et,x,y,z,rmag_rotating = read_trajectory_solution_file(file_to_plot, traj=2) # rotating frame
            add_line_trajectory_plot(p,x,y,z, scalars = rmag_rotating, show_scalar_bar = False, name='rotating')
            plot_body( p, name = 'moon', r = 1737.4, texturemap = os.path.join(PYTHON_DIR,'lroc_color_poles_1k.jpg') )
            if link_views:
                p.link_views()
        elif iplot==1:
            # inertial plot
            #p = pv.Plotter(notebook=False)
            p.subplot(0, 1)
            # read file & plot trajectory:
            et,x,y,z,rmag = read_trajectory_solution_file(file_to_plot, traj=1)
            # add_spline_trajectory_plot(p,x,y,z, name='halo', color='blue', scalars='arc_length') #, radius=100)   # doesn't show all the points ???
            add_line_trajectory_plot(p,x,y,z, scalars = rmag_rotating, show_scalar_bar = True, name='inertial')

            # https://visibleearth.nasa.gov/images/57735/the-blue-marble-land-surface-ocean-color-sea-ice-and-clouds/57737l
            # plot_body(p, name = 'earth', r = 6378.0, texturemap = 'land_ocean_ice_cloud_2048.jpg')
            # https://svs.gsfc.nasa.gov/cgi-bin/details.cgi?aid=4720
            plot_body( p, name = 'moon', r = 1737.4, texturemap = os.path.join(PYTHON_DIR,'lroc_color_poles_1k.jpg') )
        elif iplot==2:
            p.subplot(0, 2)
            file_to_plot_ = file_to_plot.replace('traj_', 'eclipse_')
            et,x,y,z,phi = read_trajectory_solution_file(file_to_plot_, traj=3)
            # this one is colored by phi ... how to get a separate colorbar?
            traj = add_line_trajectory_plot(p,x,y,z, scalars = phi, cmap='coolwarm', show_scalar_bar = False, name='rotating_se')
            # https://visibleearth.nasa.gov/images/57735/the-blue-marble-land-surface-ocean-color-sea-ice-and-clouds/57737l
            plot_body(p, name = 'earth', r = 6378.0, texturemap = os.path.join(PYTHON_DIR,'land_ocean_ice_cloud_2048.jpg') )
            # plot the sun (commented out because it makes the plot hard to deal with...)
            # plot_body(p, name = 'sun',   x=-149597871.0, y =0.0, z = 0.0, r = 696340.0,
            #           texturemap = os.path.join(PYTHON_DIR,'2k_sun.jpg'), draw_axes=False )
            # plot a basic cylinder to represent the shadow:
            # NOTE: this is just a cylinder aligned with the geometric position, it
            #       isn't a cone, isn't using the apparent position of the sun, and isn't using
            #       the eclipse bubble radius, so it's won't be quite right, but is close.  -TODO: fix this.
            cylinder = pv.Cylinder( center=[250000, 0, 0], direction=[149597871.0, 0, 0], radius=6378, height=500000.0 )
            _ = p.add_mesh(cylinder, show_edges=False, color=cylinder_color, opacity=cylinder_opacity)
            #inter = cylinder.boolean_intersection(traj)  # doesn't work ..
            # p.add_mesh(inter, show_edges=False, color='red')


        # generic star background:
        if plot_starmap:
            stars = examples.planets.download_stars_sky_background(load=False)
            p.add_background_image(stars)

        if iplot==0:
            title = os.path.basename(file_to_plot) + ' [Earth-Moon Rotating Frame]'
            axis_prefix = 'EMR'
        elif iplot==1:
            title = os.path.basename(file_to_plot) + ' [Inertial Frame]'
            axis_prefix = 'J2000'
        elif iplot==2:
            title = os.path.basename(file_to_plot) + ' [Sun-Earth Rotating Frame]'
            axis_prefix = 'SER'

        p.add_text(title, name='banner', position='upper_left', color=font_color)
        p.show_grid(color='gray', font_size=20, font_family='times',
                    xtitle = f'{axis_prefix} X (km)',
                    ytitle = f'{axis_prefix} Y (km)',
                    ztitle = f'{axis_prefix} Z (km)')
        p.set_background(color=background_color)

        p.view_isometric()

        if iplot==0:
            # hack since the original version of this script used the solution_*.txt file
            file_to_plot_ = file_to_plot.replace('traj_', 'solution_').replace('.json', '.txt')

        if iplot==0:
            ############
            # add defect plot:
            defect_filename = file_to_plot_.replace('solution', 'defects')
            defect_filename = os.path.splitext(defect_filename)[0] + '.csv'
            f = generate_defect_plot(defect_filename, save = False, ylabel='Defect', showtitle=False)
            h_chart = pv.ChartMPL(f, size=(1.0, 0.2), loc=(0.0, 0.0)) #, loc=(0.02, 0.06))
            h_chart.background_color = (1.0, 1.0, 1.0, 0.4)
            p.add_chart(h_chart)
        elif iplot==2:
            ############
            # eclipse plot
            eclipse_filename = file_to_plot_.replace('solution', 'eclipse')
            eclipse_filename = os.path.splitext(eclipse_filename)[0] + '.json'
            f = generate_eclipse_plot_from_json(eclipse_filename,
                                                save = False,
                                                ylabel='solar fraction',
                                                showtitle=False,
                                                et0 = et[0],
                                                etf = et[-1])
            h_chart = pv.ChartMPL(f, size=(1.0, 0.2), loc=(0.0, 0.0))
            h_chart.background_color = (1.0, 1.0, 1.0, 0.4)
            p.add_chart(h_chart)

        if center_camera_on_body or iplot==2:
            # recenter the camera on the Moon:
            pos = p.camera_position  # tuple: camera location, focus point, viewup vector
            p.camera_position = [pos[0], [0.0, 0.0, 0.0], pos[2]]

    if save_as_html:
        # note: background is white here, so would need to change the fonts above to black
        html_filename = os.path.splitext(file_to_plot)[0] + '.html'
        p.export_html(html_filename)

    p.show()

