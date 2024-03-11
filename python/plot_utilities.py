#
# Plot a trajectory using PyVista
#
# https://docs.pyvista.org/version/stable/
#

import os
import math
import numpy as np
import pyvista as pv
from pyvista import examples
from pyvista.core.utilities import lines_from_points
import matplotlib.pyplot as plt
import utilities

# things to test:
save_as_html = True
plot_starmap = False
subplot_test = False
# file_to_plot = '../100-rev-solution/solution_20220704131413_L2_S_NREVS=100.txt'  # file to plot
# file_to_plot = '../solution_20220704131413_L2_S_NREVS=10.txt'
# file_to_plot = '../solution_20220704131413_L2_S_NREVS=5.txt'
# file_to_plot = '../solution_20220704131415_L2_S_NREVS=100.txt'
file_to_plot = '../solution_20220627172948_L2_S_NREVS=100.txt'  # 100 rev with no eclipses
# file_to_plot = '../solution_20220704131415_L2_S_NREVS=20.txt'
# file_to_plot = '../solution_20220704131415_L2_S_NREVS=100.txt'

#TODO: make this an input to the script...

if plot_starmap:
    font_color = 'white'
    background_color = 'black'
else:
    font_color = 'black'
    background_color = 'white'


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
    p.add_mesh(spline, name=name, color=color, show_scalar_bar=show_scalar_bar, scalars=scalars)

def add_line_trajectory_plot(p,x,y,z, scalars = None, show_scalar_bar : bool = False, line_width : int = 2):
    pts = np.column_stack((x, y, z)).astype(float)
    lines = lines_from_points(pts)
    p.add_mesh(lines, show_edges=True, show_scalar_bar=False, scalars = scalars, flip_scalars = True, line_width=line_width)
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
            fmt='%10.3f' )

# this is based on a function from: https://github.com/pyvista/pyvista/blob/main/pyvista/examples/planets.py
def _sphere_with_texture_map(radius=1.0, lat_resolution=50, lon_resolution=100, lon_rotation:float = None):
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
    x = radius * np.sin(theta) * np.cos(phi + lon_rotation)
    y = radius * np.sin(theta) * np.sin(phi + lon_rotation)
    z = radius * np.cos(theta)
    sphere = pv.StructuredGrid(x, y, z)
    texture_coords = np.empty((sphere.n_points, 2))
    texture_coords[:, 0] = phi.ravel('F') / phi.max()
    texture_coords[:, 1] = theta[::-1, :].ravel('F') / theta.max()
    sphere.active_texture_coordinates = texture_coords
    return sphere.extract_surface(pass_pointid=False, pass_cellid=False)

def plot_body(p, name:str = 'earth', texturemap:str = None,
              x:float = 0.0, y:float = 0.0, z:float = 0.0,
              r:float = 6378.0, draw_axes:float = True, axis_color:str = font_color, lon_rotation : float = 0.0):

    sphere = _sphere_with_texture_map(radius=r, lat_resolution=50, lon_resolution=100, lon_rotation = lon_rotation)
    texture = pv.Texture(texturemap)
    p.add_mesh(sphere, texture=texture, smooth_shading=False)

    if draw_axes:
        x_arrow = pv.Arrow(start=np.array([x,y,z]), direction=np.array([1.0, 0.0, 0.0]), scale=r*2, tip_length=0.1, tip_radius=0.05, tip_resolution=50, shaft_radius=0.02, shaft_resolution=50)
        y_arrow = pv.Arrow(start=np.array([x,y,z]), direction=np.array([0.0, 1.0, 0.0]), scale=r*2, tip_length=0.1, tip_radius=0.05, tip_resolution=50, shaft_radius=0.02, shaft_resolution=50)
        z_arrow = pv.Arrow(start=np.array([x,y,z]), direction=np.array([0.0, 0.0, 1.0]), scale=r*2, tip_length=0.1, tip_radius=0.05, tip_resolution=50, shaft_radius=0.02, shaft_resolution=50)
        p.add_mesh(x_arrow, show_edges=False, color=axis_color)
        p.add_mesh(y_arrow, show_edges=False, color=axis_color)
        p.add_mesh(z_arrow, show_edges=False, color=axis_color)

def read_trajectory_solution_file(filename):
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
    return et,x,y,z,rmag

##########################################################################################################

if __name__ == "__main__":

    if subplot_test:
        # a test with two plots
        p = pv.Plotter(notebook=False, shape=(1, 2))
        p.subplot(0, 0)
        p.add_text('Moon Texturemap', name='banner1', position='upper_left', color=font_color)
        plot_body( p, name = 'moon', r = 1737.4, texturemap = 'lroc_color_poles_1k.jpg' )
        p.subplot(0, 1)
    else:
        p = pv.Plotter(notebook=False)

    # read file & plot trajectory:
    et,x,y,z,rmag = read_trajectory_solution_file(file_to_plot)
    # add_spline_trajectory_plot(p,x,y,z, name='halo', color='blue', scalars='arc_length') #, radius=100)   # doesn't show all the points ???
    add_line_trajectory_plot(p,x,y,z, scalars = rmag, show_scalar_bar = True)

    # https://visibleearth.nasa.gov/images/57735/the-blue-marble-land-surface-ocean-color-sea-ice-and-clouds/57737l
    # plot_body(p, name = 'earth', r = 6378.0, texturemap = 'land_ocean_ice_cloud_2048.jpg')
    # https://svs.gsfc.nasa.gov/cgi-bin/details.cgi?aid=4720
    plot_body( p, name = 'moon', r = 1737.4, texturemap = 'lroc_color_poles_1k.jpg' )

    # generic star background:
    if plot_starmap:
        stars = examples.planets.download_stars_sky_background(load=False)
        p.add_background_image(stars)

    p.add_text(file_to_plot, name='banner', position='upper_left', color=font_color)
    p.show_grid(color='gray', font_size=20, font_family='times',
                xtitle = 'J2000 X (km)',
                ytitle = 'J2000 Y (km)',
                ztitle = 'J2000 Z (km)')
    p.set_background(color=background_color)

    p.view_isometric()

    ############
    # add defect plot:
    defect_filename = file_to_plot.replace('solution', 'solution_defects')
    defect_filename = os.path.splitext(defect_filename)[0] + '.csv'
    f = utilities.generate_defect_plot(defect_filename, save = False, ylabel='Defect', showtitle=False)
    h_chart = pv.ChartMPL(f, size=(0.46, 0.25), loc=(0.02, 0.06))
    h_chart.background_color = (1.0, 1.0, 1.0, 0.4)
    p.add_chart(h_chart)

    ############
    # eclipse plot
    eclipse_filename = file_to_plot.replace('solution', 'eclipse')
    eclipse_filename = os.path.splitext(eclipse_filename)[0] + '.json'
    f = utilities.generate_eclipse_plot_from_json(eclipse_filename,
                                                  save = False,
                                                  ylabel='solar fraction',
                                                  showtitle=False,
                                                  et0 = et[0],
                                                  etf = et[-1])
    h_chart = pv.ChartMPL(f, size=(0.46, 0.25), loc=(0.02, 0.6))
    h_chart.background_color = (1.0, 1.0, 1.0, 0.4)
    p.add_chart(h_chart)

    if save_as_html:
        # note: background is white here, so would need to change the fonts above to black
        html_filename = os.path.splitext(file_to_plot)[0] + '.html'
        p.export_html(html_filename)

    p.show()

