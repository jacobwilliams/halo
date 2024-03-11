import matplotlib
import csv
import os
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter
import json

matplotlib.rcParams["font.family"] = "Serif"
# matplotlib.rcParams["font.size"] = 20
# matplotlib.rcParams["axes.labelsize"] = 25
# matplotlib.rcParams["xtick.labelsize"] = 20
# matplotlib.rcParams["ytick.labelsize"] = 20
# matplotlib.rcParams["legend.fontsize"] = 20

J2000_DATENUM = -946728000.000000  # et of python date epoch (1970-01-01)
SEC2DAY = 3600.0*24.0

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

####################################################################################
if __name__ == "__main__":

    # examples to plot
    # (to generate these files, you have to set the right flags in the config json file)

    generate_eclipse_plot('../eclipse_20220704131415_L2_S_NREVS=100.csv')

    generate_defect_plot('../solution_defects_20220704131415_L2_S_NREVS=20.csv')
    # generate_defect_plot('../solution_defects_20220704131415_L2_S_NREVS=100.csv')