import matplotlib
import csv
import os
import matplotlib
import matplotlib.pyplot as plt

matplotlib.rcParams["font.family"] = "Serif"
matplotlib.rcParams["font.size"] = 20
matplotlib.rcParams["axes.labelsize"] = 25
matplotlib.rcParams["xtick.labelsize"] = 20
matplotlib.rcParams["ytick.labelsize"] = 20
matplotlib.rcParams["legend.fontsize"] = 20

###################################################################
def generate_eclipse_plot(filename : str):

    """
    Generate a plot of the eclipses for a run.
    """

    plot_filename = os.path.splitext(filename)[0] + '.png'

    et = []
    phi = []
    with open(filename, 'r') as f:
        reader = csv.reader(f)
        for i,row in enumerate(reader):
            if i==0: continue # skip header
            et.append(float(row[0]))
            phi.append(float(row[1]))

    fig = plt.figure(figsize=(20,10),facecolor="white")
    ax = fig.add_subplot(1, 1, 1)
    ax.grid()
    ax.set_axisbelow(True)
    ax.set_xlabel("Ephemeris time (sec)")
    ax.set_ylabel("Phi (<0 means in Earth shadow)")
    ax.set_title(filename)
    ax.plot(et,phi,"b-",linewidth=2,label="phi")

    plt.savefig(plot_filename, dpi=300)

###################################################################
def generate_defect_plot(filename : str):

    """
    Generate a plot of the constraint violations for a run.
    """

    plot_filename = os.path.splitext(filename)[0] + '.png'

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

    fig = plt.figure(figsize=(20,10),facecolor="white")
    ax = fig.add_subplot(1, 1, 1)
    ax.grid()
    ax.set_axisbelow(True)
    ax.set_xlabel("Segment interface number")
    ax.set_ylabel("Constraint violation magnitude")
    ax.set_title(filename)
    ax.plot(idx,r,"b-",linewidth=2,label="r (km)")
    ax.plot(idx,v,"r-",linewidth=2,label="v (km/s)")
    ax.legend()

    plt.savefig(plot_filename, dpi=200)

####################################################################################
if __name__ == "__main__":

    # examples to plot
    # (to generate these files, you have to set the right flags in the config json file)

    generate_eclipse_plot('../eclipse_20220704131415_L2_S_NREVS=100.csv')

    generate_defect_plot('../solution_defects_20220704131415_L2_S_NREVS=20.csv')
    # generate_defect_plot('../solution_defects_20220704131415_L2_S_NREVS=100.csv')