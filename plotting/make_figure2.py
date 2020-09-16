from __future__ import division
import os
import matplotlib
import matplotlib.pyplot as plt
from processing_functions import *
matplotlib.rcParams['font.size'] = 6

def make_figure2(num_sims, base_seaspring, base_ns_seaspring, data_name, ax, m):
    "Makes a single panel of figure 2, plotting the L2 norm of results."
    times = [np.power(10.0,j)/(365*24*60*60) for j in range(0, int(num_sims))]
    folders_seaspring = [base_seaspring + str(j) for j in range(num_sims)]
    folders_ns_seaspring = [base_ns_seaspring + str(j) for j in range(num_sims)]

    if m == 0:
        values_seaspring = get_point_values(data_name, folders_seaspring)
        values_ns_seaspring = get_point_values(data_name, folders_ns_seaspring)
    else:
        values_seaspring = get_cell_values(data_name, folders_seaspring)
        values_ns_seaspring = get_cell_values(data_name, folders_ns_seaspring)

    if m==2 or m==3:
        multiplier = 1/1000
    else:
        multiplier = 1.0

    line1, = ax.plot(np.log10(times), [l2_norm(multiplier*values_seaspring[j][0]) for j in range(num_sims)], color='black', label='Seaspring', ls='--', linewidth=1)
    line2 = ax.scatter(np.log10(times), [l2_norm(multiplier*values_ns_seaspring[j][0]) for j in range(num_sims)], color='black', label='Seaspring + NS', marker='o', s=2)

    return line1, line2

if __name__ == "__main__":
    starting_directory = os.getcwd()
    os.chdir(os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', 'tests')))

    num_sims = 10
    base_seaspring = "seaspring_dt_"
    base_ns_seaspring = "ns_seaspring_dt_"
    plotName = "figure2"

    data_names = ["uz", "strainrateE", "tauMax", "sigmaP"]
    labels = ["(a)", "(b)", "(c)", "(d)"]

    fig = plt.figure(figsize=(4.72441,4.72441/16*9))
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)
    ax = [ax1, ax2, ax3, ax4]

    for j in range(len(data_names)):
        lines = (make_figure2(num_sims, base_seaspring, base_ns_seaspring, data_names[j], ax[j], j))
        ax[j].text(-0.2, 1.25, labels[j], transform=ax[j].transAxes, va='top', fontsize=8, weight='bold')
    os.chdir(starting_directory)

    plt.figlegend(lines, ['sea-spring','sea-spring + NS'], loc='lower center', borderaxespad=0.1, ncol=3)

    ax1.set(ylabel=r"$L_{2}$ Norm ($\frac{m}{a}$)", title="Vertical Velocity")
    ax2.set(ylabel=r"$L_{2}$ Norm ($\frac{1}{a}$)", title="Effective Strain Rate")
    ax3.set(ylabel=r"$L_{2}$ Norm ($kPa$)", title="Maximum Shear Stress")
    ax4.set(ylabel=r"$L_{2}$ Norm ($kPa$)", title="Greatest Principal Stress")

    ax3.set(xlabel=r"Time Step ($log_{10}(a)$)")
    ax4.set(xlabel=r"Time Step ($log_{10}(a)$)")

    plt.tight_layout(pad=1.5)

    #plt.savefig(plotName + '.png', transparent=False, dpi=300)
    plt.savefig(plotName + '.eps')