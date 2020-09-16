from __future__ import division
import os
import matplotlib
import matplotlib.pyplot as plt
from processing_functions import *
matplotlib.rcParams['font.size'] = 6

def make_figure3(num_sims, base_seaspring, base_ns_seaspring, data_name, ax):
    "Makes a single panel of figure 2, plotting the maximum value of results."
    times = [np.power(10.0,j)/(365*24*60*60) for j in range(0, int(num_sims))]
    folders_seaspring = [base_seaspring + str(j) for j in range(num_sims)]
    folders_ns_seaspring = [base_ns_seaspring + str(j) for j in range(num_sims)]

    values_seaspring = get_cell_values(data_name, folders_seaspring)
    values_ns_seaspring = get_cell_values(data_name, folders_ns_seaspring)

    multiplier = 1/1000

    line1, = ax.plot(np.log10(times), [field_max(multiplier*values_seaspring[j][0]) for j in range(num_sims)], color='black', label='Seaspring', ls='--', linewidth=1)
    line2 = ax.scatter(np.log10(times), [field_max(multiplier*values_ns_seaspring[j][0]) for j in range(num_sims)], color='black', label='Seaspring + NS', marker='o', s=2)

    return line1, line2

if __name__ == "__main__":
    starting_directory = os.getcwd()
    os.chdir(os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', 'tests')))

    num_sims = 10
    base_seaspring = "seaspring_dt_"
    base_ns_seaspring = "ns_seaspring_dt_"
    plotName = "figure3"

    data_names = ["tauMax", "sigmaP"]
    labels = ["(a)", "(b)"]

    fig = plt.figure(figsize=(3.26772,3.26772/16*9*2))
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)
    ax = [ax1, ax2]

    for j in range(len(data_names)):
        lines = (make_figure3(num_sims, base_seaspring, base_ns_seaspring, data_names[j], ax[j]))
        ax[j].text(-0.2, 1.1, labels[j], transform=ax[j].transAxes, va='top', fontsize=8, weight='bold')
    os.chdir(starting_directory)

    plt.figlegend(lines, ['sea-spring', 'sea-spring + NS'], loc='lower center', borderaxespad=0.1, ncol=3)

    ax2.set(xlabel=r"Time Step ($log_{10}(a)$)")
    ax1.set(ylabel=r"Maximum ($kPa$)", title="Maximum Shear Stress")
    ax2.set(ylabel=r"Maximum ($kPa$)", title="Greatest Principal Stress")

    plt.tight_layout(pad=1.5)

    #plt.savefig(plotName + '.png', transparent=False, dpi=300)
    plt.savefig(plotName + '.eps')