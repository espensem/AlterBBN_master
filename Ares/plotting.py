'''
* General imports.
'''
import os, sys
import imp
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter


'''
* Class containing methods to make different plots.
*
* The user is free to add methods for further plotting. Remember to call the new method
* at the end of Ares().execute, inside block commented with 'EXECUTE PLOTTING METHODS'.
* If the user wants to create a switch to enable/disable automatic plotting of the user
* defined plot(s), like in the case with the inbuilt plots, this switch must be made
* in the AlterBBN input file (follow the routine of the existing switches). If the user
* defined plotting function requires that values are printed to file from within AlterBBN,
* AlterBBN must be modified accordingly (follow the routine made for the inbuilt plots
* 'abundancesVStime' and 'NeffVSmassWIMP'), and results preferably written to a specified file
* put in the Results/ folder. Remember to read the statement(s) from the inputfile in
* AlterBBN's main-file(default="primary.c").
*
* @author Espen Sem Jenssen
'''
class MakePlots():
    
    def __init__(self):
        self.aresPath = os.path.dirname(__file__)
        self.plotSettingsPath = os.path.join(self.aresPath, "plot_settings")
        self.plotPath = os.path.join(self.aresPath, "Plots")
        # PDG2014 constraints on the observed abundances (He3 is taken from Bania et al. 2002)
        self.Yobs_min, self.Yobs_max = 0.2409, 0.2489
        self.Dobs_min, self.Dobs_max = 2.35e-5, 2.65e-5
        self.He3obs_min, self.He3obs_max = 0.9e-5, 1.3e-5
        self.Li7obs_min, self.Li7obs_max = 1.3e-10, 1.9e-10
        # CMB constraints on eta
        self.etaCMB_min = 6.06e-10
        self.etaCMB_max = 6.14e-10
        self.etaCMB = np.linspace(self.etaCMB_min, self.etaCMB_max, 10)



    '''
    * Function to plot Yp, D/H, He3/H and Li7/H, the most observationally relevant
    * elements, as a function of the parameter to vary. Saves to .eps format in folder Plots/.
    *
    * @param varyParam: the parameter that is being varied.
    * @param x: array of varyParam values.
    * @param elements: dictionary of Element objects.
    *
    * @return void: make plot.
    '''
    def abundancesVSeta(self, varyParam, x, elements):
        
        # Customize plot layout
        if varyParam == "eta":
            ltx = imp.load_source("latexify", os.path.join(self.plotSettingsPath,
                                                           "abundancesVSeta.py"))
        ltx.set_rcParams()

        # Baryon density -> Omega_b*h**2
        barDens = 1e10*(x / 273.9)

        #-------- PLOTTING --------#
        fig = plt.figure(figsize=(11,11))
        fig.subplots_adjust(hspace=0)

        # He4
        ax1 = plt.subplot2grid((3, 1), (0, 0))
        ymin1 = min(elements["He4"].cent-elements["He4"].err)-0.005
        ymax1 = max(elements["He4"].cent+elements["He4"].err)+0.005
        ax1.fill_between(x, self.Yobs_min, self.Yobs_max, facecolor="#f0e68c",
                         lw=0.0, alpha=0.95)
        ax1.fill_between(self.etaCMB, ymin1, ymax1, facecolor="#9e9e9e", lw=0.0, alpha=0.95)
        ax1.plot(x, elements["He4"].cent, color="white", lw=0.0)
        ax1.fill_between(x, elements["He4"].cent+elements["He4"].err,
                         elements["He4"].cent-elements["He4"].err,
                         facecolor='magenta', lw=0.0, alpha=0.7)
        if varyParam == "eta":
            ax11 = ax1.twiny()
            ax11.plot(barDens, elements["He4"].cent,alpha=0)
            ax11.set_xlabel(r'$Baryon \ density \ \Omega_b h^2$')
            ax11.tick_params(labelbottom=False, labeltop=True)
            ax11.set_xscale('log')
            ax1.set_xscale('log')
        ax1.tick_params(bottom=True, labelbottom=False)
        ax1.set_ylabel(r'$Y_p$')
        ax1.set_xlim(min(x), max(x))
        ax1.set_ylim(ymin1, ymax1)
        YmajorLocator1 = MultipleLocator(0.01)
        YmajorFormatter1 = FormatStrFormatter("%.2f")
        YminorLocator1 = MultipleLocator(0.001)
        ax1.yaxis.set_major_locator(YmajorLocator1)
        ax1.yaxis.set_major_formatter(YmajorFormatter1)
        ax1.yaxis.set_minor_locator(YminorLocator1)

        # D and He3
        ax2 = plt.subplot2grid((3, 1), (1, 0))
        ymin2 = 1.05e-6
        ymax2 = 1.0e-3
        ax2.fill_between(x, self.Dobs_min, self.Dobs_max, facecolor="#f0e68c",
                         lw=0.0, alpha=0.95)
        ax2.fill_between(x, self.He3obs_min, self.He3obs_max, facecolor="#f0e68c",
                         lw=0.0, alpha=0.95)
        ax2.fill_between(self.etaCMB, ymin2, ymax2, facecolor="#9e9e9e", lw=0.0, alpha=0.95)
        ax2.plot(x, elements["H2"].cent, color="white", lw=0.0)
        ax2.fill_between(x, elements["H2"].cent+elements["H2"].err,
                         elements["H2"].cent-elements["H2"].err,
                         facecolor='blue', lw=0.0, alpha=0.7)
        ax2.plot(x, elements["He3"].cent, color="white", lw=0.0)
        ax2.fill_between(x, elements["He3"].cent+elements["He3"].err,
                         elements["He3"].cent-elements["He3"].err,
                         facecolor='red', lw=0.0, alpha=0.7)
        ax2.set_xlim(min(x), max(x))
        ax2.tick_params(bottom=True, top=True, labelbottom=False, labeltop=False)
        ax2.set_ylim(ymin2, ymax2)
        if varyParam == "eta":
            ax2.set_xscale('log')
        ax2.set_yscale('log')
        ax2.set_ylabel(r'$[D/H]_p,\, [^3\!He/H]_p$')
        D_patch = mpl.patches.Patch(color="blue", label=r"$[D/H]_p$")
        He3_patch = mpl.patches.Patch(color="red", label=r"$[^3\!He/H]_p$")
        plt.legend(handles=[D_patch, He3_patch])

        # Li7
        ax3 = plt.subplot2grid((3, 1), (2, 0))
        ymin3 = 2.0e-11
        ymax3 = 1.0e-8
        ax3.fill_between(x, self.Li7obs_min, self.Li7obs_max, facecolor="#f0e68c",
                         lw=0.0, alpha=0.95)
        ax3.fill_between(self.etaCMB, ymin3, ymax3, facecolor="#9e9e9e", lw=0.0, alpha=0.95)
        ax3.plot(x, elements["Li7"].cent, color="black", lw=0.0)
        ax3.fill_between(x, elements["Li7"].cent+elements["Li7"].err,
                         elements["Li7"].cent-elements["Li7"].err,
                         facecolor='#adff2f', lw=0.0, alpha=0.7)
        ax3.set_xlim(min(x), max(x))
        ax3.tick_params(bottom=True, top=True, labelbottom=True, labeltop=False)
        ax3.set_ylim(ymin3, ymax3)
        ax3.set_yscale('log')
        ax3.set_ylabel(r'$[^7\!Li/H]_p$')
        if varyParam == "eta":
            ax3.set_xscale('log')
            ax3.set_xlabel(r'$Baryon-to-photon \ ratio, \ \eta$')
            plt.savefig(os.path.join(self.plotPath, "abundancesVSeta.eps"))
        plt.show()
