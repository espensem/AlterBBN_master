import matplotlib as mpl

def set_rcParams():
    params = {'axes.labelsize': 14, 
              'axes.titlesize': 14,
              'axes.linewidth': 0.3,
              'font.size': 14,
              'legend.fontsize': 13,
              'legend.handlelength': 3.5,
              'xtick.labelsize': 14,
              'ytick.labelsize': 14,
              'lines.linewidth': 0.3,
              'grid.linewidth': 0.4,
              'text.usetex': True,
#              'backend': 'ps',
              'font.family': 'sans-serif'
              
    }

    mpl.rcParams.update(params)
