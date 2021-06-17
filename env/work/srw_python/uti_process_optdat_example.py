import numpy as np
import sys, os
import json
import matplotlib.pyplot as plt

if __name__ == "__main__":
    f = open("data_optimize_grid_scan.txt")
    data = json.load(f)
    dat_sorted = sorted(data[:], key = lambda v:v['x']) # sorted according to 'x'
#    datax = [(np.array(v['x'], v['costf'][1]['yFWHM'], v['costf'][1]['maxI'])) for v in dat_sorted]
    datax = [(np.array(v['x']), v['costf'][1]['xFWHM'], v['costf'][1]['yFWHM'], v['costf'][1]['maxI']) for v in dat_sorted]
    datax = np.array(datax)

    plt.plot(datax[:,0]+120.6, datax[:,1]*1e6, label= 'xFWHM [um]')
    plt.plot(datax[:,0]+120.6, datax[:,2]*1e6,label = 'yFWHM [um]')
    plt.legend(loc='upper center', shadow=True, fontsize='x-large')
    plt.title("LCLS SXR Spot Size VS Sample Position")
    plt.xlabel('Sample Position [m]')
    plt.show()
