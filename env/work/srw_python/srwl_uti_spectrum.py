import numpy as np
import glob
import os
import re
from scipy import interpolate
import matplotlib.pyplot as plt
### calc fwhm by interpolation when the curve is only gaussian like.
### calc rms (root mean square )

def read_data(fname):
    # read header
    # one example of header
    #Intensity [ph/s/.1%bw/mm^2] (C-aligned, inner loop is vs Photon Energy, outer loop vs Vertical Position)
    #1000.0 #Initial Photon Energy [eV]
    #1000.0 #Final Photon Energy [eV]
    #1 #Number of points vs Photon Energy
    #-0.0001180920658348706 #Initial Horizontal Position [m]
    #0.00011735648064427847 #Final Horizontal Position [m]
    #308 #Number of points vs Horizontal Position
    #-9.013112987254392e-05 #Initial Vertical Position [m]
    #8.924050605957397e-05 #Final Vertical Position [m]
    #616 #Number of points vs Vertical Position
    #1 #Number of components

    header = {}
    for line in open(fname, 'r'):
        if not line.startswith('#'): break
        if line.startswith("#Intensity ["):
            header['label'] = line
        elif line.find("Initial Photon Energy") > 0:
            header["E0"] = float(line.split()[0][1:])
        elif line.find('#Number of points vs Horizontal Position') > 0:
            header['nx'] = int(line.split()[0][1:])
        elif line.find('#Number of points vs Vertical Position') > 0:
            header['ny'] = int(line.split()[0][1:])
        elif line.find('#Initial Horizontal Position [m]') > 0:
            header['xl'] = float(line.split()[0][1:])
        elif line.find('#Final Horizontal Position [m]') > 0:
            header['xr'] = float(line.split()[0][1:])
        elif line.find('#Initial Vertical Position [m]') > 0:
            header['yl'] = float(line.split()[0][1:])
        elif line.find('#Final Vertical Position [m]') > 0:
            header['yr'] = float(line.split()[0][1:])
        # print(line.strip())
    # print(fname, header)
    d = np.loadtxt(fname)
    # print(d[:5])
    if len(d) != header['nx'] * header['ny']:
        raise RuntimeError("invalid matrix size: {} != {}*{}".format(len(d), header['nx'], header['ny']))
    m = np.reshape(d, (header['ny'], header['nx']))
    return header, m

def Dispersion_rate_theoretical(a0, beta, f, m, E, hc):
# analytically calculate the energy to vertical displacement rate /dispersion rate
# a0: groove density, [mm^-1]
# beta: [rad]
# f: focal length [mm]
# m: the diffraction order
# E: photon energy [eV]
# hc = 1.239841984e-6 #[eV*m]
    return a0/np.cos(beta)*f*hc/E**2  # energy varies 1eV, then the spot move ocef_dE_to_dl [m]

def fwhm_interp(x, y, n=10):
    # x: array_like, A 1-D array of real values.
    # y: array_like, A 1-D array of real values.
    # n: int, optional, for the interpolate accuracy
    from scipy import interpolate
    f = interpolate.interp1d(x, y, kind='quadratic')
    x0 = np.linspace(x[0], x[-1], n*len(x))
    y0 = f(x0)
    f_interp = np.zeros([len(x0),2])
    f_interp[:,0] = x0
    f_interp[:,1] = y0
    fwhm = np.compress(f_interp[:,1] >= 0.5*np.max(f_interp[:,1]), f_interp, axis=0)
    FWHM = fwhm[-1,0] - fwhm[0,0]
    
    xc = x0-np.mean(x0)
    ynobk = np.clip(y0 - np.max(y0)/500, 0, np.max(y0)) # 1/500 of peak regarded as backgroud
    RMS = np.sqrt(np.sum(xc*xc*ynobk)/np.sum(ynobk))
    return FWHM, RMS, f_interp

def integrate_flux_trapz(fname):
    
    header, m = read_data(fname)
    x_min, x_max = header['xl'], header['xr']
    y_min, y_max = header['yl'], header['yr']
    nx, ny = np.shape(m)
    m = m - m[0,0]  # remove background
    xsum = [np.trapz(m[:,i], np.linspace(x_min, x_max, nx)) for i in range(ny)]
    integ2d = np.trapz(xsum, np.linspace(y_min, y_max, ny))
    
    return integ2d*1e6

### calculate the spectral transmission just after exit slit/ energy resolution of SXN beamline
def energy_resolution_curves_method1(s, y, m1, DR_numerical, n=10):
    # this method is used to calculate the energy resolution curve through SRW simulation code;
    # this method calculats the energy resolution curve just after the exit slit for a given slit size;
    # this method processes the intensity distribution just before the exit slit; the intensity distribution is the partially coherent intensity distribution considering the 5D electron beam distribution at fixed photon energy;
    # the steps of this method are:
    #      1, for a given exit slit, calc the flux at a fixed vertical position;
    #      2, shift the slit along the vertical position to get different amount of flux;

    # s is the exit slit size
    # y is the vertical/y-axis
    # m1 is the intensity integrate along x-axis
    # DR_numerical is the dispersion rate/ coverted coefficient
    # n is the accuracy for the interpolation
    # return FWHM_eV is the dE in the energy resolution dE/E
    cent = np.linspace(y[0], y[-1], np.shape(m1)[0])
    flux = np.zeros((len(cent), 3), 'd')
    for i,cy in enumerate(cent):
        w1, w2 = cy - s/2.0, cy+s/2.0
        f1 = np.interp(w1, y, m1)
        f2 = np.interp(w2, y, m1)
        curv = [[w1, f1],]
        for j,yi in enumerate(y):
            if yi <= w1 or yi >= w2: continue
            curv.append([yi, m1[j]])
        curv.append([w2, f2])
        curv = np.array(curv)
        flux[i,-1] = 1e6*np.trapz(curv[:,1], x=curv[:,0])
        flux[i,0] = cy
        flux[i,1] = cy/DR_numerical
    FWHM_eV, RMS_eV, spectrum_trans_curve = fwhm_interp(flux[:,1], flux[:,-1],n)
    return FWHM_eV, RMS_eV, spectrum_trans_curve

def spectral_resolution_at_exitslit(fn, DR, E):
    # fn: file name; partially coherent intensity distribution considering the 5D electron beam distribution at fixed photon energy
    # DR: dispersion rate; unit [m/eV]
    # E: photon energy; unit [eV]
    print("calc spectral resolution at exit slit....")
    header, m = read_data(fn)
    x = np.linspace(header['xl'], header['xr'], header['nx'])
    y = np.linspace(header['yl'], header['yr'], header['ny'])
    m1 = np.trapz(m - m[0,0], x, axis=1) # integrate along x-axis   # m unit is [ph/s/.1%bw/mm^2]
    
    s = np.linspace(1e-6, 50e-6,101) # different slit size
    dEE = [ ]
    for w in s:
        fwhm_eV, rms_eV, curve_data = energy_resolution_curves_method1(w, y, m1, DR, n=30)
        #dEE.append(fwhm_eV/E) # calculate deltaE by FWHM
        dEE.append(2.35482*rms_eV/E)  # calculate deltaE by RMS
    dEE = np.array(dEE)
    return s, dEE

def energy_resolution_curves_method2(file_dir):
    # this method is used to calculate the energy resolution curve through SRW simulation code;
    # this method calculats the energy resolution curve just after the exit slit for a given slit size;
    # this method processes the intensity distribution just before the exit slit; the intensity distribution is the partially coherent intensity distribution considering the 5D $
    # the steps of this method are:
    #	   1, for a given exit slit, calc the flux at a fixed vertical position;
    #	   2, shift the slit along the vertical position to get different amount of flux;
    Flux = []
    for f in glob.glob(os.path.join(file_dir, "*.dat")):
    #print(f)
        m = re.match(".+_([-\.0-9]+)eV.dat", f)
        E = float(m.group(1))
        #print(E)
        flux = integrate_flux_trapz(f)
        Flux.append([E,flux])
    Flux = sorted(Flux, key = lambda x: x[0])
    Flux = np.array(Flux)

    f = interpolate.interp1d(Flux[:,0], Flux[:,1], kind='quadratic')
    x = np.arange(-0.1,0.1,0.0001)
    y = f(x)
    FWHM_eV, RMS_eV, spectrum_trans_curve = fwhm_interp(Flux[:,0], Flux[:,1],51)
    
    plt.plot(Flux[:,0], Flux[:,1]/np.max(Flux[:,1]), '*')
    plt.plot(x, y/np.max(y), '', label = 'at sample with 0.32mm aperture, $\Delta$E_FWHM={:.3f}eV'.format(2.35482*RMS_eV))
    plt.xlabel('Energy[eV]')
    plt.ylabel('Flux/Flux_max')
    plt.title('Energy Resolution Curves')
    return FWHM_eV, RMS_eV, spectrum_trans_curve

if  __name__ == "__main__":

    f_dir = "data_spectrum"
    FWHM_eV, RMS_eV, spectrum_trans_curve = energy_resolution_curves_method2(f_dir)
    plt.show() 

    E = 1000 #[eV]
    DR = 2.059e-5/0.1 # unit [m/eV]
    fn = "res_int_pr_me_bfexitslit_testex.dat"
    s, dEE = spectral_resolution_at_exitslit(fn, DR, E)

    print("Energy Resolution of {}eV SXN Beamline at Exit Slit ".format(E))
    plt.plot(s*1e6, dEE)
    plt.xlabel('Exit slit size [$\mu$m]')#, size = 20, fontweight='bold')
    plt.ylabel('$\Delta$E/E')#, size = 20, fontweight='bold')
    plt.grid(True)#, linestyle='--')#, linewidth=1)
    plt.ticklabel_format(style='sci',scilimits=(-3,4),axis='both')
    plt.title("the energy resolution of {}eV SXN beamline ".format(E))
    plt.show()
    #plt.savefig("Spectral Resolution at 1000eV.png", bbox_inches='tight')

