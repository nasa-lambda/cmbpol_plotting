import matplotlib.pyplot as plt
from matplotlib import ticker
import numpy as np

labsize = 16
legsize = 12
texsize = 12

# This makes a plot of observed B-mode power from CMB experiments with
# significant detections. Data from BICEP2+Keck October 2015,
# BICEP2+Keck/Planck January 2015, POLARBEAR, and SPTpol are included.
# The plotted BICEP2+Keck data are the CMB component from a spectral
# decomposition into CMB, dust, and synchrotron components. The plotted
# BICEP2+Keck/Planck data have had the dust foreground subtracted based
# on the measured cross-power between Planck and BICEP2+Keck. The other
# data plotted have not had any dust foreground subtraction.
# For comparison, a theoretical curve for a LCDM model with tensor-to-scalar
# ratio r=0.1 is also plotted as a solid curve. The inflationary and
# gravitational lensing components are plotted separately as dashed
# and dotted curves, respectively.
# The data are plotted to a postscript file BB_detections.ps.

# The experimental data are read from a file bb_data_2015nov.txt, which
# should be copied to the user's local directory.
# Sources for the data and more information are given in
# http://lambda.gsfc.nasa.gov/graphics/bb_upperlimits/

# The theoretical data are read from files made by the BICEP2 team,
# from their calculations using the March 2013 version of CAMB.
# They should be copied to the user's local directory from
# http://lambda.gsfc.nasa.gov/data/suborbital/BICEP2/B2_3yr_camb_planck_lensed_uK_20140314.txt
# http://lambda.gsfc.nasa.gov/data/suborbital/BICEP2/B2_3yr_camb_planck_withB_uK_20140314.txt

def get_limit(mu, sigma, confidence):
    upper_limit = 0.
    n_iter = 10
    for i in range(n_iter):
        like = np.random.randn(10000)*sigma + mu
        like = like[like >= 0]
        n = len(like)
        frac = np.zeros(n)
        like.sort()
        for i in range(n):
            frac[i] = like[:i].sum()/like.sum()
        upper_limit += like[frac > confidence].min()/n_iter
    return upper_limit

names = ('experiment', 'l_min', 'l_center', 'l_max', 'BB', 'sigma_BB_minus', 'sigma_BB_plus')
formats = ('S20', 'i4', 'f4', 'i4', 'f4', 'f4', 'f4')
data = np.genfromtxt('bb_data_2015nov.txt', skip_header=3, max_rows=36,
                     dtype={'names': names, 'formats':formats})



plt.plot()
plt.xlabel(r'Multipole $\ell$', size=labsize)
plt.ylabel(r'$\ell(\ell+1)C_\ell^\mathrm{BB}/2\pi\ \mathrm{[\mu K^2]}$',
           size=labsize)
plt.xlim([1, 4e3])
plt.ylim([1e-3, 0.6])
plt.xscale('log')
plt.yscale('log')

experiments = [r'BICEP2+Keck', r'BICEP2+Keck/Planck', r'POLARBEAR', r'SPTpol']

colors = ['b', 'c', 'g', 'r']

for c, e in zip(colors, experiments):
    d = data[data['experiment'] == e]

    inds = (d['BB'] > 0)
    d = d[inds]
    plt.errorbar(d['l_center'], d['BB'],
                 yerr=[d['sigma_BB_minus'], d['sigma_BB_plus']],
                 xerr=[d['l_center'] - d['l_min'], d['l_max'] - d['l_center']],
                 fmt='.', label=e, color=c, capthick=0, capsize=0)

    d = data[data['experiment'] == e]
    inds = (d['BB'] > 0)
    d = d[~inds]
    ulim = get_limit(d['BB'], d['sigma_BB_plus'], 0.95)
    plt.errorbar(d['l_center'], ulim, uplims=True,
                 yerr=d['sigma_BB_plus'],
                 xerr=[d['l_center'] - d['l_min'], d['l_max'] - d['l_center']],
                 fmt='', color=c, capthick=0, capsize=0)


# oplot LCDM predictions

# read BICEP2 team results for r=0.1

names = ('l', 'TT', 'TE', 'EE', 'BB', 'TB', 'EB')
formats = ('i4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4')
theory_lensed = np.loadtxt('B2_3yr_camb_planck_lensed_uK_20140314.txt',
                           comments='#', dtype={'names':names, 'formats':formats})

theory_inflation = np.loadtxt('B2_3yr_camb_planck_withB_uK_20140314.txt',
                              comments='#', dtype={'names':names, 'formats':formats})

imax = min(len(theory_lensed), len(theory_inflation))

plt.plot(theory_lensed['l'][2:imax],
         theory_lensed['BB'][2:imax]+theory_inflation['BB'][2:imax], 'k')
plt.plot(theory_lensed['l'][2:imax], theory_lensed['BB'][2:imax], 'k:')
plt.plot(theory_lensed['l'][2:imax], theory_inflation['BB'][2:imax], 'k--')

plt.text(3., 1.3*(theory_lensed['BB'][2] + theory_inflation['BB'][2]),
         r'$r=0.1$', fontsize=texsize)

plt.legend(loc='best', fontsize=legsize)

ax = plt.gca()
for axis in [ax.xaxis, ax.yaxis]:
    axis.set_major_formatter(ticker.ScalarFormatter())

plt.savefig('BB_detections.eps', bbox_inches='tight')
plt.show()
