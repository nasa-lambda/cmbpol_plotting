import matplotlib.pyplot as plt
from matplotlib import ticker
import numpy as np

labsize = 16
legsize = 12
texsize = 12

# Makes a B-mode power spectrum plot showing 95% confidence upper limits
# for B-mode power from different experiments. Data from ACTPol, BICEP1,
# BOOMERanG, CAPMAP, CBI, DASI, MAXIPOL, QUaD, QUIET-Q, QUIET-W, and WMAP
# are included.
# For comparison, theoretical curves for a LCDM model with tensor-to-scalar
# ratio r=0.1 and r=0.01 are also plotted as solid curves. The inflationary
# and gravitational lensing components are plotted separately as dashed
# and dotted curves, respectively.
# The data are plotted to a postscript file BB_bounds.ps.

# The experimental data are read from a file bb_data_2015apr.txt, which
# should be copied to the user's local directory.
# Sources for the data and more information are given in
# http://lambda.gsfc.nasa.gov/graphics/bb_upperlimits/

# For the QUIET-W results, the larger of the upper limits from the two
# pipelines is plotted for each l-bin.

# The theoretical data are read from files made by the BICEP2 team,
# from their calculations using the March 2013 version of CAMB.
# They should be copied to the user's local directory from
# http://lambda.gsfc.nasa.gov/data/suborbital/BICEP2/B2_3yr_camb_planck_lensed_uK_20140314.txt
# http://lambda.gsfc.nasa.gov/data/suborbital/BICEP2/B2_3yr_camb_planck_withB_uK_20140314.txt

# This code calls readcol.pro and associated routines from the IDL
# Astronomy User's Library, http://idlastro.gsfc.nasa.gov/

# read and plot 95% confidence upper limits

names = ('experiment', 'l_min', 'l_max', 'BB_limit')
formats = ('S20', 'i4', 'i4', 'f4')
data = np.loadtxt('bb_data_2015apr.txt', skiprows=35,
                  dtype={'names': names, 'formats':formats})

unique_list = np.unique(data['experiment'])
n_unique = len(unique_list)

plt.plot()
plt.ylim([1e-4, 1e3])
plt.xlabel(r'Multipole $\ell$', size=labsize)
plt.ylabel(r'$\ell(\ell+1)C_\ell^\mathrm{BB}/2\pi\,\mathrm{[\mu K^2]}$',
           size=labsize)
plt.xscale('log')
plt.yscale('log')

# colors from colorbrewer
colors = ['#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99', '#e31a1c',
          '#fdbf6f', '#ff7f00', '#cab2d6', '#6a3d9a', '#ffff99'][::-1]
for i in range(n_unique):
    d = data[data['experiment'] == unique_list[i]]
    for j in range(len(d['l_min'])):
        plt.plot([d['l_min'][j], d['l_max'][j]],
                 [d['BB_limit'][j], d['BB_limit'][j]],
                 label=unique_list[i] if j == 0 else "", color=colors[i], lw=3)

names = ('l', 'TT', 'TE', 'EE', 'BB', 'TB', 'EB')
formats = ('i4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4')
theory_lensed = np.loadtxt('B2_3yr_camb_planck_lensed_uK_20140314.txt',
                           comments='#', dtype={'names':names, 'formats':formats})

theory_inflation = np.loadtxt('B2_3yr_camb_planck_withB_uK_20140314.txt',
                              comments='#', dtype={'names':names, 'formats':formats})


imax = min(len(theory_lensed), len(theory_inflation))

theory_lensed = theory_lensed[2:imax]
theory_inflation = theory_inflation[2:imax]

plt.plot(theory_lensed['l'], theory_lensed['BB'] + theory_inflation['BB'], 'k')
plt.plot(theory_lensed['l'], theory_inflation['BB'], 'k--')
plt.plot(theory_lensed['l'], theory_lensed['BB'] + 0.1*theory_inflation['BB'], 'k')
plt.plot(theory_lensed['l'], 0.1*theory_inflation['BB'], 'k--')
plt.plot(theory_lensed['l'], theory_lensed['BB'], 'k:')

plt.text(3, 1.5*(theory_lensed['BB'][0] + theory_inflation['BB'][0]),
         r'$r=0.1$', size=texsize)
plt.text(3, 1.5*(theory_lensed['BB'][0] + 0.1*theory_inflation['BB'][0]),
         r'$r=0.01$', size=texsize)
plt.legend(loc='best', fontsize=legsize)
ax = plt.gca()
ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
ax.tick_params(axis='y', which='minor', left='off')
plt.savefig('BB_bounds.eps', bbox_inches='tight')
plt.show()
