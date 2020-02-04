from __future__ import print_function

import numpy as np
import pandas
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset, inset_axes

class Plotting(object):

    def __init__(self, title=None, degreescale=False, inset=False):
        self.data = {}
        self.data['TT'] = CMBData('TT_data_2018oct_csv_format.dat', 'TT')
        self.data['EE'] = CMBData('EE_data_2018oct_csv_format.dat', 'EE')
        self.data['TE'] = CMBData('TE_data_2018oct_csv_format.dat', 'TE')
        self.data['BB'] = CMBData('BB_data_2019oct_csv_format.dat', 'BB')
        self.data['lensing'] = CMBData('lensing_data_2019dec_csv_format.dat', '')
        self.load_theory()
        
        self.degreescale = degreescale
        self.inset = inset

        #self.fig = plt.figure(tight_layout=True)
        self.fig = plt.figure()

        self.ax = self.fig.add_subplot(1, 1, 1)
        
        if self.degreescale:
            self.ax2 = self.ax.twiny()

        if title is not None:
            self.fig.suptitle(title)
        
        #if self.degreescale and title is not None:
        #    self.fig.subplots_adjust(hspace=0.3)

        if inset:
            self.axins = inset_axes(self.ax, 1.5, 1, loc=4)
            self.axins.set_xlim(2, 100)
            self.axins.set_ylim(-0.2, 1)
            self.axins.tick_params(
                axis='both',          # changes apply to the x-axis
                which='both',      # both major and minor ticks are affected
                bottom='off',      # ticks along the bottom edge are off
                top='off',         # ticks along the top edge are off
                left='off',
                right='off',
                labelbottom='off',
                labelleft='off')

        self.xlabel(r'$\ell$', string2=r'Angle ($^\circ$)')
        self.ylabel(r'$\ell (\ell+1) C_\ell / 2\pi$ ($\mu$K$^2$)')
         
        self.xlim([2, 1500])

        self.llp1 = None

    def load_theory(self):
       
        self.theory = {}
        theory_lensCl = np.loadtxt('B2_3yr_camb_planck_lensed_uK_20140314.txt')
        theory_inf = np.loadtxt('B2_3yr_camb_planck_withB_uK_20140314.txt')
        tmp = np.loadtxt('base_plikHM_TT_lowTEB.minimum.theory_cl')
        tmp2 = np.loadtxt('cl_bb_planck18_lmax4000.txt')

        self.theory['TT'] = np.array([theory_inf[:, 0], theory_inf[:, 1]])
        self.theory['EE'] = np.array([theory_inf[:, 0], theory_inf[:, 3]])
        self.theory['TE'] = np.array([theory_inf[:, 0], theory_inf[:, 2]])

        self.theory['BB-inf'] = np.array([theory_inf[:, 0], theory_inf[:, 4]])
        #self.theory['BB-lens'] = np.array([theory_lensCl[:, 0], theory_lensCl[:, 4]])
        self.theory['BB-lens'] = np.array([tmp2[:, 0], tmp2[:, 1]])

        self.theory['lensing'] = np.array([tmp[:, 0], tmp[:, 5]*1e7])

    def plot_measurement(self, experiment, cltype, color='b', bins=None, label=None, doub=False):
        '''Plot the power spectrum measurements for a given
        experiment'''

        data = self.data[cltype]

        ell_center, ell_minus, ell_plus, binval, sigma_plus, sigma_minus, upper_bound = data.get_data(experiment)

        sigmas = np.array([sigma_plus, sigma_minus])
        
        if label is None:
            label = experiment

        #if ell_minus != np.array(None):
        #    xerr = np.array([ell_minus, ell_plus])
        #else:
        #    #xerr = None
        #    xerr = np.zeros_like(sigmas)

        xerr = np.array([ell_minus, ell_plus])

        if bins is not None:
            ell_center = ell_center[bins]
            binval = binval[bins]
            upper_bound = upper_bound[bins]
            xerr = xerr[:, bins]
            sigmas = sigmas[:, bins]
            sigma_minus = sigma_minus[bins]

        #Determine points to plot as error bars and ones to plot as upperbounds
        if doub:
            i_ub = (sigma_minus > binval) | np.isnan(binval)
            i_bin = ~i_ub
        else:
            i_bin = np.isfinite(binval)
            i_ub = i_bin & False

        #Plot errorbars
        if np.any(i_bin):
            self.ax.errorbar(ell_center[i_bin], binval[i_bin], xerr=xerr[:, i_bin], yerr=sigmas[:, i_bin], color=color, fmt='o', label=label)
            label=None
            if self.inset:
                self.axins.errorbar(ell_center[i_bin], binval[i_bin], xerr=xerr[:, i_bin], yerr=sigmas[:, i_bin], color=color, fmt='o', label=label)

        if np.any(i_ub):
            self.ax.errorbar(ell_center[i_ub], upper_bound[i_ub], xerr=xerr[:, i_ub], yerr=sigmas[:, i_ub], color=color, fmt='o', label=label,
                             uplims=True)

        
        self.ax.legend(loc=0, prop={'size': 12})

    def list_experiments(self, cltype):
        '''Returns a list of the different experiments for which we
        have results given the input cltype'''

        experiments = self.data[cltype].experiments()

        print(experiments)

    def plot_theory(self, cltype, color, r=0.01, log=False, llp1=None):

        if cltype != 'BB':
            ell = self.theory[cltype][0]
        else:
            nell_a = len(self.theory['BB-inf'][0])
            nell_b = len(self.theory['BB-lens'][0])
            nell = min(nell_a, nell_b)
            ell = self.theory['BB-inf'][0]

        if llp1 is None and self.llp1 is None:
            llp1 = True
            self.llp1 = True
        elif llp1 is None:
            llp1 = self.llp1
        elif self.llp1 is None:
            self.llp1 = llp1
        elif self.llp1 != llp1:
            raise ValueError("Requested normalization does not match previous plotted lines")

        if self.llp1:
            fact = 1
        else:
            fact = ell*(ell+1) / 2*np.pi

        if 'BB' not in cltype:
            cl_theory = self.theory[cltype][1] / fact
        elif cltype == 'BB':
            cl_inf = r/0.1*self.theory['BB-inf'][1][:nell]
            cl_lens = self.theory['BB-lens'][1][:nell]
            cl_theory = cl_inf + cl_lens
            cl_theory /= fact
        elif cltype == 'BB-inf':
            cl_theory = r/0.1*self.theory[cltype][1] / fact
        elif cltype == 'BB-lens':
            cl_theory = self.theory[cltype][1] / fact
        else:
            raise ValueError('cltype is not valid')

        if log:
            self.ax.loglog(ell, cl_theory, color)
        else:
            self.ax.plot(ell, cl_theory, color)

        if self.inset:
            self.axins.semilogx(ell, cl_theory, color)

    def xlabel(self, string, string2=None):
        '''Sets the xlabel of the plot'''
        self.ax.set_xlabel(string)

        if string2 is not None and self.degreescale:
            self.ax2.set_xlabel(string2)
    
    def ylabel(self, string):
        '''Sets the xlabel of the plot'''
        self.ax.set_ylabel(string)

    def title(self, string):
        self.ax.set_title(string)

    def xlim(self, val):
        self.ax.set_xlim(val)

        if self.degreescale:
            val2 = [180.0/val[0], 180.0/val[1]]
            self.ax2.set_xlim(val2)

    def ylim(self, val):
        self.ax.set_ylim(val)

    def set_axes(self, xscale='linear', yscale='linear'):
        self.ax.set_xscale(xscale)
        self.ax.set_yscale(yscale)
        if self.degreescale:
            self.ax2.set_xscale(xscale)

    def default_BB_plot(self):
        '''Generate a default BB plot with most of the current measurements plotted'''

        self.plot_theory('BB', 'k', r=0.1)
        self.plot_theory('BB-lens', 'k')
        self.plot_theory('BB-inf', 'k', r=0.1)

        self.plot_measurement('BICEP2+Keck', 'BB', color='b')
        self.plot_measurement('BICEP2+Keck/Planck', 'BB', color='c')
        self.plot_measurement('POLARBEAR', 'BB', color='g')
        self.plot_measurement('SPTpol', 'BB', color='r')
        self.set_axes(xscale='log', yscale='log')
        self.xlim([2, 5000])
        self.ylim([1e-3, 0.6])

    def default_TT_plot(self):
        '''Generate a default TT plot'''

        self.plot_theory('TT', 'k')
        self.plot_measurement('Planck binned', 'TT', color='b', label='Planck')
        self.plot_measurement('ACTPol', 'TT', color='g')
        self.plot_measurement('SPT', 'TT', color='r')
        self.set_axes(xscale='log', yscale='log')
        self.xlim([2, 5000])
        self.ylim([0.05, 10000])

    def default_TE_plot(self):
        '''Generate a default TE plot'''

        self.plot_theory('TE', 'k')
        self.plot_measurement('ACTPol_2016', 'TE', color='g', label='ACTPol 2016')
        self.plot_measurement('BICEP2/Keck_2015', 'TE', color='c', label='BICEP2/Keck 2015')
        self.plot_measurement('Planck_2018', 'TE', color='b', label='Planck 2018')
        self.plot_measurement('SPTpol_2017', 'TE', color='r', label='SPTPol 2017')
        self.plot_measurement('WMAP_2013', 'TE', color='m', label='WMAP 2013')
        self.set_axes(xscale='log', yscale='linear')
        self.xlim([2, 5000])
        self.ylim([-200, 200])
    
    def default_EE_plot(self):
        '''Generate a default EE plot'''

        self.plot_theory('EE', 'k')
        self.plot_measurement('ACTPol_2016', 'EE', color='g', label='ACTPol 2016')
        self.plot_measurement('BICEP2/Keck_2015', 'EE', color='c', label='BICEP2/Keck 2015')
        self.plot_measurement('Planck_2018', 'EE', color='b', label='Planck 2018')
        self.plot_measurement('SPTpol_2017', 'EE', color='r', label='SPTPol 2017')
        self.plot_measurement('WMAP_2013', 'EE', color='m', label='WMAP 2013')
        self.set_axes(xscale='log', yscale='log')
        self.xlim([2, 5000])
        self.ylim([-5, 50])

    def default_lensing_plot(self):
        '''Generate a default dd plot'''

        self.plot_theory('lensing', 'k')
        self.plot_measurement('POLARBEAR_2019', 'lensing', color='g', label='POLARBEAR 2019')
        self.plot_measurement('ACTPol_2016', 'lensing', color='c', label='ACTPol 2016')
        self.plot_measurement('SPTpol_2019', 'lensing', color='r', label='SPTPol 2019')
        self.set_axes(xscale='log', yscale='log')
        self.xlim([2, 5000])
        self.ylabel(r'$10^7 \ell (\ell+1) C_\ell / 2\pi$ ($\mu$K$^2$)')
        

class CMBData(object):

    def __init__(self, filename, datatype):
        self.data = pandas.read_csv(filename, comment='#', skipinitialspace=True)

        columns = self.data.columns

        #Remove the whitespace at the end of each column name
        self.data.rename(columns=lambda x: x.rstrip(), inplace=True)
        #for i in range(len(self.data['Experiment'])):
        #    #self.data.loc[:,('Experiment', i)] = self.data['Experiment'][i].rstrip()
        #    self.data['Experiment'][i] = self.data['Experiment'][i].rstrip()
        self.data['Experiment'] = self.data['Experiment'].str.strip()
        self.data['l_min'].astype(float)
        self.data['l_center'].astype(float)
        self.data['l_max'].astype(float)
        self.data['Power'].astype(float)
        self.data['Sigma_minus'].astype(float)
        self.data['Sigma_plus'].astype(float)
        self.data['Upper Limit'].astype(float)
        
        #self.datatype = datatype

        #if datatype != 'lensing':
        #    self.sigma_plus = 'sigma_' + datatype + '_plus'
        #    self.sigma_minus = 'sigma_' + datatype + '_minus'
        #    self.upper_bound = datatype + '_limit'
        #    self.binval = datatype
        #else:
        #    self.sigma_plus = 'sigma_power'
        #    self.sigma_minus = 'sigma_power'
        #    self.binval = 'power'
        
        if datatype == 'BB':
            self._eval_ub()

    def experiments(self):
        '''Return a list of experiments'''
        return list(set(self.data['Experiment']))

    def get_data(self, experiment):

        npts = len(self.data)
       
        datatype = 'Power'
        sigplus = 'Sigma_plus'
        sigminus = 'Sigma_minus'
        ub = 'Upper Limit'

        ell_center = []

        ell_minus = None
        ell_plus = None
        binval = []
        sigma_plus = []
        sigma_minus = []
        upper_bound = []
        for i in range(npts):
            if self.data['Experiment'][i] == experiment:
                ell_center.append(self.data['l_center'][i])
                binval.append(self.data['Power'][i])
                sigma_plus.append(self.data['Sigma_plus'][i])
                sigma_minus.append(self.data['Sigma_minus'][i])
                if 'l_min' in self.data.columns and 'l_max' in self.data.columns:
                    if ell_minus is None:
                        ell_minus = []
                    if ell_plus is None:
                        ell_plus = []
                    ell_minus.append(self.data['l_center'][i] - self.data['l_min'][i])
                    ell_plus.append(self.data['l_max'][i] - self.data['l_center'][i])
                upper_bound.append(self.data['Upper Limit'][i])
        
        return np.array(ell_center), np.array(ell_minus), np.array(ell_plus), np.array(binval), np.array(sigma_plus), np.array(sigma_minus), np.array(upper_bound)
    
    def _eval_ub(self):
        '''Adding in an upper bound for bins that might need it. Some experiments report
        measurements that are not very significant and people might want to plot these as
        upper limits.'''

        #upper_bound = self.data[self.upper_bound]
        #binval = self.data[self.binval]
        #sigma_minus = self.data[self.sigma_minus]
        #sigma_plus = self.data[self.sigma_plus]

        upper_bound = self.data['Upper Limit']
        binval = self.data['Power']
        sigma_minus = self.data['Sigma_minus']
        sigma_plus = self.data['Sigma_plus']

        #Points that don't already have an upper bound and are
        #not 2 sigma measurements
        idx = np.isnan(upper_bound) & (binval - 2*sigma_minus <= 0)

        pandas.options.mode.chained_assignment = None
        upper_bound[idx] = binval[idx] + 2*sigma_plus[idx]

        #Hopefully to deal with a possible issue in pandas where 
        #upper_bound is a copy of the data and
        #not the data itself
        self.data['Upper Limit'] = upper_bound

