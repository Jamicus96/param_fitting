import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.font_manager as fm
import scipy

def setup_plot_style_Daniel():
    plt.style.use('classic')
    matplotlib.rcParams.update({'font.size': 28})
    matplotlib.rcParams.update({'font.style': "normal", 'font.family': 'serif'})
    matplotlib.rcParams['font.family'] = 'serif'
    matplotlib.rcParams['xtick.major.size'] = 10
    matplotlib.rcParams['xtick.major.width'] = 2
    matplotlib.rcParams['xtick.minor.size'] = 5
    matplotlib.rcParams['xtick.minor.width'] = 1
    matplotlib.rcParams['ytick.major.size'] = 10
    matplotlib.rcParams['ytick.major.width'] = 2
    matplotlib.rcParams['ytick.minor.size'] = 5
    matplotlib.rcParams['ytick.minor.width'] = 1
    matplotlib.rcParams['axes.linewidth'] = 2  # set the value globally
    matplotlib.rcParams['figure.facecolor'] = 'white'
    matplotlib.rcParams['figure.figsize'] = 18, 10
    matplotlib.rcParams['xtick.major.pad'] = '12'
    matplotlib.rcParams['ytick.major.pad'] = '12'
    matplotlib.rcParams['axes.unicode_minus'] = False
    font_filepath = 'LiberationSerif-Regular.ttf'
    return fm.FontProperties(fname=font_filepath, size=28)

def setup_plot_style():
    prop_font = fm.FontProperties(fname='Times_New_Roman_Normal.ttf',size = 28)

    #matplotlib.rcParams.update({'font.family': 'serif'})
    matplotlib.rcParams.update({'font.size': 28})
    matplotlib.rcParams.update({'font.style': "normal"})

    matplotlib.rcParams['xtick.major.size'] = 10
    matplotlib.rcParams['xtick.major.width'] = 2
    matplotlib.rcParams['xtick.minor.size'] = 5
    matplotlib.rcParams['xtick.minor.width'] = 1
    matplotlib.rcParams['ytick.major.size'] = 10
    matplotlib.rcParams['ytick.major.width'] = 2
    matplotlib.rcParams['ytick.minor.size'] = 5
    matplotlib.rcParams['ytick.minor.width'] = 1
    matplotlib.rcParams['axes.linewidth'] = 2 #set the value globally
    matplotlib.rcParams['figure.facecolor'] = 'white'
    matplotlib.rcParams['figure.figsize'] = 18, 10
    matplotlib.rcParams['xtick.major.pad']='12'
    matplotlib.rcParams['ytick.major.pad']='12'

    return prop_font

### COLLECT DATA ###

data_address = '/Users/jp643/Documents/Studies/PhD/Antinu/param_fitting/likelihoods/updated/param_fits_all_classCUT.txt'
# data_address = '/Users/jp643/Documents/Studies/PhD/Antinu/param_fitting/likelihoods/updated/param_fits_all.txt'
# data_address = '/Users/jp643/Documents/Studies/PhD/Antinu/param_fitting/likelihoods/replicateTony/param_fits_all.txt'

min_E = 0.9
max_E = 8.0

def read_data(data_address):
    Dm21 = []
    theta_12 = []
    minLL = []

    Ebin_centers = []
    spectra = {}

    reading_minLL = False
    reading_spectra = False
    with open(data_address, 'r') as file:
        for line in file:
            if '# Delta log-likelihood:' in line:
                first_line = True
                reading_minLL = True
                continue
            if '# Spectra:' in line:
                first_line = True
                reading_minLL = False
                reading_spectra = True
                continue

            line_splt = line.rstrip().split(' ')

            if reading_minLL:
                if first_line:
                    line_splt.pop(0)
                    Dm21 = line_splt
                    first_line = False
                else:
                    theta_12.append(line_splt[0])
                    line_splt.pop(0)
                    minLL.append(line_splt)
            
            if reading_spectra:
                if first_line:
                    line_splt.pop(0)
                    Ebin_centers = line_splt
                    first_line = False
                else:
                    spec_name = line_splt[0]
                    line_splt.pop(0)
                    spectra[spec_name] = np.asarray(line_splt, dtype=float)

    Dm21 = np.asarray(Dm21, dtype=float)[1:-1] * 1E5
    theta_12 = np.asarray(theta_12, dtype=float)[1:-1]
    minLL = np.asarray(minLL, dtype=float)[1:-1, 1:-1]
    minLL[np.where(minLL < 0)[0], np.where(minLL < 0)[1]] = 50
    minLL = minLL.transpose()

    Ebin_centers = np.asarray(Ebin_centers, dtype=float)
    print('Ebin_centers =', Ebin_centers)
    min_idx = np.where(min_E <= Ebin_centers)[0][0]
    max_idx = np.where(max_E < Ebin_centers)[0][0]
    print('min_idx =', min_idx, ', max_idx =', max_idx)
    print('Ebin_centers =', Ebin_centers[min_idx:max_idx])

    for spec in spectra:
        print(spec + ': ' + str(np.sum(spectra[spec][min_idx:max_idx])))

    return Dm21, theta_12, minLL, Ebin_centers, spectra

def ll_diff_per_nSig(n):
    '''
    Return value of -2log(L/L_max) corresponding to n-sigmas aways from max.
    ppf is the inverse of cdf. This use approximation that -2log(L/L_max) approaches
    a chi^2 distribution with 2 degrees of freedom in this case
    '''
    nDOF = 2
    return scipy.stats.chi2.ppf(2 * scipy.stats.norm.cdf(n) - 1, nDOF)

def inv_ll_diff_per_nSig(nSig):
    '''Inverse function of ll_diff_per_nSig, for plotting purposes'''
    nDOF = 2
    return scipy.stats.norm.ppf(0.5*(1 + scipy.stats.chi2.cdf(nSig, nDOF)))

def plot_LL(Dm21, theta_12, minLL):
    # Get best fit values
    min_indices = np.where(minLL == np.min(minLL))
    min_theta_idx = min_indices[1][0]
    min_Dm_idx = min_indices[0][0]

    min_theta_12 = theta_12[min_theta_idx]
    min_Dm21 = Dm21[min_Dm_idx]

    # Get uncertainties
    # Chi2Diff = ll_diff_per_nSig(1.0)
    Chi2Diff = 1.0
    for i in range(min_Dm_idx, Dm21.size):
        if minLL[i, min_theta_idx] == Chi2Diff:
            min_Dm21_max = Dm21[i] - min_Dm21
            break
        elif minLL[i, min_theta_idx] > Chi2Diff:
            min_Dm21_max = 0.5 * (Dm21[i] + Dm21[i - 1]) - min_Dm21
            break

    for i in range(min_Dm_idx):
        if minLL[min_Dm_idx - i, min_theta_idx] == Chi2Diff:
            min_Dm21_min = min_Dm21 - Dm21[min_Dm_idx - i]
            break
        elif minLL[min_Dm_idx - i, min_theta_idx] > Chi2Diff:
            min_Dm21_min = min_Dm21 - 0.5 * (Dm21[min_Dm_idx - i] + Dm21[min_Dm_idx - i - 1])
            break

    for i in range(min_theta_idx, theta_12.size):
        if minLL[min_Dm_idx, i] == Chi2Diff:
            min_theta12_max = theta_12[i] - min_theta_12
            break
        elif minLL[min_Dm_idx, i] > Chi2Diff:
            min_theta12_max = 0.5 * (theta_12[i] + theta_12[i - 1]) - min_theta_12
            break

    for i in range(min_theta_idx):
        if minLL[min_Dm_idx, min_theta_idx - i] == Chi2Diff:
            min_theta12_min = min_theta_12 - theta_12[min_theta_idx - i]
            break
        elif minLL[min_Dm_idx, min_theta_idx - i] > Chi2Diff:
            min_theta12_min = min_theta_12 - 0.5 * (theta_12[min_theta_idx - i] + theta_12[min_theta_idx - i - 1])
            break

    # Print fit results
    print('Dm21^2 = {} + {} - {} eV^2'.format(min_Dm21, min_Dm21_max, min_Dm21_min))
    print('theta_12 = {} + {} - {} degrees'.format(min_theta_12, min_theta12_max, min_theta12_min))

    ### PLOT DATA ###
    prop_font = setup_plot_style()

    # temp change
    matplotlib.rcParams['figure.figsize'] = 14, 10

    fig = plt.subplot(111)
    ax1 = plt.gca()

    # Background heat map
    im = ax1.imshow(np.sqrt(minLL), cmap=plt.cm.viridis, interpolation='none', extent=[theta_12[0], theta_12[-1], Dm21[0], Dm21[-1]],
                    aspect='auto', origin='lower', vmin=0, vmax=7.1)
    colbar = plt.colorbar(im)
    colbar.ax.set_ylabel(r'$\sqrt{-2\Delta \mathrm{log}(L)}$', rotation=90, fontproperties=prop_font)

    # Best fit point
    ax1.scatter(min_theta_12, min_Dm21, marker='+', s=100, linewidths=2, color='w')

    # Contours
    X, Y = np.meshgrid(theta_12, Dm21)
    contour_vals = [1, 2, 3, 4, 5]
    contour_cols = ['peachpuff', 'sandybrown', 'darkorange', 'red', 'brown']
    contour_styles = ['-', '--', '-.', '-.', '-.', '-.']
    contour_labels = [r'$1 \sigma$ contour', r'$2 \sigma$ contour', r'$3 \sigma$ contour', r'$4 \sigma$ contour', r'$5 \sigma$ contour']


    Chi2_contours = []
    for i in range(len(contour_vals)):
        Chi2_contours.append(ll_diff_per_nSig(contour_vals[i]))
    plot = ax1.contour(X, Y, minLL, levels=Chi2_contours, colors=contour_cols, linestyles=contour_styles, linewidths=3)
    # def fmt(x):
    #     nSig = inv_ll_diff_per_nSig(x)
    #     return ('%i' % round(nSig)) + r'$\sigma$'
    # plt.clabel(plot, inline=True, fmt=fmt, fontsize=18)

    # Legend
    legend_elems = [Line2D([0], [0],  color='w', lw=0, label='Best fit', marker='+', linestyle=None, markersize=8)]
    for i in range(len(contour_vals)):
        legend_elems.append(Line2D([0], [0],  color=contour_cols[i], linestyle=contour_styles[i], lw=3, label=contour_labels[i]))#, fontproperties=prop_font))
        
    leg = ax1.legend(handles=legend_elems,loc='upper right', fancybox=False, numpoints=1, markerscale=1.5, prop=prop_font) #, frameon=False)
    leg.get_frame().set_color('cornflowerblue')
    leg.get_frame().set_alpha(0.9)

    # Labels
    ax1.set_xlabel(r'$\theta_{12}$ $\left[{}^{\degree}\right]$', fontproperties=prop_font, x=1, ha='right')
    ax1.set_ylabel(r'$\Delta m_{12}^2$ $\left[10^{-5} \mathrm{eV}^2\right]$', fontproperties=prop_font, y=1, ha='right')

    for label in ax1.get_xticklabels():
        label.set_fontproperties(prop_font)
    for label in ax1.get_yticklabels():
        label.set_fontproperties(prop_font)
    
    ax1.minorticks_on()
    ax1.get_xaxis().set_tick_params(which='both',direction='in', width=1)
    ax1.get_yaxis().set_tick_params(which='both',direction='in', width=1)
    ax1.xaxis.set_ticks_position('both')

    ax1.set_title(r'Reactor-$\nu$ Sensitivity, from Asimov Dataset', fontproperties=prop_font)
    ax1.text(1.5, 0.3, "SNO+ Preliminary", fontproperties=prop_font)
    plt.show()

def plot_spectra(Ebin_centers, spectra):

    dE = Ebin_centers[1] - Ebin_centers[0]
    bin_edges = np.zeros(len(Ebin_centers)+1)
    bin_edges[:-1] = Ebin_centers - 0.5*dE
    bin_edges[-1] = Ebin_centers[-1] + 0.5*dE

    # ax.stackplot(bin_edges[:-1], spectra['model_BRUCE'], spectra['model_DARLINGTON'], spectra['model_PICKERING'], spectra['model_WORLD'],
    #              spectra['model_geoNu_Th'], spectra['model_geoNu_U'],
    #              spectra['model_alphaN_PR'], spectra['model_alphaN_C12'], spectra['model_alphaN_O16'],
    #              step='post',
    #              labels=[r'reactor-$\nu$ Bruce', r'reactor-$\nu$ Darlington', r'reactor-$\nu$ Pickering', r'reactor-$\nu$ rest',
    #                      r'geo-$\nu$ Th', r'geo-$\nu$ U',
    #                      r'($\alpha$,n) PR', r'($\alpha$,n) ${}^{12}$C', r'($\alpha$,n) ${}^{16}$O'])

    prop_font = setup_plot_style()
    fig = plt.subplot(111)
    ax1 = plt.gca()

    ax1.stackplot(bin_edges[:-1], spectra['Reactor::model_Esys'], spectra['geoNu::model_Esys'], spectra['alphaN::model_Esys'],
                 step='post', labels=[r'reactor-$\nu$ IBD', r'geo-$\nu$ IBD', r'($\alpha$,n)'], colors=['salmon', 'palegreen', 'cornflowerblue'], ec='grey', linewidth=1)

    # plt.stairs(spectra['Total_Model_Spectrum'], edges=bin_edges, color='k', linestyle='dashed', label='total model')
    ax1.scatter(Ebin_centers, spectra['data'], marker='+', linewidths=2, s=140, color='k', label='Asimov data')

    ax1.set_xlabel('E [MeV]', fontproperties=prop_font, size = 32, x=1, ha='right')
    ax1.set_ylabel('arbitrary', fontproperties=prop_font, size = 32, y=1, ha='right')

    for label in ax1.get_xticklabels():
        label.set_fontproperties(prop_font)
    for label in ax1.get_yticklabels():
        label.set_fontproperties(prop_font)

    ax1.set_xlim(min_E, max_E)
    ax1.set_ylim(0, 2.8)

    handles, labels = ax1.get_legend_handles_labels()
    ax1.legend(handles, labels,loc='upper right', fancybox=False, numpoints=1, prop=prop_font, frameon=False)

    ax1.minorticks_on()
    ax1.get_xaxis().set_tick_params(which='both',direction='in', width=1)
    ax1.get_yaxis().set_tick_params(which='both',direction='in', width=1)
    ax1.xaxis.set_ticks_position('both')

    ax1.set_title(r'Fit Prompt Energy Spectrum, from Asimov Dataset', fontproperties=prop_font)
    ax1.text(2.5, 0.05, "SNO+ Preliminary", fontproperties=prop_font)
    plt.show()

def testReduceLivetime(minLL):

    current_livetime = 365.25
    new_livetime = 138.906
    scale = new_livetime / current_livetime

    # ll = -norm + sum[N_i * ln(nu_i)]
    # assume norm, N_i and nu_i all get scaled the same way:
    # ll' = -scale*norm + sum[scale*N_i * ln(scale*nu_i)]
    # ll' = -scale*norm + scale * sum[N_iln(nu_i)] + scale ln(scale) * sum[N_i]
    # assume sum[N_i] = norm:
    # ll' = -scale*norm*(ln(scale) - 1) + scale * sum[N_iln(nu_i)]
    # ll' = -scale*norm*(ln(scale) - 1) + scale * [ll + norm]
    # ll' = scale * [ll - norm*ln(scale)]
    # Shift so that minimum is at 0:
    # ll' = scale * ll

    return scale * minLL


#####################
def fit_log(x, y):
    '''y = a * ln(x) + b'''
    err2 = (0.1 * y)**2
    lnx = np.log(x)

    S1 = np.sum(1 / err2)
    Sx = np.sum(lnx / err2)
    Sy = np.sum(y / err2)
    Sxx = np.sum(lnx**2 / err2)
    Sxy = np.sum(lnx*y / err2)

    Delta = S1 * Sxx - Sx**2

    m = (S1 * Sxy - Sx * Sy) / Delta
    c = (Sy * Sxx - Sx * Sxy) / Delta

    return m, c


def find_min_and_errs(minLL, Dm21, theta_12):
    # Get best fit values
    min_indices = np.where(minLL == np.min(minLL))
    min_theta_idx = min_indices[1][0]
    min_Dm_idx = min_indices[0][0]

    min_theta_12 = theta_12[min_theta_idx]
    min_Dm21 = Dm21[min_Dm_idx]

    # Get uncertainties
    # Chi2Diff = ll_diff_per_nSig(1.0)
    Chi2Diff = 1.0
    for i in range(min_Dm_idx, Dm21.size):
        if minLL[i, min_theta_idx] == Chi2Diff:
            min_Dm21_max = Dm21[i] - min_Dm21
            break
        elif minLL[i, min_theta_idx] > Chi2Diff:
            min_Dm21_max = 0.5 * (Dm21[i] + Dm21[i - 1]) - min_Dm21
            break

    for i in range(min_Dm_idx):
        if minLL[min_Dm_idx - i, min_theta_idx] == Chi2Diff:
            min_Dm21_min = min_Dm21 - Dm21[min_Dm_idx - i]
            break
        elif minLL[min_Dm_idx - i, min_theta_idx] > Chi2Diff:
            min_Dm21_min = min_Dm21 - 0.5 * (Dm21[min_Dm_idx - i] + Dm21[min_Dm_idx - i - 1])
            break

    for i in range(min_theta_idx, theta_12.size):
        if minLL[min_Dm_idx, i] == Chi2Diff:
            min_theta12_max = theta_12[i] - min_theta_12
            break
        elif minLL[min_Dm_idx, i] > Chi2Diff:
            min_theta12_max = 0.5 * (theta_12[i] + theta_12[i - 1]) - min_theta_12
            break

    for i in range(min_theta_idx):
        if minLL[min_Dm_idx, min_theta_idx - i] == Chi2Diff:
            min_theta12_min = min_theta_12 - theta_12[min_theta_idx - i]
            break
        elif minLL[min_Dm_idx, min_theta_idx - i] > Chi2Diff:
            min_theta12_min = min_theta_12 - 0.5 * (theta_12[min_theta_idx - i] + theta_12[min_theta_idx - i - 1])
            break
    
    Dm21_fit = (min_Dm21, min_Dm21_min, min_Dm21_max)
    theta12_fit = (min_theta_12, min_theta12_min, min_theta12_max)
    return Dm21_fit, theta12_fit

def sensitivity_over_time():
    Dm21, theta_12, minLL, Ebin_centers, spectra = read_data('/Users/jp643/Documents/Studies/PhD/Antinu/param_fitting/likelihoods/updated/param_fits_all.txt')
    Dm21_classCUT, theta_12_classCUT, minLL_classCUT, Ebin_centers_classCUT, spectra_classCUT = read_data('/Users/jp643/Documents/Studies/PhD/Antinu/param_fitting/likelihoods/updated/param_fits_all_classCUT.txt')
    
    livetime_init = 134.4 / 365.25
    livetimes = np.linspace(livetime_init, 4, 100)

    Dm21_errs = np.zeros(len(livetimes))
    Dm21_errs_classCUT = np.zeros(len(livetimes))

    Dm21_fit, theta12_fit = find_min_and_errs(minLL, Dm21, theta_12)
    min_Dm21, min_Dm21_min, min_Dm21_max = Dm21_fit
    min_theta_12, min_theta12_min, min_theta12_max = theta12_fit

    Dm21_fit_classCUT, theta12_fit_classCUT = find_min_and_errs(minLL_classCUT, Dm21_classCUT, theta_12_classCUT)
    min_Dm21_classCUT, min_Dm21_min_classCUT, min_Dm21_max_classCUT = Dm21_fit_classCUT
    min_theta_12_classCUT, min_theta12_min_classCUT, min_theta12_max_classCUT = theta12_fit_classCUT

    # Average negative and positive errors
    # Dm21_err = 0.5 * (min_Dm21_max + min_Dm21_min)
    # Dm21_err_classCUT = 0.5 * (min_Dm21_max_classCUT + min_Dm21_min_classCUT)
    Dm21_err = min_Dm21_max
    Dm21_err_classCUT = min_Dm21_max_classCUT

    scaled_Dm21_errs = Dm21_err / np.sqrt(livetimes / livetime_init)
    scaled_Dm21_errs_classCUT = Dm21_err_classCUT / np.sqrt(livetimes / livetime_init)

    prop_font = setup_plot_style()
    fig = plt.subplot(111)
    ax1 = plt.gca()

    ax1.plot(livetimes, scaled_Dm21_errs, linewidth=4, label='no classifier')
    ax1.plot(livetimes, scaled_Dm21_errs_classCUT, linewidth=4, label='using classifier')
    ax1.hlines(0.19, livetimes[0], livetimes[-1], linewidths=4, linestyles='dashed', colors='k', label='KamLand')
    
    ax1.set_xlabel('live-time [years]', fontproperties=prop_font, size = 32, x=1, ha='right')
    ax1.set_ylabel(r'mean $\Delta m_{12}^2$ uncertainty $\left[10^{-5} \mathrm{eV}^2\right]$', fontproperties=prop_font, size = 32, y=1, ha='right')

    for label in ax1.get_xticklabels():
        label.set_fontproperties(prop_font)
    for label in ax1.get_yticklabels():
        label.set_fontproperties(prop_font)

    handles, labels = ax1.get_legend_handles_labels()
    ax1.legend(handles, labels,loc='upper right', fancybox=False, numpoints=1, prop=prop_font, frameon=False)

    ax1.minorticks_on()
    ax1.get_xaxis().set_tick_params(which='both',direction='in', width=1)
    ax1.get_yaxis().set_tick_params(which='both',direction='in', width=1)
    ax1.xaxis.set_ticks_position('both')
    
    ax1.set_title(r'Oscillation Sensitivity Over SNO+ Live-Time, Based on Asimov Data', fontproperties=prop_font)
    ax1.text(2.5, 0.05, "SNO+ Preliminary", fontproperties=prop_font)
    plt.show()


Dm21, theta_12, minLL, Ebin_centers, spectra = read_data(data_address)
plot_LL(Dm21, theta_12, minLL)
plot_spectra(Ebin_centers, spectra)

# sensitivity_over_time()