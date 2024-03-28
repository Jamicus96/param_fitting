import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.markers as mmark
import matplotlib.font_manager as fm
import scipy

def setup_plot_style():
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


### COLLECT DATA ###

data_address = '/Users/jp643/Documents/Studies/PhD/Antinu/param_fitting/likelihoods/updated/param_fits_all_classCUT.txt'
# data_address = '/Users/jp643/Documents/Studies/PhD/Antinu/param_fitting/likelihoods/updated/param_fits_all.txt'
# data_address = '/Users/jp643/Documents/Studies/PhD/Antinu/param_fitting/likelihoods/replicateTony/param_fits_all.txt'

min_E = 0.9
max_E = 8.0

def read_data():
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

    def fmt(x):
        nSig = inv_ll_diff_per_nSig(x)
        return ('%i' % round(nSig)) + r'$\sigma$'
    
    # temp change
    matplotlib.rcParams['figure.figsize'] = 14, 10

    # Background heat map
    im = plt.imshow(np.sqrt(minLL), cmap=plt.cm.viridis, interpolation='none', extent=[theta_12[0], theta_12[-1], Dm21[0], Dm21[-1]], aspect='auto', origin="lower")
    colbar = plt.colorbar(im)
    colbar.ax.set_ylabel(r'$\sqrt{-2\Delta \mathrm{log}(L)}$', rotation=90, fontproperties=prop_font)

    # Best fit point
    plt.scatter(min_theta_12, min_Dm21, marker='+', s=100, linewidths=2, color='w')

    # Contours
    X, Y = np.meshgrid(theta_12, Dm21)
    contour_vals = [1, 2, 3, 4, 5]
    contour_cols = ['peachpuff', 'sandybrown', 'darkorange', 'orangered', 'red']
    # contour_styles = ['-', '--', ':', ':', ':', ':']
    contour_styles = ['-', '--', '-.', '-.', '-.', '-.']
    contour_labels = [r'$1 \sigma$ contour', r'$2 \sigma$ contour', r'$3 \sigma$ contour', r'$4 \sigma$ contour', r'$5 \sigma$ contour']

    # plot = ax.contour(X, Y, np.sqrt(minLL), levels=contour_vals, colors=contour_cols, linestyles=contour_styles)
    # ax.clabel(plot, inline=True, fmt=fmt, fontsize=10)

    Chi2_contours = []
    for i in range(len(contour_vals)):
        Chi2_contours.append(ll_diff_per_nSig(contour_vals[i]))
    plot = plt.contour(X, Y, minLL, levels=Chi2_contours, colors=contour_cols, linestyles=contour_styles, linewidths=2)
    # plt.clabel(plot, inline=True, fmt=fmt, fontsize=18)

    # Legend
    legend_elems = [Line2D([0], [0],  color='w', lw=0, label='Best fit', marker='+', linestyle=None, markersize=8)]
    for i in range(len(contour_vals)):
        legend_elems.append(Line2D([0], [0],  color=contour_cols[i], linestyle=contour_styles[i], lw=1.5, label=contour_labels[i]))#, fontproperties=prop_font))
    leg = plt.legend(handles=legend_elems, loc='best', markerscale=2, prop={'size': 28})
    for i in range(1, len(leg.legend_handles)):
        leg.legend_handles[i].set_linewidth(3.0)
    leg.get_frame().set_color('cornflowerblue')
    leg.get_frame().set_alpha(1)

    # Labels
    plt.xlabel(r'$\theta_{12}$ $\left[{}^{\degree}\right]$', fontproperties=prop_font)
    plt.ylabel(r'$\Delta m_{12}^2$ $\left[10^{-5} \mathrm{eV}^2\right]$', fontproperties=prop_font)
    plt.title(r'Reactor-$\nu$ Sensitivity, from Azimov Dataset (134.4 days)', fontproperties=prop_font)
    plt.text(1, Dm21[5], "SNO+ Preliminary", fontproperties=prop_font)
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

    plt.stackplot(bin_edges[:-1], spectra['Reactor::model_Esys'], spectra['geoNu::model_Esys'], spectra['alphaN::model_Esys'],
                 step='post', labels=[r'reactor-$\nu$', r'geo-$\nu$', r'($\alpha$,n)'], colors=['red', 'green', 'blue'], linewidth=0)

    # plt.stairs(spectra['Total_Model_Spectrum'], edges=bin_edges, color='k', linestyle='dashed', label='total model')
    plt.scatter(Ebin_centers, spectra['data'], marker='+', linewidths=2.5, s=180, color='k', label='Azimov data')

    plt.xlabel('E [MeV]', fontproperties=prop_font)
    plt.ylabel('arbitrary', fontproperties=prop_font)
    plt.xlim(min_E, max_E)
    plt.ylim(0, 2.8)
    plt.legend(loc='best', markerscale=1, prop={'size': 28})
    plt.title(r'Reactor-$\nu$ IBD Fit Spectrum, from Azimov Dataset', fontproperties=prop_font)
    plt.text(3, 2, "SNO+ Preliminary", fontproperties=prop_font)
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

Dm21, theta_12, minLL, Ebin_centers, spectra = read_data()
# minLL = testReduceLivetime(minLL)

# fig, ax = plt.subplots(ncols=2)
plot_LL(Dm21, theta_12, minLL)
plot_spectra(Ebin_centers, spectra)
plt.show()