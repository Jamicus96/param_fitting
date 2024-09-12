import numpy as np
import matplotlib.pyplot as plt

### TEMP ###

def plot_acc_v_ibd(Ebin_centers, spectra):

    dE = Ebin_centers[1] - Ebin_centers[0]
    bin_edges = np.zeros(len(Ebin_centers)+1)
    bin_edges[:-1] = Ebin_centers - 0.5*dE
    bin_edges[-1] = Ebin_centers[-1] + 0.5*dE

    data_bin_width = 0.1

    ibd_spec = spectra['Reactor::model_Esys']
    acc_spec = spectra['Accidentals']

    # uncorr_rate = 8.50  # Hz
    livetime = 134.5  # days
    # uncorr_norm = uncorr_rate * livetime * 24 * 3600
    uncorr_norm = 98608373.0

    uncorr_spec = acc_spec * uncorr_norm / np.sum(acc_spec)

    tag_eff_ibd = 0.78  # from geoNu
    ibd_spec /= tag_eff_ibd

    uncorr_per_year = np.sum(uncorr_spec) * 365.25 / livetime
    print('uncorr rate: {} Hz, {} per day, {} per year, {} per 134.5 days'.format(uncorr_per_year / (365.25 * 24 * 3600), uncorr_per_year / 365.25, uncorr_per_year, np.sum(uncorr_spec)))
    ibd_per_year = np.sum(ibd_spec) * 365.25 / livetime
    print('ibd rate: {} Hz, {} per day, {} per year, {} per 134.5 days'.format(ibd_per_year / (365.25 * 24 * 3600), ibd_per_year / 365.25, ibd_per_year, np.sum(ibd_spec)))
    
    fig = plt.subplot(111)
    ax1 = plt.gca()

    ax1.stairs(uncorr_spec / dE, bin_edges, color='black', label='uncorrelated')
    ax1.stairs(ibd_spec / dE, bin_edges, color='crimson', label='prompt IBD')

    ax1.plot(Ebin_centers, np.sum(ibd_spec) * scipy.stats.norm.pdf(Ebin_centers, loc=2.15, scale=0.1), color='blue', label='delayed IBD')

    ax1.vlines([0.9, 8.], 1E-3, 5E9, linestyles='dashed', colors='red', label='prompt cuts')
    ax1.vlines([1.85, 2.4], 1E-3, 5E9, linestyles='dashed', colors='royalblue', label='delayed cuts')

    ax1.set_yscale('log')
    ax1.set_ylim(1E-3, 5E9)
    ax1.set_xlim(min_E, max_E)
    ax1.legend(loc='best')
    ax1.set_xlabel('Reconstructed Energy [MeV]')
    ax1.set_ylabel(r'Events [MeV${}^{-1}$]')

    plt.show()

# plot_acc_v_ibd(Ebin_centers, spectra)