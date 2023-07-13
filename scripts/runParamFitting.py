import os
import numpy as np
import argparse
import subprocess
import json
import time

# If need to rerun jobs from an array, remember that files are labelled with integer one lower than
# the array number. For example, `qsub -t 3-3 jobScript.job` will give `outputFile_2.root`.

def argparser():
    '''Define arguments'''

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='Run parameter fitting code')

    parser.add_argument('--ntuple_repo', '-nr', type=str, dest='ntuple_repo',
                        default='/mnt/lustre/scratch/epp/jp643/antinu/param_fitting/2p2gLppo/ntuples/', help='Folder where PDFs (root hists) are saved.')
    parser.add_argument('--pdf_repo', '-pr', type=str, dest='pdf_repo',
                        default='/mnt/lustre/scratch/epp/jp643/antinu/param_fitting/2p2gLppo/PDFs/', help='Folder where param fitting results are saved (2d root hist).')
    parser.add_argument('--fit_repo', '-fr', type=str, dest='fit_repo',
                        default='/mnt/lustre/scratch/epp/jp643/antinu/param_fitting/2p2gLppo/likelihoods/', help='Folder to save recombined root files with tracking information in.')

    parser.add_argument('--Dm21_min', type=float, dest='Dm21_min', default=1E-5, help='Dm_21^2 minimum.')
    parser.add_argument('--Dm21_max', type=float, dest='Dm21_max', default=10E-5, help='Dm_21^2 maximum.')
    parser.add_argument('--theta12_min', type=float, dest='theta12_min', default=5., help='theta_12 minimum.')
    parser.add_argument('--theta12_max', type=float, dest='theta12_max', default=45., help='theta_12 maximum.')

    parser.add_argument('--classCut', '-cc', type=float, dest='classCut', default=-9999., help='Classifier cut (remove events below this)')
    parser.add_argument('--Nbins', '-N', type=int, dest='Nbins', default=200, help='Number of bins in x and y directions.')
    parser.add_argument('--bins_per_job', '-mb', type=int, dest='bins_per_job', default=200, help='Maximum number of bins looped over in one job.')

    parser.add_argument('--max_jobs', '-m', type=int, dest='max_jobs',
                        default=70, help='Max number of tasks in an array running at any one time.')
    parser.add_argument('---step', '-s', type=str, dest='step', default='all', choices=['pdf', 'fit'],
                        help='which step of the process is it in?')
    parser.add_argument('---verbose', '-v', type=bool, dest='verbose',
                        default=True, help='print and save extra info')

    args = parser.parse_args()
    return args


### Miscelaneous functions ###

def getRepoAddress():
    '''Returns the full address of the git repo containing with script'''
    repo_address = __file__[:-len('scripts/runParamFitting.py')]
    if repo_address == '':
        firt_char = None
    else:
        firt_char = repo_address[0]
    if firt_char != '/':
        repo_address = os.getcwd() + '/' + repo_address
    return repo_address

def checkRepo(repo_address, verbose=False):
    '''Check format of repo address is right to be used here. Also if repo does not exist, create it.'''

    # Format address
    new_address = repo_address
    if new_address == '.':
        new_address = ''
    elif new_address != '':
        if (new_address[int(len(new_address) - 1)] != '/'):
            new_address += '/'
        
        # If directory does not exists, make it
        if not os.path.isdir(new_address):
            os.mkdir(new_address)
            if verbose:
                print('Created new directory: ', new_address)

    return new_address


### PDF functions ###

def makePDFjobScript(example_jobScript, overall_folder, commands, verbose):
    '''Create job script to run array of rat macros'''

    new_job_address = overall_folder + 'job_scripts/'
    new_job_address = checkRepo(new_job_address, verbose)
    new_job_address += 'PDF_job.job'

    output_logFile_address = overall_folder + 'log_files/'
    output_logFile_address = checkRepo(output_logFile_address, verbose)
    output_logFile_address +=  'log_PDF.txt'

    new_jobScript = []
    for line in example_jobScript:
        # Replace placeholders in macro
        if 'output_log.txt' in line:
            new_line = line.replace('output_log.txt', output_logFile_address, 1)
        elif 'your commands' in line:
            new_line = line.replace('your commands', commands, 1)
        else:
            new_line = line

        new_jobScript.append(new_line)

    # Create job file
    with open(new_job_address, "w") as f:
        new_jobScript = "".join(new_jobScript)
        f.write(new_jobScript)

    return new_job_address

def makePDFs(args):
    '''Turn reactor IBD and alpha-n event ntuples into hist PDFs to be used for fitting'''

    # Check repos are formatted correctly and exist, otherwise make them
    ntuple_repo = checkRepo(args.ntuple_repo, args.verbose)
    pdf_repo = checkRepo(args.pdf_repo, args.verbose)

    # Get full path to this repo
    path = getRepoAddress()

    # Create command
    command = path + 'scripts/make_PDFs.exe '
    command += ntuple_repo + 'reactorIBD.ntuple.root '
    command += ntuple_repo + 'alphaN.ntuple.root '
    command += pdf_repo + 'PDFs.root '
    command += str(args.classCut)

    # Make job script
    jobScript_address = path + 'job_scripts/PDF_job.job'
    with open(jobScript_address, "r") as f:
        example_jobScript = f.readlines()
    job_address = makePDFjobScript(example_jobScript, pdf_repo, command, args.verbose)

    # Run job script!
    print('Submitting job...')
    command = 'qsub ' + job_address
    if args.verbose:
        print('Running command: ', command)
    subprocess.call(command, stdout=subprocess.PIPE, shell=True) # use subprocess to make code wait until it has finished


### Fitting functions ###

def get_param_lims_per_job(Dm21_min, Dm21_max, theta12_min, theta12_max, Nbins, bins_per_job):
    '''Create array with param lims to pass to code in each job, and the hist lims'''

    # Work out number of sections (len(idx)) the binning is split into
    n_div = Nbins // bins_per_job
    remainder = Nbins % bins_per_job

    if remainder == 0:
        idx = np.arange(0, n_div, 1, dtype=float)
    else:
        idx = np.arange(0, n_div + 1, 1, dtype=float)

    # Work out upper and lower limits of each of these section (centers of min and max bins)
    Dm21_step = (Dm21_max - Dm21_min) / float(Nbins - 1)
    theta12_step = (theta12_max - theta12_min) / float(Nbins - 1)

    Dm21_mins = Dm21_min + idx * bins_per_job * Dm21_step
    theta12_mins = theta12_min + idx * bins_per_job * theta12_step

    Dm21_maxs = Dm21_min + ((idx + 1.) * bins_per_job - 1.) * Dm21_step
    Dm21_maxs[-1] = Dm21_max
    theta12_maxs = theta12_min + ((idx + 1.) * bins_per_job - 1.) * theta12_step
    theta12_maxs[-1] = theta12_max

    # Work out overall 2d histogram edges (low edge of min bin and upper edge of max bin)
    Dm21_lower = Dm21_min - 0.5 * Dm21_step
    theta12_lower = theta12_min - 0.5 * theta12_step
    Dm21_upper = Dm21_max + 0.5 * Dm21_step
    theta12_upper = theta12_max + 0.5 * theta12_step
    hist_lims = str(Dm21_lower) + ' ' + str(Dm21_upper) + ' ' + str(theta12_lower) + ' ' + str(theta12_upper)

    # Combine x an y axis section together, to "tile" the whole 2d phase space
    job_lims = []
    for i in range(idx.size):
        for j in range(idx.size):
            job_lims.append(str(Dm21_mins[i]) + ' ' + str(Dm21_maxs[i]) + ' ' + str(theta12_mins[j]) + ' ' + str(theta12_maxs[j]))

    # Get number of bins for each job
    nBins_job = np.full((idx.size * idx.size, 2), bins_per_job)
    if remainder != 0:
        for i in range(1, idx.size + 1):
            nBins_job[i * idx.size - 1, 1] = remainder
        for i in range(idx.size * (idx.size - 1), idx.size * idx.size):
            nBins_job[i, 0] = remainder

    return job_lims, hist_lims, nBins_job

def makeJobArrayScript(example_jobScript, overall_folder, commandList_address, verbose):
    '''Create job script to run array of rat macros'''

    new_job_address = overall_folder + 'job_scripts/'
    new_job_address = checkRepo(new_job_address, verbose)
    new_job_address += 'fitting_job.job'

    output_logFile_address = overall_folder + 'log_files/'
    output_logFile_address = checkRepo(output_logFile_address, verbose)
    output_logFile_address +=  'log_fitting.txt'

    new_jobScript = []
    for line in example_jobScript:
        # Replace placeholders in macro
        if 'output_log.txt' in line:
            new_line = line.replace('output_log.txt', output_logFile_address, 1)
        elif '/Address/CommandList.txt' in line:
            new_line = line.replace('/Address/CommandList.txt', commandList_address, 1)
        else:
            new_line = line

        new_jobScript.append(new_line)

    # Create job file
    with open(new_job_address, "w") as f:
        new_jobScript = "".join(new_jobScript)
        f.write(new_jobScript)

    return new_job_address

def doFitting(args):
    '''Make 2d histogram to find best fit Dm_21^2 and theta_12 values, varying other parameters'''

    # Check repos are formatted correctly and exist, otherwise make them
    pdf_repo = checkRepo(args.pdf_repo, args.verbose)
    fit_repo = checkRepo(args.fit_repo, args.verbose)

    # Get full path to this repo
    path = getRepoAddress()

    # How to split up phase space into manageable jobs
    job_lims, hist_lims, nBins_job = get_param_lims_per_job(args.Dm21_min, args.Dm21_max, args.theta12_min, args.theta12_max, args.Nbins, args.bins_per_job)

    ### MAKE JOB SCRIPTS ###
    print('Creating split hist job scripts...')
    command_base = path + 'scripts/fit_params.exe ' + pdf_repo + 'PDFs.root ' + fit_repo + 'param_fits.root '

    # Make list of commands for job array to call
    jobScript_repo = fit_repo + 'job_scripts/'
    checkRepo(jobScript_repo, args.verbose)
    commandList_address = jobScript_repo + 'fitting_commands.txt'
    commandList_file = open(commandList_address, 'w')
    for i, job_lim in enumerate(job_lims):
        # Create all the histogram making commands
        command = command_base + hist_lims + ' ' + str(args.Nbins) + ' ' + job_lim  + ' ' + str(nBins_job[i, 0]) + ' ' + str(nBins_job[i, 1])
        commandList_file.write(command + '\n')
    commandList_file.close()

    # Create the job script to run all these macros in an array
    jobScript_address = path + 'job_scripts/fitting_job.job'
    with open(jobScript_address, "r") as f:
        example_jobScript = f.readlines()
    job_address = makeJobArrayScript(example_jobScript, fit_repo, commandList_address, args.verbose)

    ### RUN JOB SCRIPTS ###
    print('Submitting job array...')
    command = 'qsub -t 1-' + str(len(job_lims)) + ' -tc ' + str(args.max_jobs) + ' ' + job_address
    if args.verbose:
        print('Running command: ', command)
    subprocess.call(command, stdout=subprocess.PIPE, shell=True) # use subprocess to make code wait until it has finished


### MAIN ###

def main():
    # read in argument
    args = argparser()

    if args.step == 'pdf':
        makePDFs(args)
    elif args.step == 'fit':
        doFitting(args)
    else:
        print('WRONG work mode ({}). Exiting...'.format(args.step))


if __name__ == '__main__':
    main()