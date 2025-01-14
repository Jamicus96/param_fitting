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
                        default='/mnt/lustre/scratch/epp/jp643/antinu/MC_data/Geoibd_URun_709/', help='Folder where raw ntuples are saved.')
    parser.add_argument('--accidentals_ntuple_repo', '-ar', type=str, dest='accidentals_ntuple_repo',
                        default='/mnt/lustre/scratch/epp/jp643/antinu/param_fitting/replicateTony/accidentals_ntuples/', help='Folder where accidentals ntuples are saved.')
    parser.add_argument('--cut_ntuple_repo', '-cnr', type=str, dest='cut_ntuple_repo',
                        default='/mnt/lustre/scratch/epp/jp643/antinu/param_fitting/replicateTony/MC_cut_ntuples_efficiency/geoNu_U/', help='Folder where cut ntuples are saved.')
    parser.add_argument('--scaled_ntuple_repo', '-snr', type=str, dest='scaled_ntuple_repo',
                        default='/mnt/lustre/scratch/epp/jp643/antinu/param_fitting/thesis/MC_cut_ntuples/', help='Folder where re-scaled cut reactor IBD ntuples are saved.')
    parser.add_argument('--pdf_repo', '-pr', type=str, dest='pdf_repo',
                        default='/mnt/lustre/scratch/epp/jp643/antinu/param_fitting/thesis/PDFs/', help='Folder where param fitting results are saved (2d root hist).')
    parser.add_argument('--fit_repo', '-fr', type=str, dest='fit_repo',
                        default='/mnt/lustre/scratch/epp/jp643/antinu/param_fitting/thesis/data_fitting_class_varPR/', help='Folder to save recombined root files with tracking information in.')
    
    parser.add_argument('--lt_file', '-ltf', type=str, dest='lt_file',
                        default='/mnt/lustre/scratch/epp/jp643/antinu/param_fitting/replicateTony/livetime/antinuTest/livetime_per_run_antinuTest.txt',
                        help='Text file of run livetimes. Output of livetime calculator (rat-tools). Empty string or False/false means it will not be used.')
    parser.add_argument('--rl_file', '-rlf', type=str, dest='rl_file',
                        default='/mnt/lustre/scratch/epp/jp643/antinu/param_fitting/thesis/antinu_runlist_UPDATED.txt', help='Text file of runs to include (one run-number per line).')
    parser.add_argument('---use_rl', '-uRL', type=bool, dest='use_rl',
                        default=True, help='Bool to restrict only to file names including run-numbers from the run list defined in rl_file.')
    
    parser.add_argument('--Dm21_min', type=float, dest='Dm21_min', default=0.1E-5, help='Dm_21^2 minimum.')
    parser.add_argument('--Dm21_max', type=float, dest='Dm21_max', default=15.E-5, help='Dm_21^2 maximum.')
    parser.add_argument('--s12_2_min', type=float, dest='s12_2_min', default=0., help='theta_12 minimum.')
    parser.add_argument('--s12_2_max', type=float, dest='s12_2_max', default=1., help='theta_12 maximum.')

    parser.add_argument('--PDFbinWidth', '-bw', type=float, dest='PDFbinWidth', default=0.05, help='PDF bin width [MeV]')
    parser.add_argument('--classCut', '-cc', type=float, dest='classCut', default=-9999.0, help='Classifier cut (remove events below this)')
    parser.add_argument('--Nbins', '-N', type=int, dest='Nbins', default=500, help='Number of bins in x and y directions.')
    parser.add_argument('--bins_per_job', '-mb', type=int, dest='bins_per_job', default=25, help='Maximum number of bins looped over in one job.')

    parser.add_argument('---is_data', '-iD', type=bool, dest='is_data',
                        default=False, help='For energy correction: True for data, False for MC.')
    parser.add_argument('---evt_type', '-et', type=str, dest='evt_type', default='IBD',
                        choices=['any', 'IBD', 'alphaN'], help='For cut efficiency calculation (MC only).')
    parser.add_argument('---use_Azimov', '-uA', type=bool, dest='use_Azimov',
                        default=False, help='Use Azimov data set made from PDFs, instead of ntuple data from cut_ntuple_repo')

    parser.add_argument('--max_jobs', '-m', type=int, dest='max_jobs',
                        default=1000, help='Max number of tasks in an array running at any one time.')
    parser.add_argument('---step', '-s', type=str, dest='step', default='all', choices=['accidentals', 'cut', 'scale', 'pdf', 'fit', 'combi'],
                        help='which step of the process is it in?')
    parser.add_argument('---verbose', '-v', type=bool, dest='verbose',
                        default=True, help='print and save extra info')

    args = parser.parse_args()
    return args


### Miscelaneous functions ###

def getRepoAddress():
    '''Returns the full address of the git repo containing with script'''
    repo_address = __file__[:-len('runParamFitting.py')]
    if repo_address == '':
        first_char = None
    else:
        first_char = repo_address[0]
    if first_char != '/':
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


### Cutting (dat & PDF) functions ###

def selectNutples(ntuple_repo, args):

    run_numbers_str = []
    unused_runs = []
    if args.use_rl:
        f = open(args.rl_file, 'r')
        lines = f.readlines()
        f.close()
        for line in lines:
            line_strp = line.rstrip('\n')
            if line_strp != '':
                run_numbers_str.append(line_strp)
                unused_runs.append(int(line_strp))

    file_addresses = []
    runs = []
    for filename in os.listdir(ntuple_repo):
        file_address = os.path.join(ntuple_repo, filename)
        if os.path.isfile(file_address):
            if file_address[-12:] == '.ntuple.root':
                addFile = True
                run = 0
                if args.use_rl:
                    addFile = False
                    for runNum in run_numbers_str:
                        if runNum in filename:
                            addFile = True
                            if int(runNum) in unused_runs:
                                unused_runs.remove(int(runNum))
                                run = int(runNum)
                            break
                if addFile:
                    file_addresses.append(file_address)
                    runs.append(run)
    
    if args.use_rl:
        # re-order
        file_addresses = [x for _, x in sorted(zip(runs, file_addresses))]

        print('Missing {} out of {} runs'.format(len(unused_runs), len(run_numbers_str)))
        print('Missing runs:', np.array(unused_runs))

    return file_addresses

def cutNtuplesJobScript(example_jobScript, overall_folder, commandList_address, verbose):
    '''Create job script to run array of rat macros'''

    new_job_address = overall_folder + 'job_scripts/'
    new_job_address = checkRepo(new_job_address, verbose)
    new_job_address += 'cutting_job.job'

    output_logFile_address = overall_folder + 'log_files/'
    output_logFile_address = checkRepo(output_logFile_address, verbose)
    # output_logFile_address +=  'log_cutting.txt'  # remove to make each subjob have its own log file

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

def cutNtuples(args):
    ''''''

    # Check repos are formatted correctly and exist, otherwise make them
    ntuple_repo = checkRepo(args.ntuple_repo, args.verbose)
    cut_ntuple_repo = checkRepo(args.cut_ntuple_repo, args.verbose)
    job_script_repo = checkRepo(cut_ntuple_repo + 'job_scripts/', args.verbose)
    commandList_address = job_script_repo + 'commandList.txt'

    # Get full path to this repo
    path = getRepoAddress()

    # Get ntuple file addresses
    ntuple_addresses = selectNutples(ntuple_repo, args)

    # Create command
    commandBase = path + 'cutting/cut_data.exe '
    
    commandList_file = open(commandList_address, 'w')
    for i, file_address in enumerate(ntuple_addresses):
        outNtuple_address = 'CUT_' + file_address.split('/')[-1][:-12] + '.ntuple.root'
        if i != 0 and args.use_rl:
            command = commandBase + file_address + ' ' + ntuple_addresses[i-1] + ' '
        else:
            command = commandBase + file_address + ' 0 '
        command += cut_ntuple_repo + outNtuple_address
        command += ' ' + str(args.classCut) + ' ' + str(int(args.is_data)) + ' ' + args.evt_type
        commandList_file.write(command + '\n')
    commandList_file.close()

    # Create new job file
    jobScript_address = path + 'job_scripts/fitting_job.job'
    with open(jobScript_address, "r") as f:
        example_jobScript = f.readlines()
    job_address = cutNtuplesJobScript(example_jobScript, cut_ntuple_repo, commandList_address, args.verbose)

    # Run job script!
    print('Submitting job...')
    command = 'qsub -t 1-' + str(len(ntuple_addresses)) + ' -tc ' + str(args.max_jobs) + ' ' + job_address
    if args.verbose:
        print('Running command: ', command)
    subprocess.call(command, stdout=subprocess.PIPE, shell=True) # use subprocess to make code wait until it has finished


def get_accidentals_and_Vetos(args):
    ''''''

    # Check repos are formatted correctly and exist, otherwise make them
    ntuple_repo = checkRepo(args.ntuple_repo, args.verbose)
    accidentals_ntuple_repo = checkRepo(args.accidentals_ntuple_repo, args.verbose)
    job_repo = accidentals_ntuple_repo + 'job_scripts/'
    job_repo = checkRepo(job_repo, args.verbose)
    commandList_address = job_repo + 'commandList.txt'

    # Get full path to this repo
    path = getRepoAddress()

    # Get ntuple file addresses
    ntuple_addresses = selectNutples(ntuple_repo, args)

    # Create command
    commandBase = path + 'cutting/get_accidentals_and_vetos.exe '
    
    commandList_file = open(commandList_address, 'w')
    for i, file_address in enumerate(ntuple_addresses):
        outNtuple_address = 'CUT_' + file_address.split('/')[-1][:-12] + '.ntuple.root'
        outTxt_address = 'Veto_' + file_address.split('/')[-1][:-12] + '_vetooutput.txt'  # livetime calculator uses `*vetooutput.txt` format
        if i != 0 and args.use_rl:
            command = commandBase + file_address + ' ' + ntuple_addresses[i-1] + ' '
        else:
            command = commandBase + file_address + ' 0 '
        command += accidentals_ntuple_repo + outNtuple_address + ' ' + accidentals_ntuple_repo + outTxt_address
        command += ' ' + str(args.classCut) + ' ' + str(int(args.is_data))
        commandList_file.write(command + '\n')
    commandList_file.close()

    # Create new job file
    jobScript_address = path + 'job_scripts/fitting_job.job'
    with open(jobScript_address, "r") as f:
        example_jobScript = f.readlines()
    job_address = cutNtuplesJobScript(example_jobScript, accidentals_ntuple_repo, commandList_address, args.verbose)

    # Run job script!
    print('Submitting job...')
    command = 'qsub -l m_mem_free=4G -t 1-' + str(len(ntuple_addresses)) + ' -tc ' + str(args.max_jobs) + ' ' + job_address
    if args.verbose:
        print('Running command: ', command)
    subprocess.call(command, stdout=subprocess.PIPE, shell=True) # use subprocess to make code wait until it has finished


def makePDFsjobScript(example_jobScript, overall_folder, commands, args):
    '''Create job script to run array of rat macros'''

    new_job_address = overall_folder + 'job_scripts/'
    new_job_address = checkRepo(new_job_address, args.verbose)
    new_job_address += 'PDF_job.job'

    output_logFile_address = overall_folder + 'log_files/'
    output_logFile_address = checkRepo(output_logFile_address, args.verbose)
    output_logFile_address +=  'log_PDF_' + str(args.PDFbinWidth) + '_cut' + str(args.classCut) + '.txt'

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
    cut_ntuple_repo = checkRepo(args.cut_ntuple_repo, args.verbose)
    accidentals_ntuple_repo = checkRepo(args.accidentals_ntuple_repo, args.verbose)
    pdf_repo = checkRepo(args.pdf_repo, args.verbose)

    # Get full path to this repo
    path = getRepoAddress()

    # Create command
    command = path + 'cutting/make_PDFs.exe '
    command += cut_ntuple_repo + 'scaled_CUT_reactorIBD.ntuple.root '
    command += cut_ntuple_repo + 'CUT_alphaN.ntuple.root '
    command += cut_ntuple_repo + 'CUT_geoNu_Th.ntuple.root '
    command += cut_ntuple_repo + 'CUT_geoNu_U.ntuple.root '
    command += accidentals_ntuple_repo + 'CUT_accidentals.ntuple.root '
    command += pdf_repo + 'PDFs_' + str(args.PDFbinWidth) + '_cut' + str(args.classCut) + '.root '
    command += str(args.PDFbinWidth) + ' ' + str(args.classCut)

    # Make job script
    jobScript_address = path + 'job_scripts/PDF_job.job'
    with open(jobScript_address, "r") as f:
        example_jobScript = f.readlines()
    job_address = makePDFsjobScript(example_jobScript, pdf_repo, command, args)

    # Run job script!
    print('Submitting job...')
    command = 'qsub -l m_mem_free=4G ' + job_address
    if args.verbose:
        print('Running command: ', command)
    subprocess.call(command, stdout=subprocess.PIPE, shell=True) # use subprocess to make code wait until it has finished


def scale_reactorNtuples(args):
    '''Re-scale reactor IBD ntuples to true reactor power'''

    # Check repos are formatted correctly and exist, otherwise make them
    cut_ntuple_repo = checkRepo(args.cut_ntuple_repo, args.verbose)
    scaled_cut_ntuple_repo = checkRepo(args.scaled_ntuple_repo, args.verbose)
    commandList_address = scaled_cut_ntuple_repo + 'commandList.txt'

    # Get full path to this repo
    path = getRepoAddress()

    # Get ntuple file addresses
    ntuple_addresses = selectNutples(cut_ntuple_repo, args)

    # Create command
    commandBase = path + 'cutting/reScaleReactorIBD.exe '
    
    commandList_file = open(commandList_address, 'w')
    for file_address in ntuple_addresses:
        outNtuple_address = 'scaled_' + file_address.split('/')[-1][:-12] + '.ntuple.root'
        command = commandBase + file_address + ' ' + args.lt_file + ' '
        command += scaled_cut_ntuple_repo + outNtuple_address
        commandList_file.write(command + '\n')
    commandList_file.close()

    # Create new job file
    jobScript_address = path + 'job_scripts/fitting_job.job'
    with open(jobScript_address, "r") as f:
        example_jobScript = f.readlines()
    job_address = cutNtuplesJobScript(example_jobScript, scaled_cut_ntuple_repo, commandList_address, args.verbose)

    # Run job script!
    print('Submitting job...')
    command = 'qsub -t 1-' + str(len(ntuple_addresses)) + ' ' + job_address
    if args.verbose:
        print('Running command: ', command)
    subprocess.call(command, stdout=subprocess.PIPE, shell=True) # use subprocess to make code wait until it has finished

### Fitting functions ###

def get_param_lims_per_job(Dm21_min, Dm21_max, s12_2_min, s12_2_max, Nbins, bins_per_job):
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
    theta12_step = (s12_2_max - s12_2_min) / float(Nbins - 1)

    Dm21_mins = Dm21_min + idx * bins_per_job * Dm21_step
    s12_2_mins = s12_2_min + idx * bins_per_job * theta12_step

    Dm21_maxs = Dm21_min + ((idx + 1.) * bins_per_job - 1.) * Dm21_step
    Dm21_maxs[-1] = Dm21_max
    s12_2_maxs = s12_2_min + ((idx + 1.) * bins_per_job - 1.) * theta12_step
    s12_2_maxs[-1] = s12_2_max

    # Work out lower bin number edges of these sections
    Dm21_indices = np.arange(0, Nbins, bins_per_job)
    theta12_indices = Dm21_indices

    # Work out overall 2d histogram edges (low edge of min bin and upper edge of max bin)
    Dm21_lower = Dm21_min - 0.5 * Dm21_step
    theta12_lower = s12_2_min - 0.5 * theta12_step
    Dm21_upper = Dm21_max + 0.5 * Dm21_step
    theta12_upper = s12_2_max + 0.5 * theta12_step
    hist_lims = str(Dm21_lower) + ' ' + str(Dm21_upper) + ' ' + str(theta12_lower) + ' ' + str(theta12_upper)

    # Combine x an y axis section together, to "tile" the whole 2d phase space
    job_lims = []
    for i in range(idx.size):
        for j in range(idx.size):
            job_lims.append(str(Dm21_mins[i]) + ' ' + str(Dm21_maxs[i]) + ' ' + str(s12_2_mins[j]) + ' ' + str(s12_2_maxs[j])\
                             + ' ' + str(Dm21_indices[i]) + ' ' + str(theta12_indices[j]))

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
    ntuple_repo = checkRepo(args.cut_ntuple_repo, args.verbose)

    # Get full path to this repo
    path = getRepoAddress()

    # How to split up phase space into manageable jobs
    job_lims, hist_lims, nBins_job = get_param_lims_per_job(args.Dm21_min, args.Dm21_max, args.s12_2_min, args.s12_2_max, args.Nbins, args.bins_per_job)

    ### MAKE JOB SCRIPTS ###
    print('Creating split hist job scripts...')
    command_base = path + 'fitting/fit_params.exe ' + pdf_repo + 'PDFs_' + str(args.PDFbinWidth) + '_cut' + str(args.classCut) + '.root ' + ntuple_repo + 'CUT_data.ntuple.root ' + fit_repo + 'param_fits_'

    # Make list of commands for job array to call
    jobScript_repo = fit_repo + 'job_scripts/'
    checkRepo(jobScript_repo, args.verbose)
    commandList_address = jobScript_repo + 'fitting_commands.txt'
    commandList_file = open(commandList_address, 'w')
    for i, job_lim in enumerate(job_lims):
        # Create all the histogram making commands
        command = command_base + str(i) + '.root ' + str(int(args.use_Azimov)) + ' ' + hist_lims + ' ' + str(args.Nbins) + ' '\
                               + job_lim  + ' ' + str(nBins_job[i, 0]) + ' ' + str(nBins_job[i, 1]) + ' ' + str(int(args.verbose))
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


def makeCombiJobScript(example_jobScript, overall_folder, commands, verbose):
    '''Create job script to run array of rat macros'''

    new_job_address = overall_folder + 'job_scripts/'
    new_job_address = checkRepo(new_job_address, verbose)
    new_job_address += 'combi_job.job'

    output_logFile_address = overall_folder + 'log_files/'
    output_logFile_address = checkRepo(output_logFile_address, verbose)
    output_logFile_address +=  'log_combi.txt'

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

def combi_fits(args):

    # Check fitting repo
    pdf_repo = checkRepo(args.pdf_repo, args.verbose)
    fit_repo = checkRepo(args.fit_repo, args.verbose)
    ntuple_repo = checkRepo(args.cut_ntuple_repo, args.verbose)

    # Get full path to this repo
    path = getRepoAddress()

    # Work out lower bin number edges of these sections
    Dm21_start_indices = np.arange(0, args.Nbins, args.bins_per_job, dtype=int)
    Dm21_end_indices = np.zeros(Dm21_start_indices.size, dtype=int)
    Dm21_end_indices[:-1] = Dm21_start_indices[1:] - 1
    Dm21_end_indices[-1] = args.Nbins - 1

    theta12_start_indices = Dm21_start_indices
    theta12_end_indices = Dm21_end_indices

    # make job command
    command = path + 'fitting/re_combine_fits.exe ' + pdf_repo + 'PDFs_' + str(args.PDFbinWidth) + '_cut' + str(args.classCut) + '.root ' + ntuple_repo + 'CUT_data.ntuple.root ' + fit_repo + 'param_fits_all.root ' + str(int(args.use_Azimov))
    k = 0
    for i in range(Dm21_start_indices.size):
        for j in range(theta12_start_indices.size):
            command += ' ' + fit_repo + 'param_fits_' + str(k) + '.root '
            command += str(Dm21_start_indices[i]) + ' ' + str(Dm21_end_indices[i]) + ' '
            command += str(theta12_start_indices[j]) + ' ' + str(theta12_end_indices[j])
            k += 1

    # Make job script
    jobScript_address = path + 'job_scripts/PDF_job.job'
    with open(jobScript_address, "r") as f:
        example_jobScript = f.readlines()
    job_address = makeCombiJobScript(example_jobScript, fit_repo, command, args.verbose)

    # Run job script!
    command = 'qsub -l m_mem_free=4G ' + job_address
    if args.verbose:
        print('Running command: ', command)
    subprocess.call(command, stdout=subprocess.PIPE, shell=True) # use subprocess to make code wait until it has finished
    

### MAIN ###

def main():
    # read in argument
    args = argparser()

    if args.step == 'accidentals':
        get_accidentals_and_Vetos(args)
    elif args.step == 'cut':
        cutNtuples(args)
    elif args.step == 'scale':
        scale_reactorNtuples(args)
    elif args.step == 'pdf':
        makePDFs(args)
    elif args.step == 'fit':
        doFitting(args)
    elif args.step == 'combi':
        combi_fits(args)
    else:
        print('WRONG work mode ({}). Exiting...'.format(args.step))


if __name__ == '__main__':
    main()