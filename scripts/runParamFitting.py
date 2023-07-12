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
                        default='/mnt/lustre/scratch/epp/jp643/antinu/param_fitting/data/', help='Folder where PDFs (root hists) are saved.')
    parser.add_argument('--pdf_repo', '-pr', type=str, dest='pdf_repo',
                        default='/mnt/lustre/scratch/epp/jp643/antinu/param_fitting/PDFs/', help='Folder where param fitting results are saved (2d root hist).')
    parser.add_argument('--fit_repo', '-fr', type=str, dest='fit_repo',
                        default='/mnt/lustre/scratch/epp/jp643/antinu/param_fitting/likelihoods/', help='Folder to save recombined root files with tracking information in.')

    parser.add_argument('--Dm21_min', type=float, dest='Dm21_min', default=1E-5, help='Dm_21^2 minimum.')
    parser.add_argument('--Dm21_max', type=float, dest='Dm21_max', default=10E-5, help='Dm_21^2 maximum.')
    parser.add_argument('--theta12_min', type=float, dest='theta12_min', default=5., help='theta_12 minimum.')
    parser.add_argument('--theta12_max', type=float, dest='theta12_max', default=45., help='theta_12 maximum.')

    parser.add_argument('--Nbins', '-N', type=int, dest='Nbins', default=200, help='Number of bins in x and y directions.')
    parser.add_argument('--max_loops', '-ml', type=int, dest='max_loops', default=200, help='Maximum number of loops in one job.')

    parser.add_argument('---step', '-s', type=str, dest='step', default='all', choices=['pdf', 'fitting'],
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

def makeJobArrayScript(jobName_str, example_jobScript, overall_folder, commandList_address, info, verbose):
    '''Create job script to run array of rat macros'''

    new_job_address = overall_folder + 'job_scripts/'
    new_job_address = checkRepo(new_job_address, verbose)
    new_job_address += job_str_map(jobName_str, info) + '.job'

    output_logFile_address = overall_folder + 'log_files/'
    output_logFile_address = checkRepo(output_logFile_address, verbose)
    output_logFile_address +=  'log_' + jobName_str + filename_format(info) + '.txt'

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

def makeJobSingleScript(jobName_str, example_jobScript, overall_folder, commands, info, verbose):
    '''Create job script to run array of rat macros'''

    new_job_address = overall_folder + 'job_scripts/'
    new_job_address = checkRepo(new_job_address, verbose)
    new_job_address += job_str_map(jobName_str, info) + '.job'

    output_logFile_address = overall_folder + 'log_files/'
    output_logFile_address = checkRepo(output_logFile_address, verbose)
    output_logFile_address +=  'log_' + jobName_str + filename_format(info) + '.txt'

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

def getNevtsPerMacro(nevts_total, nevts_persim):
    '''Create array with number of events to simulate per macro'''

    n_macros = nevts_total // nevts_persim
    remainder = nevts_total % nevts_persim
    if remainder == 0:
        n_evts = np.an_array = np.full(n_macros, nevts_persim)
    else:
        n_evts = np.an_array = np.full(n_macros + 1, nevts_persim)
        n_evts[n_macros] = remainder

    return n_evts

def checkJobsDone(jobName_substr, input_info, wait_time, verbose):
    '''Wait until submitted jobs of certain forma are finished. Wait time in seconds.'''

    # Turns out the name of the job is only the 10 first characters of the job file name

    running = True
    while running:
        running = False
        output = subprocess.Popen('qstat -u $USER', stdout=subprocess.PIPE, shell=True).communicate()[0]
        lines = output.decode("utf-8").split('\n')
        for line in lines:
            if running:
                break
            else:
                for info in input_info:
                    map_str = job_str_map(jobName_substr, info)
                    if len(map_str) > 10:
                        map_str = map_str[:9]
                    if map_str in line:
                        running = True
                        if verbose:
                            print('Waiting for jobs to finish...')
                        break
        time.sleep(wait_time)

    return True


### PDF functions ###

def getHists(args, input_info):
    '''Get histograms from split up simulation outputs, and recombine them'''
    print('Running getHists().')

    # Read in example macro and job script + info
    repo_address = getRepoAddress()

    jobArrayScript_address = repo_address + 'job_scripts/jobArray.job'
    with open(jobArrayScript_address, "r") as f:
        example_jobArrayScript = f.readlines()

    # Make sure folders are of the correct format to  use later
    save_sims_folder = checkRepo(args.sim_repo, args.verbose)
    save_splithists_folder = checkRepo(args.splithist_repo, args.verbose)

    # Folder for job scripts and command lists they use
    jobScript_repo = save_sims_folder + 'job_scripts/'
    jobScript_repo = checkRepo(jobScript_repo, args.verbose)

    # How to split up sims into manageable macros
    n_evts = getNevtsPerMacro(args.nevts_total, args.nevts_persim)

    ### Check that all simulations have finished running ###
    checkJobsDone('sims_', input_info, 10, args.verbose)

    ### MAKE JOB SCRIPTS TO CREATE HISTOGRAMS ###
    print('Creating split hist job scripts...')
    hist_command_base = repo_address + 'scripts/GetHists.exe '
    job_addresses = []
    for info in input_info:
        if args.verbose:
            print('geo_file=', info[0], ', inner_av_material=', info[1], ', LED=', info[2], ', fibre=', info[3], ', reemis=', info[4], ', abs=', info[5])
        # Make list of commands for job array to call
        commandList_address = jobScript_repo + 'hist_commandList_' + filename_format(info) + '.txt'
        commandList_file = open(commandList_address, 'w')
        wavelength = info[2][3:]
        for i in range(len(n_evts)):
            # Create all the histogram making commands
            hist_command = hist_command_base + save_sims_folder + 'simOut_' + filename_format(info) + '_' + str(i) + '.root '\
                                             + save_splithists_folder + 'splitHist_' + filename_format(info) + '_' + str(i) + '.root '\
                                             + info[3] + ' ' + wavelength + ' ' + str(int(args.verbose))
            commandList_file.write(hist_command + '\n')
        commandList_file.close()

        # Create the job script to run all these macros in an array
        new_job_address = makeJobArrayScript('splitHist_', example_jobArrayScript, save_sims_folder, commandList_address, info, args.verbose)
        job_addresses.append(new_job_address)

    ### RUN JOB SCRIPTS ###
    print('Submitting jobs...')
    for job_address in job_addresses:
        command = 'qsub -t 1-' + str(len(n_evts)) + ' -tc ' + str(args.max_jobs) + ' ' + job_address
        if args.verbose:
            print('Running command: ', command)
        subprocess.call(command, stdout=subprocess.PIPE, shell=True) # use subprocess to make code wait until it has finished

    return True

def combiHists(args, input_info):
    '''Combine split histograms together'''
    print('Running combiHists().')

    # Read in example macro and job script + info
    repo_address = getRepoAddress()

    # Make sure folders are of the correct format to  use later
    save_splithists_folder = checkRepo(args.splithist_repo, args.verbose)
    save_tothists_folder = checkRepo(args.tothist_repo, args.verbose)

    jobSingleScript_address = repo_address + 'job_scripts/jobSingle.job'
    with open(jobSingleScript_address, "r") as f:
        example_jobSingleScript = f.readlines()

    # Making job scripts
    hist_command_base = 'hadd ' + save_tothists_folder + 'tot_hists_'
    job_addresses = []
    for info in input_info:
        if args.verbose:
            print('geo_file=', info[0], ', inner_av_material=', info[1], ', LED=', info[2], ', fibre=', info[3], ', reemis=', info[4], ', abs=', info[5])
        # Make list of command for job to call
        command = hist_command_base + filename_format(info) + '.root ' + save_splithists_folder + 'splitHist_' + filename_format(info) + '_*.root'

        # Create the job script to run all these macros in an array
        new_job_address = makeJobSingleScript('tot_hists_', example_jobSingleScript, save_splithists_folder, command, info, args.verbose)
        job_addresses.append(new_job_address)

    # Wait until these job arrays are done
    checkJobsDone('splitHist_', input_info, 10, args.verbose)
    
    ### RUN JOB SCRIPTS ###
    print('Submitting jobs...')
    for job_address in job_addresses:
        command = 'qsub ' + job_address
        if args.verbose:
            print('Running command: ', command)
        subprocess.call(command, stdout=subprocess.PIPE, shell=True) # use subprocess to make code wait until it has finished

### Fitting functions ###

def makeStatsJobScript(hist_files, info, outer_lims_list, repo_address, save_stats_folder, save_tothists_folder, verbose):
    '''Create job script and associated command list file to run Analysis on array of same simulations
    (except absorption factor), with different region limits.'''

    # Folder for job scripts and command lists they use
    jobScript_repo = save_tothists_folder + 'job_scripts/'
    jobScript_repo = checkRepo(jobScript_repo, verbose)

    # For each set of absorptions, apply all the region limits
    commandList_address = jobScript_repo + 'stats_commandList_' + filename_format(info, True) + '.txt'
    commandList_file = open(commandList_address, 'w')
    output_stats_files = []
    for i in range(len(outer_lims_list)):
        outer_lims = outer_lims_list[i]
        if verbose:
            print('outer_lims=', outer_lims)

        output_stats_file = save_stats_folder + 'stats_' + filename_format(info, True) + '_' + outer_lims.replace(' ', '_') + '.txt'
        slope_command = repo_address + 'scripts/GetStats.exe ' + output_stats_file + ' ' + outer_lims + ' ' + str(int(verbose)) + ' ' + hist_files
        commandList_file.write(slope_command)
        output_stats_files.append(output_stats_file)
        if i != len(outer_lims_list) - 1:
            commandList_file.write('\n')
    commandList_file.close()

    jobArrayScript_address = repo_address + 'job_scripts/jobArray.job'
    with open(jobArrayScript_address, "r") as f:
        example_jobArrayScript = f.readlines()

    job_address = makeJobArrayScript('stats_', example_jobArrayScript, save_tothists_folder, commandList_address, info, verbose)

    return job_address, output_stats_files

def writeJsonStatsFile(output_info, region_lims_list, json_stats_folder):
    '''Write all stats to json files'''

    for output_stats_file_list, info, abs_factors in output_info:
        # Write sim info
        table = {}
        table['geo_file'] = info[0]
        table['inner_av_material'] = info[1]
        table['LED'] = info[2]
        table['fibre'] = info[3]
        table['reemis'] = info[4]

        # Write results for each region applied
        table['regions'] = []
        for i in range(len(output_stats_file_list)):
            output_stats_file = output_stats_file_list[i]
            region_lims = region_lims_list[i]
            table['regions'].append({})

            # Read in results to make json table for easier use
            with open(output_stats_file, "r") as f:
                stats = f.readlines()

            table['regions'][i]['region_lims'] = {}
            table['regions'][i]['region_lims']['direct_x_max'] = region_lims[0]
            table['regions'][i]['region_lims']['reflected_x_min'] = region_lims[1]
            table['regions'][i]['region_lims']['direct_y_centre'] = region_lims[2]
            table['regions'][i]['region_lims']['direct_dy'] = region_lims[3]
            table['regions'][i]['region_lims']['reflected_y_centre'] = region_lims[4]
            table['regions'][i]['region_lims']['reflected_dy'] = region_lims[5]

            # Get stats for each individual absorption factor
            table['regions'][i]['stats'] = {}
            for j in range(len(stats)):
                absorb = abs_factors[j]
                stats_info = stats[j].split(' ')  # [tot_hits, direct_hits, reflected_hits]
                table['regions'][i]['stats'][absorb] = {}
                table['regions'][i]['stats'][absorb]['tot_hits'] = int(stats_info[0])
                table['regions'][i]['stats'][absorb]['direct_hits'] = int(stats_info[1])
                table['regions'][i]['stats'][absorb]['reflected_hits'] = int(stats_info[2])

        # Save table to json file
        save_file = json_stats_folder + 'FinalStats_' + filename_format(info, True) + '.json'
        with open(save_file, 'w') as f:
            json.dump(table, f)

def getStats(args, input_info):
    '''Compute stats (number of hits in direct and reflected regions, and total hits) from different absorption scalings.
    If verbose flag is true, also create 2-D hists with regions overlain.'''
    print('Running getStats().')

    # Read in example macro and job script + info
    repo_address = getRepoAddress()

    # Make sure folders are of the correct format to  use later
    save_tothists_folder = checkRepo(args.tothist_repo, args.verbose)
    save_stats_folder = checkRepo(args.stats_repo, args.verbose)
    json_stats_folder = checkRepo(args.json_repo, args.verbose)

    # Wait until previous jobs are done
    checkJobsDone('tot_hists_', input_info, 10, args.verbose)

    # Apply region limits from list
    if args.verbose:
        print('Reading in region limits...')
    textFile = open(args.region_lims, "r")
    lines = [line.split(', ') for line in textFile]
    textFile.close()
    region_lims_list = np.asarray(lines).astype(np.float)

    # Get outer region limit strings
    outer_lims_list = getOuterLims(region_lims_list, args.verbose)


    ### MAKE JOB SCRIPTS TO RUN ANALYSIS ###
    print('Creating analysis job scripts...')
    jobScript_addresses = []
    output_info = []
    abs_factors = []
    info_1 = input_info[0]
    hist_files = ''
    for i in range(len(input_info)):
        info = input_info[i]
        if args.verbose:
            print('geo_file=', info[0], ', inner_av_material=', info[1], ', LED=', info[2], ', fibre=', info[3], ', reemis=', info[4], ', abs=', info[5])
        # Check whether to add arguments to current command, or start new one
        is_same = info_1[0] == info[0] and info_1[1] == info[1] and info_1[2] == info[2] and info_1[3] and info_1[4] == info[4]
        tot_hist_file_address = save_tothists_folder + 'tot_hists_' + filename_format(info) + '.root'

        # Group up files which are the same except for absorption scaling
        if is_same:
            hist_files += ' ' + tot_hist_file_address
            abs_factors.append(info[5])
        else:
            jobScript_address, output_stats_files = makeStatsJobScript(hist_files, info, outer_lims_list, repo_address, save_stats_folder, save_tothists_folder, args.verbose)
            jobScript_addresses.append(jobScript_address)
            output_info.append((output_stats_files, info_1, abs_factors))
            abs_factors = []
            hist_files = tot_hist_file_address
            info_1 = info
            abs_factors.append(info[5])
        if i == len(input_info) - 1:
            jobScript_address, output_stats_files = makeStatsJobScript(hist_files, info, outer_lims_list, repo_address, save_stats_folder, save_tothists_folder, args.verbose)
            jobScript_addresses.append(jobScript_address)
            output_info.append((output_stats_files, info_1, abs_factors))


    ### RUN JOB SCRIPTS ###
    print('Submitting job(s)...')
    for job_address in jobScript_addresses:
        command = 'qsub -t 1-' + str(len(outer_lims_list)) + ' -tc ' + str(args.max_jobs) + ' ' + job_address
        if args.verbose:
            print('Running command: ', command)
        subprocess.call(command, stdout=subprocess.PIPE, shell=True) # use subprocess to make code wait until it has finished

    # Wait until these job arrays are done
    checkJobsDone('stats_', input_info, 10, args.verbose)

    ### READ STATS AND WRITE TO JSON FILE ###
    print('Writing stats to json file...')
    writeJsonStatsFile(output_info, region_lims_list, json_stats_folder)

    return True


### MAIN ###

def main():
    # read in argument
    args = argparser()

    # read in file info
    textFile = open(args.list_file, "r")
    lines = [line.replace('\n', '').split(', ') for line in textFile]
    textFile.close()
    input_info = np.asarray(lines)

    # geo_files = input_info[:, 0]
    # inner_av_material = input_info[:, 1]
    # wavelengths = input_info[:, 2]
    # fibres = input_info[:, 3]
    # reemissions = input_info[:, 4]
    # abs_factors = input_info[:, 5]

    work_modes = {
        'sim': runSims,
        'hist': getHists,
        'combi': combiHists,
        'stats': getStats,

        'sim-hist': runSims_getHists,
        'hist-combi': getHists_combine,
        'combi-stats': combine_getStats,
        'hist-combi-stats': getHists_Combine_getStats,
        'all': runAll
    }

    result = work_modes[args.step](args, input_info)


if __name__ == '__main__':
    main()