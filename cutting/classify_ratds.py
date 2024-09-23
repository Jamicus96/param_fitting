import os

def getRepoAddress():
    '''Returns the full address of the git repo containing with script'''
    repo_address = __file__[:-len('cutting/classify_ratds.py')]
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

run_GTID_map = {
    '300840' : 16172306,
    '300901' : 8076817,
    '300991' : 13283754,
    '300996' : 8933258,
    '301016' : 6319826,
    '301018' : 9587233,
    '301020' : 8512433,
    '301071' : 510346,
    '301096' : 15373145,
    '301213' : 6603410,
    '301356' : 15094346,
    '301357' : 13236217,
    '301476' : 3223113,
    '301653' : 1814583,
    '301664' : 3695230,
    '301708' : 7245924,
    '301892' : 7616837,
    '301947' : 7872258,
    '301988' : 9445857,
    '301994' : 10756052,
    '302005' : 6336779,
    '302895' : 4803480,
    '303471' : 1333077,
    '303738' : 13173840,
    '303767' : 14585202,
    '304039' : 12707480,
    '304120' : 14539193,
    '304928' : 11611175,
    '305236' : 10232334,
    '305536' : 6599503,
    '305564' : 5834383,
    '305601' : 16160926,
    '305665' : 15419458,
    '305688' : 8780761,
    '305960' : 8656140,
    '306357' : 8397987,
    '306376' : 2229033,
    '306756' : 5389833,
    '306806' : 4749727,
    '307191' : 15357177,
    '307210' : 13146082,
    '307315' : 10764972,
    '307455' : 14047256,
    '307468' : 3953659,
    '307580' : 6423278,
    '307604' : 4791488,
    '307982' : 11426723,
    '308131' : 3046897,
    '308132' : 8667209,
    '308165' : 2666701,
    '308174' : 182501,
    '308406' : 13106488,
    '308562' : 14264382,
    '308632' : 7415327,
    '308693' : 15247845,
    '308948' : 7864297,
    '308956' : 11630946
}

ratds_repo = '/mnt/lustre/scratch/epp/jp643/antinu/param_fitting/thesis/data_ratds/'
commandList_address = ratds_repo + 'job_scripts/commands.txt'


new_job_address = ratds_repo + 'job_scripts/'
new_job_address = checkRepo(new_job_address, True)
new_job_address += 'cutting_job.job'

output_logFile_address = ratds_repo + 'log_files/'
output_logFile_address = checkRepo(output_logFile_address, True)

path = getRepoAddress()
commandBase = path + 'cutting/cut_ratds.exe '
commandList_file = open(commandList_address, 'w')
for filename in os.listdir(ratds_repo):
    file_address = os.path.join(ratds_repo, filename)
    if os.path.isfile(file_address):
        if '.root' in filename:
            run = int(filename.split('_')[1][5:])
            GTID = run_GTID_map[str(run)]

            command = commandBase + file_address + ' ' + str(-8.81) + ' ' + str(run)\
                      + ' ' + str(GTID) + ' ' + str(int(True)) + ' ' + str(int(True))
            commandList_file.write(command + '\n')
commandList_file.close()


jobScript_address = path + 'job_scripts/fitting_job.job'
with open(jobScript_address, "r") as f:
    example_jobScript = f.readlines()

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