import shlex
import subprocess


def callSubprocess(cmdString, debugString):
    cmdargs = shlex.split(cmdString)
    #print(cmdString)
    #print(cmdargs)
    debFile = open('debug_file.txt', 'w')
    debFile.write('Starting %s \n' % debugString)
    retValue = subprocess.call(cmdString,stdout=debFile,shell=True)               # use shell=True with a single string of commands; shell=False with list of strings with first element the executable
    if (retValue==0):
        debFile.write('%s Successful\n' % debugString)
        debFile.close()
    else:
        debFile.write('There was error in %s\n' % debugString)
        debFile.close()
