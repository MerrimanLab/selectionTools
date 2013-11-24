import os
import subprocess
import sys
import logging

#queue for threads
#regex for hash at start of line


import Queue
from threading import Thread
logger = logging.getLogger(__name__)
SUBPROCESS_FAILED_EXIT=10
MISSING_EXECUTABLE_ERROR=5

# Get end of haps file
def haps_start_and_end(file):
    return 0
# Get start of vcf file
def get_vcf_line_count(input_file):
    with open(input_file,'r') as vcf:
        line_count = 0
        for line in vcf:
            if re.match("^#",line) is not None:
                line_count = 1
                break
        for line in vcf:
            line_count += 1
        return(line_count)
        
def vcf_start_and_end(file):
    return 0
#

def __is_script__(fpath):
        return os.path.isfile(fpath)
def __is_exe__(fpath):
    return os.path.isfile(fpath) and os.access(fpath, os.X_OK)
    #Stolen code from 
    #http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python
def which(program,program_name):
    fpath, fname = os.path.split(program)
    if fpath:
        if __is_exe__(program):
            return program
        elif (__is_script__(program)):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if __is_exe__(exe_file):
                return exe_file
    logger.error(program_name +" path = " + fpath+" not locatable path or in the directory specified in your config file ")
    return None



def run_subprocess(command,tool,stdout=None):
        try:
            if(stdout is None):
                exit_code = subprocess.Popen(command,stderr=subprocess.PIPE,stdout=subprocess.PIPE)
            else:
                exit_code = subprocess.Popen(command,stdout=stdout,stderr=subprocess.PIPE)
        except:
            logger.error(tool + " failed to run " + ' '.join(command))
            sys.exit(SUBPROCESS_FAILED_EXIT)
        exit_code.wait()
        if(exit_code.returncode != 0):
            logger.error(tool + "failed to run " +  ' '.join(command))
            while True:
                line = exit_code.stderr.readline()
                if not line:
                    break
                logger.info(tool +" STDERR: " + line.strip())
            sys.exit(SUBPROCESS_FAILED_EXIT)
        if(stdout is None):
            while True:
                line = exit_code.stdout.readline()
                if not line:
                    break
                logger.info(tool + " STDOUT: " +line.strip())
        while True:
            line = exit_code.stderr.readline()
            if not line:
                break
            logger.info(tool +" STDERR: " + line.strip())
        logger.error("Finished tool " + tool)

def __queue_worker__(q):
    while True:
        cmd=q.get()
        run_subprocess(cmd,'impute2')
        q.task_done()

def queue_jobs(commands,threads):
    q = queue.Queue()
    for i in range(int(self.threads)):
        t = Thread(target=self.queue_worker,args=[commands])
        t.daemon = True
        t.start()
    for command in commands:
        q.put(command)  
    q.join()

def clean_folder(folder):
    for the_file in os.listdir(folder):
        file_path = os.path.join(folder, the_file)
        try:
            if os.path.isfile(file_path):
                os.unlink(file_path)
        except Exception, e:
            logger.error(e)
