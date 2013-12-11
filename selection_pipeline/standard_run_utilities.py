import os
import subprocess
import sys
import logging
import re
import time
#queue for threads
#regex for hash at start of line

try:
    import Queue
except ImportError:
    import queue
from threading import Thread
logger = logging.getLogger(__name__)
SUBPROCESS_FAILED_EXIT = 10
MISSING_EXECUTABLE_ERROR = 5

# returns Split VCF files that can be used
# by the vcf-subset function to take advantage of the cores avaliable
#
#
# split positions looks like an array
# 0,100,200 etc. The line ranges not
# including the header are seperated into
# seperate files which the function returns
#
# potential edge cases with really really small vcf files these should not be
# in use
#


def split_vcf(input_file, split_positions):
    """ Split a vcf file by input_positions

    """
    header = ''
    output_vcfs = []
    file_id = 0
    line_count = 1
    # get splits files positions 0 and 1
    # for a 1 core setup these will be
    # the start and the end of the file
    # and so the file will not change
    i = 0
    pos1 = split_positions[i]
    output_vcf = open(os.path.basename(input_file)+str(file_id), 'w')
    output_vcfs.append(os.path.basename(input_file)+str(file_id))
    with open(input_file, 'r') as vcf:
        for line in vcf:
            if re.match("^#", line) is not None:
                header += line
            else:
                output_vcf.write(header)
                output_vcf.write(line)
                break
        for line in vcf:
            if(line_count < pos1):
                output_vcf.write(line)
            else:
                i = i + 1
                pos1 = split_positions[i]
                file_id += 1
                out_name = input_file + str(file_id)
                output_vcfs.append(out_name)
                output_vcf = open(input_file+str(file_id), 'w')
                output_vcf.write(header)
                output_vcf.write(line)
            line_count += 1
    return(output_vcfs)


def get_vcf_line_count(input_file):
    """ Return the line count of a vcf file

    """
    with open(input_file, 'r') as vcf:
        line_count = 0
        for line in vcf:
            if re.match("^#", line) is not None:
                line_count = 1
            else:
                break
        for line in vcf:
            line_count += 1
        return(line_count)


def __is_script__(fpath):
    """ Return true if the path is a file

    """
    return os.path.isfile(fpath)


def __is_exe__(fpath):
    """ Return true if the path is a file and the executable bit is set

    """
    return os.path.isfile(fpath) and os.access(fpath, os.X_OK)


def which(program, program_name):
    """ Checks whether the file exists on the path or the system path

    """
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
    logger.error(program_name + " path = " + fpath +
                 " not locatable path or in the directory specified \
                 in your config file ")
    return None


def run_subprocess(
    command,tool,stdout=None,
    stderr=None,stdoutlog=False,
    working_dir = None):
    """ Runs a command on the system shell and forks a new process
   
	also creates a file for stderr and stdout if needed
	to avoid deadlock. 
    """
    # Very dirty hack 
    if(tool== 'selection_pipeline'):
        stderr = 'selection_stderr.tmp'
        stdout = 'selection_stdout.tmp'
    if (working_dir is not None):
        orig_dir = os.getcwd()
        os.chdir(working_dir)
    if(stderr is None):
        stderr = 'stderr.tmp'
        standard_err = open(stderr,'w')
    else:
        standard_err = open(stderr,'w')  
    try:
        if(stdout is None):
            standard_out = open('stdout.tmp','w')
            exit_code = subprocess.Popen(
                command,stderr=standard_out,stdout=standard_err)
        else:
        # find out what kind of exception to try here
            if(hasattr(stdout,'read')):
                exit_code = subprocess.Popen(
                    command,stdout=stdout,stderr=standard_err)
            else:
                stdout=open(stdout,'w')
                exit_code = subprocess.Popen(
                    command,stdout=stdout,stderr=standard_err)
            standard_out = stdout
    except:
        logger.error(tool + " failed to run " + ' '.join(command))
        standard_err = open(stderr,'r')
        while True:
            line = standard_err.readline() 
            if not line:
                break
            logger.info(tool +" STDERR: " + line.strip())
        standard_err.close()
        sys.exit(SUBPROCESS_FAILED_EXIT)
    exit_code.wait()
    standard_err.close()
    standard_out.close()
    standard_err = open(stderr,'r')
    if(exit_code.returncode != 0):
        logger.error(tool + " failed to run " +  ' '.join(command))
        while True:
            line = standard_err.readline() 
            if not line:
                break
            logger.info(tool +" STDERR: " + line.strip())
        sys.exit(SUBPROCESS_FAILED_EXIT)
    stdout_log = False
    if(stdout is None):
        standard_out = open('stdout.tmp', 'r')  
        stdout_log = True
    elif(stdoutlog):
        if(hasattr(stdout,'write')):   
            standard_out = open(stdout.name,'r')
        else: 
            standard_out = open(stdout, 'r')  
        stdout_log = True
    if(stdout_log):
        while True:
            line = standard_out.readline()
            if not line:
                break
            logger.info(tool + " STDOUT: " +line.strip())
        standard_out.close()
    while True:
        line = standard_err.readline()
        if not line:
            break
        logger.info(tool +" STDERR: " + line.strip())
    logger.info("Finished tool " + tool)
    standard_err.close()
    standard_out.close()
    # Removed stdout if it either was not specified
    # or the log was specified.
    if(stdout is None):
        os.remove('stdout.tmp')
    elif(stdoutlog):
        os.remove(standard_out.name)
    os.remove(stderr)
    if (working_dir is not None):
        os.chdir(orig_dir)


def __queue_worker__(q,tool_name):
    while True:
        queue_item = q.get()
        try:
            cmd = queue_item[0]
            stdout = queue_item[1]
            stdoutlog = queue_item[2]
            stderr = queue_item[3]
            folder_names = queue_item[4]
        except IndexError:
            cmd = queue_item[0]
            stdout = queue_item[1] 
            stdoutlog = False
            stderr = None
        try:
            run_subprocess(
                cmd, tool_name, stdout=stdout,
                stdoutlog=stdoutlog,stderr=stderr)
        except SystemExit:
            logger.error(tool_name + ": Failed to run in thread")
            sys.exit(SUBPROCESS_FAILED_EXIT)
        q.task_done()


def queue_jobs(commands, tool_name, threads, stdouts=None,folder_names=None):
    """ Creates a queue for running jobs

        Using a synchronized queue to spawn jobs equal
        to the number of cores specified to the user.
        The method blocks until all tasks are complete
    """
    q = Queue.Queue()
    for i in range(int(threads)):
        t = Thread(target=__queue_worker__, args=[q, tool_name])
        t.daemon = True
        t.start()
    for i, cmd in enumerate(commands):
        stderr = 'stderr' + str(i) + '.tmp'
        if(folder_names is None):
            folder_name = '.'
        else:
            folder_name = folder_names[i]
        if (stdouts is not None):
            q.put([cmd, stdouts[i], False, stderr,folder_name])
        else:
            stdout = 'stdout' + str(i) + '.tmp' 
            q.put([cmd, stdout, True, stderr,folder_name])
    q.join()

# clean folder expecting a list containing
# files to keep from that folder
# only required if the user
# runs the analysis from their root directory


def clean_folder(folder, keep=None):
    """ Cleans the working directory

        Takes as a parameter a list of the files
        to not delete.
    """
    for the_file in os.listdir(folder):
        file_path = os.path.join(folder, the_file)
        print(keep)
        if keep is not None:
            if (file_path in [os.path.join(folder, x) for x in keep]):
                continue
        try:
            if os.path.isfile(file_path):
                os.unlink(file_path)
        except Exception as e:
            logger.error(e)
