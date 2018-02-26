#!/usr/bin/python

import os, sys, subprocess
from subprocess import Popen, PIPE

tool = "MSGFPlus.jar"
script_home = os.path.dirname(os.path.realpath(__file__))
tool_path = os.path.join(script_home, tool)

def main(argv):
    command = ["java"]
    # extract Xmx, since it's a special parameter that
    # must go at the beginning of the command line
    for i, arg in enumerate(argv):
        if arg.lower().startswith("-xmx"):
            command.append(arg)
            del argv[i]
            break
    # build the rest of the command line
    command += ["-cp", tool_path, "edu.ucsd.msjava.ui.MSGFPlus"]
    for arg in argv:
        command.append(arg)
    print ' '.join(command) + "\n----------"
    # run the command, capture its console output
    process = Popen(command, stdin=PIPE, stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()
    print stdout + "\n----------\n" + stderr
    # check the output to determine if a silent error occurred
    if process.returncode == 0 and ("exception" in stdout.lower() or "exception" in stderr.lower()):
        print "Found exception, exiting with status [1]."
        sys.exit(1);

if __name__ == "__main__":
   main(sys.argv[1:])
