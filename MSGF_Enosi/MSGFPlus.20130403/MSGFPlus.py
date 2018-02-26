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
    if process.returncode == 0 and (has_bad_exception(stdout) or has_bad_exception(stderr)):
        print "Found exception, exiting with status [1]."
        sys.exit(1);
    else:
        sys.exit(process.returncode)

def has_bad_exception(content):
    lines = content.splitlines()
    for i, line in enumerate(lines):
        if "exception" in line.lower():
            if "javax.xml.stream.XMLStreamException: ParseError" in line and len(lines) > (i+1) and "xsi:nil" in lines[i+1]:
                continue
            if 'com.ctc.wstx.exc.WstxParsingException: Undeclared namespace prefix "xsi" (for attribute "nil")' in line:
                return False
            else:
                return True
    return False



if __name__ == "__main__":
   main(sys.argv[1:])
