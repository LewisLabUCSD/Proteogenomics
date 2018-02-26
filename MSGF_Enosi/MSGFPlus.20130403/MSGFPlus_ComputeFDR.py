#!/usr/bin/python

import os, sys, subprocess
from subprocess import Popen, PIPE

# child process tool parameters
tool = "MSGFPlus.jar"
script_home = os.path.dirname(os.path.realpath(__file__))
tool_path = os.path.join(script_home, tool)

# MS-GF+ TSV output column names
spec_file_column = "#SpecFile"
spec_id_column = "SpecID"
peptide_column = "Peptide"
protein_column = "Protein"
score_column = "SpecEValue"

def main(argv):
    command = ["java"]
    # get TSV file parameters from command line arguments
    input_file = None
    xmx_command_index = None
    protein_column_command_index = None
    score_column_command_index = None
    for i, arg in enumerate(argv):
        # extract Xmx, since it's a special parameter that
        # must go at the beginning of the command line
        if arg.lower().startswith("-xmx"):
            command.append(arg)
            xmx_command_index = i
        elif arg == "-f":
            if len(argv) <= i + 2:
                print "The \"-f\" parameter must be followed by two arguments."
                sys.exit(1);
            input_file = argv[i + 1]
            protein_column_command_index = i + 2
        elif arg == "-s":
            if len(argv) <= i + 1:
                print "The \"-s\" parameter must be followed by one argument."
                sys.exit(1);
            score_column_command_index = i + 1
    # be sure all relevant arguments were found
    if input_file is None or protein_column_command_index is None:
        print "A \"-f\" parameter, followed by a valid input TSV filename, must be provided."
        sys.exit(1)
    elif score_column_command_index is None:
        print "A \"-s\" parameter, followed by a valid score sort argument (0/1), must be provided."
        sys.exit(1)
    # parse input TSV file to gather relevant column indices from its header
    spec_file_column_index = None
    spec_id_column_index = None
    peptide_column_index = None
    protein_column_index = None
    score_column_index = None
    with open(input_file, 'r') as tsv:
        header = tsv.readline().rstrip().split()
        for i, column in enumerate(header):
            if column == spec_file_column:
                spec_file_column_index = i
            elif column == spec_id_column:
                spec_id_column_index = i
            elif column == peptide_column:
                peptide_column_index = i
            elif column == protein_column:
                protein_column_index = i
            elif column == score_column:
                score_column_index = i
    # be sure all relevant column headers were found
    if spec_file_column_index is None:
        print "A \"" + spec_file_column + "\" column was not found in TSV file [" + input_file + "]"
        sys.exit(1)
    elif spec_id_column_index is None:
        print "A \"" + spec_id_column + "\" column was not found in TSV file [" + input_file + "]"
        sys.exit(1)
    elif peptide_column_index is None:
        print "A \"" + peptide_column + "\" column was not found in TSV file [" + input_file + "]"
        sys.exit(1)
    elif protein_column_index is None:
        print "A \"" + protein_column + "\" column was not found in TSV file [" + input_file + "]"
        sys.exit(1)
    elif score_column_index is None:
        print "A \"" + score_column + "\" column was not found in TSV file [" + input_file + "]"
        sys.exit(1)
    # append column header arguments
    argv.append("-i")
    argv.append(`spec_file_column_index`)
    argv.append("-n")
    argv.append(`spec_id_column_index`)
    argv.append("-p")
    argv.append(`peptide_column_index`)
    # clean arguments to be added to child process command
    if xmx_command_index is not None:
        del argv[xmx_command_index]
        if protein_column_command_index > xmx_command_index:
            protein_column_command_index -= 1
        if score_column_command_index > xmx_command_index:
            score_column_command_index -= 1
    argv.insert(protein_column_command_index, `protein_column_index`)
    if score_column_command_index > protein_column_command_index:
        score_column_command_index += 1
    argv.insert(score_column_command_index, `score_column_index`)
    # build the rest of the command line
    command += ["-cp", tool_path, "edu.ucsd.msjava.fdr.ComputeFDR"]
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
        sys.exit(1)

if __name__ == "__main__":
   main(sys.argv[1:])
