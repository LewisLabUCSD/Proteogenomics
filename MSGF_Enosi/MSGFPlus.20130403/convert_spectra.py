#!/usr/bin/python


import sys
import getopt
import os
import shutil

def usage():
    print "<input folder> <output folder> <path to msconvert> <further mappings of conversion eg. mzXML:mzML >"

def get_file_stats(file_path):
    return {}

def convert_single_file(parameters):
    input_file = parameters["input_file"]
    extension_conversion_mapping = parameters["extension_conversion_mapping"]
    output_folder = parameters["output_folder"]
    allowed_passthrough_extensions = parameters["allowed_passthrough_extensions"]
    path_to_msconvert = parameters["path_to_msconvert"]
    clean_mzXML = parameters["clean_mzXML"]

    input_extension = os.path.splitext(input_file)[1][1:]

    if input_extension in allowed_passthrough_extensions:
        output_filename = os.path.join(output_folder, os.path.basename(input_file))
        shutil.copyfile(input_file, output_filename)
        if input_extension == "mzXML" and clean_mzXML == 1:
            #We need to do the removal for xsi:nil
            sed_removal_cmd = "sed -i -e 's/%s/%s/g' %s" % ("xsi:nil=\"true\"", "", output_filename)
            #sed_removal_cmd = "sed -i -e 's/%s/%s/g' %s" % ("xsi:nil", "", output_filename)
            print sed_removal_cmd
            os.system(sed_removal_cmd)
        return {"input_filename" : os.path.basename(input_file), "status" : "Copied"}
    elif input_extension in extension_conversion_mapping:
        output_filename = os.path.splitext(os.path.basename(input_file))[0] + "." + extension_conversion_mapping[input_extension]
        full_output_filename = os.path.join(output_folder, output_filename)

        cmd = path_to_msconvert + " " + input_file + " "
        cmd += "--32 "
        cmd += "--" + extension_conversion_mapping[input_extension] + " "
        cmd += "-o " + output_folder
        cmd += " --outfile " + output_filename
        os.system(cmd)

        #Take a look at the file and if the target was an mzXML, remove xsi:nil string
        if extension_conversion_mapping[input_extension] == "mzXML" and clean_mzXML == 1:
            #We need to do the removal for xsi:nil
            sed_removal_cmd = "sed -i -e 's/%s/%s/g' %s" % ("xsi:nil=\"true\"", "", full_output_filename)
            print sed_removal_cmd
            os.system(sed_removal_cmd)

        return {"input_filename" : os.path.basename(input_file), "status" : "Converted to " + extension_conversion_mapping[input_extension]}
    else:
        return {"input_filename" : os.path.basename(input_file), "status" : "Extension not allowed"}

def convert_input_folder(input_folder, output_folder, path_to_msconvert, allowed_passthrough_extensions, extension_conversion_mapping, clean_mzXML):
    all_input_files = [ os.path.join(input_folder,f) for f in os.listdir(input_folder) if os.path.isfile(os.path.join(input_folder,f)) ]

    output_status = []

    all_parameters = []

    for input_file in all_input_files:
        parameters = {}
        parameters["input_file"] = input_file
        parameters["extension_conversion_mapping"] = extension_conversion_mapping
        parameters["output_folder"] = output_folder
        parameters["allowed_passthrough_extensions"] = allowed_passthrough_extensions
        parameters["path_to_msconvert"] = path_to_msconvert
        parameters["clean_mzXML"] = clean_mzXML

        all_parameters.append(parameters)

    output_status = []
    for parameter in all_parameters:
        convert_single_file(parameter)

    return output_status

def main():
    input_folder = sys.argv[1]
    output_folder = sys.argv[2]
    path_to_msconvert = sys.argv[3]
    clean_mzXML = int(sys.argv[4])

    allowed_passthrough_extensions = []
    extension_conversion_mapping = {}

    for i in range(5, len(sys.argv)):
        conversion_parameter = sys.argv[i]
        from_extension = conversion_parameter.split(":")[0]
        to_extension = conversion_parameter.split(":")[1]
        extension_conversion_mapping[from_extension] = to_extension

        if from_extension == to_extension:
            allowed_passthrough_extensions.append(from_extension)

    output_status = convert_input_folder(input_folder, output_folder, path_to_msconvert, allowed_passthrough_extensions, extension_conversion_mapping, clean_mzXML)




if __name__ == "__main__":
    main()
