#!/usr/bin/python


import sys
import getopt
import os
import shutil
import ming_fileio_library

def usage():
    print "<input folder> <output folder> <path to msconvert> <further mappings of conversion eg. mzXML:mzML >"


def main():
    input_folder = sys.argv[1]
    input_tsvfile = sys.argv[2]
    output_tsvfile = sys.argv[3]

    allowed_passthrough_extensions = []
    extension_conversion_mapping = {}

    for i in range(4, len(sys.argv)):
        print(i)
        conversion_parameter = sys.argv[i]
        print(conversion_parameter)
        from_extension = conversion_parameter.split(":")[0]
        to_extension = conversion_parameter.split(":")[1]
        extension_conversion_mapping[from_extension] = to_extension

        if from_extension == to_extension:
            allowed_passthrough_extensions.append(from_extension)

    file_renaming_reverse_mapping = {}

    all_input_files = [ os.path.join(input_folder,f) for f in os.listdir(input_folder) if os.path.isfile(os.path.join(input_folder,f)) ]
    for input_file in all_input_files:
        input_extension = os.path.splitext(input_file)[1][1:]
        if input_extension in extension_conversion_mapping:
            renamed = os.path.splitext(os.path.basename(input_file))[0] + "." + extension_conversion_mapping[input_extension]
            file_renaming_reverse_mapping[renamed] = os.path.basename(input_file)


    row_count, table_data = ming_fileio_library.parse_table_with_headers(input_tsvfile)

    for header in table_data:
        for i in range(row_count):
            for find_to_replace in file_renaming_reverse_mapping:
                table_data[header][i] = table_data[header][i].replace(find_to_replace, file_renaming_reverse_mapping[find_to_replace])



    ming_fileio_library.write_dictionary_table_data(table_data, output_tsvfile)



if __name__ == "__main__":
    main()
