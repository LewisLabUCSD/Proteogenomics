#!/usr/bin/python


import sys
import getopt
import os
import shutil
import math
import ming_fileio_library
import ming_proteosafe_library
import ming_protein_library
from collections import defaultdict

def usage():
    print("<param.xml> <sequence> <input filename> <output_filename>")

def add_decoy_to_results(table_data, row_count, decoy_marker):
    table_data["decoy"] = []
    table_data["sorting_score"] = []
    for i in range(row_count):
        if table_data["Protein"][i].find(decoy_marker) != -1:
            table_data["decoy"].append(1)
        else:
            table_data["decoy"].append(0)
        table_data["sorting_score"].append(-math.log10(float(table_data["EValue"][i])))

def reconcile_protein_repeats(psm_list):
    #Group by specID + peptide
    unique_psm_map = defaultdict(list)
    for psm in psm_list:
        key = psm["SpecID"] + ":" + psm["Peptide"]
        unique_psm_map[key].append(psm)

    output_list = []
    for key in unique_psm_map:
        if len(unique_psm_map[key]) == 0:
            output_list += unique_psm_map[key]
        else:
            protein_list = []
            for psm in unique_psm_map[key]:
                protein_list.append(psm["Protein"])
            first_psm = unique_psm_map[key][0]
            first_psm["Protein"] = ";".join(protein_list)
            output_list.append(first_psm)
    return output_list

def add_fdr_to_results(table_data, row_count):
    psm_list = []
    for i in range(row_count):
        psm_obj = {}
        for header in table_data:
            psm_obj[header] = table_data[header][i]

        psm_list.append(psm_obj)

    #Filter out repeats for proteins
    psm_list = reconcile_protein_repeats(psm_list)

    #Sort this B
    sorted_psm_list = compute_fdr_to_list(psm_list)

    peptide_dict = {}
    peptide_set = set()

    for psm in sorted_psm_list:
        if not psm["Peptide"] in peptide_dict:
            peptide_dict[psm["Peptide"]] = {}
            peptide_dict[psm["Peptide"]]["Peptide"] = psm["Peptide"]
            peptide_dict[psm["Peptide"]]["sorting_score"] = psm["sorting_score"]
            peptide_dict[psm["Peptide"]]["decoy"] = psm["decoy"]
        peptide_dict[psm["Peptide"]]["sorting_score"] = max(peptide_dict[psm["Peptide"]]["sorting_score"], psm["sorting_score"])


    #Do Peptide Level FDR
    peptide_object_list = []
    for peptide in peptide_dict:
        peptide_object_list.append(peptide_dict[peptide])
    sorted_peptides = compute_fdr_to_list(peptide_object_list)
    peptide_FDR_map = {}
    for peptide_obj in sorted_peptides:
        peptide_FDR_map[peptide_obj["Peptide"]] = peptide_obj["QValue"]

    output_psm_list = []
    #filter out decoys
    for psm in sorted_psm_list:
        psm["PepQValue"] = peptide_FDR_map[psm["Peptide"]]
        if psm["decoy"] == 0:
            output_psm_list.append(psm)



    return output_psm_list


def compute_fdr_to_list(input_list):
    running_target_count = 0
    running_decoy_count = 0
    sorted_psm_list = sorted(input_list, key=lambda psm: psm["sorting_score"], reverse=True)
    for psm in sorted_psm_list:
        if psm["decoy"] == 1:
            running_decoy_count += 1
        else:
            running_target_count += 1

        current_fdr = 1.0
        try:
            current_fdr = float(running_decoy_count) / float(running_target_count)
        except:
            continue
        psm["QValue"] = current_fdr

    sorted_psm_list.reverse()
    min_fdr = 1.0
    for psm in sorted_psm_list:
        min_fdr = min(min_fdr, psm["QValue"])
        psm["QValue"] = min_fdr

    return sorted_psm_list

def main():
    params = ming_proteosafe_library.parse_xml_file(open(sys.argv[1]))
    proteome = ming_protein_library.parse_fasta_proteome_file(sys.argv[2])

    row_count, table_data = ming_fileio_library.parse_table_with_headers(sys.argv[3])
    decoy_marker = sys.argv[5]

    add_decoy_to_results(table_data, row_count, decoy_marker)
    psm_results = add_fdr_to_results(table_data, row_count)

    output_table = defaultdict(list)

    #Performing filters
    filter_type = params["filter.filter"][0]
    if filter_type == "FDR":
        fdr_threshold = float(params["FDR.FDR"][0])
        for psm in psm_results:
            if psm["QValue"] < fdr_threshold:
                for key in psm:
                    output_table[key].append(psm[key])
    if filter_type == "PepFDR":
        fdr_threshold = float(params["PepFDR.PepFDR"][0])
        for psm in psm_results:
            if psm["PepQValue"] < fdr_threshold and psm["QValue"] < fdr_threshold:
                for key in psm:
                    output_table[key].append(psm[key])
    if filter_type == "FPR":
        print("Lets do nothing, don't know what this is")

    ming_fileio_library.write_dictionary_table_data(output_table, sys.argv[4])





if __name__ == "__main__":
    main()
