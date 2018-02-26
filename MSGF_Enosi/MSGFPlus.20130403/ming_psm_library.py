#!/usr/bin/python

"""

PSM Utilities to read psms

"""

import ming_fileio_library
import math
import re
from pyteomics import mass

class PSM:
    def __init__(self, filename, scan, annotation, score, decoy, protein, charge):
        self.filename = filename
        self.scan = scan
        self.annotation = annotation
        self.score = score
        self.decoy = decoy
        self.protein = protein
        self.charge = charge
        self.fdr = -1.0
        self.ppm_error = -1.0
        self.extra_metadata = {}

    @staticmethod
    def output_header():
        return_headers = "sequence\tscore\tdecoy\tFDR\tfilename\tscan\tcharge\tppm_error"
        return return_headers

    def __str__(self):
        return "%s\t%s\t%s\t%s\t%s\t%s\t%d\t%f" % (self.annotation, str(self.score), str(self.decoy), str(self.fdr), self.filename, str(self.scan), self.charge, self.ppm_error)

    def __repr__(self):
        return str(self)

    def get_extra_metadata_headers(self):
        return "\t".join(self.extra_metadata.keys())

    def get_extra_metadata_values(self):
        return "\t".join(self.extra_metadata.values())

    def is_decoy(self):
        return self.decoy

    def sorting_value(self):
        return self.score

    def get_stripped_sequence(self):
        sequence = self.annotation
        p = re.compile('\W|\d')
        sequence = p.sub("", sequence)
        return sequence

    def get_annotation_without_charge(self):
        if self.annotation[-2] == ".":
            return self.annotation[:-2]
        return self.annotation


class PSMset:
    def __init__(self, name):
        self.name = name
        self.psms = []

    def __len__(self):
        return len(self.psms)



    #Loading a TSV File from MSGFDB
    def load_MSGF_tsvfile(self, filename):
        self.psms += parse_MSGF_tsvfile(filename)

    def load_PSM_tsvfile(self, filename, load_extra_metadata=False):
        self.psms = parse_psm_file(filename, load_extra_metadata)

    #Filter PSMs to given FDR
    def filter_to_fdr(self, fdr):
        filtered_psms = filter_psm_fdr(self.psms, fdr)
        print("Filtered " + str(len(self.psms)) + " to " + str(len(filtered_psms)))
        self.psms = filtered_psms

    def filter_to_fdr_by_length(self, fdr):
        output_psms = []
        peptide_length_map = {}
        for psm in self.psms:
            peptide_length = len(psm.get_stripped_sequence())
            if not peptide_length in peptide_length_map:
                peptide_length_map[peptide_length] = []
            peptide_length_map[peptide_length].append(psm)

        for peptide_length in peptide_length_map:
            filtered_psms = filter_psm_fdr(peptide_length_map[peptide_length], fdr)
            print("Filtered Length " + str(peptide_length) + " " + str(len(peptide_length_map[peptide_length])) + " to " + str(len(filtered_psms)))
            output_psms += filtered_psms
        self.psms = output_psms

    #Calculate FDR of PSM Set
    def calculate_fdr(self):
        running_target_count = 0
        running_decoy_count = 0
        for psm in self.psms:
            if psm.is_decoy() == 0:
                running_target_count += 1
            else:
                running_decoy_count += 1

        current_fdr = float(running_decoy_count) / float(running_target_count)
        return current_fdr

    def write_output(self, output_file, write_extra_metadata=False):
        if write_extra_metadata:
            if len(self.psms) > 0:
                output_headers = PSM.output_header() + "\t" + self.psms[0].get_extra_metadata_headers() + "\n"
                output_file.write(output_headers)
                for psm in self.psms:
                    output_file.write(str(psm) + "\t" + psm.get_extra_metadata_values() + "\n")

            else:
                output_file.write(PSM.output_header() + "\n")
        else:
            output_file.write(PSM.output_header() + "\n")
            for psm in self.psms:
                output_file.write(str(psm) + "\n")


class PeptideVariant:
    def __init__(self, variant_sequence):
        self.variant_sequence = variant_sequence
        self.psms = []
        self.fdr = -1

    @staticmethod
    def output_header():
        return "variant_sequence\tscore\tdecoy\tFDR\tfilename\tscan\tcharge\tppm_error"

    def __str__(self):
        max_psm = self.get_best_psm()
        return str(max_psm)

    def add_psm(self, psm_object):
        self.psms.append(psm_object)

    def is_decoy(self):
        return self.psms[0].is_decoy()

    def sorting_value(self):
        max_score = 0
        for psm in self.psms:
            max_score = max(psm.sorting_value(), max_score)
        return max_score

    def get_best_psm(self):
        max_score = -1
        max_psm = None
        for psm in self.psms:
            if psm.sorting_value() > max_score:
                max_score = psm.sorting_value()
                max_psm = psm
        return max_psm

    def get_spectrum_count(self):
        return len(self.psms)

    def get_stripped_sequence(self):
        sequence = self.variant_sequence
        p = re.compile('\W|\d')
        sequence = p.sub("", sequence)
        return sequence

    def sequence_length(self):
        return len(self.get_stripped_sequence())

    #Research-y Portion

###
 # Class to hold a set of library peptides
###
class PeptideVariantSet:
    def __init__(self, name):
        self.name = name
        self.peptide_list = []
        self.peptide_map = {}

    def __len__(self):
        return len(self.peptide_list)

    def get_total_spectra_count(self):
        return sum(self.get_spectra_count_list())
        total_ms_ms_count = 0
        for variant in self.peptide_list:
            total_ms_ms_count += len(variant.psms)
        return total_ms_ms_count

    #Get total peptides regardless of modifications
    def get_total_unique_sequence_count(self):
        return len(self.get_unique_sequences())

    def get_unique_sequences_spectrum_count_map(self):
        sequence_map = {}
        for variant in self.peptide_list:
            sequence = variant.variant_sequence
            p = re.compile('\W|\d')
            sequence = p.sub("", sequence)
            if not sequence in sequence_map:
                sequence_map[sequence] = 0
            sequence_map[sequence] += variant.get_spectrum_count()

        return sequence_map

    def get_unique_sequences(self):
        sequence_map = {}
        for variant in self.peptide_list:
            sequence = variant.variant_sequence
            p = re.compile('\W|\d')
            sequence = p.sub("", sequence)
            sequence_map[sequence] = 1

        return sequence_map.keys()

    #Returns a list of spectral counts for each variant
    def get_spectra_count_list(self):
        ms_ms_count_list = []
        for variant in self.peptide_list:
            ms_ms_count_list.append(len(variant.psms))
        return ms_ms_count_list

    def add_psms_set(self, psm_set):
        self.add_psms_list(psm_set.psms)

    def add_psms_list(self, psm_list):
        for psm in psm_list:
            if not(psm.annotation in self.peptide_map):
                peptide_variant = PeptideVariant(psm.annotation)
                self.peptide_list.append(peptide_variant)
                self.peptide_map[psm.annotation] = peptide_variant
            self.peptide_map[psm.annotation].add_psm(psm)

    #Appending a variant set
    def add_variant_set(self, variant_set):
        for variant in variant_set.peptide_list:
            if not(variant.variant_sequence in self.peptide_map):
                self.peptide_list.append(variant)
                self.peptide_map[variant.variant_sequence] = variant
            else:
                self.add_psms_list(variant.psms)

    def add_variant(self, variant_obj):
        if not(variant_obj.variant_sequence in self.peptide_map):
            self.peptide_list.append(variant_obj)
            self.peptide_map[variant_obj.variant_sequence] = variant_obj
        else:
            self.add_psms_list(variant_obj.psms)

    def remove_variant(self, variant_obj):
        self.peptide_list.remove(variant_obj)
        del self.peptide_map[variant_obj.variant_sequence]

    def filter_to_fdr(self, fdr):
        filtered_peptides = filter_psm_fdr(self.peptide_list, fdr)
        print("Filtered " + str(len(self.peptide_list)) + " to " + str(len(filtered_peptides)))
        self.peptide_list = filtered_peptides
        self.peptide_map = {}
        for variant in self.peptide_list:
            self.peptide_map[variant.variant_sequence] = variant

    def filter_to_fdr_by_length(self, fdr):
        output_peptides = []
        peptide_length_map = {}
        for peptide_obj in self.peptide_list:
            peptide_length = len(peptide_obj.get_stripped_sequence())
            if not peptide_length in peptide_length_map:
                peptide_length_map[peptide_length] = []
            peptide_length_map[peptide_length].append(peptide_obj)

        for peptide_length in peptide_length_map:
            filtered_peptides = filter_psm_fdr(peptide_length_map[peptide_length], fdr)
            print("Filtered Length " + str(peptide_length) + " " + str(len(peptide_length_map[peptide_length])) + " to " + str(len(filtered_peptides)))
            output_peptides += filtered_peptides
        self.peptide_list = output_peptides
        self.peptide_map = {}
        for variant in self.peptide_list:
            self.peptide_map[variant.variant_sequence] = variant

    def calculate_fdr(self):
        running_target_count = 0
        running_decoy_count = 0
        for psm in self.peptide_list:
            if psm.is_decoy() == 0:
                running_target_count += 1
            else:
                running_decoy_count += 1

        current_fdr = float(running_decoy_count) / float(running_target_count)
        return current_fdr

    def write_output(self, output_file):
        output_file.write(PeptideVariant.output_header() + "\n")
        for variant in self.peptide_list:
            output_file.write(str(variant) + "\n")


###
 # Parsing peptide string into a list of literals including mods in the string
###
def get_peptide_modification_list_inspect_format(peptide):
    return re.findall('[^A-Z]*[A-Z][^A-Z]*', peptide)

known_modification_masses = mass.std_aa_comp

def create_theoretical_peak_map(peptide, ion_type_list, charge_set=[1]):
    amino_acid_list = get_peptide_modification_list_inspect_format(peptide)

    only_letters_list = [letter for letter in peptide if letter.isalpha()]

    only_mods_mass_add_list = []
    for amino_acid in amino_acid_list:
        mod_mass_to_add = 0.0
        mod_strings_tokenized = re.findall('[+-][1-9]*', re.sub("[A-Z]", "", amino_acid))
        for mod_tokenized in mod_strings_tokenized:
            mod_mass_to_add += float(mod_tokenized)
        only_mods_mass_add_list.append(mod_mass_to_add)

    ion_to_mass_mapping = {}
    for charge in charge_set:
        for ion_type in ion_type_list:
            for i in range(len(amino_acid_list)):
                peak_mass = 0.0
                if ion_type in "abc":
                    peak_annotation = ion_type + ":" + str(i+1) + ":" + str(charge)
                    peak_mass = mass.fast_mass("".join(only_letters_list[:i+1]), ion_type=ion_type, charge=charge) + sum(only_mods_mass_add_list[:i+1])
                else:
                    peak_annotation = ion_type + ":" + str(len(amino_acid_list) - i) + ":" + str(charge)
                    peak_mass = mass.fast_mass("".join(only_letters_list[i:]), ion_type=ion_type, charge=charge) + sum(only_mods_mass_add_list[i:])
                ion_to_mass_mapping[peak_annotation] = peak_mass

    return ion_to_mass_mapping

#Returns both annotated and unannotated peaks
def extract_annotated_peaks(ion_peak_mapping, peak_list, tolerance):
    extracted_peaks = []
    unannotated_peaks = []
    for peak in peak_list:
        mass = peak[0]
        isAnnotated = False
        for ion_peak in ion_peak_mapping:
            if abs(mass - ion_peak_mapping[ion_peak]) < tolerance:
                #extracted_peaks.append(peak)
                isAnnotated = True
                break
        if isAnnotated:
            extracted_peaks.append(peak)
        else:
            unannotated_peaks.append(peak)
    return extracted_peaks, unannotated_peaks



###
 # Takes as input a filename
 # Returns a list of PSM
###
def parse_MSGF_tsvfile(filename):
    rows, table_data = ming_fileio_library.parse_table_with_headers(filename)

    scan_header = "Scan#"
    peptide_header = "Peptide"
    protein_header = "Protein"
    score_header = "P-value"
    filename_header = "#SpecFile"
    charge_header = "Charge"
    ppm_error_header = "PMError(ppm)"
    da_pm_error_header = "PMError(Da)"
    precursor_header = "Precursor"

    parse_da_error = False
    if not ppm_error_header in table_data:
        parse_da_error = True


    decoy_indicator = "REV_"

    psm_list = []

    for i in range(rows):
        scan = table_data[scan_header][i]
        peptide = table_data[peptide_header][i]
        protein = table_data[protein_header][i]
        score = -math.log10(float(table_data[score_header][i]))
        #print table_data[score_header][i] + "\t" + str(score)
        filename = table_data[filename_header][i]
        charge = int(table_data[charge_header][i])
        if parse_da_error:
            ppm_error = float(table_data[da_pm_error_header][i])/float(table_data[precursor_header][i]) * 1000000
        else:
            ppm_error = float(table_data[ppm_error_header][i])
        decoy = 0

        #Stripping peptide dots
        if peptide[1] == "." and peptide[-2] == ".":
            peptide = peptide[2:-2]


        if protein.find(decoy_indicator) != -1:
            decoy = 1

        #Adding charge state to peptide name
        peptide += "." + str(charge)

        new_psm = PSM(filename, scan, peptide, score, decoy, protein, charge)
        new_psm.ppm_error = ppm_error
        psm_list.append(new_psm)

    return psm_list

def parse_MSGFPlus_tsvfile(filename):
    rows, table_data = ming_fileio_library.parse_table_with_headers(filename)

    scan_header = "ScanNum"
    peptide_header = "Peptide"
    protein_header = "Protein"
    score_header = "EValue"
    filename_header = "#SpecFile"
    charge_header = "Charge"
    ppm_error_header = "PrecursorError(ppm)"
    da_pm_error_header = "PrecursorError(Da)"
    precursor_header = "Precursor"

    parse_da_error = False
    if not ppm_error_header in table_data:
        parse_da_error = True


    decoy_indicator = "REV_"

    psm_list = []

    for i in range(rows):
        scan = table_data[scan_header][i]
        peptide = table_data[peptide_header][i]
        protein = table_data[protein_header][i]
        score = -math.log10(float(table_data[score_header][i]))
        #print table_data[score_header][i] + "\t" + str(score)
        filename = table_data[filename_header][i]
        charge = int(table_data[charge_header][i])
        if parse_da_error:
            ppm_error = float(table_data[da_pm_error_header][i])/float(table_data[precursor_header][i]) * 1000000
        else:
            ppm_error = float(table_data[ppm_error_header][i])
        decoy = 0

        #Stripping peptide dots
        if peptide[1] == "." and peptide[-2] == ".":
            peptide = peptide[2:-2]


        if protein.find(decoy_indicator) != -1:
            decoy = 1

        #Adding charge state to peptide name
        peptide += "." + str(charge)

        new_psm = PSM(filename, scan, peptide, score, decoy, protein, charge)
        new_psm.ppm_error = ppm_error
        psm_list.append(new_psm)

    return psm_list

###
 # Takes as input a filename for a variant file output by this code
###
def parse_variant_file(filename):
    rows, table_data = ming_fileio_library.parse_table_with_headers(filename)

    psm_list = []
    for i in range(rows):
        filename = table_data["filename"][i]
        scan = int(table_data["scan"][i])
        score = float(table_data["score"][i])
        decoy = int(table_data["decoy"][i])
        variant_sequence = table_data["variant_sequence"][i]
        charge = 0
        if "charge" in table_data:
            charge = int(table_data["charge"][i])
        else:
            charge = int(variant_sequence.split(".")[-1])
        protein = "NONE"

        if "unmangled_name" in table_data:
            filename = table_data["unmangled_name"][i]

        new_psm = PSM(filename, scan, variant_sequence, score, decoy, protein, charge)
        psm_list.append(new_psm)

    return psm_list

###
 # Takes as input a filename for a variant file output by this code
###
def parse_psm_file(filename, load_extra_metadata=False):
    rows, table_data = ming_fileio_library.parse_table_with_headers(filename)

    known_headers = ["filename", "scan", "score", "decoy", "sequence", "charge", "ppm_error", "unmangled_name", "FDR"]
    extra_metadata_headers = set(table_data.keys()).difference(set(known_headers))

    psm_list = []
    for i in range(rows):
        filename = table_data["filename"][i]
        scan = int(table_data["scan"][i])
        score = float(table_data["score"][i])
        decoy = int(table_data["decoy"][i])
        variant_sequence = table_data["sequence"][i]
        charge = int(table_data["charge"][i])
        ppm_error = float(table_data["ppm_error"][i])
        fdr = float(table_data["FDR"][i])
        protein = "NONE"

        if "unmangled_name" in table_data:
            filename = table_data["unmangled_name"][i]

        new_psm = PSM(filename, scan, variant_sequence, score, decoy, protein, charge)
        new_psm.ppm_error = ppm_error
        new_psm.fdr = fdr

        if load_extra_metadata:
            extra_metadata = {}
            for header in extra_metadata_headers:
                extra_metadata[header] = table_data[header][i]
            new_psm.extra_metadata = extra_metadata

        psm_list.append(new_psm)

    return psm_list

###
 # Filtering PSM results so that the returned set is at the
 # prescribed FDR
###
def filter_psm_fdr(input_psms, fdr_percentage):
    #TODO Add Sorting
    #Sorting
    #print "Sorting shit to " + str(fdr_percentage*100) + "%"
    #input_psms = sorted(input_psms, key=PSM.sorting_value)
    input_psms = sorted(input_psms, key=lambda psm: psm.sorting_value(), reverse=True)

    running_target_count = 1
    running_decoy_count = 0

    output_psms = []

    for psm in input_psms:
        if psm.is_decoy() == 0:
            running_target_count += 1
        else:
            running_decoy_count += 1

        current_fdr = float(running_decoy_count) / float(running_target_count)

        psm.fdr = current_fdr

    #Properly finding the min q value for PSM
    min_fdr = 1
    for i in range(len(input_psms)):
        index_to_check = len(input_psms) - i - 1
        min_fdr = min(min_fdr, input_psms[index_to_check].fdr)
        input_psms[index_to_check].fdr = min_fdr

    for psm in input_psms:
        if psm.fdr < fdr_percentage:
            output_psms.append(psm)

    return output_psms
