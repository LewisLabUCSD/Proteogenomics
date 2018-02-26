#!/usr/bin/python


import sys
import getopt
import os
import re
from collections import defaultdict

if (sys.version_info > (3, 0)):
    # Python 3 code in this block
    import _pickle as pickle
else:
    # Python 2 code in this block
    import cPickle as pickle


class ProteinAnnotation:
    def __init__(self, annotation_class, start_site, end_site, annotation_text):
        self.annotation_class = annotation_class
        self.start_site = start_site
        self.end_site = end_site
        self.annotation_text = annotation_text


def parse_fasta_proteome_file(filename):
    protein_sequence = ""
    protein_name = ""
    gene_name = ""
    proteome = Proteome()

    for line in open(filename, "r"):
        line = line.strip()
        if len(line) == 0:
            continue
        if line[0] == ">":
            if not len(protein_sequence) == 0:
                new_protein = Protein(protein_name, protein_sequence, gene_name)
                proteome.add_protein(new_protein)

            protein_name = line.split(" ")[0][1:]
            protein_sequence = ""
            #Try to get gene name
            name_splits = line.split(" ")
            for split in name_splits:
                if split.find("GN=") != -1:
                    gene_name = split[3:]
        else:
            protein_sequence += line.rstrip()

    #Get last protein
    new_protein = Protein(protein_name, protein_sequence, gene_name)
    proteome.add_protein(new_protein)

    return proteome

class Proteome:
    def __init__(self):
        self.protein_list = []
        self.protein_map = {}

    def add_protein(self, new_protein):
        self.protein_list.append(new_protein)
        self.protein_map[new_protein.protein] = new_protein

    def get_protein(self, protein_name):
        if protein_name in self.protein_map:
            return self.protein_map[protein_name]
        return None
        #for protein in self.protein_list:
        #    if protein.protein == protein_name:
        #        return protein
        #return None

    def get_peptides_mapped_to_proteins(self, peptide_list):
        peptide_mapping = defaultdict(list)
        for peptide in peptide_list:
            for protein in self.protein_list:
                if protein.contain_peptide(peptide) == True:
                    peptide_mapping[peptide].append(protein.protein)
        return peptide_mapping

    def get_peptides_mapped_to_proteins_efficient(self, peptide_list):
        peptide_mapping = defaultdict(list)

        #hash 4-mers
        fourmer_hash_to_proteins = defaultdict(list)
        for protein in self.protein_list:
            substring_set = set()
            for i in range(len(protein.sequence)):
                substring = protein.sequence[i:i+4]
                substring_set.add(substring)
            for substring in substring_set:
                fourmer_hash_to_proteins[substring].append(protein.protein)

        for peptide in peptide_list:
            protein_candidates = fourmer_hash_to_proteins[peptide[:4]]
            for protein_name in protein_candidates:
                if self.protein_map[protein_name].contain_peptide(peptide) == True:
                    peptide_mapping[peptide].append(protein_name)

        return peptide_mapping

    #Return a list of proteins that are covered by a set of peptides
    def get_proteins_covered_by_peptides_hash(self, peptide_list, do_reverse=False, protein_pkl_filename = None):
        peptide_to_protein_hash = {}
        create_hash = False

        if protein_pkl_filename == None:
            create_hash = True
        else:
            if os.path.exists(protein_pkl_filename):
                #Load that shit
                create_hash = False
                peptide_to_protein_hash = pickle.load(open(protein_pkl_filename, 'rb'))
            else:
                create_hash = True


        if create_hash == True:
            print("Creating Hash")
            protein_count = 0
            for protein in self.protein_list:
                protein_count += 1
                print(protein.protein + "\t" + str(protein_count) + " of " + str(len(self.protein_list)))
                #print peptide_to_protein_hash
                for i in range(len(protein.sequence)):
                    for j in range(len(protein.sequence)):
                        if j < i or j - i > 20 or j - i < 4:
                            continue
                        substring = protein.sequence[i:j+1]
                        if not substring in peptide_to_protein_hash:
                            peptide_to_protein_hash[substring] = set()
                        peptide_to_protein_hash[substring].add(protein.protein)

            if protein_pkl_filename != None:
                pickle.dump(peptide_to_protein_hash, open(protein_pkl_filename, 'rb'))

        #Doing Actual Search
        search_peptide_list = []

        if do_reverse == True:
            for peptide in peptide_list:
                rev_peptide = peptide[::-1]
                search_peptide_list.append(peptide)
                search_peptide_list.append(rev_peptide)
        else:
            search_peptide_list = peptide_list

        #Do shit now with it
        return_proteins = []
        for peptide in search_peptide_list:
            if peptide in peptide_to_protein_hash:
                return_proteins += list(peptide_to_protein_hash[peptide])

        return return_proteins


    #Return a list of proteins that are covered by a set of peptides
    def get_proteins_covered_by_peptides(self, peptide_list, do_reverse=False):
        return_proteins = []

        search_peptide_list = []

        if do_reverse == True:
            for peptide in peptide_list:
                rev_peptide = peptide[::-1]
                search_peptide_list.append(peptide)
                search_peptide_list.append(rev_peptide)
        else:
            search_peptide_list = peptide_list

        #Lets make one big string and see how it goes
        protein_sequence_list = []
        protein_position_list = []
        protein_name_list = []
        running_position = 0
        for protein in self.protein_list:
            protein_sequence_list.append(protein.sequence)
            protein_position_list.append(running_position)
            protein_name_list.append(protein)
            running_position += len(protein.sequence) + 1

        #make big string
        proteome_string = ";".join(protein_sequence_list)

        peptide_count = 0
        for peptide in search_peptide_list:
            peptide_count += 1
            print(str(peptide_count) + " of " + str(len(search_peptide_list)))

            find_position = 0
            while find_position != -1:
                find_position = proteome_string.find(peptide, find_position + 1)

                if find_position == -1:
                    break

                #Determine which protein is that position
                for i in range(len(protein_position_list)):
                    if find_position < protein_position_list[i]:
                        protein_found = protein_name_list[i-1]
                        return_proteins.append(protein_found)
                        break

        return return_proteins

    #Return a list of proteins that are covered by a set of peptides with at least K per protein
    def get_proteins_covered_by_k_peptides(self, peptide_list, minimum_peptides=1, do_reverse=False):
        return_proteins = []
        search_peptide_list = []

        if do_reverse == True:
            for peptide in peptide_list:
                rev_peptide = peptide[::-1]
                search_peptide_list.append(peptide)
                search_peptide_list.append(rev_peptide)
        else:
            search_peptide_list = peptide_list

        search_peptide_list = list(set(search_peptide_list))

        protein_name_to_number_of_peptides_count = self.get_proteins_with_number_of_peptides_covered_map(search_peptide_list)

        for protein in self.protein_list:
            number_of_peptides = protein_name_to_number_of_peptides_count[protein.protein]
            if number_of_peptides >= minimum_peptides:
                return_proteins.append(protein)

        return return_proteins

    #Given a list of peptides, returns a map of protein names to the number of peptides that cover it
    def get_proteins_with_number_of_peptides_covered_map(self, peptide_list):
        return_proteins = {}
        for protein in self.protein_list:
            number_of_peptides = protein.number_of_peptides_on_protein(peptide_list)
            return_proteins[protein.protein] = number_of_peptides
        return return_proteins

    #List of peptide sequences
    def calculate_peptide_coverage_of_all_proteins(self, peptides):
        coverage_map = {}
        count = 0
        for protein in self.protein_list:
            count += 1
            coverage_number = protein.coverage(peptides)
            coverage_map[protein.protein] = coverage_number

            #print str(count) + " of " + str(len(self.protein_list))
            #if coverage_number > 0:
            #    print protein.protein + "\t" + str(coverage_number)
        return coverage_map

    #Prints out the percetnage of tryptic peptides covered per protein
    def calculate_tryptic_peptide_coverage_of_all_proteins(self, peptides):
        coverage_map = {}
        for protein in self.protein_list:
            coverage_number = protein.get_tryptic_peptide_coverage(peptides)
            coverage_map[protein.protein] = coverage_number
            #if coverage_number > 0:
            #    print protein.protein + "\t" + str(coverage_number)
        return coverage_map

    def calculate_unique_trypic_peptide_coverage_of_all_proteins(self, peptides):
        coverage_map = {}
        count = 0
        unique_peptides = self.get_unique_peptides()
        for protein in self.protein_list:
            count += 1
            #print str(count) + " of " + str(len(self.protein_list))
            coverage_number = protein.get_unique_tryptic_peptide_coverage(peptides, unique_peptides)
            coverage_map[protein.protein] = coverage_number
            #if coverage_number > 0:
            #    print protein.protein + "\t" + str(coverage_number)
        return coverage_map

    #Returns a list of protein objects for where the peptide matches
    def get_proteins_with_sequence(self, peptide):
        protein_list = []
        for protein in self.protein_list:
            if protein.contain_peptide(peptide):
                protein_list.append(protein)
        return protein_list

    #Returns a list of all trypic peptides
    def get_tryptic_peptides(self, min_length=0):
        tryptic_peptides = []
        for protein in self.protein_list:
            tryptic_peptides += protein.get_tryptic_peptides(min_length)

        return list(set(tryptic_peptides))

    #Given a list of peptides that were identified
    #Returns the percentage covered of proteome, percentage of input peptides in proteome
    def calculate_tryptic_peptides_covered(self, covered_peptides):
        proteome_peptides = self.get_tryptic_peptides()
        #proteome_peptides_map = {}
        #output_proteome_file = open("output_proteome.out", "w")
        #for peptide in proteome_peptides:
        #    proteome_peptides_map[peptide] = 1
        #    output_proteome_file.write(peptide + "\n")

        #intersection_set = []
        #for peptide in covered_peptides:
        #    if peptide in proteome_peptides_map:
        #        intersection_set.append(peptide)

        intersection_set = (set(covered_peptides)).intersection(set(proteome_peptides))
        return float(len(intersection_set))/float(len(proteome_peptides)), float(len(intersection_set))/float(len(covered_peptides))

    def set_coverage_with_proteome(self, covered_peptides):
        proteome_peptides = self.get_tryptic_peptides()
        intersection_set = list((set(covered_peptides)).intersection(set(proteome_peptides)))
        cover_exclusive = list((set(covered_peptides)).difference(set(proteome_peptides)))
        proteome_exclusive = list((set(proteome_peptides)).difference(set(covered_peptides)))
        return intersection_set, proteome_exclusive, cover_exclusive

    #Give set of unique peptides that are unique to proteins
    def get_unique_peptides(self, min_length=0):
        peptide_sharing_count = {}
        pattern = re.compile('[KR][^P]')
        for protein in self.protein_list:
            tryptic_peptides = list(set(protein.get_tryptic_peptides(min_length)))
            for peptide in tryptic_peptides:
                if not peptide in peptide_sharing_count:
                    peptide_sharing_count[peptide] = 0
                peptide_sharing_count[peptide]+= 1

        unique_peptides = []
        for peptide in peptide_sharing_count:
            if peptide_sharing_count[peptide] == 1:
                unique_peptides.append(peptide)

        return unique_peptides

    def calculate_unique_peptide_coverage(self, covered_peptides):
        unique_peptides = self.get_unique_peptides()
        intersection_set = (set(covered_peptides)).intersection(set(unique_peptides))
        return float(len(intersection_set))/float(len(unique_peptides))

    #Get Genes Covered
    def get_genes_covered_by_peptides(self, peptides):
        genes_covered = []

        for protein in self.protein_list:
            coverage_number = protein.coverage(peptides)
            if coverage_number > 0:
                genes_covered.append(protein.gene_name)
        return list(set(genes_covered))

    #Get All Genes
    def get_all_genes(self):
        all_genes = []
        for protein in self.protein_list:
            all_genes.append(protein.gene_name)
        return list(set(all_genes))

    #Returns a list of peptides that uniquely identify a gene
    def get_unique_gene_peptides(self):
        gene_to_peptide = {} #Map of all peptides per gene
        for protein in self.protein_list:
            if not protein.gene_name in gene_to_peptide:
                gene_to_peptide[protein.gene_name] = []
            gene_to_peptide[protein.gene_name] += protein.get_tryptic_peptides()

        #Making the peptides unique
        for gene in gene_to_peptide:
            gene_to_peptide[gene] = list(set(gene_to_peptide[gene]))

        #Count Uniqueness of peptides
        peptide_sharing_count = {}
        for gene in gene_to_peptide:
            for peptide in gene_to_peptide[gene]:
                if not peptide in peptide_sharing_count:
                    peptide_sharing_count[peptide] = 0
                peptide_sharing_count[peptide] += 1

        #Find Unique Peptides
        unique_peptides = []
        for peptide in peptide_sharing_count:
            if peptide_sharing_count[peptide] == 1:
                unique_peptides.append(peptide)

        return unique_peptides

class Protein:
    def __init__(self, protein_name, sequence, gene_name):
        self.protein = protein_name
        self.sequence = sequence
        self.peptides = []
        self.gene_name = gene_name


    #Returns the coverage of protein, given peptides
    def coverage(self, peptides):
        coverage_list = self.coverage_by_amino_acids(peptides)
        return float(sum(coverage_list))/float(len(coverage_list))

    def coverage_by_amino_acids(self, peptides):
        coverage_list = [0] * len(self.sequence)

        if len(coverage_list) == 0:
            print(self.protein)
            return 0

        substring_map = {}
        #print "Creating Map " + self.protein + " " + str(len(self.sequence))


        #for peptide in tryptic_peptides:

        #for i in range(len(self.sequence)):
        #    for j in range(len(self.sequence)):
        #        if j < i or j - i > 20:
        #            continue
        #        substring = self.sequence[i:j+1]
        #        if not substring in substring_map:
        #            substring_map[substring] = []
        #        substring_map[substring].append(i)

        #for peptide in peptides:
        #    if peptide in substring_map:
        #        found_locations = substring_map[peptide]
                #for found_location in found_locations:
                #    for i in range(found_location, found_location + len(peptide)):
                #        coverage_list[i] = 1


        tryptic_peptides = set(self.get_tryptic_peptides())
        for peptide in peptides:
            #if not peptide in tryptic_peptides:
            #    continue

            starting_location = 0
            found_location = self.sequence.find(peptide, starting_location)
            #print self.sequence
            while found_location != -1:
                #coverage_list[found_location:found_location + len(self.sequence)] = 1
                for i in range(found_location, found_location + len(peptide)):
                    coverage_list[i] = 1

                starting_location = found_location + 1
                found_location = self.sequence.find(peptide, starting_location)

        return coverage_list


    #returns percentage of tryptic peptides covered in protein
    def get_tryptic_peptide_coverage(self, peptides):
        tryptic_peptides = set(self.get_tryptic_peptides())
        found_tryptic_peptides = set()

        for peptide in peptides:
            if peptide in tryptic_peptides:
                found_tryptic_peptides.add(peptide)

        return float(len(found_tryptic_peptides))/float(len(tryptic_peptides))

    def get_unique_tryptic_peptide_coverage(self, peptides, unique_peptides):
        tryptic_peptides = set(self.get_tryptic_peptides())
        unique_tryptic_peptides = set(unique_peptides).intersection(tryptic_peptides)
        found_tryptic_peptides = set()

        if(len(unique_tryptic_peptides) == 0):
            return -1

        for peptide in peptides:
            if peptide in unique_tryptic_peptides:
                found_tryptic_peptides.add(peptide)

        return float(len(found_tryptic_peptides))/float(len(unique_tryptic_peptides))

    def contain_peptide(self, peptide):
        if peptide in self.sequence:
            return True
        else:
            return False

    def number_of_peptides_on_protein(self, peptides):
        peptide_map_count = 0
        for peptide in peptides:
            if self.contain_peptide(peptide):
                peptide_map_count += 1
        return peptide_map_count

    #Tryptic Peptides for protein
    def get_tryptic_peptides(self, min_length=0):
        tryptic_peptides = re.findall(r".(?:(?<![KR](?!P)).)*", self.sequence)
        return_list = []
        for peptide in tryptic_peptides:
            if len(peptide) > min_length:
                return_list.append(peptide)
        return return_list

    def to_fasta(self):
        output_string = ">" + self.protein +  " " + "GN=" + self.gene_name + "\n"
        output_string += self.sequence.strip()
        return output_string
