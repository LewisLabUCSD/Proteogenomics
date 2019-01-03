* MSGF_Enosi: This folder has the pipeline developed by Dr. Vineet Bafna's lab. It processes raw mass spectrum and generates matched peptides.
The rests are custom codes of processing the data.
* Annotation_stringtie_pasa.ipynb: generate a draft annotation using stringtie and PASA.
* Refseq_picr_pr_cmp.ipynb: compare the pr sequence between refseq and maddy's annotation.
* Ribotaper_result_analysis.ipynb: analyze CDS predicted by ribotaper.
* Update_other_events.ipynb: analyze the novel translational events.
* Update_transcript_CDS_event.ipynb: integrate novel transcript_CDS_evnet to draft annotation.
* Virus_analysis.ipynb: Detect retrovirus in CHO cells using proteomics.
* db_01_copydb.py: generate commands to run MSGF in parallel.
* db_01_splicedb.py: create peptide database using spliced RNA-Seq reads.
* db_01_snpdb.py: create peptide database using snps called from RNA-Seq reads.
* draft_refseq_cmp.ipynb: compare draft annotation with refseq protein sequences.
* enosi_local.py: run the pipeline part (FDR correction, peptide event calling in local machine)
* get_perfect_pr_map.py: compare protein sequence from different genome assemblies.
* gff_statistics.ipynb: some statistics function for parsing gff file.

p01~p08: code to run the proteogenomics pipeline.
m01_parse_event.py: code to process the event results of the pipeline. Here each event represents a call (eg:new splice, new gene) based on the identified peptides.
