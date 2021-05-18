import gzip
import os
import sys
import shutil
import yaml

from Bio import SeqIO


"""
Takes as input:
    1. A receipt obtained from ENA submission tool. 
    A txt file that includes a YAML section with 

    2. A fasta file with fasta entries ids defined after the files used for the raw submission.

    3. Path to write generated manifests
    4. Path to write generated fasta files
    5. manifest template path: the manifest with the global values set (e.g COVERAGE, MINGAPLENGHT..)
"""

def get_section_string(f, start_line, end_line):
    # consume starting lines
    start_string = iter(f.readline, start_line)
    start_string = ''.join(line for line in start_string)
    # read YAML lines
    yaml_string = iter(f.readline, end_line)
    return ''.join(x for x in yaml_string)


def main():
    input_file = open(sys.argv[1])
    fasta_in = open(sys.argv[2])
    out_manifest_base = sys.argv[3] 
    out_fasta_base = sys.argv[4] 
    manifest_template = sys.argv[5]
    yaml_delimiter = 'YAML -------------\n'
    yaml_only = yaml.safe_load(get_section_string(input_file, start_line=yaml_delimiter, end_line=yaml_delimiter))
    # print(yaml_only)
    submission_tuples_list = []
    # parse the sequence IDs
    for record in SeqIO.parse(fasta_in, "fasta"):
        seq_id = record.id
        # need to map the seq ID to 1 or more seq files submitted: for single files these is the exact file name? 
        #  ... but for paired it may not be the exact same
        # ...initially I should attempt to find a 
        #  .. if this 
        # in any case, if I cant make the right match then I should do a kind of substring match

        # find the exp_alias associated with the file
        exp_alias = None
        for index,run in yaml_only['ENA_run'].items():
            if run['file_name'] == seq_id:
                ## TODO: match also cases when the seq entry name is == entry_[1|2].fastq.gz or something
                exp_alias = run['experiment_alias']
                break
        if not exp_alias:
            raise Exception("No run files match for the sequence entry {seq_id}")
        # find the sample and study for that experiment
        sample_alias = None
        study_alias = None
        for index,exp in yaml_only['ENA_experiment'].items():
            if exp['alias'] == exp_alias:
                sample_alias = exp['sample_alias']
                study_alias = exp['study_alias']
                platform = exp['platform']
                break
        if not sample_alias:
            raise Exception("No sample associated with experiment {exp_alias}")
        if not study_alias:
            raise Exception("No study associated with experiment {exp_alias}")


        # and finally create a fasta file for each sequence (e.g named with the seq id or the run ID)
        fasta_path = os.path.join(out_fasta_base, seq_id + '.fasta')
        with open(fasta_path, "w") as output_handle:
            SeqIO.write([record], output_handle, "fasta")
        #gzip the file (required by ENA upload tool)
        fasta_path_gz = fasta_path + '.gz'
        with open(fasta_path, 'rb') as f_in:
            with gzip.open(fasta_path_gz, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        # create the manifest
        # add to the manifest the: 
        # 
        manifest_path = os.path.join(out_manifest_base, seq_id + '.manifest.txt')
        with open(manifest_path, "w") as output_handle:
            # first dump the contents of manifest template
            # containing the global vars
            with open(manifest_template) as m_template:
                output_handle.write(m_template.read())
            output_handle.write("ASSEMBLYNAME\tconsensus_" + seq_id + "\n")
            output_handle.write("PLATFORM\t" + platform + "\n")
            output_handle.write("STUDY\t" + study_alias + "\n")
            output_handle.write("SAMPLE\t" + sample_alias + "\n")
            output_handle.write("FASTA\t" + fasta_path_gz + "\n")

        # ... and a dict  (or tuple list???) that contains for each study - sample  the name of the file that has the consensus sequence
        # ****  is it ok to use the unique ids of the study and sample in the manifest?? or should I use the accessions??
        # in the latest case then I also need to parse the  Study accession details: and Sample accession details: entries
        # samples_dir[study][sample] = seq_id + '.fasta'
        submission_tuples_list.append((manifest_path, fasta_path))

    with open('submit_list.tab', "w") as output_handle:
        for submit_tuple in submission_tuples_list:
            output_handle.write('\t'.join(submit_tuple) + '\n')
    ## DEBUG CASE
    #study details
    # start_study = 'Study accession details:\n'
    # empty_end = '\n'
    # study_data = get_section_string(input_file, start_line=start_study, end_line=empty_end)
    # if len(study_data.split('\n')) > 2:
        # # more than 1 study accession
        # raise Exception("Multiple study accessions found")
    # out_manifest.write(f'STUDY\t{study_data.split()[1]}\n')
    # start_sample = 'Sample accession details:\n'
    # sample_data = get_section_string(input_file, start_line=start_sample, end_line=empty_end)
    # if len(sample_data.split('\n')) > 2:
        # # more than 1 study accession
        # raise Exception("Multiple sample accessions found")
    # out_manifest.write(f'SAMPLE\t{sample_data.split()[1]}\n')
    # platform = 'Ion Torrent'
    # out_manifest.write(f"PLATFORM\t{platform}\n")
    # out_manifest.close()

if __name__ == '__main__':
    main()
