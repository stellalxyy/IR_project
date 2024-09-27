#!/usr/bin/python
# -*- encoding: utf-8 -*-
"""
Created by Stella Li on 2024/08/07 18:49
"""
import os


class Config:
    def __init__(self, task_id: str):
        self.TASK_ID = task_id

        # data
        self.GENOME_DIR = '/root/genome_dir'
        self.GENOME_FASTA = '/root/data/original/genome/Homo/Homo_sapiens.GRCh38.dna.primary_assembly.fa'
        self.INTRON_GTF_FILE_PATH = '/root/data/original/genome/Homo/Homo_sapiens.GRCh38.107.intron.all.nd.id.gtf'
        self.INTRON_GTF_FILE_HOST_PATH = '/data/youpu/Homo_sapiens.GRCh38.107.intron.all.nd.id.gtf'
        
        # software
        self.SAMTOOLS = '/usr/local/bin/samtools'
        self.STAR = '/root/data/youpu/software/STAR-2.7.11a/bin/Linux_x86_64/STAR'
        self.STRINGTIE = '/root/software/stringtie/stringtie'
        self.TRANSDECODER = '/root/software/TransDecoder-v5.7.1/util/gtf_genome_to_cdna_fasta.pl'
        self.TRANSDECODER_LONGORFS = '/root/software/TransDecoder-v5.7.1/TransDecoder.LongOrfs'
        self.TRANSDECODER_PREDICT = '/root/software/TransDecoder-v5.7.1/TransDecoder.Predict'
        self.MAXQUANT = '/root/software/MaxQuant_v_2.4.9.0/bin/MaxQuantCmd.exe'
        self.FEATURECOUNTS = '/data/youpu/miniconda3/bin/featureCounts'

        # running dir
        self.BASE_DOCKER_DIR = '/root/data/youpu'
        self.BASE_HOST_DIR = '/data/youpu'
        # docker folders
        self.DOCKER_MAIN_FOLDER = os.path.join(self.BASE_DOCKER_DIR, 'Stella_IR_Project')
        self.DOCKER_DATA_SUBFOLDER = os.path.join(self.DOCKER_MAIN_FOLDER, 'data')
        self.DOCKER_TEMP_SUBFOLDER = os.path.join(self.DOCKER_MAIN_FOLDER, 'temp')
        self.DOCKER_OUTPUT_SUBFOLDER = os.path.join(self.DOCKER_MAIN_FOLDER, 'output')
        self.DOCKER_SCRIPT_SUBFOLDER = os.path.join(self.DOCKER_MAIN_FOLDER, 'scripts')
        self.DOCKER_DATA_TASK_SUBFOLDER = os.path.join(self.DOCKER_DATA_SUBFOLDER, self.TASK_ID)
        self.DOCKER_TEMP_TASK_SUBFOLDER = os.path.join(self.DOCKER_TEMP_SUBFOLDER, self.TASK_ID)
        self.DOCKER_OUTPUT_TASK_SUBFOLDER = os.path.join(self.DOCKER_OUTPUT_SUBFOLDER, self.TASK_ID)
        # host folders
        self.HOST_MAIN_FOLDER = os.path.join(self.BASE_HOST_DIR, 'Stella_IR_Project')
        self.HOST_DATA_SUBFOLDER = os.path.join(self.HOST_MAIN_FOLDER, 'data')
        self.HOST_TEMP_SUBFOLDER = os.path.join(self.HOST_MAIN_FOLDER, 'temp')
        self.HOST_OUTPUT_SUBFOLDER = os.path.join(self.HOST_MAIN_FOLDER, 'output')
        self.HOST_SCRIPT_SUBFOLDER = os.path.join(self.HOST_MAIN_FOLDER, 'scripts')
        self.HOST_DATA_TASK_SUBFOLDER = os.path.join(self.HOST_DATA_SUBFOLDER, self.TASK_ID)
        self.HOST_TEMP_TASK_SUBFOLDER = os.path.join(self.HOST_TEMP_SUBFOLDER, self.TASK_ID)
        self.HOST_OUTPUT_TASK_SUBFOLDER = os.path.join(self.HOST_OUTPUT_SUBFOLDER, self.TASK_ID)
        self.HOST_TEMP_TRANSCRIPTS_DIR = os.path.join(self.HOST_TEMP_TASK_SUBFOLDER, 'transcripts')
        self.HOST_TEMP_SORTED_BAM_DIR = os.path.join(self.HOST_TEMP_TASK_SUBFOLDER, 'sorted_bam_files')
        self.HOST_TEMP_EXPRESSION_MATRIX = os.path.join(self.HOST_TEMP_TASK_SUBFOLDER, 'expression_matrix.csv')
        self.HOST_TEMP_GTF_GROUP_INFO = os.path.join(self.HOST_TEMP_TASK_SUBFOLDER, 'transcript_group_info.csv')
        self.HOST_OUTPUT_DIFFERENTIAL_ANALYSIS_DIR = os.path.join(self.HOST_OUTPUT_TASK_SUBFOLDER, 'differential_analysis')
        self.HOST_ALL_BAM_FILES_DIR = os.path.join(self.HOST_MAIN_FOLDER, 'all_bam_files')
        
        # self.TEMP_IR_DIR = os.path.join(self.BASE_DOCKER_DIR, self.TASK_ID)
        self.TEMP_BAM_STORAGE_DIR = os.path.join(self.DOCKER_MAIN_FOLDER, 'all_bam_files')
        self.TEMP_MAPPING_DIR = os.path.join(self.DOCKER_TEMP_TASK_SUBFOLDER, 'mapping')  # '/root/data/youpu/temp/mapping'
        self.TEMP_TRANSCRIPTS_DIR = os.path.join(self.DOCKER_TEMP_TASK_SUBFOLDER, 'transcripts')  # '/root/data/TEMP/transcripts'
        self.TEMP_GTF_LIST_DIR = os.path.join(self.DOCKER_TEMP_TASK_SUBFOLDER, 'gtf_list')  # '/root/data/TEMP/gtf_list'
        self.TEMP_OUTPUT_TRANSCRIPT_FASTA = os.path.join(self.DOCKER_TEMP_TASK_SUBFOLDER, 'fasta')  # '/root/data/TEMP/fasta'
        ## shell
        self.TEMP_MERGE_BASH_SCRIPT = os.path.join(self.DOCKER_SCRIPT_SUBFOLDER, 'merge.sh')  # '/root/temp/merge.sh'
        self.TEMP_PREDICT_BASH_SCRIPT = os.path.join(self.DOCKER_SCRIPT_SUBFOLDER, 'predict.sh')  # '/root/temp/predict.sh'
        self.TEMP_ASSEMBLE_BASH_SCRIPT = os.path.join(self.DOCKER_SCRIPT_SUBFOLDER, 'assemble.sh')  #  '/root/temp/assemble.sh'
        self.TEMP_SORT_BASH_SCRIPT = os.path.join(self.DOCKER_SCRIPT_SUBFOLDER, 'sort_bam_files.sh')  # '/root/temp/sort_bam_files.sh'
        self.TEMP_MAPPING_BASH_SCRIPT = os.path.join(self.DOCKER_SCRIPT_SUBFOLDER, 'mapping.sh')  # '/root/data/youpu/temp/mapping.sh'
        self.TEMP_MAXQUANT_BASH_SCRIPT = os.path.join(self.DOCKER_SCRIPT_SUBFOLDER, 'exec_maxquant.sh')

        ## MONO
        # self.TEMP_MONO_DIR = os.path.join(self.BASE_DOCKER_DIR, self.TASK_ID)
        self.TEMP_MONO_OUTPUT_TRANSCRIPT_FASTA = os.path.join(self.TEMP_OUTPUT_TRANSCRIPT_FASTA, 'transcript.fasta.transdecoder.pep')  # '/root/data/TEMP/fasta'
        # self.TEMP_MONO_RAW_FILE_DIR = os.path.join(self.BASE_DOCKER_DIR, 'raw_file_dir', 'DDA') # '/root/raw_file_dir/DDA/DDA'
        self.TEMP_MONO_NEW_XML_FILE = os.path.join(self.DOCKER_TEMP_TASK_SUBFOLDER, 'new_mqpar.xml')  # '/root/data/TEMP/new_mqpar.xml'

        self.ORIGIN_XML_FILE = os.path.join(self.HOST_MAIN_FOLDER, 'mqpar.xml')
        self.NEW_XML_FILE = os.path.join(self.HOST_TEMP_TASK_SUBFOLDER, 'new_mqpar.xml')

        # docker id
        self.YOU_IR_CONTAINER = 'dac5aa0958e2'
        self.MONO_CONTAINER = '8692c1d0c266'

        self.MONO = '/usr/bin/mono'


    def change(self):
        pass

    def __repr__(self):
        return f"TASKID: {self.TASK_ID}"



if __name__ == "__main__":
    c = Config("111")
    print(c)
