#!/usr/bin/python
# -*- encoding: utf-8 -*-
'''
Created by Stella Li on 2024/09/26 17:21
'''

import os
import logging
import subprocess
import pandas as pd

from stella_config import Config
from RNA_seq_pipeline import TASK_ID

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

config = Config(TASK_ID)

class RunningDeseq2Analysis:
    def __init__(self):
        self.group_info_file_path = config.HOST_TEMP_GTF_GROUP_INFO
        self.annotation_file_path = config.INTRON_GTF_FILE_HOST_PATH
        self.output_counts_file_path = os.path.join(config.HOST_TEMP_TASK_SUBFOLDER, 'counts.txt')
        self.all_bam_files_path = config.HOST_ALL_BAM_FILES_DIR
        self.featurecounts = config.FEATURECOUNTS
        self.temp_task_sorted_bam_dir = config.HOST_TEMP_SORTED_BAM_DIR
        self.r_script_path = '/data/youpu/Stella_IR_Project/utils/deseq2_differential_analysis.r'
        self.deseq2_output_file_path = os.path.join(config.HOST_OUTPUT_DIFFERENTIAL_ANALYSIS_DIR, 'deseq2_analysis.csv')
        self.counts_expression_matrix_path = os.path.join(config.HOST_TEMP_TASK_SUBFOLDER, 'counts_expression_matrix.csv')
    
    def get_sorted_bam(self):
        df = pd.read_csv(self.group_info_file_path, sep='\t', header=0, index_col=None)
        sample_name_list = df['sample'].tolist()
        sorted_bam_files_list = []
        for i in os.listdir(self.all_bam_files_path):
            sample_name = i.split("Align_sort.bam")[0]
            if sample_name in sample_name_list:
                sorted_bam_files_list.append(i)
        for j in sorted_bam_files_list:
            sample = j.split("Align_sort.bam")[0]
            origin_bam_path = os.path.join(self.all_bam_files_path, j)
            new_bam_path = os.path.join(self.temp_task_sorted_bam_dir, sample+'.bam')
            command = f"cp {origin_bam_path} {new_bam_path}"
            subprocess.run(command, shell=True, check=True)
        return sample_name_list
    
    def generate_counts(self, sample_name_list: list):
        bam_file_list = []
        for i in sample_name_list:
            bam_file_path = os.path.join(self.temp_task_sorted_bam_dir, i+'.bam')
            bam_file_list.append(bam_file_path)
        bam_files_str = ' '.join(bam_file_list)
        command = f'{self.featurecounts} -a {self.annotation_file_path} -o {self.output_counts_file_path} {bam_files_str}'
        subprocess.run(command, shell=True, check=True)
        
    def process_counts_file(self, sample_name_list: list):
        counts_df = pd.read_csv(self.output_counts_file_path, sep='\t', header=1)
        sample_list_length = len(sample_name_list)
        selected_columns = pd.concat([counts_df.iloc[:, [0]], counts_df.iloc[:, -sample_list_length:]], axis=1)
        selected_columns.set_index(selected_columns.columns[0], inplace=True)
        selected_columns.index.name = None
        selected_columns.columns = sample_name_list
        selected_columns.to_csv(self.counts_expression_matrix_path, sep='\t', header=True, index=True)
        
    def main(self):
        logging.info(f"Start running deseq2 differential analysis with TASK ID {TASK_ID}")
        
        sample_name_list = self.get_sorted_bam()
        self.generate_counts(sample_name_list)
        self.process_counts_file(sample_name_list)
        cmd = [
                'Rscript',
                self.r_script_path,
                self.counts_expression_matrix_path,
                self.group_info_file_path,
                self.deseq2_output_file_path
            ]
        result = subprocess.run(cmd)
        if result.returncode == 0:
            logging.info(f"Command succeeded with output: {result.stdout}")
        else:
            logging.error(f"Command failed with error: {result.stderr}")
        
        logging.info(f"Finish running deseq2 differential analysis with TASK ID {TASK_ID}")
        
        
if __name__ == '__main__':
    deseq_analysis = RunningDeseq2Analysis()
    deseq_analysis.main()
        