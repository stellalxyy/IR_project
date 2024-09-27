#!/usr/bin/python
# -*- encoding: utf-8 -*-
"""
Created by Stella Li on 2024/08/08 13:59
"""

import datetime


def generate_task_id():
    dt = datetime.datetime.now().strftime('%Y-%m-%d_%H:%M:%S')
    return dt

TASK_ID = generate_task_id()

import os
import merge
import shutil
import logging
import predict
import mapping
import argparse
import assemble
import exec_maxquant
import sort_bam_files
import pandas as pd

from tqdm import tqdm
from pathlib import Path
from stella_config import Config
from XML_processor import XMLProcessor
from construct_group_info_matrix import ConstructGroupInfoMatrix
from construct_expression_matrix import ConstructExpressionMatrix
from running_differential_analysis import RunningDifferentialAnalysis
from running_deseq2_differential_analysis import RunningDeseq2Analysis

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

config = Config(TASK_ID)


class RNASeqPipeline:
    def __init__(self):
        self.parse_args()
        self.fastaq_dir = self.args.fastaq_dir
        self.raw_file_path = self.args.raw_file_path
        self.gtf_group_info_dir = self.args.gtf_group_info_dir
        self.raw_group_info_dir = self.args.raw_group_info_dir
        self.main_dir = config.HOST_MAIN_FOLDER
        self.data_dir = config.HOST_DATA_SUBFOLDER
        self.temp_dir = config.HOST_TEMP_SUBFOLDER
        self.output_dir = config.HOST_OUTPUT_SUBFOLDER
        self.scripts_dir = config.HOST_SCRIPT_SUBFOLDER
        self.temp_task_dir = os.path.join(self.temp_dir, TASK_ID)
        self.output_task_dir = os.path.join(self.output_dir, TASK_ID)
        self.data_task_dir = os.path.join(self.data_dir, TASK_ID)
        self.data_task_fastq_dir = os.path.join(self.data_task_dir, 'fastq')
        self.data_task_raw_files_dir = os.path.join(self.data_task_dir, 'raw_files')
        self.data_task_group_info_dir = os.path.join(self.data_task_dir, 'group_info')
        self.data_task_gtf_group_info = os.path.join(self.data_task_group_info_dir, 'gtf_groups.csv')
        self.data_task_raw_group_info = os.path.join(self.data_task_group_info_dir, 'raw_groups.csv')
        self.temp_task_fasta_dir = os.path.join(self.temp_task_dir, 'fasta')
        self.temp_task_mapping_dir = os.path.join(self.temp_task_dir, 'mapping')
        self.temp_task_sorted_bam_dir = os.path.join(self.temp_task_dir, 'sorted_bam_files')
        self.bam_file_storage = os.path.join(self.main_dir, 'all_bam_files')
        self.docker_fastq_dir = '/root' + self.data_task_fastq_dir

    def parse_args(self):
        parser = argparse.ArgumentParser(description='Stella IR PROJECT RNA-seq pipeline')

        parser.add_argument('--fastaq_dir', '-f',
                            type=str,
                            help='Path to the directory containing the fastq files',
                            required=True)
        parser.add_argument('--raw_file_path', '-r',
                            type=str,
                            help='Path to the directory containing the raw files',
                            required=True)
        parser.add_argument('--gtf_group_info_dir', '-gg',
                            type=str,
                            help='Path to the directory containing the gtf group information',
                            required=True)
        parser.add_argument('--raw_group_info_dir', '-rg',
                            type=str,
                            help='Path to the directory containing the raw group information',
                            required=True)
        self.args = parser.parse_args()

    def copy_file(self, file: str, data_type: str):
        _, filename = os.path.split(file)
        if data_type == 'fastq':
            logging.info(f'Copying {file} to {os.path.join(self.data_task_fastq_dir, filename)}')
            shutil.copy(file, os.path.join(self.data_task_fastq_dir, filename))
        elif data_type == 'raw_files':
            logging.info(f'Copying {file} to {os.path.join(self.data_task_raw_files_dir, filename)}')
            shutil.copy(file, os.path.join(self.data_task_raw_files_dir, filename))
        elif data_type == 'gtf':
            logging.info(f'Copying {file} to {os.path.join(self.data_task_group_info_dir, 'gtf_groups.csv')}')
            shutil.copy(file, os.path.join(self.data_task_group_info_dir, 'gtf_groups.csv'))
        else:
            logging.info(f'Copying {file} to {os.path.join(self.data_task_group_info_dir, 'raw_groups.csv')}')
            shutil.copy(file, os.path.join(self.data_task_group_info_dir, 'raw_groups.csv'))
            
    def _copy(self, source_dir: str, char: str, data_type: str):
        filename_list = os.listdir(source_dir)
        for f in filename_list:
            if f.endswith(char):
                before_file = os.path.join(source_dir, f)
                self.copy_file(before_file, data_type)
            pass
        
    def catch_input_idx(self):
        df = pd.read_csv(self.data_task_gtf_group_info, sep = ',', header = None)
        sample_list = df.iloc[:,0].tolist()
        if not sample_list[0].endswith('.raw'):
            df.drop([0], axis=0, inplace=True)
        pass
        df.columns = ['sample', 'group']
        df.reset_index(inplace=True)
        df.drop(columns = ['index'], inplace=True)
        sample_list = df['sample'].tolist()
        idx_list = [i.split('.')[0] for i in sample_list]
        df['sample'] = idx_list
        return idx_list
    
    def match_sample_id(self, idx_list):
        saved_sample_list = [i.split('_')[0].replace('Align', '') for i in os.listdir(self.bam_file_storage)]
        new_idx_list = [i for i in idx_list if i not in saved_sample_list]
        return new_idx_list
    
    def single_pair_judgement(self):
        fastq_file_name_list = os.listdir(self.data_task_fastq_dir)
        has_single = False
        has_paired = False

        for f in fastq_file_name_list:
            if f.endswith('_1.fastq'):
                has_paired = True
            elif f.endswith('_2.fastq'):
                has_paired = True
            else:
                has_single = True

        if has_single and has_paired:
            return 'mixed'
        elif has_single:
            return 'single'
        elif has_paired:
            return 'paired'
        else:
            return 'none'
        
    def mapping(self, new_idx_list, form):
        logging.info(f"Task ID {TASK_ID}: Genome mapping")
        cmd_cf = f'docker exec -i {config.YOU_IR_CONTAINER} mkdir -p {config.TEMP_MAPPING_DIR}'
        os.system(cmd_cf)
        num_new_idx = len(new_idx_list)
        for i, idx in enumerate(tqdm(new_idx_list, desc="Mapping", unit="sample")):
            if num_new_idx == 1:
                mapping.main(idx, self.docker_fastq_dir, form, "last")
            else:
                if i == 0:
                    mapping.main(idx, self.docker_fastq_dir, form, "begin")
                elif i == len(new_idx_list) - 1:
                    mapping.main(idx, self.docker_fastq_dir, form, "last")
                else:
                    mapping.main(idx, self.docker_fastq_dir, form, 'stella')
        logging.info(f"Task ID {TASK_ID}: Genome mapping completed")

    @staticmethod
    def sort_bam_files(new_idx_list):
        logging.info(f"Task ID {TASK_ID}: Sorting BAM files")
        num_new_idx = len(new_idx_list)
        for i, idx in enumerate(tqdm(new_idx_list, desc="Sorting BAM files", unit="file")):
            if num_new_idx == 1:
                sort_bam_files.main(idx, "last")
            else:
                if i == 0:
                    sort_bam_files.main(idx, 'begin')
                elif i == len(new_idx_list) - 1:
                    sort_bam_files.main(idx, 'last')
                else:
                    sort_bam_files.main(idx, 'stella')
        logging.info(f"Task ID {TASK_ID}: Sorting BAM files completed")
        
    def save_bam_files(self):
        logging.info(f"TASK ID {TASK_ID}: Saving BAM files")
        file_name_list = os.listdir(self.temp_task_mapping_dir)
        bam_file_list = [i for i in file_name_list if i.endswith("_sort.bam")]
        for i in bam_file_list:
            shutil.copy(os.path.join(self.temp_task_mapping_dir, i), os.path.join(self.bam_file_storage, i))
        logging.info(f"TASK ID {TASK_ID}: Saving BAM files completed")

    def find_bad_input_sample(self, idx_list):
        saved_sample_list = [i.split('_')[0].replace('Align', '') for i in os.listdir(self.bam_file_storage)]
        bad_input = [i for i in idx_list if i not in saved_sample_list]
        final_idx_list = [i for i in idx_list if i not in bad_input]
        return final_idx_list, bad_input
        
    @staticmethod
    def assemble(final_idx_list):
        logging.info(f"Task ID {TASK_ID}: Assembling transcripts")
        cmd_cf = f'docker exec -i {config.YOU_IR_CONTAINER} mkdir -p {config.TEMP_TRANSCRIPTS_DIR}'
        os.system(cmd_cf)
        for i, idx in enumerate(tqdm(final_idx_list, desc="Assembling", unit="sample")):
            if i == 0:
                assemble.main(idx, 'begin')
            elif i == len(final_idx_list) - 1:
                assemble.main(idx, 'last')
            else:
                assemble.main(idx, 'stella')
        logging.info(f"Task ID {TASK_ID}: Assembling transcripts completed")

    def differential_analysis(self):
        logging.info(f"Task ID {TASK_ID}: Starting differential analysis")

        logging.info("Start updating transcripts samples group info...")
        group_info = ConstructGroupInfoMatrix(self.data_task_gtf_group_info, 
                                              config.HOST_TEMP_TRANSCRIPTS_DIR,
                                              self.temp_task_dir)
        group_info.main()
        
        logging.info("Start constructing expression matrix...")
        matrix_constructor = ConstructExpressionMatrix(os.path.join(self.temp_task_dir, 'transcript_group_info.csv'),
                                                       config.HOST_TEMP_TRANSCRIPTS_DIR, 
                                                       self.temp_task_dir)
        matrix_constructor.main()

        logging.info("Start running limma differential analysis...")
        analysis = RunningDifferentialAnalysis()
        analysis.run()
        
        logging.info("Start running deseq2 differential analysis...")
        deseq_analysis = RunningDeseq2Analysis()
        deseq_analysis.main()

        logging.info(f"Task ID {TASK_ID}: Finish differential analysis")

    @staticmethod
    def merge():
        logging.info(f"Task ID {TASK_ID}: Merging transcripts")
        merge.main()
        logging.info(f"Task ID {TASK_ID}: Merging transcripts completed")

    def predict(self):
        logging.info(f"Task ID {TASK_ID}: Predicting disease-specific coding sequences")
        predict.main()
        shutil.copy(os.path.join(self.temp_task_fasta_dir, 'transcript.fasta.transdecoder.pep'),
                    os.path.join(self.output_task_dir, 'transcript.fasta.transdecoder.pep'))
        logging.info(f"Task ID {TASK_ID}: Predicting disease-specific coding sequences completed")

    def exec_maxquant(self):
        logging.info(f"Task ID {TASK_ID}: Executing MaxQuant")

        logging.info('Start writing XML file...')
        processor = XMLProcessor(config.TEMP_MONO_OUTPUT_TRANSCRIPT_FASTA,
                                 self.data_task_raw_files_dir,
                                 self.data_task_raw_group_info)
        processor.run()

        try:
            exec_maxquant.main()
        except Exception as e:
            print(f"Error: {e}")

        logging.info(f"Task ID {TASK_ID}: Executing MaxQuant completed")

    def run(self):
        os.makedirs(Path(self.bam_file_storage), exist_ok=True)
        os.makedirs(Path(self.temp_task_dir), exist_ok=True)
        os.makedirs(Path(self.data_task_dir), exist_ok=True)
        os.makedirs(Path(self.output_task_dir), exist_ok=True)
        os.makedirs(Path(self.data_task_fastq_dir), exist_ok=True)
        os.makedirs(Path(self.data_task_raw_files_dir), exist_ok=True)
        os.makedirs(Path(self.data_task_group_info_dir), exist_ok=True)
        os.makedirs(Path(self.temp_task_sorted_bam_dir), exist_ok=True)
        
        logging.info(f'Copying fastq files to data folder...')
        self._copy(self.fastaq_dir, '.fastq', 'fastq')
        logging.info(f'Copying fastq files done.')

        logging.info(f'Copying raw files to data folder...')
        self._copy(self.raw_file_path, '.raw', 'raw_files')
        logging.info(f'Copying raw files done.')

        logging.info(f'Copying gtf group info to data folder...')
        self._copy(self.gtf_group_info_dir, '.csv', 'gtf')
        logging.info(f'Copying gtf group info done.')
        
        logging.info(f'Copying raw group info to data folder...')
        self._copy(self.raw_group_info_dir, '.csv', 'raw')
        logging.info(f'Copying raw group info done.')

        logging.info(f'Start of program with Task ID {TASK_ID}')
        idx_list = self.catch_input_idx()
        logging.info(f'Input idx list: {idx_list}')
        logging.info(f'length of input idx list: {len(idx_list)}')

        logging.info('Start sample id matching...')
        new_idx_list = self.match_sample_id(idx_list)
        if len(new_idx_list) == 0:
            logging.info('New sample id not find')
        else:
            logging.info(f'New sample id: {new_idx_list}')
        
        logging.info('Start judging single-end mapping or paired-end mapping...')
        form = self.single_pair_judgement()
        logging.info(f'FORM: {form}')
        
        logging.info('Start genome mapping...')
        if len(new_idx_list) != 0:
            self.mapping(new_idx_list, form)
        else:
            logging.info('New sample id not find. Skip mapping process.')

        logging.info('Start sorting bam files...')
        if len(new_idx_list) != 0:
            self.sort_bam_files(new_idx_list)
        else:
            logging.info('New sample id not find. Skip sorting process.')

        logging.info('Start saving bam files...')
        if len(new_idx_list) != 0:
            self.save_bam_files()
        else:
            logging.info('New sample id not find. Skip saving bam files process.')
        
        logging.info('Start updating sample idx list...')
        final_idx_list, bad_input = self.find_bad_input_sample(idx_list)
        logging.info(f'Final sample idx list: {final_idx_list}')
        logging.info(f'Invalid input sample id list: {bad_input}. Please check your input data.')
            
        logging.info('Start assembling transcripts...')
        self.assemble(final_idx_list)

        logging.info('Start differential analysis...')
        self.differential_analysis()

        logging.info('Start merging transcripts...')
        with tqdm(total=1, desc="Merging transcripts") as pbar:
            self.merge()
            pbar.update(1)

        logging.info('Start predicting disease-specific coding sequences...')
        with tqdm(total=10, desc="Predicting coding sequences") as pbar:
            self.predict()
            pbar.update(10)

        logging.info('Start executing Maxquant...')
        with tqdm(total=10, desc="Executing Maxquant") as pbar:
            self.exec_maxquant()
            pbar.update(10)

        logging.info('Copying results...')
        cmd = "cp -r " + os.path.join(self.data_task_raw_files_dir, 'combined') + ' ' + self.output_task_dir
        logging.info(cmd)
        os.system(cmd)

        logging.info(f'End of program with Task ID {TASK_ID}')


if __name__ == '__main__':
    pipeline = RNASeqPipeline()
    pipeline.run()
