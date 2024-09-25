#!/usr/bin/python
# -*- encoding: utf-8 -*-
'''
Created by Stella Li on 2024/09/25 14:04
'''

import os
import logging
import pandas as pd

from stella_config import Config
from RNA_seq_pipeline import TASK_ID

config = Config(TASK_ID)

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class ConstructGroupInfoMatrix:
    def __init__(self, group_file_path: str, trans_folder_path: str, output_dir: str):
        self.group_file_path = group_file_path
        self.trans_folder_path = trans_folder_path
        self.output_dir = output_dir
        self.df = self.load_group_info()
        self.trans_name_list = self.load_transcript_names()
        
    def load_group_info(self):
        df = pd.read_csv(self.group_file_path, header=None, names=['sample', 'group'])
        df['sample'] = df['sample'].apply(lambda x: x.split('.')[0])
        return df

    def load_transcript_names(self):
        return [n.split('.')[0] for n in os.listdir(self.trans_folder_path) if n.startswith('SRR') and n.endswith('.gtf')]

    def filter_invalid_samples(self):
        valid_samples = self.df['sample'].isin(self.trans_name_list)
        return self.df[valid_samples]

    def main(self):
        logging.info(f"Start updating input group information matrix with TASK ID {TASK_ID}")
        self.df = self.filter_invalid_samples()
        output_path = os.path.join(self.output_dir, 'transcript_group_info.csv')
        self.df.to_csv(output_path, sep='\t', index=False, header=True)
        logging.info(f"Finish updating input group information matrix with TASK ID {TASK_ID}")

if __name__ == '__main__':
    group_file_path = '/data/youpu/Stella_IR_Project/data/2024-09-12_01:57:43/group_info/gtf_groups.csv'
    trans_folder_path = '/data/youpu/Stella_IR_Project/temp/2024-09-12_01:57:43/transcripts'
    output_dir = '/data/youpu/Stella_IR_Project/utils'
    
    matrix_constructor = ConstructGroupInfoMatrix(group_file_path, trans_folder_path, output_dir)
    matrix_constructor.main()