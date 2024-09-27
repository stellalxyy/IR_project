#!/usr/bin/python
# -*- encoding: utf-8 -*-
"""
Created by Stella Li on 2024/08/30 12:00
"""

import os
import logging
import itertools
import subprocess
import pandas as pd

from pathlib import Path
from stella_config import Config
from RNA_seq_pipeline import TASK_ID

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

config = Config(TASK_ID)

class RunningDifferentialAnalysis:
    def __init__(self):
        self.group_info_file_path = config.HOST_TEMP_GTF_GROUP_INFO
        self.expression_matrix_path_param = config.HOST_TEMP_EXPRESSION_MATRIX
        self.output_path_param = config.HOST_OUTPUT_DIFFERENTIAL_ANALYSIS_DIR
        self.group_info_params = self.produce_group_params()
        self.r_script_path = '/data/youpu/Stella_IR_Project/utils/differential_analysis.r'

    def produce_group_params(self):
        df = pd.read_csv(self.group_info_file_path, sep='\t', header=0, index_col=0)
        groups = df['group'].unique().tolist()

        comb_list = list(itertools.combinations(groups, 2))
        group_info_params = [f"{i[0]}-{i[1]}" for i in comb_list]
        return group_info_params

    def run(self):
        os.makedirs(Path(self.output_path_param), exist_ok=True)

        for group_pair in self.group_info_params:
            logging.info(f"Start running differential analysis with TASK ID {TASK_ID} and GROUP INFO {group_pair}")
            output_file_name = group_pair + '.csv'
            cmd = [
                'Rscript',
                self.r_script_path,
                self.expression_matrix_path_param,
                self.group_info_file_path,
                os.path.join(self.output_path_param, output_file_name),
                group_pair
            ]
            logging.info(f"Running command: {' '.join(cmd)}")
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode == 0:
                logging.info(f"Command succeeded with output: {result.stdout}")
            else:
                logging.error(f"Command failed with error: {result.stderr}")

        logging.info(f"Finish running differential analysis with TASK ID {TASK_ID}")


if __name__ == '__main__':
    analysis = RunningDifferentialAnalysis()
    analysis.run()
