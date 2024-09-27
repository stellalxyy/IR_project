#!/usr/bin/python
# -*- encoding: utf-8 -*-
"""
Created by Stella Li on 2024/08/28 10:48
"""

import os
import logging
import pandas as pd

from gtfparse import read_gtf

from stella_config import Config
from RNA_seq_pipeline import TASK_ID

config = Config(TASK_ID)

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class ConstructExpressionMatrix:
    def __init__(self, group_info_path: str, trans_folder_path: str, output_dir: str):
        self.group_info_path = group_info_path
        self.trans_folder_path = trans_folder_path
        self.sample_name_list = self.extract_sample_name()
        self.output_dir = output_dir

    def extract_sample_name(self):
        group_info_df = pd.read_csv(self.group_info_path, sep='\t', header=0, index_col=None)
        sample_name_list = group_info_df['sample'].tolist()
        return sample_name_list

    def create_trans_df(self, sample_file: str):
        gtf_file_path = os.path.join(self.trans_folder_path, sample_file)
        old_df = read_gtf(gtf_file_path, column_converters={"TPM": float})
        new_df = pd.DataFrame(old_df, columns=old_df.columns)
        return new_df

    @staticmethod
    def extract_trans_id(new_df: pd.DataFrame):
        transcript_id = new_df['transcript_id'].to_list()
        return transcript_id

    @staticmethod
    def remove_duplicates(lst: list):
        seen = set()
        result = []
        for item in lst:
            if item not in seen:
                seen.add(item)
                result.append(item)
        return result

    def extract_unique_trans_id(self, new_df: pd.DataFrame):
        transcript_id = self.extract_trans_id(new_df)
        unique_trans_id = self.remove_duplicates(transcript_id)
        return unique_trans_id

    @staticmethod
    def extract_tpm(new_df: pd.DataFrame):
        tpm = list(new_df['TPM'].dropna())
        return tpm

    @staticmethod
    def create_trans_tpm_dict(transcript_id: list, tpm: list):
        di = {}
        for i in zip(transcript_id, tpm):
            k = i[0]
            v = i[1]
            di[k] = v
        return di

    def create_combined_list(self):
        combined_list = []
        for sample_name in self.sample_name_list:
            sample_file = sample_name + '.transcripts.gtf'
            new_df = self.create_trans_df(sample_file)
            transcript_id = self.extract_unique_trans_id(new_df)
            combined_list.append(transcript_id)
        return combined_list

    @staticmethod
    def extract_common_trans(combined_list: list):
        common_elms = set(combined_list[0])
        for lst in combined_list[1:]:
            common_elms.intersection_update(lst)
        common_elms = list(common_elms)
        return common_elms

    @staticmethod
    def create_sample_tpm_pairs(common_elms: list, trans_tpm_dict: dict):
        new_tpm = []
        for elm in common_elms:
            tpm = float(trans_tpm_dict[elm])
            new_tpm.append(tpm)
        return new_tpm

    def create_expression_matrix(self, t_list: list, common_elms: list):
        tpm_dict = {}
        for t in t_list:
            sample_name = t[0]
            tpm = t[1]
            tpm_dict[sample_name] = tpm
        new_tpm_dict = {key: tpm_dict[key] for key in self.sample_name_list}
        expression_matrix = pd.DataFrame(new_tpm_dict, index=common_elms)
        return expression_matrix

    def main(self):
        logging.info(f"Start constructing expression matrix with Task ID {TASK_ID}")

        combined_trans_list = self.create_combined_list()
        common_elms = self.extract_common_trans(combined_trans_list)

        t_list = []
        for sample_name in self.sample_name_list:
            sample_file = sample_name + '.transcripts.gtf'
            new_df = self.create_trans_df(sample_file)
            transcript_id = self.extract_unique_trans_id(new_df)
            tpm = self.extract_tpm(new_df)
            trans_tpm_dict = self.create_trans_tpm_dict(transcript_id, tpm)
            new_tpm = self.create_sample_tpm_pairs(common_elms, trans_tpm_dict)
            sample_tpm_pair = (sample_name, new_tpm)
            t_list.append(sample_tpm_pair)
       
        expression_matrix = self.create_expression_matrix(t_list, common_elms)
        output_path = os.path.join(self.output_dir, "expression_matrix.csv")
        expression_matrix.to_csv(output_path, sep='\t', index=True, header=True)

        logging.info(f"Finish constructing expression matrix with Task ID {TASK_ID}")


if __name__ == '__main__':
    group_info_path = '/data/youpu/Stella_IR_Project/utils/transcript_group_info.csv'
    trans_folder_path = '/data/youpu/Stella_IR_Project/temp/2024-09-12_01:57:43/transcripts'
    output_dir = '/data/youpu/Stella_IR_Project/utils'

    matrix_constructor = ConstructExpressionMatrix(group_info_path, trans_folder_path, output_dir)
    matrix_constructor.main()
