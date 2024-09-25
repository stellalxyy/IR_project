#!/usr/bin/python
# -*- encoding: utf-8 -*-
"""
Created by Stella Li on 2024/08/07 16:27
"""

import os
import pandas as pd

from stella_config import Config
from RNA_seq_pipeline import TASK_ID

config = Config(TASK_ID)


class XMLProcessor:
    def __init__(self,
                 ms_input_fasta_path: str,
                 raw_file_path: str,
                 raw_group_info_path: str):
        
        self.ms_input_fasta_path = ms_input_fasta_path
        self.raw_file_path = raw_file_path
        self.raw_group_info_path = raw_group_info_path
        self.raw_file_list, self.experiments_name = self.generate_experiments_name()
        self.num_raw_files = len(self.raw_file_list)
    
    def generate_experiments_name(self):
        df = pd.read_csv(self.raw_group_info_path, sep=',', header=None)
        sample_list = df.iloc[:,0].tolist()
        if not sample_list[0].endswith('.raw'):
            df.drop([0], axis=0, inplace=True)
        pass
        df.columns = ['sample', 'group']
        df.reset_index(inplace=True)
        df.drop(columns = ['index'], inplace=True)
        raw_files = df['sample'].tolist()
        groups = df['group'].tolist()
        return raw_files, groups

    def generate_tags(self):
        tag_raw = ''
        tag_experiments = ''
        tag_fractions = ''
        tag_ptms = ''
        tag_paramGroupIndices = ''
        tag_referenceChannel = ''

        for i in self.raw_file_list:
            bb = '\t<string>' + '/root' + self.raw_file_path + '/' + str(i) + '</string>\n\t'
            tag_raw += bb

        for n in self.experiments_name:
            bb = '\t<string>' + str(n) + '</string>\n\t'
            tag_experiments += bb

        for i in range(self.num_raw_files):
            bb1 = '\t<short>32767</short>\n\t'
            tag_fractions += bb1
            bb2 = '\t<boolean>False</boolean>\n\t'
            tag_ptms += bb2
            bb3 = '\t<int>0</int>\n\t'
            tag_paramGroupIndices += bb3
            bb4 = '\t<string></string>\n\t'
            tag_referenceChannel += bb4

        return tag_raw, tag_experiments, tag_fractions, tag_ptms, tag_paramGroupIndices, tag_referenceChannel

    def write_new_xml_file(self, tag_raw, tag_experiments, tag_fractions, tag_ptms, tag_paramGroupIndices,
                           tag_referenceChannel):
        if not os.path.isfile(config.ORIGIN_XML_FILE):
            raise FileNotFoundError(f"Origin XML file {config.ORIGIN_XML_FILE} does not exist.")

        new_xml_list = []
        with open(config.ORIGIN_XML_FILE, 'r') as f:
            bb = f.readlines()
            for l in bb:
                if '<fastaFilePath>' in l:
                    new_l = '         <fastaFilePath>' + self.ms_input_fasta_path + '</fastaFilePath>\n'
                elif '<filePaths>' in l:
                    new_l = l.replace('<filePaths></filePaths>', '<filePaths>\n\t' + tag_raw + '</filePaths>')
                elif '<experiments>' in l:
                    new_l = l.replace('<experiments></experiments>',
                                      '<experiments>\n\t' + tag_experiments + '</experiments>')
                elif '<fractions>' in l:
                    new_l = l.replace('<fractions></fractions>', '<fractions>\n\t' + tag_fractions + '</fractions>')
                elif '<ptms>' in l:
                    new_l = l.replace('<ptms></ptms>', '<ptms>\n\t' + tag_ptms + '</ptms>')
                elif '<paramGroupIndices>' in l:
                    new_l = l.replace('<paramGroupIndices></paramGroupIndices>',
                                      '<paramGroupIndices>\n\t' + tag_paramGroupIndices + '</paramGroupIndices>')
                elif '<referenceChannel>' in l:
                    new_l = l.replace('<referenceChannel></referenceChannel>',
                                      '<referenceChannel>\n\t' + tag_referenceChannel + '</referenceChannel>')
                else:
                    new_l = l
                new_xml_list.append(new_l)

        with open(config.NEW_XML_FILE, 'w') as f:
            f.writelines(new_xml_list)

    def run(self):
        try:
            tags = self.generate_tags()
            self.write_new_xml_file(*tags)
        except Exception as e:
            print(f"Error: {e}")


if __name__ == '__main__':
    _ms_input_fasta_path = '/root/data/original/genome/Homo/AD_ir_orf001.fasta'
    _raw_file_path = '/root/data/original/ms/PD/PXD034120/brain/DDA'
    _new_xml_file = '/root/new_mqpar.xml'

    # processor = XMLProcessor(ms_input_fasta_path, raw_file_path, new_xml_file)
    # processor.run()
