#!/usr/bin/python
# -*- encoding: utf-8 -*-
"""
Created by Stella Li on 2024/08/07 17:01
"""

import os

from stella_config import Config
from RNA_seq_pipeline import TASK_ID


def main(idx, fastaq_dir, form, state):
    config = Config(TASK_ID)
    
    if form == 'paired':
        if state == 'begin':
            cmd1 = ("echo '#!/bin/sh\n"
                    f"{config.STAR} "
                    f"--runThreadN 48 "
                    f"--genomeDir {config.GENOME_DIR} "
                    f"--readFilesIn {os.path.join(fastaq_dir, idx + '_1.fastq')} {os.path.join(fastaq_dir, idx + '_2.fastq')} "
                    f"--limitGenomeGenerateRAM 300000000000 "
                    f"--outSAMtype BAM Unsorted "
                    f"--outFileNamePrefix {os.path.join(config.TEMP_MAPPING_DIR, idx)} "
                    f"--sjdbOverhang 649\n' "
                    f"| docker exec -i {config.YOU_IR_CONTAINER} sh -c 'cat > {config.TEMP_MAPPING_BASH_SCRIPT}'")
            os.system(cmd1)
        else:
            cmd1 = (f"echo '{config.STAR} "
                    f"--runThreadN 48 "
                    f"--genomeDir {config.GENOME_DIR} "
                    f"--readFilesIn {os.path.join(fastaq_dir, idx + '_1.fastq')} {os.path.join(fastaq_dir, idx + '_2.fastq')} "
                    f"--limitGenomeGenerateRAM 300000000000 "
                    f"--outSAMtype BAM Unsorted "
                    f"--outFileNamePrefix {os.path.join(config.TEMP_MAPPING_DIR, idx)} "
                    f"--sjdbOverhang 649\n' "
                    f"| docker exec -i {config.YOU_IR_CONTAINER} sh -c 'cat >> {config.TEMP_MAPPING_BASH_SCRIPT}'")
            os.system(cmd1)
        if state == 'last':
            cmd2 = f"docker exec {config.YOU_IR_CONTAINER} sh -c 'chmod 777 {config.TEMP_MAPPING_BASH_SCRIPT} && sh {config.TEMP_MAPPING_BASH_SCRIPT}'"
            os.system(cmd2)
    else:
        if state == 'begin':
            cmd1 = ("echo '#!/bin/sh\n"
                    f"{config.STAR} "
                    f"--runThreadN 48 "
                    f"--genomeDir {config.GENOME_DIR} "
                    f"--readFilesIn {os.path.join(fastaq_dir, idx + '.fastq')} "
                    f"--limitGenomeGenerateRAM 300000000000 "
                    f"--outSAMtype BAM Unsorted "
                    f"--outFileNamePrefix {os.path.join(config.TEMP_MAPPING_DIR, idx)} "
                    f"--sjdbOverhang 649\n' "
                    f"| docker exec -i {config.YOU_IR_CONTAINER} sh -c 'cat > {config.TEMP_MAPPING_BASH_SCRIPT}'")
            os.system(cmd1)
        else:
            cmd1 = (f"echo '{config.STAR} "
                    f"--runThreadN 48 "
                    f"--genomeDir {config.GENOME_DIR} "
                    f"--readFilesIn {os.path.join(fastaq_dir, idx + '.fastq')} "
                    f"--limitGenomeGenerateRAM 300000000000 "
                    f"--outSAMtype BAM Unsorted "
                    f"--outFileNamePrefix {os.path.join(config.TEMP_MAPPING_DIR, idx)} "
                    f"--sjdbOverhang 649\n' "
                    f"| docker exec -i {config.YOU_IR_CONTAINER} sh -c 'cat >> {config.TEMP_MAPPING_BASH_SCRIPT}'")
            os.system(cmd1)
        if state == 'last':
            cmd2 = f"docker exec {config.YOU_IR_CONTAINER} sh -c 'chmod 777 {config.TEMP_MAPPING_BASH_SCRIPT} && sh {config.TEMP_MAPPING_BASH_SCRIPT}'"
            os.system(cmd2)

if __name__ == '__main__':
    main()
