#!/usr/bin/python
# -*- encoding: utf-8 -*-
"""
Created by Stella Li on 2024/08/07 16:31
"""

import os

from RNA_seq_pipeline import TASK_ID
from stella_config import Config


def main(idx, i):
    config = Config(TASK_ID)
    
    if i == 'begin':
        cmd1 = ("echo '#!/bin/sh\n"
                f"{config.STRINGTIE} "
                f"{os.path.join(config.TEMP_BAM_STORAGE_DIR, idx + 'Align_sort.bam')} "
                f"-p 16 "
                f"-G {config.INTRON_GTF_FILE_PATH} "
                f"-B "
                f"-o {os.path.join(config.TEMP_TRANSCRIPTS_DIR, idx + ".transcripts.gtf")} "
                f"-e\n' "
                f"| docker exec -i {config.YOU_IR_CONTAINER} sh -c 'cat > {config.TEMP_ASSEMBLE_BASH_SCRIPT}'")
        os.system(cmd1)
    else:
        cmd1 = (f"echo '{config.STRINGTIE} "
                f"{os.path.join(config.TEMP_BAM_STORAGE_DIR, idx + 'Align_sort.bam')} "
                f"-p 16 "
                f"-G {config.INTRON_GTF_FILE_PATH} "
                f"-B "
                f"-o {os.path.join(config.TEMP_TRANSCRIPTS_DIR, idx + ".transcripts.gtf")} "
                f"-e\n' "
                f"| docker exec -i {config.YOU_IR_CONTAINER} sh -c 'cat >> {config.TEMP_ASSEMBLE_BASH_SCRIPT}'")
        os.system(cmd1)
    if i == 'last':
        cmd2 = f"docker exec {config.YOU_IR_CONTAINER} sh -c 'chmod 777 {config.TEMP_ASSEMBLE_BASH_SCRIPT} && sh {config.TEMP_ASSEMBLE_BASH_SCRIPT}'"
        os.system(cmd2)

    
if __name__ == '__main__':
    main()
