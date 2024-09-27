#!/usr/bin/python
# -*- encoding: utf-8 -*-
"""
Created by Stella Li on 2024/08/07 16:46
"""

import os

from stella_config import Config
from RNA_seq_pipeline import TASK_ID


def main(idx, i):
    config = Config(TASK_ID)
    
    if i == 'begin':
        cmd1 = ("echo '#!/bin/sh\n"
                f"{config.SAMTOOLS} "
                f"sort -@16 "
                f"{os.path.join(config.TEMP_MAPPING_DIR, idx + 'Aligned.out.bam')} "
                f"-o {os.path.join(config.TEMP_MAPPING_DIR, idx + 'Align_sort.bam')} && "
                f"rm {os.path.join(config.TEMP_MAPPING_DIR, idx + 'Aligned.out.bam')}\n' "
                f"| docker exec -i {config.YOU_IR_CONTAINER} sh -c 'cat > {config.TEMP_SORT_BASH_SCRIPT}'")
        print(cmd1)
        os.system(cmd1)
    else:
        cmd1 = (f"echo '{config.SAMTOOLS} "
                f"sort -@16 "
                f"{os.path.join(config.TEMP_MAPPING_DIR, idx + 'Aligned.out.bam')} "
                f"-o {os.path.join(config.TEMP_MAPPING_DIR, idx + 'Align_sort.bam')} && "
                f"rm {os.path.join(config.TEMP_MAPPING_DIR, idx + 'Aligned.out.bam')}\n' "
                f"| docker exec -i {config.YOU_IR_CONTAINER} sh -c 'cat >> {config.TEMP_SORT_BASH_SCRIPT}'")
        print(cmd1)
        os.system(cmd1)
    if i == 'last':
        cmd2 = f"docker exec {config.YOU_IR_CONTAINER} sh -c 'chmod 777 {config.TEMP_SORT_BASH_SCRIPT} && sh {config.TEMP_SORT_BASH_SCRIPT}'"
        os.system(cmd2)


if __name__ == '__main__':
    main()
