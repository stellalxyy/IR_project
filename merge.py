#!/usr/bin/python
# -*- encoding: utf-8 -*-
"""
Created by Stella Li on 2024/08/07 15:22
"""

import os
from stella_config import Config
from RNA_seq_pipeline import TASK_ID


def main():
    config = Config(TASK_ID)

    cmd1 = ("echo '#!/bin/sh\n"
            f"mkdir -p {config.TEMP_GTF_LIST_DIR} && "
            f"find {config.TEMP_TRANSCRIPTS_DIR} -type f -name '*.gtf' > "
            f"{os.path.join(config.TEMP_GTF_LIST_DIR, 'gtf_list.txt')} && "
            f"{config.STRINGTIE} "
            f"--merge -G {config.INTRON_GTF_FILE_PATH} "
            f"-p 16 "
            f"-o {os.path.join(config.TEMP_TRANSCRIPTS_DIR, 'merged.transcripts.gtf')} {os.path.join(config.TEMP_GTF_LIST_DIR, 'gtf_list.txt')} "
            f"-i' "
            f"| docker exec -i {config.YOU_IR_CONTAINER} sh -c 'cat > {config.TEMP_MERGE_BASH_SCRIPT}'")
    # print(cmd1)
    os.system(cmd1)
    cmd2 = f"docker exec {config.YOU_IR_CONTAINER} sh -c 'chmod 777 {config.TEMP_MERGE_BASH_SCRIPT} && sh {config.TEMP_MERGE_BASH_SCRIPT}'"
    # print(cmd2)
    os.system(cmd2)


if __name__ == "__main__":
    main()
