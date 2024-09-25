#!/usr/bin/python
# -*- encoding: utf-8 -*-
"""
Created by Stella Li on 2024/08/07 15:26
"""

import os

from stella_config import Config
from RNA_seq_pipeline import TASK_ID


def main():
    config = Config(TASK_ID)

    cmd1 = ("echo '#!/bin/sh\n"
            f"mkdir -p {config.TEMP_OUTPUT_TRANSCRIPT_FASTA} && "
            f"{config.TRANSDECODER} "
            f"{os.path.join(config.TEMP_TRANSCRIPTS_DIR, 'merged.transcripts.gtf')} "
            f"{config.GENOME_FASTA} > "
            f"{os.path.join(config.TEMP_OUTPUT_TRANSCRIPT_FASTA, 'transcript.fasta')} && "
            f"sleep 5 && "
            f"{config.TRANSDECODER_LONGORFS} "
            f"-t {os.path.join(config.TEMP_OUTPUT_TRANSCRIPT_FASTA, 'transcript.fasta')} "
            f"-O {config.TEMP_OUTPUT_TRANSCRIPT_FASTA} && "
            f"sleep 5 && "
            f"{config.TRANSDECODER_PREDICT} "
            f"-t {os.path.join(config.TEMP_OUTPUT_TRANSCRIPT_FASTA, 'transcript.fasta')} "
            f"-O {config.TEMP_OUTPUT_TRANSCRIPT_FASTA}' "
            f"| docker exec -i {config.YOU_IR_CONTAINER} sh -c 'cat > {config.TEMP_PREDICT_BASH_SCRIPT}'")
    os.system(cmd1)

    cmd2 = f"docker exec {config.YOU_IR_CONTAINER} sh -c 'chmod 777 {config.TEMP_PREDICT_BASH_SCRIPT} && sh {config.TEMP_PREDICT_BASH_SCRIPT}'"
    os.system(cmd2)


if __name__ == '__main__':
    main()
