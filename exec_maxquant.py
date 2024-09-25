#!/usr/bin/python
# -*- encoding: utf-8 -*-
"""
Created by Stella Li on 2024/08/09 11:07
"""

import os

from stella_config import Config
from RNA_seq_pipeline import TASK_ID


def main():
    config = Config(TASK_ID)

    cmd1 = ("echo '#!/bin/sh\n"
            f"{config.MONO} "
            f"{config.MAXQUANT} "
            f"{config.TEMP_MONO_NEW_XML_FILE}' "
            f"| docker exec -i {config.MONO_CONTAINER} sh -c 'cat > {config.TEMP_MAXQUANT_BASH_SCRIPT}'")
    os.system(cmd1)

    cmd2 = f"docker exec {config.MONO_CONTAINER} sh -c 'chmod 777 {config.TEMP_MAXQUANT_BASH_SCRIPT} && sh {config.TEMP_MAXQUANT_BASH_SCRIPT}'"
    os.system(cmd2)


if __name__ == '__main__':
    main()
