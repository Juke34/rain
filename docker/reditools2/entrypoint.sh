#!/bin/bash
source /opt/conda/etc/profile.d/conda.sh
conda activate myenv
exec "$@"
