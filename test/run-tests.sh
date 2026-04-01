 #!/usr/bin/env bash

set -euox pipefail

script_dir="$(dirname $0)"
DOCKER=${1:-no}

if [ "$DOCKER" == "gh" ]
then
    export NXF_CONTAINER_ENGINE=docker
    docker_flag='-profile gh -stub'
else
    export SINGULARITY_FAKEROOT=1
    docker_flag='-stub '
fi

nextflow run "$script_dir"/.. \
    -resume $docker_flag \
    --sample_sheet "$script_dir"/sra/inputs/sample-sheet.csv \
    --inputs "$script_dir"/sra/inputs \
    --outputs "$script_dir"/sra/outputs \
    -c "$script_dir"/sra/nextflow.config
