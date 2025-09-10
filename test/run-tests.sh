 #!/usr/bin/env bash

set -euox pipefail

script_dir="$(dirname $0)"
DOCKER=${1:-no}

if [ "$DOCKER" == "gh" ]
then
    export NXF_CONTAINER_ENGINE=docker
    docker_flag='-profile gh'
else
    export SINGULARITY_FAKEROOT=1
    docker_flag=''
fi

nextflow run "$script_dir"/.. \
    -resume $docker_flag \
    --sample_sheet "$script_dir"/sra/inputs/sample-sheet.csv \
    --inputs "$script_dir"/sra/inputs \
    --outputs "$script_dir"/sra/outputs
