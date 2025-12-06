#!/bin/bash
set -ex

export OMPI_ALLOW_RUN_AS_ROOT=1
export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

while [[ $# -gt 0 ]]; do
  key="$1"
  case $key in
    --working-dir)
      working_dir="$2"
      shift; shift
      ;;
    --setup-file)
      setup_file="$2"
      shift; shift
      ;;
    --case-name)
      case_name="$2"
      shift; shift
      ;;
    *)
      shift
      ;;
  esac
done

rm -rf "$working_dir/$case_name"

poetry run python3 ../preparation/prepare_case.py \
    --working-dir "$working_dir" \
    --setup-file "$setup_file" \
    --case-name "$case_name"

cp *.sh "$working_dir/$case_name/"

cd "$working_dir/$case_name"

# Pass the setup_file path directly
bash ./run_meshing.sh "$setup_file"
bash ./run_simulation.sh "$setup_file"
# bash ./run_postprocess.sh
