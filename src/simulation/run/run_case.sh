#!/bin/bash
set -e

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
    --case-name "$case_name" 2>&1 | tee prepare_case.log

cp *.sh "$working_dir/$case_name/"
cp ${setup_file} "$working_dir/$case_name/${case_name}.json"

cd "$working_dir/$case_name"

bash ./run_meshing.sh "$setup_file" 2>&1 | tee run_meshing.log
bash ./run_simulation.sh "$setup_file" 2>&1 | tee run_simulation.log
bash ./run_postprocess.sh "$setup_file" 2>&1 | tee run_postprocess.log
