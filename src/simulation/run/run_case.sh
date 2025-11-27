#!/bin/bash
# run_case.sh

set -e

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
    --resolution)
      resolution="$2"
      shift; shift
      ;;
    *)
      shift
      ;;
  esac
done

n_proc=$(jq -r '.Simulation.Decomposition.NumberOfSubdomains' "$setup_file")

python3 ../preparation/prepare_case.py \ 
    --working-dir "$working_dir" \
    --setup-file "$setup_file" \
    --case-name "$case_name" \
    --resolution "$resolution"

cp *.sh $working_dir/$case_name/

cd $working_dir/$case_name 

./run_meshing.sh "$n_proc"

./run_simulation.sh "$n_proc"

./run_postprocess.sh



