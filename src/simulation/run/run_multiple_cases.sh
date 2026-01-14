#!/bin/bash

# set -e

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 WORKING_DIR"
    exit 1
fi

WORKING_DIR="$1"

SETUP_FILES=(
    "/app/examples/naca0012_study/aoa_0.json"
    "/app/examples/naca0012_study/aoa_5.json"
    "/app/examples/naca0012_study/aoa_10.json"
    "/app/examples/naca0012_study/aoa_15.json"
)

CASE_NAMES=(
    "case_aoa_0"
    "case_aoa_5"
    "case_aoa_10"
    "case_aoa_15"
)

for i in "${!SETUP_FILES[@]}"; do
    setup_file="${SETUP_FILES[$i]}"
    case_name="${CASE_NAMES[$i]}"
    echo "Launching $case_name with $setup_file"
    ./run_case.sh \
        --working-dir "$WORKING_DIR" \
        --setup-file "$setup_file" \
        --case-name "$case_name"
done

poetry run python3 /app/src/postprocess/aggregate_summaries.py \
    --working-path "$WORKING_DIR" \
    --output "$WORKING_DIR/summary_naca0012.csv"

poetry run python3 /app/src/postprocess/plotting/generate_polar_plots.py \
    --summary-csvs "$WORKING_DIR/summary_naca0012.csv" \
    --reference-csv "/app/examples/naca0012_study/nasa_landson_experiment.csv"