#!/bin/bash

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 WORKING_DIR"
    exit 1
fi

WORKING_DIR="$1"

SETUP_FILES=(
    "/app/test_input/aoa_5.json"
    "/app/test_input/aoa_10.json"
    "/app/test_input/aoa_15.json"
)

CASE_NAMES=(
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