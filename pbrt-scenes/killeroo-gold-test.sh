#!/usr/bin/bash

NAME="killeroo-gold-test"
PBRT_FILE="${NAME}.pbrt"
EXR_FILE="killeroo-gold.exr"
OUTPUT_DIR="${NAME}"

mkdir -p "${OUTPUT_DIR}"

# Search for unused sequence number
for i in $(seq -f '%03g' 0 999); do
    TIF_FILE="${NAME}.${i}.tif"
    if [ ! -e "${OUTPUT_DIR}/$TIF_FILE" ]; then
        break
    fi
done

if [ -e "${OUTPUT_DIR}/$TIF_FILE" ]; then
    echo >&2 "All sequence numbers are used. Abort"
    exit -1
fi

echo >&2 "Write to ${OUTPUT_DIR}/${TIF_FILE}"

pbrt "$PBRT_FILE"
exrtotiff "$EXR_FILE" "${OUTPUT_DIR}/${TIF_FILE}"
i_view32 "${OUTPUT_DIR}/${TIF_FILE}"

