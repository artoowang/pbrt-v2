#!/usr/bin/bash

if [ $# -lt 1 ]; then
    echo >&2 "Usage: $0 <pbrt_file_name>"
    exit -1
fi

NAME=$(echo "$1" | sed -e 's/\.pbrt$//')
PBRT_FILE="${NAME}.pbrt"
EXR_FILE="${NAME}.exr"
OUTPUT_DIR="${NAME}"

if [ ! -r "$PBRT_FILE" ]; then
    echo >&2 "Cannot read $PBRT_FILE"
    exit -1
fi

DISPLAY_CMD="i_view32"

# Determine if display command is available
if ! which "$DISPLAY_CMD"; then
    DISPLAY_CMD="eog"
    if ! which "$DISPLAY_CMD"; then
        DISPLAY_CMD="display"
    fi
fi

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
exrtotiff -gamma 2.2 "$EXR_FILE" "${OUTPUT_DIR}/${TIF_FILE}"
#exrdisplay "$EXR_FILE" &
"$DISPLAY_CMD" "${OUTPUT_DIR}/${TIF_FILE}" &

