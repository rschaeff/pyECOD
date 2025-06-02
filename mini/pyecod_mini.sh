#!/bin/bash
# pyecod_mini - Executable wrapper script

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Execute the Python script with all arguments
exec python3 "$SCRIPT_DIR/pyecod_mini.py" "$@"
