#!/bin/bash
# submit_repopulate_batches.sh - Submit SLURM jobs to repopulate created batches

# Default configuration
CONFIG_PATH="config/config.yml"
SCRIPT_PATH="scripts/core/repopulate_ecod_schema.py"
THREADS=4
MEMORY="8G"
TIME="4:00:00"
LOG_DIR="logs/repopulate_jobs"
VERBOSE=false
DRY_RUN=false

# Parse command line arguments
while [[ $# -gt 0 ]]; do
  case $1 in
    --config)
      CONFIG_PATH="$2"
      shift 2
      ;;
    --script)
      SCRIPT_PATH="$2"
      shift 2
      ;;
    --threads)
      THREADS="$2"
      shift 2
      ;;
    --memory)
      MEMORY="$2"
      shift 2
      ;;
    --time)
      TIME="$2"
      shift 2
      ;;
    --dry-run)
      DRY_RUN=true
      shift
      ;;
    -v|--verbose)
      VERBOSE=true
      shift
      ;;
    *)
      echo "Unknown option: $1"
      echo "Usage: $0 [--config CONFIG_PATH] [--script SCRIPT_PATH] [--threads N] [--memory XG] [--time HH:MM:SS] [--dry-run] [-v|--verbose]"
      exit 1
      ;;
  esac
done

# Create log directory
mkdir -p "$LOG_DIR"
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_FILE="${LOG_DIR}/repopulate_submission_${TIMESTAMP}.log"

# Build command
CMD="python ecod/scripts/submit_repopulate_jobs.py --config ${CONFIG_PATH} --script ${SCRIPT_PATH} --threads ${THREADS} --memory ${MEMORY} --time ${TIME} --log-file ${LOG_FILE}"

# Add optional flags
if [ "$VERBOSE" = true ]; then
  CMD="${CMD} --verbose"
fi

if [ "$DRY_RUN" = true ]; then
  CMD="${CMD} --dry-run"
fi

# Echo command
echo "Executing: $CMD"

# Execute command
eval "$CMD"

# Check exit code
EXIT_CODE=$?
if [ $EXIT_CODE -eq 0 ]; then
  echo "Job submission completed successfully. Log saved to: ${LOG_FILE}"
else
  echo "Job submission failed with exit code ${EXIT_CODE}. Check log file: ${LOG_FILE}"
fi

exit $EXIT_CODE