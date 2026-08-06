#!/usr/bin/env bash
# Project-wide defaults. Override by exporting variables before running a script
# or by editing config/local_resources.sh after server inspection.

MIN_LENGTH="${MIN_LENGTH:-1200}"
MAX_LENGTH="${MAX_LENGTH:-1800}"
MIN_QUALITY="${MIN_QUALITY:-10}"
THREADS="${THREADS:-8}"
EMU_MIN_PID="${EMU_MIN_PID:-80}"

MINIMAP_MIN_QUERY_COVERAGE="${MINIMAP_MIN_QUERY_COVERAGE:-0.80}"
MINIMAP_MIN_IDENTITY="${MINIMAP_MIN_IDENTITY:-0.85}"
MINIMAP_MIN_MAPQ="${MINIMAP_MIN_MAPQ:-10}"
MINIMAP_MIN_AS_DELTA="${MINIMAP_MIN_AS_DELTA:-20}"
MINIMAP_MIN_ALIGNMENT_LENGTH="${MINIMAP_MIN_ALIGNMENT_LENGTH:-1200}"
MINIMAP_MAX_ALIGNMENT_LENGTH="${MINIMAP_MAX_ALIGNMENT_LENGTH:-1800}"

EXPECTED_FASTQ_COUNT="${EXPECTED_FASTQ_COUNT:-20}"
ALLOW_NON_20_FASTQ="${ALLOW_NON_20_FASTQ:-0}"

OLD_PROJECT_ROOT="${OLD_PROJECT_ROOT:-}"
BROAD_EMU_DB="${BROAD_EMU_DB:-}"
NCBI_TAXONOMY_DIR="${NCBI_TAXONOMY_DIR:-}"
TAXONOMY_LIST="${TAXONOMY_LIST:-}"
VALIDATED_METADATA_FLAG="${VALIDATED_METADATA_FLAG:-metadata/.metadata_reviewed}"
VALIDATED_REFERENCES_FLAG="${VALIDATED_REFERENCES_FLAG:-references/combined/.references_validated}"
VALIDATED_CUSTOM_DB_FLAG="${VALIDATED_CUSTOM_DB_FLAG:-databases/custom_emu/.custom_emu_db_validated}"

# The Plasmidsaurus identifier/prefix column should be edited on the server if
# uploaded original filenames do not begin with these species-style prefixes.
# Format: clean_species_name|display species|original filename prefix|clean FASTA filename|clean FASTA record ID
REFERENCE_SPECS=(
  "Bacteroides_thetaiotaomicron|Bacteroides thetaiotaomicron|Bacteroides_thetaiotaomicron|Bacteroides_thetaiotaomicron.fasta|Bacteroides_thetaiotaomicron_ref"
  "Cutibacterium_acnes|Cutibacterium acnes|Cutibacterium_acnes|Cutibacterium_acnes.fasta|Cutibacterium_acnes_ref"
  "Fusobacterium_nucleatum|Fusobacterium nucleatum|Fusobacterium_nucleatum|Fusobacterium_nucleatum.fasta|Fusobacterium_nucleatum_ref"
  "Prevotella_melaninogenica|Prevotella melaninogenica|Prevotella_melaninogenica|Prevotella_melaninogenica.fasta|Prevotella_melaninogenica_ref"
  "Veillonella_parvula|Veillonella parvula|Veillonella_parvula|Veillonella_parvula.fasta|Veillonella_parvula_ref"
)

EXPECTED_SPECIES_ORDER=(
  "Bacteroides thetaiotaomicron"
  "Cutibacterium acnes"
  "Fusobacterium nucleatum"
  "Prevotella melaninogenica"
  "Veillonella parvula"
)

if [[ -f "$(dirname "${BASH_SOURCE[0]}")/local_resources.sh" ]]; then
  # shellcheck source=/dev/null
  source "$(dirname "${BASH_SOURCE[0]}")/local_resources.sh"
fi
