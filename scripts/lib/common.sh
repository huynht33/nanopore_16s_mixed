#!/usr/bin/env bash

set -euo pipefail

die() {
  printf 'ERROR: %s\n' "$*" >&2
  exit 1
}

warn() {
  printf 'WARNING: %s\n' "$*" >&2
}

project_root_from_script() {
  local script_dir
  if [[ -n "${SCRIPT_DIR:-}" ]]; then
    script_dir="${SCRIPT_DIR}"
  else
    local idx
    idx=$((${#BASH_SOURCE[@]} - 1))
    script_dir="$(cd "$(dirname "${BASH_SOURCE[${idx}]}")" && pwd)"
  fi
  cd "${script_dir}/.." && pwd
}

load_project_config() {
  PROJECT_ROOT="${PROJECT_ROOT:-$(project_root_from_script)}"
  [[ -d "${PROJECT_ROOT}" ]] || die "PROJECT_ROOT does not exist: ${PROJECT_ROOT}"
  # shellcheck source=/dev/null
  source "${PROJECT_ROOT}/config/project_config.sh"
}

timestamp() {
  date '+%Y%m%d_%H%M%S'
}

make_log() {
  local name="$1"
  mkdir -p "${PROJECT_ROOT}/logs"
  LOG_FILE="${PROJECT_ROOT}/logs/$(timestamp)_${name}.log"
  exec >> "${LOG_FILE}" 2>&1
  printf 'Log: %s\n' "${LOG_FILE}"
  printf 'Project root: %s\n' "${PROJECT_ROOT}"
}

have_cmd() {
  command -v "$1" >/dev/null 2>&1
}

md5_file() {
  local file="$1"
  if have_cmd md5sum; then
    md5sum "${file}" | awk '{print $1}'
  elif have_cmd md5; then
    md5 -q "${file}"
  else
    die "Neither md5sum nor md5 is available"
  fi
}

abs_path() {
  local path="$1"
  if [[ -d "${path}" ]]; then
    (cd "${path}" && pwd)
  else
    local dir base
    dir="$(dirname "${path}")"
    base="$(basename "${path}")"
    printf '%s/%s\n' "$(cd "${dir}" && pwd)" "${base}"
  fi
}

record_version() {
  local label="$1"
  shift
  printf '\n## %s\n' "${label}"
  if "$@" >/tmp/nano16s_version_stdout.$$ 2>/tmp/nano16s_version_stderr.$$; then
    cat /tmp/nano16s_version_stdout.$$
    cat /tmp/nano16s_version_stderr.$$ >&2
  else
    printf 'Command failed: %s\n' "$*"
    cat /tmp/nano16s_version_stdout.$$ || true
    cat /tmp/nano16s_version_stderr.$$ || true
  fi
  rm -f /tmp/nano16s_version_stdout.$$ /tmp/nano16s_version_stderr.$$
}

fastq_list_file() {
  printf '%s/metadata/fastq_files.list\n' "${PROJECT_ROOT}"
}
