#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
LOG_DIR="${PROJECT_ROOT}/08_logs"
TIMESTAMP="$(date +"%Y%m%d_%H%M%S")"
LOG_FILE="${LOG_DIR}/environment_audit_${TIMESTAMP}.log"

mkdir -p "${LOG_DIR}"

run_section() {
  local title="$1"
  shift
  {
    printf '\n## %s\n' "${title}"
    "$@" 2>&1 || printf 'Command failed: %s\n' "$*"
  } >> "${LOG_FILE}"
}

record_version() {
  local tool="$1"
  local exe
  exe="$(command -v "${tool}" 2>/dev/null || true)"
  {
    printf '\n### %s\n' "${tool}"
    if [[ -z "${exe}" ]]; then
      printf 'status: unavailable\n'
      return 0
    fi
    printf 'path: %s\n' "${exe}"
    case "${tool}" in
      esearch|efetch)
        "${tool}" -help 2>&1 | head -n 3 || true
        ;;
      snakemake)
        "${tool}" --version 2>&1 || true
        ;;
      *)
        "${tool}" --version 2>&1 | head -n 3 || true
        ;;
    esac
  } >> "${LOG_FILE}"
}

{
  printf '# Environment Audit\n'
  printf 'timestamp: %s\n' "${TIMESTAMP}"
  printf 'project_root: %s\n' "${PROJECT_ROOT}"
  printf 'user: redacted\n'
  printf 'credential_policy: environment variables and credential files were not printed\n'
  printf 'conda_default_env: %s\n' "${CONDA_DEFAULT_ENV:-}"
  printf 'virtual_env_active: %s\n' "${VIRTUAL_ENV:+yes}"
} > "${LOG_FILE}"

run_section "Operating System" uname -a

if command -v sw_vers >/dev/null 2>&1; then
  run_section "macOS Version" sw_vers
fi

if command -v sysctl >/dev/null 2>&1; then
  {
    printf '\n## CPU\n'
    cpu_brand="$(sysctl -n machdep.cpu.brand_string 2>/dev/null || true)"
    cpu_count="$(sysctl -n hw.ncpu 2>/dev/null || true)"
    if [[ -n "${cpu_brand}${cpu_count}" ]]; then
      printf 'brand: %s\n' "${cpu_brand:-unknown}"
      printf 'logical_cpus: %s\n' "${cpu_count:-unknown}"
    elif command -v system_profiler >/dev/null 2>&1; then
      system_profiler SPHardwareDataType 2>/dev/null \
        | awk '
            /Chip:/ || /Total Number of Cores:/ || /Processor Name:/ || /Number of Processors:/ || /Total Number of Cores:/ {
              sub(/^[[:space:]]+/, "")
              print
            }
          '
    else
      printf 'brand: unknown\nlogical_cpus: unknown\n'
    fi

    printf '\n## Memory\n'
    mem_bytes="$(sysctl -n hw.memsize 2>/dev/null || true)"
    if [[ -n "${mem_bytes}" ]]; then
      printf 'bytes: %s\n' "${mem_bytes}"
      awk -v bytes="${mem_bytes}" 'BEGIN {printf "gib: %.2f\n", bytes/1024/1024/1024}'
    elif command -v system_profiler >/dev/null 2>&1; then
      system_profiler SPHardwareDataType 2>/dev/null \
        | awk '/Memory:/ {sub(/^[[:space:]]+/, ""); print}'
    else
      printf 'bytes: unknown\n'
    fi
  } >> "${LOG_FILE}"
else
  run_section "CPU" sh -c 'lscpu 2>/dev/null || true'
  run_section "Memory" sh -c 'free -h 2>/dev/null || true'
fi

run_section "Disk" df -h "${PROJECT_ROOT}"

{
  printf '\n## Software Versions\n'
} >> "${LOG_FILE}"

for tool in git python pip conda mamba R Rscript curl wget aria2c fasterq-dump prefetch esearch efetch git-lfs snakemake nextflow; do
  record_version "${tool}"
done

printf 'Environment audit written to %s\n' "${LOG_FILE}"
