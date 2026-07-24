#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TETCAL_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
WORK_DIR="$(cd "${TETCAL_DIR}/.." && pwd)"
cd "${WORK_DIR}"

EXE="./build/TetCal"
OUT_DIR="${OUT_DIR:-addedDoseCases}"
NPS="${NPS:-100000}"
THREADS="${THREADS:-4}"
MACRO_DIR="${OUT_DIR}/macros"
mkdir -p "${OUT_DIR}" "${MACRO_DIR}"

if [[ ! -x "${EXE}" ]]; then
  echo "Executable not found: ${WORK_DIR}/${EXE}" >&2
  exit 1
fi

make_run_macro() {
  local src_macro="$1"
  local dst_macro="$2"

  awk -v nps="${NPS}" -v threads="${THREADS}" '
    $1 == "/run/numberOfThreads" {
      print "/run/numberOfThreads " threads
      next
    }
    $1 == "/run/beamOn" {
      print "/run/beamOn " nps
      next
    }
    { print }
  ' "${src_macro}" > "${dst_macro}"
}

run_case() {
  local case_name="$1"
  local src_macro="TetCal/batch_cases/${case_name}_dose_TetCal_NM.in"
  local run_macro="${MACRO_DIR}/${case_name}_dose_TetCal_NM.in"

  if [[ ! -f "${src_macro}" ]]; then
    echo "Macro not found: ${WORK_DIR}/${src_macro}" >&2
    exit 1
  fi

  make_run_macro "${src_macro}" "${run_macro}"
  echo "=== ${case_name} ==="
  "${EXE}" -m "${run_macro}" -o "${OUT_DIR}/${case_name}.out"
}

cases=(
  "InjectionArmVein_I-131_1"
  "InjectionArmVein_I-131_2"
  "InjectionArmVein_Tc-99m_1"
  "InjectionArmVein_Tc-99m_2"
  "InjectionArmVein_F-18_1"
  "InjectionArmVein_F-18_2"
  "PostinjectionAllBlood_I-131_1"
  "PostinjectionAllBlood_I-131_2"
  "PostinjectionAllBlood_Tc-99m_1"
  "PostinjectionAllBlood_Tc-99m_2"
  "PostinjectionAllBlood_F-18_1"
  "PostinjectionAllBlood_F-18_2"
)

for case_name in "${cases[@]}"; do
  run_case "${case_name}"
done
