#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TETCAL_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
WORK_DIR="$(cd "${TETCAL_DIR}/.." && pwd)"
cd "${WORK_DIR}"

EXE="./build/TetCal"
OUT_DIR="importanceTest"
mkdir -p "${OUT_DIR}"

if [[ ! -x "${EXE}" ]]; then
  echo "Executable not found: ${WORK_DIR}/${EXE}" >&2
  exit 1
fi

run_case() {
  local macro="$1"
  local stem
  stem="$(basename "${macro}" .in)"
  echo "=== ${stem} ==="
  "${EXE}" -m "${macro}" -o "${OUT_DIR}/${stem}.out"
}

run_case "TetCal/batch_cases/Preinjection_I-131_1_eye_tune1_TetCal_NM.in"
run_case "TetCal/batch_cases/Preinjection_I-131_2_eye_tune1_TetCal_NM.in"
run_case "TetCal/batch_cases/Preinjection_I-131_3_eye_tune1_TetCal_NM.in"
run_case "TetCal/batch_cases/Preinjection_Tc-99m_1_eye_tune1_TetCal_NM.in"
run_case "TetCal/batch_cases/Preinjection_Tc-99m_2_eye_tune1_TetCal_NM.in"
run_case "TetCal/batch_cases/Preinjection_Tc-99m_3_eye_tune1_TetCal_NM.in"
run_case "TetCal/batch_cases/Preinjection_F-18_1_eye_tune1_TetCal_NM.in"
run_case "TetCal/batch_cases/Preinjection_F-18_2_eye_tune1_TetCal_NM.in"
run_case "TetCal/batch_cases/Preinjection_F-18_3_eye_tune1_TetCal_NM.in"
run_case "TetCal/batch_cases/Injection_I-131_1_eye_tune1_TetCal_NM.in"
run_case "TetCal/batch_cases/Injection_I-131_2_eye_tune1_TetCal_NM.in"
run_case "TetCal/batch_cases/Injection_I-131_3_eye_tune1_TetCal_NM.in"
run_case "TetCal/batch_cases/Injection_Tc-99m_1_eye_tune1_TetCal_NM.in"
run_case "TetCal/batch_cases/Injection_Tc-99m_2_eye_tune1_TetCal_NM.in"
run_case "TetCal/batch_cases/Injection_F-18_1_eye_tune1_TetCal_NM.in"
run_case "TetCal/batch_cases/Injection_F-18_2_eye_tune1_TetCal_NM.in"
