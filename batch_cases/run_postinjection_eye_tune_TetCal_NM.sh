#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TETCAL_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
WORK_DIR="$(cd "${TETCAL_DIR}/.." && pwd)"
cd "${WORK_DIR}"

EXE="./build/TetCal"
OUT_DIR="importanceTest"
EYE_MACRO_DIR="${OUT_DIR}/postEyeMacros"
mkdir -p "${OUT_DIR}" "${EYE_MACRO_DIR}"

if [[ ! -x "${EXE}" ]]; then
  echo "Executable not found: ${WORK_DIR}/${EXE}" >&2
  exit 1
fi

run_macro() {
  local macro="$1"
  local out_stem="$2"
  echo "=== ${out_stem} ==="
  "${EXE}" -m "${macro}" -o "${OUT_DIR}/${out_stem}.out"
}

make_eye_macro() {
  local eye_macro="$1"
  local out_macro="$2"

  awk '
    $1 == "/tetcal/eyeImportanceTuneMaxIterations" {
      print "/tetcal/eyeImportanceTuneMaxIterations 1"
      next
    }
    $1 == "/run/beamOn" {
      print "/run/beamOn 100000"
      next
    }
    { print }
  ' "${eye_macro}" > "${out_macro}"
}

eye_tune_converged() {
  local tune_data="$1"
  [[ -f "${tune_data}" ]] || return 1
  awk '
    NR > 1 { status = $NF }
    END { exit(status == "yes" ? 0 : 1) }
  ' "${tune_data}"
}

run_eye_tune() {
  local case_name="$1"
  local max_iterations=10
  local base_macro="TetCal/batch_cases/${case_name}_eye_tune1_TetCal_NM.in"
  local current_macro="${EYE_MACRO_DIR}/${case_name}_eye_iter1_TetCal_NM.in"

  make_eye_macro "${base_macro}" "${current_macro}"

  for (( iteration = 1; iteration <= max_iterations; iteration++ )); do
    local out_stem="${case_name}_eye_tune_iter${iteration}"
    run_macro "${current_macro}" "${out_stem}"

    local tune_data="${current_macro}.tuneData"
    local tuned_macro="${current_macro%.in}_tuned.in"
    if eye_tune_converged "${tune_data}"; then
      echo "Eye tune converged for ${case_name} at iteration ${iteration}: ${tuned_macro}"
      break
    fi

    if [[ ! -f "${tuned_macro}" ]]; then
      echo "Expected tuned macro not found: ${WORK_DIR}/${tuned_macro}" >&2
      exit 1
    fi

    if (( iteration == max_iterations )); then
      echo "Eye tune reached ${max_iterations} iterations for ${case_name}: ${tuned_macro}"
      break
    fi

    local next_macro="${EYE_MACRO_DIR}/${case_name}_eye_iter$((iteration + 1))_TetCal_NM.in"
    cp "${tuned_macro}" "${next_macro}"
    current_macro="${next_macro}"
  done
}

run_eye_tune "Postinjection_I-131_1"
run_eye_tune "Postinjection_I-131_2"
run_eye_tune "Postinjection_Tc-99m_1"
run_eye_tune "Postinjection_Tc-99m_2"
run_eye_tune "Postinjection_F-18_1"
run_eye_tune "Postinjection_F-18_2"
