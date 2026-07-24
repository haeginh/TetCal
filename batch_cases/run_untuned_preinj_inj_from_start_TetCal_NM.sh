#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TETCAL_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
WORK_DIR="$(cd "${TETCAL_DIR}/.." && pwd)"
cd "${WORK_DIR}"

EXE="./build/TetCal"
OUT_DIR="importanceTest"
EYE_MACRO_DIR="${OUT_DIR}/retuneEyeMacros"
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

shield_value_from_tuned_macro() {
  local tuned_macro="$1"
  awk '$1 == "/tetcal/shieldImportanceValue" { value = $2 } END { if (value != "") print value; else exit 1 }' "${tuned_macro}"
}

make_initial_eye_macro_with_shield_value() {
  local eye_macro="$1"
  local shield_tuned_macro="$2"
  local out_macro="$3"
  local shield_value

  shield_value="$(shield_value_from_tuned_macro "${shield_tuned_macro}")"

  awk -v shield_value="${shield_value}" '
    $1 == "/tetcal/eyeImportanceValues" {
      for (i = 2; i <= NF; ++i) {
        $i = 1
      }
      print
      next
    }
    $1 == "/tetcal/shieldImportanceValue" {
      print "/tetcal/shieldImportanceValue " shield_value
      next
    }
    $1 == "/tetcal/shieldImportanceTune" {
      print "/tetcal/shieldImportanceTune true"
      next
    }
    $1 == "/tetcal/eyeImportanceTuneMaxIterations" {
      print "/tetcal/eyeImportanceTuneMaxIterations 1"
      next
    }
    $1 == "/tetcal/shieldImportanceTuneMaxIterations" {
      print "/tetcal/shieldImportanceTuneMaxIterations 1"
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

run_eye_tune_from_start() {
  local case_name="$1"
  local eye_macro="TetCal/batch_cases/${case_name}_eye_tune1_TetCal_NM.in"
  local shield_tuned="TetCal/batch_cases/${case_name}_shield_tune_TetCal_NM_tuned.in"
  local max_iterations=10
  local current_macro="${EYE_MACRO_DIR}/${case_name}_eye_iter1_TetCal_NM.in"

  if [[ ! -f "${eye_macro}" ]]; then
    echo "Eye tune macro not found: ${WORK_DIR}/${eye_macro}" >&2
    exit 1
  fi
  if [[ ! -f "${shield_tuned}" ]]; then
    echo "Shield tuned macro not found: ${WORK_DIR}/${shield_tuned}" >&2
    exit 1
  fi

  make_initial_eye_macro_with_shield_value "${eye_macro}" "${shield_tuned}" "${current_macro}"

  for (( iteration = 1; iteration <= max_iterations; iteration++ )); do
    local out_stem="${case_name}_retune_eye_after_shield_iter${iteration}"
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

case_name="Preinjection_I-131_1"

echo "Retuning currently unconverged preinjection/injection case from initial macros: ${case_name}"
run_macro "TetCal/batch_cases/${case_name}_shield_tune_TetCal_NM.in" "${case_name}_retune_shield"
run_eye_tune_from_start "${case_name}"
