#!/usr/bin/env python3
import csv
import io
import math
import os
import urllib.parse
import urllib.request
from collections import defaultdict


BASE_URL = "https://www-nds.iaea.org/relnsd/v1/data"
OUT_DIR = "iaea_spectra"

NUCLIDES = [
    {
        "label": "I-131",
        "api_nuclide": "131i",
        "parent_energy_keV": 0.0,
        "rad_types": ["g", "x", "e", "bm", "bp", "a"],
        "beta_types": ["bm", "bp"],
    },
    {
        "label": "Tc-99m",
        "api_nuclide": "99tc",
        "parent_energy_keV": 142.6836,
        "rad_types": ["g", "x", "e", "bm", "bp", "a"],
        "beta_types": ["bm", "bp"],
    },
    {
        "label": "F-18",
        "api_nuclide": "18f",
        "parent_energy_keV": 0.0,
        "rad_types": ["g", "x", "e", "bm", "bp", "a"],
        "beta_types": ["bm", "bp"],
    },
]


def fetch_csv(params):
    query = urllib.parse.urlencode(params)
    req = urllib.request.Request(
        f"{BASE_URL}?{query}",
        headers={"User-Agent": "Livechart/1.0"},
    )
    with urllib.request.urlopen(req, timeout=60) as response:
        return response.read().decode("ISO-8859-1")


def rows_from_csv(text):
    if not text.strip() or text.strip() == "0":
        return []
    return list(csv.DictReader(io.StringIO(text)))


def as_float(value, default=math.nan):
    if value is None:
        return default
    value = value.strip()
    if not value:
        return default
    try:
        return float(value)
    except ValueError:
        return default


def same_parent_energy(row, wanted_keV):
    return abs(as_float(row.get("p_energy"), 0.0) - wanted_keV) < 1.0e-3


def write_tsv(path, fieldnames, rows):
    with open(path, "w", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def radiation_particle(rad_type, row):
    if rad_type in ("g", "x"):
        return "gamma" if rad_type == "g" else "xray"
    if rad_type == "e":
        return "electron"
    if rad_type == "a":
        return "alpha"
    if rad_type == "bm":
        return "beta-"
    if rad_type == "bp":
        return "beta+"
    return rad_type


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    raw_dir = os.path.join(OUT_DIR, "raw")
    os.makedirs(raw_dir, exist_ok=True)

    summary_rows = []

    for nuclide in NUCLIDES:
        label = nuclide["label"]
        api_nuclide = nuclide["api_nuclide"]
        parent_energy = nuclide["parent_energy_keV"]

        discrete_rows = []
        beta_branch_rows = []
        beta_bin_rows = []

        for rad_type in nuclide["rad_types"]:
            text = fetch_csv(
                {
                    "fields": "decay_rads",
                    "nuclides": api_nuclide,
                    "rad_types": rad_type,
                }
            )
            raw_path = os.path.join(raw_dir, f"{label}_decay_rads_{rad_type}.csv")
            with open(raw_path, "w", newline="") as output:
                output.write(text)

            rows = [r for r in rows_from_csv(text) if same_parent_energy(r, parent_energy)]
            if rad_type in ("g", "x", "e", "a"):
                for row in rows:
                    energy_keV = as_float(row.get("energy"))
                    intensity_percent = as_float(row.get("intensity"))
                    if math.isnan(energy_keV) or math.isnan(intensity_percent):
                        continue
                    discrete_rows.append(
                        {
                            "nuclide": label,
                            "particle": radiation_particle(rad_type, row),
                            "radiation_code": row.get("type") or rad_type.upper(),
                            "shell": row.get("shell", ""),
                            "energy_keV": f"{energy_keV:.12g}",
                            "energy_MeV": f"{energy_keV / 1000.0:.12g}",
                            "yield_per_decay": f"{intensity_percent / 100.0:.12g}",
                            "intensity_percent": f"{intensity_percent:.12g}",
                            "parent_energy_keV": row.get("p_energy", ""),
                            "decay": row.get("decay", ""),
                            "decay_percent": row.get("decay_%", ""),
                            "daughter": row.get("d_symbol", ""),
                            "unc_energy_keV": row.get("unc_en", ""),
                            "unc_intensity_percent": row.get("unc_i", ""),
                            "source": "IAEA LiveChart decay_rads",
                        }
                    )
            elif rad_type in ("bm", "bp"):
                for row in rows:
                    beta_intensity = as_float(row.get("intensity_beta"))
                    if not math.isnan(beta_intensity):
                        beta_branch_rows.append(
                            {
                                "nuclide": label,
                                "particle": radiation_particle(rad_type, row),
                                "mean_energy_keV": row.get("mean_energy", ""),
                                "mean_energy_MeV": format_optional_mev(row.get("mean_energy")),
                                "endpoint_energy_keV": row.get("max_energy", ""),
                                "endpoint_energy_MeV": format_optional_mev(row.get("max_energy")),
                                "yield_per_decay": f"{beta_intensity / 100.0:.12g}",
                                "intensity_percent": f"{beta_intensity:.12g}",
                                "daughter_level_energy_keV": row.get("daughter_level_energy", ""),
                                "transition_type": row.get("transition_type", ""),
                                "log_ft": row.get("log_ft", ""),
                                "parent_energy_keV": row.get("p_energy", ""),
                                "decay": row.get("decay", ""),
                                "decay_percent": row.get("decay_%", ""),
                                "daughter": row.get("d_symbol", ""),
                                "source": "IAEA LiveChart decay_rads",
                            }
                        )
                    ec_intensity = as_float(row.get("intensity_ec"))
                    if not math.isnan(ec_intensity):
                        beta_branch_rows.append(
                            {
                                "nuclide": label,
                                "particle": "electron-capture",
                                "mean_energy_keV": "",
                                "mean_energy_MeV": "",
                                "endpoint_energy_keV": row.get("energy_ec", ""),
                                "endpoint_energy_MeV": format_optional_mev(row.get("energy_ec")),
                                "yield_per_decay": f"{ec_intensity / 100.0:.12g}",
                                "intensity_percent": f"{ec_intensity:.12g}",
                                "daughter_level_energy_keV": row.get("daughter_level_energy", ""),
                                "transition_type": row.get("transition_type", ""),
                                "log_ft": row.get("log_ft", ""),
                                "parent_energy_keV": row.get("p_energy", ""),
                                "decay": row.get("decay", ""),
                                "decay_percent": row.get("decay_%", ""),
                                "daughter": row.get("d_symbol", ""),
                                "source": "IAEA LiveChart decay_rads",
                            }
                        )

        for beta_type in nuclide["beta_types"]:
            text = fetch_csv(
                {
                    "fields": "bin_beta",
                    "nuclides": api_nuclide,
                    "rad_types": beta_type,
                }
            )
            raw_path = os.path.join(raw_dir, f"{label}_bin_beta_{beta_type}.csv")
            with open(raw_path, "w", newline="") as output:
                output.write(text)

            rows = [r for r in rows_from_csv(text) if same_parent_energy(r, parent_energy)]
            for row in rows:
                bin_keV = as_float(row.get("bin_en"))
                dn_de_keV = as_float(row.get("dn_de"))
                dn_de_nu_keV = as_float(row.get("dn_de_nu"))
                if math.isnan(bin_keV) or math.isnan(dn_de_keV):
                    continue
                beta_bin_rows.append(
                    {
                        "nuclide": label,
                        "particle": radiation_particle(beta_type, row),
                        "bin_energy_keV": f"{bin_keV:.12g}",
                        "bin_energy_MeV": f"{bin_keV / 1000.0:.12g}",
                        "dn_dE_per_decay_per_keV": f"{dn_de_keV:.12g}",
                        "dn_dE_per_decay_per_MeV": f"{dn_de_keV * 1000.0:.12g}",
                        "unc_dn_dE_per_decay_per_keV": row.get("unc_dn_de", ""),
                        "nu_dn_dE_per_decay_per_keV": "" if math.isnan(dn_de_nu_keV) else f"{dn_de_nu_keV:.12g}",
                        "nu_unc_dn_dE_per_decay_per_keV": row.get("unc_dn_de_nu", ""),
                        "parent_energy_keV": row.get("p_energy", ""),
                        "daughter": row.get("d_symbol", ""),
                        "source": "IAEA LiveChart bin_beta Betashape",
                    }
                )

        agg = defaultdict(float)
        agg_meta = {}
        for row in discrete_rows:
            key = (
                row["particle"],
                row["radiation_code"],
                row["shell"],
                row["energy_keV"],
            )
            agg[key] += float(row["yield_per_decay"])
            agg_meta[key] = row
        agg_rows = []
        for key, total_yield in sorted(agg.items(), key=lambda item: (item[0][0], float(item[0][3]))):
            meta = agg_meta[key]
            agg_rows.append(
                {
                    "nuclide": label,
                    "particle": key[0],
                    "radiation_code": key[1],
                    "shell": key[2],
                    "energy_keV": key[3],
                    "energy_MeV": meta["energy_MeV"],
                    "yield_per_decay": f"{total_yield:.12g}",
                    "intensity_percent": f"{total_yield * 100.0:.12g}",
                    "source": "IAEA LiveChart decay_rads, grouped by particle/code/shell/energy",
                }
            )

        prefix = os.path.join(OUT_DIR, label)
        write_tsv(
            f"{prefix}_per_decay_discrete.tsv",
            [
                "nuclide",
                "particle",
                "radiation_code",
                "shell",
                "energy_keV",
                "energy_MeV",
                "yield_per_decay",
                "intensity_percent",
                "parent_energy_keV",
                "decay",
                "decay_percent",
                "daughter",
                "unc_energy_keV",
                "unc_intensity_percent",
                "source",
            ],
            discrete_rows,
        )
        write_tsv(
            f"{prefix}_per_decay_discrete_aggregated.tsv",
            [
                "nuclide",
                "particle",
                "radiation_code",
                "shell",
                "energy_keV",
                "energy_MeV",
                "yield_per_decay",
                "intensity_percent",
                "source",
            ],
            agg_rows,
        )
        write_tsv(
            f"{prefix}_beta_branches.tsv",
            [
                "nuclide",
                "particle",
                "mean_energy_keV",
                "mean_energy_MeV",
                "endpoint_energy_keV",
                "endpoint_energy_MeV",
                "yield_per_decay",
                "intensity_percent",
                "daughter_level_energy_keV",
                "transition_type",
                "log_ft",
                "parent_energy_keV",
                "decay",
                "decay_percent",
                "daughter",
                "source",
            ],
            beta_branch_rows,
        )
        write_tsv(
            f"{prefix}_beta_continuum.tsv",
            [
                "nuclide",
                "particle",
                "bin_energy_keV",
                "bin_energy_MeV",
                "dn_dE_per_decay_per_keV",
                "dn_dE_per_decay_per_MeV",
                "unc_dn_dE_per_decay_per_keV",
                "nu_dn_dE_per_decay_per_keV",
                "nu_unc_dn_dE_per_decay_per_keV",
                "parent_energy_keV",
                "daughter",
                "source",
            ],
            beta_bin_rows,
        )

        discrete_yield_by_particle = defaultdict(float)
        for row in discrete_rows:
            discrete_yield_by_particle[row["particle"]] += float(row["yield_per_decay"])
        beta_branch_yield_by_particle = defaultdict(float)
        for row in beta_branch_rows:
            beta_branch_yield_by_particle[row["particle"]] += float(row["yield_per_decay"])
        beta_integral_by_particle = integrate_beta(beta_bin_rows)

        particles = sorted(
            set(discrete_yield_by_particle)
            | set(beta_branch_yield_by_particle)
            | set(beta_integral_by_particle)
        )
        for particle in particles:
            summary_rows.append(
                {
                    "nuclide": label,
                    "particle": particle,
                    "discrete_yield_sum_per_decay": f"{discrete_yield_by_particle.get(particle, 0.0):.12g}",
                    "beta_branch_yield_sum_per_decay": f"{beta_branch_yield_by_particle.get(particle, 0.0):.12g}",
                    "beta_continuum_integral_per_decay": f"{beta_integral_by_particle.get(particle, 0.0):.12g}",
                    "n_discrete_rows": sum(1 for r in discrete_rows if r["particle"] == particle),
                    "n_beta_branch_rows": sum(1 for r in beta_branch_rows if r["particle"] == particle),
                    "n_beta_bins": sum(1 for r in beta_bin_rows if r["particle"] == particle),
                }
            )

    write_tsv(
        os.path.join(OUT_DIR, "summary.tsv"),
        [
            "nuclide",
            "particle",
            "discrete_yield_sum_per_decay",
            "beta_branch_yield_sum_per_decay",
            "beta_continuum_integral_per_decay",
            "n_discrete_rows",
            "n_beta_branch_rows",
            "n_beta_bins",
        ],
        summary_rows,
    )

    with open(os.path.join(OUT_DIR, "README.md"), "w") as output:
        output.write(
            "# IAEA LiveChart per-decay spectra\n\n"
            "Generated from IAEA LiveChart API v1 on demand for I-131, Tc-99m, and F-18.\n\n"
            "Files:\n"
            "- `*_per_decay_discrete.tsv`: ungrouped discrete gamma, x-ray, Auger/conversion electron, and alpha lines.\n"
            "- `*_per_decay_discrete_aggregated.tsv`: discrete lines grouped by particle, radiation code, shell, and energy.\n"
            "- `*_beta_branches.tsv`: beta-/beta+/electron-capture branches from `decay_rads`.\n"
            "- `*_beta_continuum.tsv`: binned beta spectra from `bin_beta`; `dn_dE` is already per decay.\n"
            "- `raw/*.csv`: raw API responses.\n"
            "- `summary.tsv`: yield/integral checks by nuclide and particle.\n\n"
            "IAEA units are keV and percent. The generated files also include MeV and yield per decay.\n"
            "For beta continuum, `dn_dE_per_decay_per_keV` integrates over keV to the beta yield per decay.\n"
            "Tc-99m is selected from the `99tc` API response by `p_energy = 142.6836 keV`.\n"
        )


def format_optional_mev(value):
    number = as_float(value)
    if math.isnan(number):
        return ""
    return f"{number / 1000.0:.12g}"


def integrate_beta(rows):
    points_by_particle = defaultdict(list)
    for row in rows:
        points_by_particle[row["particle"]].append(
            (
                float(row["bin_energy_keV"]),
                float(row["dn_dE_per_decay_per_keV"]),
            )
        )
    integrals = {}
    for particle, points in points_by_particle.items():
        points.sort()
        total = 0.0
        for idx in range(len(points) - 1):
            x0, y0 = points[idx]
            x1, y1 = points[idx + 1]
            total += (x1 - x0) * (y0 + y1) / 2.0
        integrals[particle] = total
    return integrals


if __name__ == "__main__":
    main()
