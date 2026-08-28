"""Analyse the OPC3 settings sweep with the diagnostics machinery, unchanged.

Usage: analyse_sweep.py TAG [TAG ...]   -> writes <tag>.json beside the run
       analyse_sweep.py --table          -> prints the comparison from the .json files
"""
import json
import re
import sys
from pathlib import Path

from water_benchmarking import diagnostics, experiment
from water_benchmarking.analysis import density_hov

RUN = Path("/ssd1_nas_md/water_benchmarking/opc3/settings")
REFERENCE = Path("/ssd1_nas_md/water_benchmarking/opc3/gromacs")

CONDITIONS = {
    "o3rf78": ("RF eps_rf = 78.4", 1.4, "OPC3's own eps as eps_rf"),
    "o3rf61": ("RF eps_rf = 61", 1.4, "the SPC value the main protocol uses"),
    "o3pme": ("PME", 1.4, "what OPC3 was parameterised under"),
}


def ns_per_day(log: Path) -> float:
    m = re.search(r"Performance:\s+([\d.]+)", log.read_text())
    return float(m.group(1)) if m else float("nan")


def analyse(tag: str) -> None:
    label, cutoff, changed = CONDITIONS[tag]
    d = diagnostics.Diagnostic(tag=tag, label=label, n_molecules=2048, cutoff=cutoff,
                               electrostatics=label, changed=changed)
    diagnostics.analyse(d, RUN, model="opc3")
    out = {"values": d.values, "uncertainties": d.uncertainties,
           "ns_per_day": ns_per_day(RUN / f"{tag}.log")}
    (RUN / f"{tag}.json").write_text(json.dumps(out, indent=1))
    print(tag, json.dumps({k: float(f"{v:.4g}") if isinstance(v, float) else v
                           for k, v in d.values.items()}))


def table() -> None:
    corr = density_hov.POLARISATION_CORRECTION["opc3"]
    # The 1.8 nm / 1 fs / RF61 main run, for the reference row.
    ref_perf = [float(x) for x in re.findall(r"Performance:\s+([\d.]+)",
                "".join(p.read_text() for p in REFERENCE.glob("w_opc3-*.out")))]
    ref = {"values": dict(density=994.399, hov=51.6563, diffusion=2.36989e-09,
                          tau2_HH=2.27839, dielectric=78.2918),
           "uncertainties": {},
           "ns_per_day": sum(ref_perf[-10:]) / len(ref_perf[-10:])}
    rows = [("1.8 nm, 1 fs, RF 61 (main)", ref)]
    for tag in CONDITIONS:
        f = RUN / f"{tag}.json"
        if f.exists():
            rows.append((f"1.4 nm, 2 fs, {CONDITIONS[tag][0]}", json.loads(f.read_text())))

    cols = [("density", "rho", 1, "{:.1f}"), ("hov", "dHvap-pol", 1, "{:.2f}"),
            ("diffusion", "D", 1e9, "{:.2f}"), ("tau2_HH", "tau2(HH)", 1, "{:.2f}"),
            ("dielectric", "eps", 1, "{:.1f}")]
    hdr = f"{'condition':30}" + "".join(f"{c[1]:>11}" for c in cols) + f"{'ns/day':>9}{'rel':>6}"
    print(hdr)
    base = rows[0][1]["ns_per_day"]
    for name, r in rows:
        v = r["values"]
        line = f"{name:30}"
        for key, _, scale, fmt in cols:
            x = v[key] - corr if key == "hov" else v[key]
            line += f"{fmt.format(x * scale):>11}"
        line += f"{r['ns_per_day']:>9.0f}{r['ns_per_day'] / base:>6.2f}"
        print(line)
    e = experiment.EXPERIMENT
    print(f"{'experiment':30}{e['density'].value:>11.1f}{e['hov'].value:>11.2f}"
          f"{e['diffusion'].value * 1e9:>11.2f}{e['tau2_HH'].value:>11.2f}{e['dielectric'].value:>11.1f}")


if __name__ == "__main__":
    if sys.argv[1] == "--table":
        table()
    else:
        for tag in sys.argv[1:]:
            analyse(tag)
