"""MS1 stick spectra (negative mode) for every labelled LC peak in Figures 2-4 and S20-S24.

Peaks are located from the EIC of the label's [M-H]- near the figure's retention time, and the
spectrum is the sum of the MS1 scans across the UV peak's half-height width (UV time + offset).
Reads the OneDrive files only; writes to OUT.
"""
import csv
import re
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.signal import find_peaks, peak_widths

BACKEND = Path(r"C:\dev\MFP_Analysis_app\MFP_analysis_app\web\backend")
sys.path[:0] = [str(BACKEND), str(BACKEND.parents[1])]
from app.services.lcms_service import (  # noqa: E402
    registry, attach_uv_from_csv, summed_spectrum_in_rt_range, extracted_ion_chromatogram,
)
from app.services.lcms_pipeline import align_uv_to_ms  # noqa: E402

SRC = Path(r"C:\Users\owner\OneDrive - huji.ac.il\שולחן העבודה\Omer daniil paper results\files for review")
OUT = Path(r"C:\Users\owner\Documents\LCMS peak spectra 2026-11-10")

H2O, H, PROTON = 18.010565, 1.007825, 1.007276
RES = {"His": 155.069477, "tg": 91.993200, "Ala": 89.047678}

HIS_FULL = [  # Figure 4, His+tg no metal (superset of the His labels in Fig 3A-C and S20-S24)
    ("1His", 2.2), ("2His1tg", 4.4), ("1His1tg", 4.8), ("1tg", 7.2), ("1His1tg cyclic", 8.4),
    ("2His2tg_S-S", 11.0), ("2His2tg", 11.3), ("1His2tg_S-S", 11.8), ("2His2tg", 12.0),
    ("1His2tg", 12.3), ("2His4tg", 13.3), ("1His3tg", 14.0), ("2tg_S-S", 14.05), ("1His3tg", 14.6),
    ("2tg", 15.7), ("3tg_S-S", 16.9), ("4tg", 18.6),
]
HIS_BASIC = [("1His", 2.2), ("1His1tg", 4.8), ("1tg", 7.3), ("1His1tg cyclic", 8.8), ("1His2tg", 12.3), ("2tg", 15.7)]
ALA = [("1tg", 7.3), ("1tg1Ala", 12.6), ("1tg2Ala", 12.9), ("1tg2Ala", 13.2), ("2tg3Ala", 14.3),
       ("2tg", 15.2), ("4tg1Ala", 18.3), ("4tg", 18.6)]

SAMPLES = [  # (folder, file stem, tag used in file names, figures, labels)
    ("His+tg", "D15 His+tg (11.09.24)06", "His+tg_no-metal", "Fig3A-C_Fig4_FigS20-S24",
     HIS_FULL),
    ("His+Zn+tg", "D1(0.2)5 His+Zn+tg (11.09.24)16", "His+tg+Zn_0.2eq", "Fig4",
     [l for l in HIS_FULL if l[0] != "2His1tg"] + [("2His1tg", 4.4)]),
    ("His+Zn+tg", "D115 His+Zn+tg (11.09.24)10", "His+tg+Zn_1eq", "Fig3A_Fig4",
     [("1His", 2.2), ("1His1tg", 5.0), ("1tg", 7.2), ("1His1tg cyclic", 8.7), ("2His2tg_S-S", 11.0),
      ("1His2tg_S-S", 11.8), ("1His2tg", 12.4), ("1His3tg", 14.0), ("2tg_S-S", 14.05),
      ("1His3tg", 14.6), ("2tg", 15.7), ("3tg_S-S", 16.9), ("4tg", 18.6)]),
    ("His+Zn+tg", "D155 His+Zn+tg (11.09.24)22", "His+tg+Zn_5eq", "Fig4",
     [("1His", 2.2), ("1His1tg", 2.7), ("1His1tg", 5.1), ("1tg", 7.1), ("1His2tg_S-S", 11.8),
      ("2tg_S-S", 14.2), ("2tg", 15.7), ("3tg_S-S", 16.9), ("4tg", 18.6), ("5tg", 20.1)]),
    ("His+Co+tg", "D115 His+Co+tg (11.09.24)08", "His+tg+Co", "Fig3B", HIS_BASIC),
    ("His+Na+tg", "D115 His+Na+tg (11.09.24)07", "His+tg+Na", "Fig3C", HIS_BASIC),
    ("His+Li+tg", "DD115 HIs+Li+tg 10mM63", "His+tg+Li", "FigS20", HIS_BASIC),
    ("His+K+tg", "DD115 His+K+tg 10mM61", "His+tg+K", "FigS21", HIS_BASIC),
    ("His+Mg+tg", "DD115 HIs+Mg+tg 10mM64", "His+tg+Mg", "FigS22", HIS_BASIC),
    ("His+Ni+tg", "DD115 His+Ni+tg 10mM66", "His+tg+Ni", "FigS23", HIS_BASIC),
    ("His+Cu+tg", "D115 His+Cu+tg (11.09.24)09", "His+tg+Cu", "FigS24", HIS_BASIC),
    ("Ala+tg", "DD15 tg+Ala 1mM05", "Ala+tg_no-metal", "Fig2C_Fig3D", ALA),
    ("Ala+Zn+tg", "DD115 Ala+Zn+tg 1mM11", "Ala+tg+Zn", "Fig3D", ALA),
]


def neutral_mass(label: str) -> float:
    counts = {k: int(n) for n, k in re.findall(r"(\d+)(His|tg|Ala)", label)}
    n = sum(counts.values())
    total = sum(RES[k] * c for k, c in counts.items())
    if "S-S" in label:  # two thiol chains joined by a disulfide
        return total - (n - 2) * H2O - 2 * H
    m = total - (n - 1) * H2O
    return m - H2O if "cyclic" in label else m


def half_height_window(rt, y, apex_rt):
    i = int(np.argmin(np.abs(rt - apex_rt)))
    lo, hi = max(0, i - 40), min(len(y), i + 41)
    j = lo + int(np.argmax(y[lo:hi]))
    w = peak_widths(y, [j], rel_height=0.5)
    return float(np.interp(w[2][0], np.arange(len(rt)), rt)), float(np.interp(w[3][0], np.arange(len(rt)), rt)), float(rt[j])


def centroid(mz, inten, half=3):
    """One stick per local maximum on the 0.1 Da grid: m/z of the maximum, intensity summed over +-0.3 Da."""
    if mz.size == 0:
        return mz, inten
    k = np.rint(mz / 0.1).astype(np.int64)
    k0 = k.min()
    dense = np.zeros(k.max() - k0 + 1)
    np.add.at(dense, k - k0, inten)
    pk, _ = find_peaks(np.concatenate([[0], dense, [0]]), distance=5)
    pk -= 1
    c_int = np.array([dense[max(0, i - half): i + half + 1].sum() for i in pk])
    return (pk + k0) * 0.1, c_int


def plot(mz, inten, title, path_stem, label_min=5.0, max_labels=15):
    mz, inten = centroid(mz, inten)
    rel = 100 * inten / inten.max()
    keep = rel >= 1.0
    mz, rel = mz[keep], rel[keep]
    fig, ax = plt.subplots(figsize=(7.09, 3.5))  # 180 mm wide
    ax.vlines(mz, 0, rel, color="#1f3b73", lw=1.0)
    order = np.argsort(rel)[::-1]
    labelled = []
    for k in order:
        if rel[k] < label_min or len(labelled) >= max_labels:
            break
        if any(abs(mz[k] - m) < 2.5 for m in labelled):  # isotope / neighbour of a bigger stick
            continue
        labelled.append(mz[k])
        ax.text(mz[k], rel[k] + 1.5, f"{mz[k]:.1f}", ha="center", va="bottom", fontsize=7, rotation=90)
    ax.set_xlim(max(0, mz.min() - 20), mz.max() + 20)
    ax.set_ylim(0, 125)
    ax.set_xlabel("m/z")
    ax.set_ylabel("Relative intensity (%)")
    ax.set_title(title, fontsize=9)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(f"{path_stem}.svg")
    fig.savefig(f"{path_stem}.png", dpi=600)
    plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for folder, stem, tag, figs, labels in SAMPLES:
        d = SRC / folder
        state = registry.restore_from_path(f"s-{tag}", d / f"{stem}.mzML")
        csvs = sorted(d.glob(f"{stem}*.CSV"))
        uv = attach_uv_from_csv(state, csvs[0], filename=csvs[0].name)
        uv_rt, uv_y = np.asarray(uv.rt_min, float), np.asarray(uv.signal, float)
        metas = [m for m in state.index.ms1 if m.polarity == "negative"]
        ms_rt = np.array([m.rt_min for m in metas]); ms_tic = np.array([m.tic for m in metas])
        offset = align_uv_to_ms(ms_rt, ms_tic, uv_rt, uv_y, max_lag_min=1.0).offset_min
        sub = OUT / f"{figs}_{tag}"
        sub.mkdir(exist_ok=True)
        print(f"== {tag}: {len(metas)} neg scans, UV->MS offset {offset:+.3f} min")
        for n, (label, fig_rt) in enumerate(labels, 1):
            ion = neutral_mass(label) - PROTON
            eic = extracted_ion_chromatogram(state, ion, tolerance=0.5, polarity="negative")
            ert, ey = np.asarray(eic["rt_min"], float), np.asarray(eic["intensity"], float)
            near = (ert > fig_rt + offset - 0.6) & (ert < fig_rt + offset + 0.6)
            note = ""
            if near.any() and ey[near].max() > 0:
                pk, _ = find_peaks(np.where(near, ey, 0))
                pk = pk[ey[pk] > 0.2 * ey[near].max()] if pk.size else np.array([int(np.argmax(np.where(near, ey, 0)))])
                apex_ms = float(ert[pk[np.argmin(np.abs(ert[pk] - fig_rt - offset))]])
            else:
                apex_ms, note = fig_rt + offset, "no EIC signal; used figure RT"
            # window = the EIC peak at half height (0.03-0.3 min each side), so shoulders and small
            # peaks next to big ones get their own scans
            i = int(np.argmin(np.abs(ert - apex_ms)))
            half = ey[i] / 2
            a = i
            while a > 0 and ey[a - 1] > half and apex_ms - ert[a - 1] <= 0.3:
                a -= 1
            b = i
            while b < len(ey) - 1 and ey[b + 1] > half and ert[b + 1] - apex_ms <= 0.3:
                b += 1
            lo_ms, hi_ms = min(ert[a], apex_ms - 0.03), max(ert[b], apex_ms + 0.03)
            lo, hi, uv_apex = lo_ms - offset, hi_ms - offset, apex_ms - offset
            spec = summed_spectrum_in_rt_range(state, rt_min=lo_ms, rt_max=hi_ms, polarity="negative", bin_width=0.1)
            mz, inten = np.asarray(spec["mz"], float), np.asarray(spec["intensity"], float)
            # is the label's ion in the spectrum?
            cmz, cint = centroid(mz, inten)
            hit = cint[np.abs(cmz - ion) <= 0.5]
            rel_ion = 100 * hit.max() / cint.max() if hit.size else 0.0
            if any(r[1] == tag and r[3] == label and r[5] == round(uv_apex, 2) for r in rows):
                print(f"  {n:02d} {label}: same peak as an earlier label, skipped")
                continue
            name = f"{figs}_{tag}_peak{n:02d}_RT{uv_apex:.2f}_{label.replace(' ', '-')}"
            title = f"{tag.replace('_', ' ')}  |  UV peak {uv_apex:.2f} min  (MS {lo + offset:.2f}-{hi + offset:.2f} min)"
            plot(mz, inten, title, sub / name)
            rows.append([figs, tag, n, label, fig_rt, round(uv_apex, 2), round(lo + offset, 2), round(hi + offset, 2),
                         round(ion, 2), round(rel_ion, 1), int(spec["n_scans"]), note, name])
            print(f"  {n:02d} {label:<16} fig {fig_rt:5.2f}  UV {uv_apex:5.2f}  [M-H]- {ion:7.2f}  {rel_ion:5.1f}%  {note}")
        state.close_reader()
    with open(OUT / "index.csv", "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow(["figures", "sample", "peak", "figure label", "figure RT (min)", "UV apex (min)",
                    "MS window start", "MS window end", "expected [M-H]- m/z", "label ion % of base peak",
                    "scans summed", "note", "file name"])
        w.writerows(rows)


if __name__ == "__main__":
    main()
