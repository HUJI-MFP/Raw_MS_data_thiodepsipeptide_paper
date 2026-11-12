# LC-MS raw data: His / Ala + thioglycolic acid dry-down products

Raw LC-MS data and per-peak mass spectra behind the LC-MS figures of the thiodepsipeptide–metal manuscript: main text Figures 2–4 and SI Figures S20–S24.

## Contents

| Folder | What it holds |
| --- | --- |
| `raw_data/` | One folder per sample: the LC-MS run as **mzML** (converted from the Agilent .D files) and the UV (DAD) chromatogram as **CSV** |
| `peak_spectra/` | One MS spectrum (stick plot, PNG + SVG) for every peak labelled in the figures, grouped by figure and sample |
| `peak_spectra/index.csv` | One row per spectrum: figure, sample, peak number, label, retention times, summed MS window, expected *m/z*, file name |
| `peak_spectra/peak_spectra.py` | The script that produced the spectra |

## Samples

| Sample | Raw data folder | mzML file | Figures |
| --- | --- | --- | --- |
| His + tg, no metal | `His+tg` | `D15 His+tg (11.09.24)06` | 3A–C, 4, S20–S24 |
| His + tg + Zn 0.2 eq | `His+Zn+tg` | `D1(0.2)5 His+Zn+tg (11.09.24)16` | 4 |
| His + tg + Zn 1 eq | `His+Zn+tg` | `D115 His+Zn+tg (11.09.24)10` | 3A, 4 |
| His + tg + Zn 5 eq | `His+Zn+tg` | `D155 His+Zn+tg (11.09.24)22` | 4 |
| His + tg + Co | `His+Co+tg` | `D115 His+Co+tg (11.09.24)08` | 3B |
| His + tg + Na | `His+Na+tg` | `D115 His+Na+tg (11.09.24)07` | 3C |
| His + tg + Li | `His+Li+tg` | `DD115 HIs+Li+tg 10mM63` | S20 |
| His + tg + K | `His+K+tg` | `DD115 His+K+tg 10mM61` | S21 |
| His + tg + Mg | `His+Mg+tg` | `DD115 HIs+Mg+tg 10mM64` | S22 |
| His + tg + Ni | `His+Ni+tg` | `DD115 His+Ni+tg 10mM66` | S23 |
| His + tg + Cu | `His+Cu+tg` | `D115 His+Cu+tg (11.09.24)09` | S24 |
| Ala + tg, no metal | `Ala+tg` | `DD15 tg+Ala 1mM05` | 2C, 3D |
| Ala + tg + Zn | `Ala+Zn+tg` | `DD115 Ala+Zn+tg 1mM11` | 3D |

Each mzML file has a UV CSV with the same name (`….dx_DAD1A.CSV`).

## How the spectra were made

- Negative-mode ESI. For each labelled peak, the extracted-ion chromatogram (EIC) of the expected [M−H]⁻ ion (±0.5 Da) was taken near the figure's retention time, after aligning the UV and MS time axes.
- All MS scans across the EIC peak at half height were summed, binned at 0.1 Da and centroided. Peaks of 5% or more of the base peak are labelled with their *m/z*.
- Spectra are numbered to match the peak numbers in the figures. A peak that appears in several figures (the no-metal His + tg trace) is shown once.
- File names read `{figures}_{sample}_peak{NN}_RT{UV retention time}_{label}`.

Many products also appear as [M+Cl]⁻ (with its ³⁷Cl isotope), [M+HCOO]⁻ and non-covalent dimers ([2M−H]⁻, [2M+Cl]⁻), not only as [M−H]⁻.

The mzML files open in any open MS viewer, for example MZmine, OpenMS/TOPPView or pyteomics.
