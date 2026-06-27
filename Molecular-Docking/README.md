# 🔬 Molecular Docking Pipeline

Automated virtual screening pipeline using **AutoDock Vina**. Designed for high-throughput docking of natural compound libraries against prepared protein targets.

## Features
- Batch docking of entire compound libraries (PDBQT format)
- Auto-parsing of binding affinities from Vina output logs
- Ranked results CSV with all binding modes
- Hit flagging based on configurable affinity cutoff
- Detailed run log with per-compound status

## Quick Start
```bash
pip install -r requirements.txt

python autodock_vina_pipeline.py \
    --receptor protein.pdbqt \
    --ligands_dir ligands_pdbqt/ \
    --center_x -26.5 --center_y 15.2 --center_z -16.8 \
    --size_x 20 --size_y 20 --size_z 20 \
    --exhaustiveness 8 \
    --output results/ \
    --top_n 20
```

## Output
| File | Description |
|------|-------------|
| `results/screening_results_all.csv` | Full ranked compound list |
| `results/screening_results_top20.csv` | Top-20 hits |
| `results/docked_poses/` | Best pose PDBQT for each compound |
| `results/vina_logs/` | Raw Vina output for each run |

## See Also
- [SARS-CoV-2 Mpro Study](SARS-CoV2-Mpro-Study/README.md) — real application of this pipeline

---
**Author:** Sana Begum | Research Analyst, Era's Lucknow Medical College | [GitHub](https://github.com/Sanabegum09)
