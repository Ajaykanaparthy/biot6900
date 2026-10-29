# Module 4: Structural Druggability of LUAD Candidates

This assignment examines **SCN7A (Q01118; Module 2 rank 5)** and **MASP1 (P48740; rank 11)** using the existing AlphaFold models, PAE matrices, experimental structures, and retained fpocket outputs. These are computational candidates, not experimentally established LUAD therapeutic targets.

## Verified results

| Metric | SCN7A | MASP1 |
|---|---:|---:|
| Model residues | 1,682 | 699 |
| Mean pLDDT | 73.74 | 90.46 |
| Residues below pLDDT 50 | 14.74% | 2.72% |
| Experimental comparison | 7TJ8 chain A | 4IGD chain A |
| Identical, confident C-alpha pairs | 954 | 377 |
| RMSD over those pairs | 5.35 Angstrom | 1.68 Angstrom |
| Highest-druggability kept pocket | 98 | 1 |
| fpocket **Score** | -0.081 | 0.341 |
| fpocket **Druggability Score** | 0.982 | 0.687 |
| Pocket mean lining pLDDT | 85.06 | 70.91 |
| Pocket mean lining PAE | 3.88 Angstrom | 5.89 Angstrom |
| Confidence verdict | keep | keep |
| Final structural call | Structure supports a small-molecule pocket | Possible but uncertain |

The original `structure_assessment_luad.csv` values were independently reproduced and remain byte-for-byte unchanged. Its `best_pocket_score` means **Druggability Score**, which is distinct from fpocket's separate `Score` descriptor. The more explicit [verified summary](results/verified_structure_summary.csv) contains both quantities.

## Pocket confidence and provenance

The professor's rules are applied to every pocket in this order:

1. **Reject** when mean lining-residue pLDDT is below 70.
2. Otherwise, mark **suspect** when mean PAE among lining residues exceeds 10 Angstrom.
3. Otherwise, **keep**.

pLDDT is averaged once per unique lining residue. Mean PAE uses every ordered pair of those residues, including the diagonal, following the lab helper. A kept pocket can still include individual low-confidence residues; MASP1's best pocket only narrowly passes the mean-pLDDT threshold.

The retained results contain 165 SCN7A pockets (109 keep, 56 reject) and 50 MASP1 pockets (49 keep, 1 reject). No pocket meets the suspect condition after applying the pLDDT rule. All five highest-druggability pockets per target are shown, including SCN7A pocket 149, which has high druggability score 0.926 but is rejected for mean pLDDT 64.20.

These fpocket-format files were attributed to Northeastern Explorer in the existing coursework. **The original cluster job log, command, and fpocket version are unavailable.** Reuse is supported by exact residue/atom identity and coordinate agreement: every one of the 5,105 pocket-lining atom entries matches the preserved input model. This confirms the structural inputs represented by the cached pockets; it does not independently establish the historical execution host or version. No new local fpocket execution is claimed. Source hashes and checks are in [structural_input_provenance.json](results/structural_input_provenance.json).

## Interpretation and limitations

- **SCN7A:** mean inter-component PAE is 20.72 Angstrom, making relative placement uncertain. The provided `find_domains` helper identifies connected confidence components with broad sequence spans, not curated biological domains. A domain-only fit for the first component reproduces RMSD 3.37 Angstrom; the second has fewer than ten identical, confident experimental matches, so a reliable separate fit is unavailable. The global RMSD is sensitive to this placement and alignment coverage. The geometric pocket is not evidence of drug binding, selectivity, or LUAD dependency.
- **MASP1:** [4IGD](https://www.rcsb.org/structure/4IGD) is an engineered **R448Q zymogen catalytic-region fragment, UniProt residues 298-699**, not the full-length wild-type protein. Its RMSD validates only the aligned region. One confidence component does not imply one biological domain. The marginal pocket confidence supports a cautious structural call.
- **Ligand context:** [7TJ8](https://www.rcsb.org/structure/7TJ8) represents the human NaX/SCN7A channel in a membrane-associated complex. The provided ligand-site check finds lipids and glycans; these are not validated therapeutic inhibitors. Pocket 98 covers 42% of the mapped PEV A2401 contact residues. The 4IGD ligand check excludes glycerol as an additive and supplies no functional inhibitor-site validation.
- **p53 anchor:** retained Part 2 structures reproduce RMSD 0.68 Angstrom over 192 confident residues for AlphaFold p53 versus 1TUP chain B. C176, H179, C238, and C242 are zinc ligands under the distance rule; R248 and R273 contact DNA. This helps interpret language-model intolerance without proving mutation-specific therapeutic effects. See [anchor measurements](results/p53_anchor_verified.csv).
- Structural confidence, pocket geometry, and sequence similarity do not establish disease causality or clinical utility. Experimental binding, functional modulation, target selectivity, and LUAD-context validation remain necessary.

## Reproduce the existing notebook

Use the course `biot6900` environment with NumPy, pandas, matplotlib, nbformat, nbclient, and ipykernel. The original assignment notebook remains [BIOT6900_Module4_Assignment.ipynb](BIOT6900_Module4_Assignment.ipynb).

```powershell
conda activate biot6900
cd C:\Users\kanap\OneDrive\Desktop\biot6900\module_4
jupyter notebook BIOT6900_Module4_Assignment.ipynb
```

Set the working directory to `module_4`, then **Kernel -> Restart & Run All**. For automated fresh-kernel execution and comparison with the preserved original summary:

```powershell
python -B scripts/execute_assignment.py
```

The notebook automatically restores only missing files from [the retained-input ZIP](reproducibility/luad_structure_inputs.zip), checking any existing files instead of overwriting them. The ZIP preserves 452 existing source files, including both AlphaFold models/PAE matrices, experimental structures, pocket atom/vertex files, and p53 anchor inputs. `structures/` remains ignored as required by the professor; the portable ZIP supports reproduction without relying on new model downloads. A new fpocket run requires Linux/Explorer and fpocket, but is unnecessary for the verified retained-output analysis. A missing or incompatible cache raises a real error.

Fresh-kernel execution passed all 11 code cells with zero error outputs; four notebook plots are embedded. The execution record is [notebook_validation.json](results/notebook_validation.json). The original notebook, including its fpocket error and unfinished cells, is preserved under `archive/portfolio_deliverable_1_originals/`; no assignment questions or instructions were removed.

## Figures and supporting tables

- [SCN7A structure, pLDDT, PAE, and top pockets](Figures/scn7a_structure_confidence_pockets.png)
- [MASP1 structure, pLDDT, PAE, and top pockets](Figures/masp1_structure_confidence_pockets.png)
- [Figure captions](Figures/captions.md)
- [All 215 pocket assessments](results/pocket_quality_luad.csv)
- [Top five pockets per target](results/top_five_pockets_luad.csv)
- [Verified target-level summary](results/verified_structure_summary.csv)

The structure panels are actual C-alpha backbone traces colored by pLDDT, with retained pocket-lining residues highlighted. Existing AlphaFold structures and PyMOL/ChimeraX scripts were preserved. These artifacts support synthesis into the separate portfolio report; no final portfolio PDF was created.

Primary resources: [AlphaFold Q01118](https://alphafold.ebi.ac.uk/entry/Q01118), [AlphaFold P48740](https://alphafold.ebi.ac.uk/entry/P48740), [RCSB 7TJ8](https://www.rcsb.org/structure/7TJ8), [RCSB 4IGD](https://www.rcsb.org/structure/4IGD), and [fpocket](https://github.com/Discngine/fpocket).
