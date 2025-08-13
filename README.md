# intern_hot5f3_gem
This repository contains the code and results from my MS internship in the Segre Lab: comparative genomics and genome-scale metabolic modeling (GEM) of the marine bacterium Marinovum sp. Hot5_F3. The goal is to predict carbon utilization capabilities and, by comparing to the closely related strain Marinovum algicola DG898, generate experimentally testable hypotheses.

Core components

Two modeling pipelines: bottom‑up (KBase/OMEGGA) and top‑down (eggNOG + CarveMe).
Dynamic growth simulations in COMETS with qualitative comparisons to experimental growth curves.
Comparative genomics at the KO and pathway levels, including KEGG Mapper two‑organism overlays.
Reproducible workflows: scripted model building, gap‑fill auditing, FBA/COMETS simulations, and figure generation.
Key takeaways

After media‑specific, targeted gap‑filling, both pipelines produce consistent grow/no‑grow predictions on single‑carbon‑source media (glucose, lactate, pyruvate); COMETS dynamics reproduce the trend “glucose fastest, lactate slowest” observed experimentally.
The KO sets of Hot5_F3 and DG898 largely overlap, but differences are enriched in aromatic compound degradation, nitrogen metabolism/two‑component systems, pentose and glucuronate interconversions, nucleotide‑sugar pathways, and aromatic amino acid metabolism—suggesting testable hypotheses about substrate preferences and regulatory differences.