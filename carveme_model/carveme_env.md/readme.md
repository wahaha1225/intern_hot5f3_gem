Introduction

This document records a reproducible workflow to install CarveMe on macOS using a conda environment, configure SCIP/PySCIPOpt as the MILP solver, and run CarveMe to build a draft model from protein sequences followed by iterative gapfilling across multiple media.
Key point: SCIP and PySCIPOpt must be version-compatible, and the environment variable SCIPOPTDIR must point to a prefix that contains include/scip and lib (for example, /opt/anaconda3 or the current environment’s $CONDA_PREFIX).
Environment prerequisites

Have Anaconda/Miniconda (or Miniforge) installed and your target conda environment activated.
Optional: You may have a Homebrew SCIP (check with brew list scip), but it’s recommended to use the conda-forge build of SCIP and PySCIPOpt to avoid mixing toolchains.
Installation steps (based on your actual session)

1. Try installing CarveMe (may fail before the solver is configured)
bash:
$pip install carveme
# error (typically due to missing MILP solver dependencies, especially SCIP/PySCIPOpt)

2. Install dependencies (diamond, SCIP, PySCIPOpt) via conda-forge/bioconda
bash:
$conda install -c bioconda -c conda-forge diamond
$conda install -c conda-forge pyscipopt
$conda install -c conda-forge scip=9.2.1

3. If you installed the wrong versions or hit conflicts, remove and reinstall
bash：
$conda remove scip pyscipopt
conda install -c conda-forge scip=9.2.1 pyscipopt

Configure SCIPOPTDIR

Confirm SCIP path:
bash:
$which scip   # my output was: /opt/anaconda3/bin/scip

Check for headers:
bash:
$ls /opt/anaconda3/include/scip

Set the environment variable (use your actual prefix, my actual prefix is  /opt/anaconda3/bin/scip ):
bash:
$export SCIPOPTDIR=/opt/anaconda3

If using the current conda environment, you can also set:
bash:
$export SCIPOPTDIR="$CONDA_PREFIX"

(Optional) I Reinstall PySCIPOpt to match the current SCIP
bash:
$pip uninstall pyscipopt -y
$pip install pyscipopt --no-cache-dir --verbose

Install CarveMe (after dependencies are in place)
bash:
$pip install carveme

Verify the setup
bash:
which scip
scip -v
ls /opt/anaconda3/include/scip
python -c "import pyscipopt as p; print('PySCIPOpt OK')"
carve --help
gapfill --help
diamond -v

How to call SCIP (CarveMe and gapfill)

Specify the solver on the command line with --solver scip. Ensure SCIPOPTDIR is set correctly and both scip and pyscipopt work.
Example:
$carve NDPXLJ_6_hot5f3.faa -o output.xml --solver scip