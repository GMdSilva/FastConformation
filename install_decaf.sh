#!/bin/bash -e

type wget 2>/dev/null || { echo "wget is not installed. Please install it using apt or yum." ; exit 1 ; }

CURRENTPATH=`pwd`
DECAFDIR="${CURRENTPATH}/decaf_e_dev"

mkdir -p "${DECAFDIR}"
cd "${DECAFDIR}"
wget -q -P . https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
bash ./Mambaforge-Linux-x86_64.sh -b -p "${DECAFDIR}/conda"
rm Mambaforge-Linux-x86_64.sh

source "${DECAFDIR}/conda/etc/profile.d/conda.sh"
export PATH="${DECAFDIR}/conda/condabin:${PATH}"
conda update -n base conda -y
conda create -p "$DECAFDIR/decaf_e_dev-conda" -c conda-forge -c bioconda -c biocore \
    git python=3.10 openmm==7.7.0 pdbfixer \
    kalign2=2.04 hhsuite=3.3.0 mmseqs2=15.6f452 hmmer scikit-learn mdanalysis seaborn scipy -y
conda activate "$DECAFDIR/decaf_e_dev-conda"

# Install additional Python packages for the GUI
"$DECAFDIR/decaf_e_dev-conda/bin/pip" install PyQt5 pandas matplotlib silence_tensorflow pyqtgraph

# Install ColabFold and Jaxlib
"$DECAFDIR/decaf_e_dev-conda/bin/pip" install --no-warn-conflicts \
    "colabfold[alphafold-minus-jax] @ git+https://github.com/GMdSilva/ColabFold"
"$DECAFDIR/decaf_e_dev-conda/bin/pip" install "colabfold[alphafold]"
"$DECAFDIR/decaf_e_dev-conda/bin/pip" install --upgrade "jax[cuda12]"==0.4.28
"$DECAFDIR/decaf_e_dev-conda/bin/pip" install --upgrade tensorflow
"$DECAFDIR/decaf_e_dev-conda/bin/pip" install pdb-tools

"$DECAFDIR/decaf_e_dev-conda/bin/pip" install --no-warn-conflicts \
    "decaf_e_dev @ git+https://github.com/GMdSilva/decaf_e_dev"

# use 'agg' for non-GUI backend
cd "${DECAFDIR}/decaf_e_dev-conda/lib/python3.10/site-packages/colabfold"
sed -i -e "s#from matplotlib import pyplot as plt#import matplotlib\nmatplotlib.use('Agg')\nimport matplotlib.pyplot as plt#g" plot.py
# modify the default params directory
sed -i -e "s#appdirs.user_cache_dir(__package__ or \"colabfold\")#\"${DECAFDIR}/colabfold\"#g" download.py
# suppress warnings related to tensorflow
sed -i -e "s#from io import StringIO#from io import StringIO\nfrom silence_tensorflow import silence_tensorflow\nsilence_tensorflow()#g" batch.py
# remove cache directory
rm -rf __pycache__

sudo apt install tm-align

# Download weights
"$DECAFDIR/decaf_e_dev-conda/bin/python3" -m colabfold.download
echo "Download of alphafold2 weights finished."
echo "-----------------------------------------"
echo "Installation of decaf_e_dev finished."
echo "Add ${DECAFDIR}/decaf_e_dev-conda/bin to your PATH environment variable to run 'decaf_e_dev'."
echo -e "i.e. for Bash:\n\texport PATH=\"${DECAFDIR}/decaf_e_dev-conda/bin:\$PATH\""
echo "For more details, please run 'decaf_e_dev --help'."

# Run the GUI
"$DECAFDIR/decaf_e_dev-conda/bin/python" -m decaf_e_dev.gui.run_gui
