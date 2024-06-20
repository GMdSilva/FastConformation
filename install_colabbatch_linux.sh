#!/bin/bash -e

type wget 2>/dev/null || { echo "wget is not installed. Please install it using apt or yum." ; exit 1 ; }

CURRENTPATH=`pwd`
COLABFOLDDIR="${CURRENTPATH}/decaf_e_dev"

mkdir -p "${COLABFOLDDIR}"
cd "${COLABFOLDDIR}"
wget -q -P . https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
bash ./Mambaforge-Linux-x86_64.sh -b -p "${COLABFOLDDIR}/conda"
rm Mambaforge-Linux-x86_64.sh

source "${COLABFOLDDIR}/conda/etc/profile.d/conda.sh"
export PATH="${COLABFOLDDIR}/conda/condabin:${PATH}"
conda update -n base conda -y
conda create -p "$COLABFOLDDIR/decaf_e_dev" -c conda-forge -c bioconda -c biocore \
    git python=3.10 openmm==7.7.0 pdbfixer \
    kalign2=2.04 hhsuite=3.3.0 mmseqs2=15.6f452 hmmer scikit-learn mdanalysis tqdm seaborn pandas scipy -y
conda activate "$COLABFOLDDIR/decaf_e_dev-conda"

# install ColabFold and Jaxlib
"$COLABFOLDDIR/decaf_e_dev-conda/bin/pip" install --no-warn-conflicts \
    "colabfold[alphafold-minus-jax] @ git+https://github.com/sokrypton/ColabFold"
"$COLABFOLDDIR/decaf_e_dev-conda/bin/pip" install "colabfold[alphafold]"
"$COLABFOLDDIR/decaf_e_dev-conda/bin/pip" install --upgrade "jax[cuda12]"==0.4.28
"$COLABFOLDDIR/decaf_e_dev-conda/bin/pip" install --upgrade tensorflow
"$COLABFOLDDIR/decaf_e_dev-conda/bin/pip" install silence_tensorflow
"$COLABFOLDDIR/decaf_e_dev-conda/bin/pip" install absl-py
"$COLABFOLDDIR/decaf_e_dev-conda/bin/pip" install pdb-tools

"$COLABFOLDDIR/decaf_e_dev-conda/bin/pip" install --no-warn-conflicts \
    "decaf_e_dev @ git+https://github.com/GMdSilva/decaf_e_dev"


# use 'agg' for non-GUI backend
cd "${COLABFOLDDIR}/decaf_e_dev-conda/lib/python3.10/site-packages/colabfold"
sed -i -e "s#from matplotlib import pyplot as plt#import matplotlib\nmatplotlib.use('Agg')\nimport matplotlib.pyplot as plt#g" plot.py
# modify the default params directory
sed -i -e "s#appdirs.user_cache_dir(__package__ or \"colabfold\")#\"${COLABFOLDDIR}/colabfold\"#g" download.py
# suppress warnings related to tensorflow
sed -i -e "s#from io import StringIO#from io import StringIO\nfrom silence_tensorflow import silence_tensorflow\nsilence_tensorflow()#g" batch.py
# remove cache directory
rm -rf __pycache__

# Download weights
"$COLABFOLDDIR/decaf_e_dev-conda/bin/python3" -m colabfold.download
echo "Download of alphafold2 weights finished."
echo "-----------------------------------------"
echo "Installation of ColabFold finished."
echo "Add ${COLABFOLDDIR}/colabfold-conda/bin to your PATH environment variable to run 'colabfold_batch'."
echo -e "i.e. for Bash:\n\texport PATH=\"${COLABFOLDDIR}/colabfold-conda/bin:\$PATH\""
echo "For more details, please run 'colabfold_batch --help'."

sudo apt install tm-align