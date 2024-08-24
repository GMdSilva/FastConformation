CURRENTPATH=`pwd`
FENSEMBLEDIR="${CURRENTPATH}/fast_ensemble"
conda activate $FENSEMBLEDIR/fast_ensemble-conda
python -m fast_ensemble.mmseqs2_msa "$@"