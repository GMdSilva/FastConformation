CURRENTPATH=`pwd`
FENSEMBLEDIR="${CURRENTPATH}/fast_ensemble"
conda activate $FENSEMBLEDIR/fast_ensemble-conda
python -m fast_ensemble.jackhmmer_msa "$@"