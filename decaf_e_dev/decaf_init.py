import json
import os

# Configuration dictionary
config = {
    "_comment": "GENERAL CONFIGURATIONS",
    "output_path": "sample_predictions",
    "jobname": "abl_wt",

    "_comment2": "CONFIGURATIONS TO BUILD MSA",
    "sequence_path": "abl1.fasta",
    "tmp_dir": "tmp",
    "homooligomers": 1,
    "use_ramdisk": True,

    "_comment3": "CONFIGURATIONS TO MAKE PREDICTIONS",
    "engine": "alphafold2",
    "msa_path": None,
    "msa_from": "mmseqs2",
    "seq_pairs": [
        [64, 128],
        [128, 256],
        [256, 512],
        [512, 1024]
    ],
    "seeds": 10,
    "platform": "cpu",
    "save_all": False,
    "models": [1, 2, 3, 4, 5],
    "recycles": 4,
    "subset_msa_to": None,

    "_comment4": "GENERAL ANALYSIS CONFIGURATIONS",
    "align_range": "backbone",
    "ref1d": None,
    "ref2d1": None,
    "ref2d2": None,

    "_comment5": "RMSF ANALYSIS CONFIGURATIONS",
    "detect_mobile": True,
    "peak_width": 3,
    "peak_prominence": 1,
    "peak_height": 2,
    "starting_residue": 200,

    "_comment6": "RMSD ANALYSIS CONFIGURATIONS",
    "analysis_range": "backbone and resid 358-365",
    "analysis_range_name": "aloop",

    "_comment7": "RMSD_2D OR TMSCORE_2D CONFIGURATIONS",
    "mode_results": None,
    "n_stdevs": 5,

    "_comment8": "TMSCORE CONFIGURATIONS",
    "slice_predictions": "backbone and resid 210-459",

    "_comment10": "PCA CONFIGURATIONS",
    "n_pca_clusters": 3,

    "_comment11": "TRAJECTORY SAVING CONFIGURATIONS",
    "reorder": "rmsd_1d",
    "traj_format": "pdb",

    "n_clusters": None
}

def main():
    # Extract output path
    output_path = config.get("output_path")

    # Create directory if it doesn't exist
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    # Save the config dictionary to config.json
    with open("config.json", "w") as config_file:
        json.dump(config, config_file, indent=4)

if __name__ == "__main__":
    main()
