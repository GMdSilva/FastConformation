import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import os
from urllib import request
from concurrent import futures
import pickle
import re

from tqdm import tqdm

from fast_ensemble.msa_generation import jackhmmer
from fast_ensemble.msa_generation import parsers
from fast_ensemble.msa_generation import colabfold as cf
from fast_ensemble.msa_generation import pairmsa

TQDM_BAR_FORMAT = '{l_bar}{bar}| {n_fmt}/{total_fmt} [elapsed: {elapsed} remaining: {remaining}]'

#######################################################################################################################################
# prep_inputs
#######################################################################################################################################

def prep_inputs(sequence, jobname="test", homooligomer="1", output_dir=None, clean=False, verbose=True):
    """
    Prepares the input sequence and parameters for MSA generation.

    Args:
        sequence (str): The protein sequence to be processed.
        jobname (str): The name of the job. Default is "test".
        homooligomer (str): A string specifying the number of homooligomers for each sequence segment.
                            Default is "1".
        output_dir (str): The directory where output files will be saved. If None, a default directory
                          is created based on the jobname and sequence hash.
        clean (bool): If True, cleans the output directory by removing existing files. Default is False.
        verbose (bool): If True, prints warnings and information during execution. Default is True.

    Returns:
        dict: A dictionary containing the processed inputs, including sequences, homooligomer information,
              and output directory.
    """
    # process inputs
    sequence = str(sequence)
    sequence = re.sub("[^A-Z:/]", "", sequence.upper())
    sequence = re.sub(":+", ":", sequence)
    sequence = re.sub("/+", "/", sequence)
    sequence = re.sub("^[:/]+", "", sequence)
    sequence = re.sub("[:/]+$", "", sequence)
    jobname = re.sub(r'\W+', '', jobname)
    homooligomer = str(homooligomer)
    homooligomer = re.sub("[:/]+", ":", homooligomer)
    homooligomer = re.sub("^[:/]+", "", homooligomer)
    homooligomer = re.sub("[:/]+$", "", homooligomer)

    if len(homooligomer) == 0: homooligomer = "1"
    homooligomer = re.sub("[^0-9:]", "", homooligomer)

    # define inputs
    I = {"ori_sequence": sequence,
         "sequence": sequence.replace("/", "").replace(":", ""),
         "seqs": sequence.replace("/", "").split(":"),
         "homooligomer": homooligomer,
         "homooligomers": [int(h) for h in homooligomer.split(":")],
         "msas": [], "deletion_matrices": []}

    # adjust homooligomer option
    if len(I["seqs"]) != len(I["homooligomers"]):
        if len(I["homooligomers"]) == 1:
            I["homooligomers"] = [I["homooligomers"][0]] * len(I["seqs"])
        else:
            if verbose:
                print("WARNING: Mismatch between number of breaks ':' in 'sequence' and 'homooligomer' definition")
            while len(I["seqs"]) > len(I["homooligomers"]):
                I["homooligomers"].append(1)
            I["homooligomers"] = I["homooligomers"][:len(I["seqs"])]
        I["homooligomer"] = ":".join([str(h) for h in I["homooligomers"]])

    # define full sequence being modelled
    I["full_sequence"] = ''.join([s * h for s, h in zip(I["seqs"], I["homooligomers"])])
    I["lengths"] = [len(seq) for seq in I["seqs"]]

    # prediction directory
    if output_dir is None:
        I["output_dir"] = 'prediction_' + jobname + '_' + cf.get_hash(I["full_sequence"])[:5]
    else:
        I["output_dir"] = output_dir
    os.makedirs(I["output_dir"], exist_ok=True)

    # delete existing files in working directory
    if clean:
        for f in os.listdir(I["output_dir"]):
            os.remove(os.path.join(I["output_dir"], f))

    if verbose and len(I["full_sequence"]) > 1400:
        print(
            f"WARNING: For a typical Google-Colab-GPU (16G) session, the max total length is ~1400 residues. You are at {len(I['full_sequence'])}!")
        print(
            f"Run Alphafold may crash, unless you trim to the protein(s) to a short length. (See trim options below).")

    if verbose:
        print(f"\nhomooligomer: {I['homooligomer']}")
        print(f"total_length: {len(I['full_sequence'])}")
        print(f"output_dir: {I['output_dir']}\n")

    return I


#######################################################################################################################################
# prep_msa
#######################################################################################################################################

def run_jackhmmer(sequence, prefix, jackhmmer_binary_path='jackhmmer', verbose=True, use_ramdisk=False):
    """
    Runs the jackhmmer tool to search for homologous sequences in a protein sequence database.

    Args:
        sequence (str): The query protein sequence.
        prefix (str): The prefix for output files.
        jackhmmer_binary_path (str): Path to the jackhmmer binary executable. Default is 'jackhmmer'.
        verbose (bool): If True, prints progress and information during execution. Default is True.
        use_ramdisk (bool): If True, uses a RAM disk for temporary storage. Default is False.

    Returns:
        tuple: A tuple containing the MSAs, deletion matrices, and sequence names.
    """
    fasta_path = f"{prefix}.fasta"
    with open(fasta_path, 'wt') as f:
        f.write(f'>query\n{sequence}')

    pickled_msa_path = f"{prefix}.jackhmmer.pickle"
    if os.path.isfile(pickled_msa_path):
        msas_dict = pickle.load(open(pickled_msa_path, "rb"))
        msas, deletion_matrices, names = (msas_dict[k] for k in ['msas', 'deletion_matrices', 'names'])
        full_msa = []
        for msa in msas:
            full_msa += msa
    else:
        # --- Find the closest source ---
        test_url_pattern = 'https://storage.googleapis.com/alphafold-colab{:s}/latest/uniref90_2021_03.fasta.1'
        ex = futures.ThreadPoolExecutor(3)

        def fetch(source):
            request.urlretrieve(test_url_pattern.format(source))
            return source

        fs = [ex.submit(fetch, source) for source in ['', '-europe', '-asia']]
        source = None
        for f in futures.as_completed(fs):
            source = f.result()
            ex.shutdown()
            break

        dbs = []

        num_jackhmmer_chunks = {'uniref90': 59, 'smallbfd': 17, 'mgnify': 71}
        total_jackhmmer_chunks = sum(num_jackhmmer_chunks.values())
        disable_tqdm = not verbose
        with tqdm(total=total_jackhmmer_chunks, bar_format=TQDM_BAR_FORMAT, disable=disable_tqdm) as pbar:
            def jackhmmer_chunk_callback(i):
                pbar.update(n=1)

            pbar.set_description('Searching uniref90')
            jackhmmer_uniref90_runner = jackhmmer.Jackhmmer(
                binary_path=jackhmmer_binary_path,
                database_path=f'https://storage.googleapis.com/alphafold-colab{source}/latest/uniref90_2021_03.fasta',
                get_tblout=True,
                num_streamed_chunks=num_jackhmmer_chunks['uniref90'],
                streaming_callback=jackhmmer_chunk_callback,
                z_value=135301051,
                use_ramdisk=use_ramdisk)
            dbs.append(('uniref90', jackhmmer_uniref90_runner.query(fasta_path)))

            pbar.set_description('Searching smallbfd')
            jackhmmer_smallbfd_runner = jackhmmer.Jackhmmer(
                binary_path=jackhmmer_binary_path,
                database_path=f'https://storage.googleapis.com/alphafold-colab{source}/latest/bfd-first_non_consensus_sequences.fasta',
                get_tblout=True,
                num_streamed_chunks=num_jackhmmer_chunks['smallbfd'],
                streaming_callback=jackhmmer_chunk_callback,
                z_value=65984053,
                use_ramdisk=use_ramdisk)
            dbs.append(('smallbfd', jackhmmer_smallbfd_runner.query(fasta_path)))

            pbar.set_description('Searching mgnify')
            jackhmmer_mgnify_runner = jackhmmer.Jackhmmer(
                binary_path=jackhmmer_binary_path,
                database_path=f'https://storage.googleapis.com/alphafold-colab{source}/latest/mgy_clusters_2019_05.fasta',
                get_tblout=True,
                num_streamed_chunks=num_jackhmmer_chunks['mgnify'],
                streaming_callback=jackhmmer_chunk_callback,
                z_value=304820129,
                use_ramdisk=use_ramdisk)
            dbs.append(('mgnify', jackhmmer_mgnify_runner.query(fasta_path)))

        # --- Extract the MSAs and visualize ---
        # Extract the MSAs from the Stockholm files.
        # NB: deduplication happens later in pipeline.make_msa_features.

        mgnify_max_hits = 501
        msas = []
        deletion_matrices = []
        names = []
        for db_name, db_results in dbs:
            unsorted_results = []
            for i, result in enumerate(db_results):
                msa, deletion_matrix, target_names = parsers.parse_stockholm(result['sto'])
                e_values_dict = parsers.parse_e_values_from_tblout(result['tbl'])
                e_values = [e_values_dict[t.split('/')[0]] for t in target_names]
                zipped_results = zip(msa, deletion_matrix, target_names, e_values)
                if i != 0:
                    # Only take query from the first chunk
                    zipped_results = [x for x in zipped_results if x[2] != 'query']
                unsorted_results.extend(zipped_results)
            sorted_by_evalue = sorted(unsorted_results, key=lambda x: x[3])
            db_msas, db_deletion_matrices, db_names, _ = zip(*sorted_by_evalue)
            if db_msas:
                if db_name == 'mgnify':
                    db_msas = db_msas[:mgnify_max_hits]
                    db_deletion_matrices = db_deletion_matrices[:mgnify_max_hits]
                    db_names = db_names[:mgnify_max_hits]
                msas.append(db_msas)
                deletion_matrices.append(db_deletion_matrices)
                names.append(db_names)
                msa_size = len(set(db_msas))
                print(f'{msa_size} Sequences Found in {db_name}')

            pickle.dump({"msas": msas,
                         "deletion_matrices": deletion_matrices,
                         "names": names}, open(pickled_msa_path, "wb"))
    return msas, deletion_matrices, names


def prep_msa(I, msa_method="mmseqs2", add_custom_msa=False, msa_format="fas",
             pair_mode="unpaired", pair_cov=50, pair_qid=20,
             hhfilter_loc="hhfilter", reformat_loc="reformat.pl", TMP_DIR="tmp",
             custom_msa=None, precomputed=None,
             mmseqs_host_url="https://a3m.mmseqs.com",
             verbose=True, use_ramdisk=False):
    """
    Prepares and processes MSAs for the given sequences using the specified MSA method.

    Args:
        I (dict): Dictionary containing input sequences and other parameters.
        msa_method (str): Method used to generate MSAs. Default is "mmseqs2".
        add_custom_msa (bool): Whether to add a custom MSA. Default is False.
        msa_format (str): The format of the MSA file. Default is "fas".
        pair_mode (str): Pairing mode for sequences. Can be "unpaired", "paired", or "unpaired+paired".
                         Default is "unpaired".
        pair_cov (int): Coverage threshold for pairing sequences. Default is 50.
        pair_qid (int): Identity threshold for pairing sequences. Default is 20.
        hhfilter_loc (str): Path to the hhfilter binary. Default is "hhfilter".
        reformat_loc (str): Path to the reformat.pl script. Default is "reformat.pl".
        TMP_DIR (str): Path to the temporary directory. Default is "tmp".
        custom_msa (str): Path to a custom MSA file (optional).
        precomputed (str): Path to a precomputed MSA file (optional).
        mmseqs_host_url (str): URL of the MMseqs2 server. Default is "https://a3m.mmseqs.com".
        verbose (bool): If True, prints progress and information during execution. Default is True.
        use_ramdisk (bool): If True, uses a RAM disk for temporary storage. Default is False.

    Returns:
        dict: The updated dictionary I containing the generated MSAs and deletion matrices.
    """
    # make temp directory
    os.makedirs(TMP_DIR, exist_ok=True)

    # clear previous inputs
    I["msas"] = []
    I["deletion_matrices"] = []

    _blank_seq = ["-" * L for L in I["lengths"]]
    _blank_mtx = [[0] * L for L in I["lengths"]]

    def _pad(ns, vals, mode):
        if mode == "seq": _blank = _blank_seq.copy()
        if mode == "mtx": _blank = _blank_mtx.copy()
        if isinstance(ns, list):
            for n, val in zip(ns, vals): _blank[n] = val
        else:
            _blank[ns] = vals
        if mode == "seq": return "".join(_blank)
        if mode == "mtx": return sum(_blank, [])

    if len(I["seqs"]) == 1 or "unpaired" in pair_mode:
        # gather msas
        for n, seq in enumerate(I["seqs"]):
            # tmp directory
            prefix = cf.get_hash(seq)
            prefix = os.path.join(TMP_DIR, prefix)

            print(f"Running jackhmmer")
            # run jackhmmer
            msas_, mtxs_, names_ = ([sum(x, ())] for x in run_jackhmmer(seq, prefix, use_ramdisk=use_ramdisk))

            # pad sequences
            for msa_, mtx_ in zip(msas_, mtxs_):
                msa, mtx = [I["sequence"]], [[0] * len(I["sequence"])]
                for s, m in zip(msa_, mtx_):
                    msa.append(_pad(n, s, "seq"))
                    mtx.append(_pad(n, m, "mtx"))

                I["msas"].append(msa)
                I["deletion_matrices"].append(mtx)

        # PAIR_MSA
        if len(I["seqs"]) > 1 and (pair_mode == "paired" or pair_mode == "unpaired+paired"):
            print("attempting to pair some sequences...")

            _data = []
            for a in range(len(I["seqs"])):
                print(f"prepping seq_{a}")
                _seq = I["seqs"][a]
                _prefix = os.path.join(TMP_DIR, cf.get_hash(_seq))

                _msas, _mtxs, _names = run_jackhmmer(_seq, _prefix, use_ramdisk=use_ramdisk)
                _msa, _mtx, _lab = pairmsa.get_uni_jackhmmer(_msas[0], _mtxs[0], _names[0],
                                                             filter_qid=pair_qid / 100,
                                                             filter_cov=pair_cov / 100)

                if len(_msa) > 1:
                    _data.append(pairmsa.hash_it(_msa, _lab, _mtx, call_uniprot=False))
                else:
                    _data.append(None)

            Ln = len(I["seqs"])
            O = [[None for _ in I["seqs"]] for _ in I["seqs"]]
            for a in range(Ln):
                if _data[a] is not None:
                    for b in range(a + 1, Ln):
                        if _data[b] is not None:
                            print(f"attempting pairwise stitch for {a} {b}")
                            O[a][b] = pairmsa._stitch(_data[a], _data[b])
                            _seq_a, _seq_b, _mtx_a, _mtx_b = (*O[a][b]["seq"], *O[a][b]["mtx"])

                            # filter to remove redundant sequences
                            ok = []
                            with open(f"{TMP_DIR}/tmp.fas", "w") as fas_file:
                                fas_file.writelines(
                                    [f">{n}\n{a + b}\n" for n, (a, b) in enumerate(zip(_seq_a, _seq_b))])
                            os.system(
                                f"{hhfilter_loc} -maxseq 1000000 -i {TMP_DIR}/tmp.fas -o {TMP_DIR}/tmp.id90.fas -id 90")
                            for line in open(f"{TMP_DIR}/tmp.id90.fas", "r"):
                                if line.startswith(">"): ok.append(int(line[1:]))

                            if verbose:
                                print(f"found {len(_seq_a)} pairs ({len(ok)} after filtering)")

                            if len(_seq_a) > 0:
                                msa, mtx = [I["sequence"]], [[0] * len(I["sequence"])]
                                for s_a, s_b, m_a, m_b in zip(_seq_a, _seq_b, _mtx_a, _mtx_b):
                                    msa.append(_pad([a, b], [s_a, s_b], "seq"))
                                    mtx.append(_pad([a, b], [m_a, m_b], "mtx"))
                                I["msas"].append(msa)
                                I["deletion_matrices"].append(mtx)

    # save MSA as pickle
    pickle.dump({"msas": I["msas"], "deletion_matrices": I["deletion_matrices"]},
                open(os.path.join(I["output_dir"], "msa.pickle"), "wb"))
    return I