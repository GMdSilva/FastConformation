import numpy as np
from string import ascii_uppercase, ascii_lowercase
import urllib.parse
import urllib.request
import time

def parse_a3m(a3m_lines=None, a3m_file=None, filter_qid=0.15, filter_cov=0.5, N=100000):
    """
    Parses an A3M file or lines and filters sequences based on sequence identity and coverage.

    Args:
        a3m_lines (list of str): Lines from an A3M file (optional).
        a3m_file (str): Path to an A3M file (optional).
        filter_qid (float): Minimum sequence identity to retain a sequence. Default is 0.15.
        filter_cov (float): Minimum coverage to retain a sequence. Default is 0.5.
        N (int): Maximum number of sequences to retain. Default is 100000.

    Returns:
        tuple: (sequences, deletion_matrices, names) where each is a list.
    """
    def seqid(a, b):
        return sum(c1 == c2 for c1, c2 in zip(a, b))

    def nongaps(a):
        return sum(c != "-" for c in a)

    def chk(seq, ref_seq):
        rL = len(ref_seq)
        L = nongaps(seq)
        return not (L > filter_cov * rL and seqid(seq, ref_seq) > filter_qid * L)

    rm_lower = str.maketrans('', '', ascii_lowercase)

    # prep inputs
    if a3m_lines is None:
        a3m_lines = open(a3m_file, "r")
    else:
        a3m_lines = a3m_lines.splitlines()

    # parse inputs
    n, nams, seqs, mtx = 0, [], [], []

    def do_filter():
        seq = seqs[-1].translate(rm_lower)
        if "_UPI" in nams[-1] or chk(seq, ref_seq):
            nams.pop()
            seqs.pop()
        else:
            # deletion matrix
            deletion_vec = []
            deletion_count = 0
            for j in seqs[-1]:
                if j.islower():
                    deletion_count += 1
                else:
                    deletion_vec.append(deletion_count)
                    deletion_count = 0
            mtx.append(deletion_vec)
            seqs[-1] = seq

    for line in a3m_lines:
        line = line.rstrip()
        if line.startswith(">"):
            if n == 1:
                ref_seq = seqs[0].translate(rm_lower)
            if n >= 1:
                # filter previous entry
                do_filter()
            # start new sequence entry
            nam = line.split()[0][1:]
            if "_" not in nam: nam = f"X_{nam}"
            nams.append(nam)
            seqs.append("")
            n += 1
        else:
            seqs[-1] += line

    # filter last entry
    do_filter()

    if len(seqs) > N + 1:
        print(f"found too many sequences ({len(seqs)}), taking the top{N} (sorted by qid)")
        sid = np.argsort([seqid(seq, ref_seq) for seq in seqs])[::-1][:N + 1]
        seqs = [seqs[i] for i in sid]
        mtx = [mtx[i] for i in sid]
        nams = [nams[i] for i in sid]
    return seqs[1:], mtx[1:], nams[1:]


def get_uni_jackhmmer(msa, mtx, lab, filter_qid=0.15, filter_cov=0.5):
    """
    Filters sequences to retain only UniProt entries from a multiple sequence alignment (MSA).

    Args:
        msa (list of str): List of sequences in the MSA.
        mtx (list of list of int): List of deletion matrices corresponding to the MSA.
        lab (list of str): List of labels corresponding to the MSA.
        filter_qid (float): Minimum sequence identity to retain a sequence. Default is 0.15.
        filter_cov (float): Minimum coverage to retain a sequence. Default is 0.5.

    Returns:
        tuple: Filtered (msa, mtx, lab) where each is a list.
    """
    lab_, msa_, mtx_ = [], [], []
    ref_seq = np.array(list(msa[0]))
    rL = len(ref_seq)
    for l, s, x in zip(lab[1:], msa[1:], mtx[1:]):
        if l.startswith("UniRef"):
            l = l.split("/")[0]
            if "_UPI" not in l:
                tar_seq = np.array(list(s))
                L = (tar_seq != "-").sum()
                qid = (ref_seq == tar_seq).sum()
                if L > filter_cov * rL and qid > filter_qid * L:
                    lab_.append(l)
                    msa_.append(s)
                    mtx_.append(x)
    return msa_, mtx_, lab_


def uni_num(ids):
    """
    Converts UniProt IDs to numerical representations.

    Args:
        ids (list of str): List of UniProt IDs.

    Returns:
        list of int: Numerical representations of the UniProt IDs.
    """
    pa = {a: 0 for a in ascii_uppercase}
    for a in ["O", "P", "Q"]: pa[a] = 1
    ma = [[{} for k in range(6)], [{} for k in range(6)]]
    for n, t in enumerate(range(10)):
        for i in [0, 1]:
            for j in [0, 4]: ma[i][j][str(t)] = n
    for n, t in enumerate(list(ascii_uppercase) + list(range(10))):
        for i in [0, 1]:
            for j in [1, 2]: ma[i][j][str(t)] = n
        ma[1][3][str(t)] = n
    for n, t in enumerate(ascii_uppercase):
        ma[0][3][str(t)] = n
        for i in [0, 1]: ma[i][5][str(t)] = n

    nums = []
    for uni in ids:
        p = pa[uni[0]]
        tot, num = 1, 0
        if len(uni) == 10:
            for n, u in enumerate(reversed(uni[-4:])):
                num += ma[p][n][u] * tot
                tot *= len(ma[p][n].keys())
        for n, u in enumerate(reversed(uni[:6])):
            num += ma[p][n][u] * tot
            tot *= len(ma[p][n].keys())
        nums.append(num)
    return nums


def map_retrieve(ids, call_uniprot=False):
    """
    Maps UniRef IDs to UniProt accession numbers.

    Args:
        ids (list of str): List of UniRef IDs.
        call_uniprot (bool): Whether to query UniProt for mapping information. Default is False.

    Returns:
        dict: Mapping from UniRef IDs to UniProt accession numbers.
    """
    if call_uniprot:
        mode = "NF100" if "UniRef100" in ids[0] else "NF90"
        url = 'https://www.uniprot.org/uploadlists/'
        out = []
        for i in range(0, len(ids), 5000):
            params = {
                'from': mode,
                'to': 'ACC',
                'format': 'tab',
                'query': " ".join(ids[i:i + 5000])
            }
            data = urllib.parse.urlencode(params)
            data = data.encode('utf-8')
            req = urllib.request.Request(url, data)
            with urllib.request.urlopen(req) as f:
                response = f.read()
            out += [line.split() for line in response.decode('utf-8').splitlines()]
            time.sleep(5)

        # combine mapping
        mapping = {}
        for i, j in out:
            if i != "From":
                if i not in mapping:
                    mapping[i] = [j]
                else:
                    mapping[i].append(j)
    else:
        mapping = {}

    for i in ids:
        if i not in mapping:
            mapping[i] = [i.split("_")[1]]

    return mapping


def hash_it(_seq, _lab, _mtx, call_uniprot=False):
    """
    Generates a hash for a given sequence and label.

    Args:
        _seq (list of str): List of sequences.
        _lab (list of str): List of labels corresponding to the sequences.
        _mtx (list of list of int): List of deletion matrices corresponding to the sequences.
        call_uniprot (bool): Whether to query UniProt for mapping information. Default is False.

    Returns:
        dict: Contains mappings of sequences, labels, and hashes.
    """
    if _seq is None or _lab is None:
        _seq, _lab = parse_a3m(a3m_lines)

    _lab_to_seq = {L: S for L, S in zip(_lab, _seq)}
    _lab_to_mtx = {L: M for L, M in zip(_lab, _mtx)}

    # call uniprot
    _lab_to_uni = map_retrieve(_lab, call_uniprot=call_uniprot)

    _uni_to_lab = {}
    for L, U in _lab_to_uni.items():
        for u in U: _uni_to_lab[u] = L

    _uni, __lab = [], []
    for U, L in _uni_to_lab.items():
        _uni.append(U)
        __lab.append(L)

    _hash = uni_num(_uni)
    _uni_to_hash = {u: h for u, h in zip(_uni, _hash)}
    _hash_to_lab = {h: l for h, l in zip(_hash, __lab)}

    _lab_to_hash = {}
    for L, U in _lab_to_uni.items():
        _lab_to_hash[L] = []
        for u in U: _lab_to_hash[L].append(_uni_to_hash[u])

    return {"_lab_to_seq": _lab_to_seq,
            "_lab_to_mtx": _lab_to_mtx,
            "_lab_to_hash": _lab_to_hash,
            "_hash_to_lab": _hash_to_lab}


def stitch(_hash_a, _hash_b, stitch_min=1, stitch_max=20, filter_id=None):
    """
    Stitches two hashed sequences together based on their alignment.

    Args:
        _hash_a (dict): First sequence hash information.
        _hash_b (dict): Second sequence hash information.
        stitch_min (int): Minimum allowed distance between aligned sequences. Default is 1.
        stitch_max (int): Maximum allowed distance between aligned sequences. Default is 20.
        filter_id (None): Placeholder for a potential filtering ID (not used).

    Returns:
        tuple: (sequences, deletion matrices) for the stitched sequences.
    """
    o = _stitch(_hash_a, _hash_b, stitch_min, stitch_max)
    return (*o["seq"], *o["mtx"])


def _stitch(_hash_a, _hash_b, stitch_min=1, stitch_max=20):
    """
    Internal function to stitch two sequences based on their hashes.

    Args:
        _hash_a (dict): First sequence hash information.
        _hash_b (dict): Second sequence hash information.
        stitch_min (int): Minimum allowed distance between aligned sequences. Default is 1.
        stitch_max (int): Maximum allowed distance between aligned sequences. Default is 20.

    Returns:
        dict: Contains sequences, deletion matrices, labels, and delta gene information.
    """
    _seq, _mtx, _lab, _delta_gene = [[], []], [[], []], [[], []], []
    TOTAL = len(_hash_a["_lab_to_hash"])
    H_A = np.asarray(list(_hash_a["_hash_to_lab"].keys()))
    H_B = np.asarray(list(_hash_b["_hash_to_lab"].keys()))

    def hit(h, H):
        h = np.asarray(h)
        match = np.abs(h[:, None] - H[None, :]).min(0)
        match_min = match.min()
        if match_min >= stitch_min and match_min <= stitch_max:
            return True, H[match.argmin()], match_min
        else:
            return False, None, None

    for n, (l_a, h_a) in enumerate(_hash_a["_lab_to_hash"].items()):
        chk_b, h_b, dg = hit(h_a, H_B)
        if chk_b:
            l_b = _hash_b["_hash_to_lab"][h_b]
            h_b = _hash_b["_lab_to_hash"][l_b]
            chk_c, h_c, _ = hit(h_b, H_A)
            if chk_c and _hash_a["_hash_to_lab"][h_c] == l_a:
                _seq[0].append(_hash_a["_lab_to_seq"][l_a])
                _mtx[0].append(_hash_a["_lab_to_mtx"][l_a])
                _lab[0].append(l_a)
                _seq[1].append(_hash_b["_lab_to_seq"][l_b])
                _mtx[1].append(_hash_b["_lab_to_mtx"][l_b])
                _lab[1].append(l_b)
                _delta_gene.append(dg)

    return {"seq": _seq,
            "mtx": _mtx,
            "lab": _lab,
            "delta_gene": _delta_gene}
