# fmt: off

############################################
# imports
###########################################
import hashlib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects

from string import ascii_uppercase, ascii_lowercase
import hashlib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects
from string import ascii_uppercase, ascii_lowercase

# List of colors used in PyMOL
pymol_color_list = ["#33ff33", "#00ffff", "#ff33cc", "#ffff00", "#ff9999", "#e5e5e5", "#7f7fff", "#ff7f00",
                    "#7fff7f", "#199999", "#ff007f", "#ffdd5e", "#8c3f99", "#b2b2b2", "#007fff", "#c4b200",
                    "#8cb266", "#00bfbf", "#b27f7f", "#fcd1a5", "#ff7f7f", "#ffbfdd", "#7fffff", "#ffff7f",
                    "#00ff7f", "#337fcc", "#d8337f", "#bfff3f", "#ff7fff", "#d8d8ff", "#3fffbf", "#b78c4c",
                    "#339933", "#66b2b2", "#ba8c84", "#84bf00", "#b24c66", "#7f7f7f", "#3f3fa5", "#a5512b"]

pymol_cmap = matplotlib.colors.ListedColormap(pymol_color_list)
alphabet_list = list(ascii_uppercase + ascii_lowercase)

aatypes = set('ACDEFGHIKLMNPQRSTVWY')

def get_hash(x):
    """
    Generate a SHA-1 hash for a given string.

    Args:
        x (str): The input string to be hashed.

    Returns:
        str: The SHA-1 hash of the input string.
    """
    return hashlib.sha1(x.encode()).hexdigest()

def homooligomerize(msas, deletion_matrices, homooligomer=1):
    """
    Homooligomerizes the input MSAs (Multiple Sequence Alignments) and deletion matrices.

    Args:
        msas (list of lists): A list of MSAs.
        deletion_matrices (list of lists): A list of deletion matrices corresponding to the MSAs.
        homooligomer (int): The number of homooligomeric copies. Default is 1 (no homooligomerization).

    Returns:
        tuple: A tuple containing the homooligomerized MSAs and deletion matrices.
    """
    if homooligomer == 1:
        return msas, deletion_matrices
    else:
        new_msas = []
        new_mtxs = []
        for o in range(homooligomer):
            for msa, mtx in zip(msas, deletion_matrices):
                num_res = len(msa[0])
                L = num_res * o
                R = num_res * (homooligomer - (o + 1))
                new_msas.append(["-" * L + s + "-" * R for s in msa])
                new_mtxs.append([[0] * L + m + [0] * R for m in mtx])
        return new_msas, new_mtxs

# Keeping the typo for cross-compatibility
def homooliomerize(msas, deletion_matrices, homooligomer=1):
    """
    Homooligomerizes the input MSAs and deletion matrices.
    This function is a typo version of `homooligomerize` for cross-compatibility.

    Args:
        msas (list of lists): A list of MSAs.
        deletion_matrices (list of lists): A list of deletion matrices corresponding to the MSAs.
        homooligomer (int): The number of homooligomeric copies. Default is 1 (no homooligomerization).

    Returns:
        tuple: A tuple containing the homooligomerized MSAs and deletion matrices.
    """
    return homooligomerize(msas, deletion_matrices, homooligomer=homooligomer)

def homooligomerize_heterooligomer(msas, deletion_matrices, lengths, homooligomers):
    """
    Homooligomerizes the input MSAs and deletion matrices for heterooligomeric complexes.

    Args:
        msas (list of lists): A list of MSAs.
        deletion_matrices (list of lists): A list of deletion matrices corresponding to the MSAs.
        lengths (list of int): A list of lengths for each component in the complex.
        homooligomers (list of int): A list of homooligomeric copies for each component.

    Returns:
        tuple: A tuple containing the homooligomerized MSAs and deletion matrices.
    """
    if max(homooligomers) == 1:
        return msas, deletion_matrices
    elif len(homooligomers) == 1:
        return homooligomerize(msas, deletion_matrices, homooligomers[0])
    else:
        frag_ij = [[0, lengths[0]]]
        for length in lengths[1:]:
            j = frag_ij[-1][-1]
            frag_ij.append([j, j + length])

        mod_msas, mod_mtxs = [], []
        for msa, mtx in zip(msas, deletion_matrices):
            mod_msa, mod_mtx = [], []
            for n, (s, m) in enumerate(zip(msa, mtx)):
                _s, _m, _ok = [], [], []
                for i, j in frag_ij:
                    _s.append(s[i:j])
                    _m.append(m[i:j])
                    _ok.append(max([o != "-" for o in _s[-1]]))

                if n == 0:
                    mod_msa.append("".join([x * h for x, h in zip(_s, homooligomers)]))
                    mod_mtx.append(sum([x * h for x, h in zip(_m, homooligomers)], []))
                elif sum(_ok) == 1:
                    a = _ok.index(True)
                    for h_a in range(homooligomers[a]):
                        _blank_seq = [["-" * l] * h for l, h in zip(lengths, homooligomers)]
                        _blank_mtx = [[[0] * l] * h for l, h in zip(lengths, homooligomers)]
                        _blank_seq[a][h_a] = _s[a]
                        _blank_mtx[a][h_a] = _m[a]
                        mod_msa.append("".join(["".join(x) for x in _blank_seq]))
                        mod_mtx.append(sum([sum(x, []) for x in _blank_mtx], []))
                else:
                    for a in range(len(lengths) - 1):
                        if _ok[a]:
                            for b in range(a + 1, len(lengths)):
                                if _ok[b]:
                                    for h_a in range(homooligomers[a]):
                                        for h_b in range(homooligomers[b]):
                                            _blank_seq = [["-" * l] * h for l, h in zip(lengths, homooligomers)]
                                            _blank_mtx = [[[0] * l] * h for l, h in zip(lengths, homooligomers)]
                                            for c, h_c in zip([a, b], [h_a, h_b]):
                                                _blank_seq[c][h_c] = _s[c]
                                                _blank_mtx[c][h_c] = _m[c]
                                            mod_msa.append("".join(["".join(x) for x in _blank_seq]))
                                            mod_mtx.append(sum([sum(x, []) for x in _blank_mtx], []))
            mod_msas.append(mod_msa)
            mod_mtxs.append(mod_mtx)
        return mod_msas, mod_mtxs

def chain_break(idx_res, Ls, length=200):
    """
    Adds a large number to residue indices to indicate chain breaks in a sequence.

    Args:
        idx_res (ndarray): The array of residue indices.
        Ls (list of int): The lengths of different segments in the sequence.
        length (int): The value to add to the residue index at chain breaks.

    Returns:
        ndarray: The updated array of residue indices with chain breaks.
    """
    L_prev = 0
    for L_i in Ls[:-1]:
        idx_res[L_prev + L_i:] += length
        L_prev += L_i
    return idx_res

##################################################
# plotting functions
##################################################

def plot_plddt_legend(dpi=100):
    """
    Plots a legend for pLDDT (predicted Local Distance Difference Test) scores.

    Args:
        dpi (int): Dots per inch setting for the plot.

    Returns:
        matplotlib.pyplot: The plot object with the pLDDT legend.
    """
    thresh = ['plDDT:', 'Very low (<50)', 'Low (60)', 'OK (70)', 'Confident (80)', 'Very high (>90)']
    plt.figure(figsize=(1, 0.1), dpi=dpi)
    for c in ["#FFFFFF", "#FF0000", "#FFFF00", "#00FF00", "#00FFFF", "#0000FF"]:
        plt.bar(0, 0, color=c)
    plt.legend(thresh, frameon=False, loc='center', ncol=6, handletextpad=1, columnspacing=1, markerscale=0.5)
    plt.axis(False)
    return plt

def plot_ticks(Ls):
    """
    Plots tick marks indicating segment boundaries on a plot.

    Args:
        Ls (list of int): The lengths of different segments in the sequence.
    """
    Ln = sum(Ls)
    L_prev = 0
    for L_i in Ls[:-1]:
        L = L_prev + L_i
        L_prev += L_i
        plt.plot([0, Ln], [L, L], color="black")
        plt.plot([L, L], [0, Ln], color="black")
    ticks = np.cumsum([0] + Ls)
    ticks = (ticks[1:] + ticks[:-1]) / 2
    plt.yticks(ticks, alphabet_list[:len(ticks)])

def plot_confidence(plddt, pae=None, Ls=None, dpi=100):
    """
    Plots predicted confidence metrics (pLDDT and PAE) for a protein structure.

    Args:
        plddt (ndarray): Array of predicted Local Distance Difference Test (pLDDT) scores.
        pae (ndarray): Array of Predicted Aligned Error (PAE) scores (optional).
        Ls (list of int): The lengths of different segments in the sequence (optional).
        dpi (int): Dots per inch setting for the plot.

    Returns:
        matplotlib.pyplot: The plot object displaying the confidence metrics.
    """
    use_ptm = False if pae is None else True
    if use_ptm:
        plt.figure(figsize=(10, 3), dpi=dpi)
        plt.subplot(1, 2, 1)
    else:
        plt.figure(figsize=(5, 3), dpi=dpi)
    plt.title('Predicted lDDT')
    plt.plot(plddt)
    if Ls is not None:
        L_prev = 0
        for L_i in Ls[:-1]:
            L = L_prev + L_i
            L_prev += L_i
            plt.plot([L, L], [0, 100], color="black")
    plt.ylim(0, 100)
    plt.ylabel('plDDT')
    plt.xlabel('position')
    if use_ptm:
        plt.subplot(1, 2, 2)
        plt.title('Predicted Aligned Error')
        Ln = pae.shape[0]
        plt.imshow(pae, cmap="bwr", vmin=0, vmax=30, extent=(0, Ln, Ln, 0))
        if Ls is not None and len(Ls) > 1:
            plot_ticks(Ls)
        plt.colorbar()
        plt.xlabel('Scored residue')
        plt.ylabel('Aligned residue')
    return plt

def plot_msas(msas, ori_seq=None, sort_by_seqid=True, deduplicate=True, dpi=100, return_plt=True):
    """
    Plots Multiple Sequence Alignments (MSAs).

    Args:
        msas (list of lists): A list of MSAs to be plotted.
        ori_seq (str): The original sequence (optional).
        sort_by_seqid (bool): Whether to sort sequences by sequence identity (default: True).
        deduplicate (bool): Whether to remove duplicate sequences (default: True).
        dpi (int): Dots per inch setting for the plot.
        return_plt (bool): Whether to return the plot object (default: True).

    Returns:
        matplotlib.pyplot: The plot object displaying the MSAs, if return_plt is True.
    """
    if ori_seq is None:
        ori_seq = msas[0][0]
    seqs = ori_seq.replace("/", "").split(":")
    seqs_dash = ori_seq.replace(":", "").split("/")

    Ln = np.cumsum(np.append(0, [len(seq) for seq in seqs]))
    Ln_dash = np.cumsum(np.append(0, [len(seq) for seq in seqs_dash]))
    Nn, lines = [], []
    for msa in msas:
        msa_ = set(msa) if deduplicate else msa
        if len(msa_) > 0:
            Nn.append(len(msa_))
            msa_ = np.asarray([list(seq) for seq in msa_])
            gap_ = msa_ != "-"
            qid_ = msa_ == np.array(list("".join(seqs)))
            gapid = np.stack([gap_[:, Ln[i]:Ln[i + 1]].max(-1) for i in range(len(seqs))], -1)
            seqid = np.stack([qid_[:, Ln[i]:Ln[i + 1]].mean(-1) for i in range(len(seqs))], -1).sum(-1) / (
                        gapid.sum(-1) + 1e-8)
            non_gaps = gap_.astype(float)
            non_gaps[non_gaps == 0] = np.nan
            if sort_by_seqid:
                lines.append(non_gaps[seqid.argsort()] * seqid[seqid.argsort(), None])
            else:
                lines.append(non_gaps[::-1] * seqid[::-1, None])

    Nn = np.cumsum(np.append(0, Nn))
    lines = np.concatenate(lines, 0)

    if return_plt:
        plt.figure(figsize=(8, 5), dpi=dpi)
        plt.title("Sequence coverage")
    plt.imshow(lines,
               interpolation='nearest', aspect='auto',
               cmap="rainbow_r", vmin=0, vmax=1, origin='lower',
               extent=(0, lines.shape[1], 0, lines.shape[0]))
    for i in Ln[1:-1]:
        plt.plot([i, i], [0, lines.shape[0]], color="black")
    for i in Ln_dash[1:-1]:
        plt.plot([i, i], [0, lines.shape[0]], "--", color="black")
    for j in Nn[1:-1]:
        plt.plot([0, lines.shape[1]], [j, j], color="black")

    plt.plot((np.isnan(lines) == False).sum(0), color='black')
    plt.xlim(0, lines.shape[1])
    plt.ylim(0, lines.shape[0])
    plt.colorbar(label="Sequence identity to query")
    plt.xlabel("Positions")
    plt.ylabel("Sequences")
    if return_plt:
        return plt
