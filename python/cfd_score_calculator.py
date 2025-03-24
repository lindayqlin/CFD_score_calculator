"""
Cutting frequency determination (CFD) score calculator
Written by Linda Lin 6/14/2022 (last modified on 3/23/2025)

Command line usage:
    python3 cfd_score_calculator.py <20nt_aligned_spacer_sequence>
    <20nt_aligned_protospacer_sequence> <last_2nt_of_PAM>
Example:
    python3 cfd_score_calculator.py CTAACAGTTGCTTTTATCAC tT-ACAGcTGCaTTTATCAC GG

Module usage:
    from cfd_score_calculator import calculate_cfd

Input: mismatch_scores.pkl & pam_scores.pkl
Output: CFD score
"""

from sys import argv
from pickle import load


def calculate_cfd(spacer, protospacer, pam):
    """Calculate CFD score
    Args:
        spacer [str]: 20 nt aligned spacer sequence (gRNA)
        protospacer [str]: 20 nt aligned protospacer sequence (DNA)
        pam [str]: last 2nt of PAM sequence
    Returns:
        score [float]: CFD score
    """
    # Check for expected input lengths
    if len(spacer) != 20 or len(protospacer) != 20 or len(pam) != 2:
        raise ValueError("Incorrect input sequence length\n" \
                         "Usage: python3 cfd_score_calculator.py " \
                         "<20nt_aligned_spacer_sequence> " \
                         "<20nt_aligned_protospacer_sequence> " \
                         "<last_2nt_of_PAM>")

    # Get empirical scoring matrices and pre-process sequences
    mm_scores, pam_scores = get_mm_pam_scores()
    spacer_list = list(spacer.upper().replace("T","U"))
    protospacer_list = list(protospacer.upper().replace("T","U"))

    # Calculate CFD score for alignment by nucleotide
    score = 1
    for i, nt in enumerate(protospacer_list):
        if spacer_list[i] == nt:
            # No penalty for perfect match
            score *= 1
        elif i == 0 and (spacer_list[i] == "-" or nt == "-"):
            # No penalty for gap at most PAM-distal nucleotide
            # Conservative choice to maximize sensitivity given lack of empirical data
            score *= 1
        else:
            # Incorporate score for given RNA-DNA basepair at this position
            key = f"r{spacer_list[i]}:d{reverse_complement_nt(nt)},{i + 1}"
            score *= mm_scores[key]
    # Incorporate PAM score
    score *= pam_scores[pam.upper()]

    return score


def get_mm_pam_scores():
    """Parse mismatch and PAM scoring matrices from pickled files
    Returns:
        mm_scores [dict]: mismatch scores
        pam_scores [dict]: PAM scores
    """
    try:
        mm_scores = load(open("mismatch_scores.pkl", "rb"))
        pam_scores = load(open("pam_scores.pkl", "rb"))
        return mm_scores, pam_scores
    except FileNotFoundError as e:
        raise FileNotFoundError("Requires mismatch_scores.pkl and " \
                                "pam_scores.pkl in directory") from e


def reverse_complement_nt(seq):
    """Get reverse complement of nucleotide (supports bulges)
    Args:
        seq [str]: nucleotide in caps
    Returns:
        rc [str]: reverse complement of nucleotide
    """
    return seq.translate(str.maketrans("ACTUG-", "TGAAC-"))


def main():
    # Get args
    try:
        spacer, protospacer, pam = argv[1:4]
    except ValueError as e:
        raise ValueError("Too few arguments\n" \
                         "Usage: python3 cfd_score_calculator.py " \
                         "<20nt_aligned_spacer_sequence> " \
                         "<20nt_aligned_protospacer_sequence> " \
                         "<last_2nt_of_PAM>") from e

    # Calculate CFD score
    print(calculate_cfd(spacer, protospacer, pam))


if __name__ == "__main__":
    main()
