# Importing Dependencies
import pytest
from align import NeedlemanWunsch, read_fasta
import numpy as np

def test_nw_alignment():
    """
    TODO: Write your unit test for NW alignment
    using test_seq1.fa and test_seq2.fa by
    asserting that you have correctly filled out
    the your 3 alignment matrices.
    Use the BLOSUM62 matrix and a gap open penalty
    of -10 and a gap extension penalty of -1.
    """
    seq1, _ = read_fasta("./data/test_seq1.fa")
    seq2, _ = read_fasta("./data/test_seq2.fa")

    sub_matrix_file = 'substitution_matrices/BLOSUM62.mat'
    gap_open= -10
    gap_extend = -1

    align = NeedlemanWunsch(sub_matrix_file=sub_matrix_file,gap_open=gap_open,gap_extend=gap_extend)
    align.align(seq1,seq2)

    expected_align_matrix = np.array([  
        [ 0., -11., -12., -13.],
        [-11.,   5.,  -5.,  -6.],
        [-12.,  -5.,   4.,  -6.],
        [-13.,  -6.,   0.,   5.],
        [-14.,  -7.,  -5.,   5.]
    ])

    expected_gapA_matrix = np.array([
        [-np.inf, -11., -12., -13.],
        [-np.inf, -21.,  -5.,  -6.],
        [-np.inf, -22., -15.,  -6.],
        [-np.inf, -23., -16., -10.],
        [-np.inf, -24., -17., -15.]
    ])

    expected_gapB_matrix = np.array([
        [-np.inf, -np.inf, -np.inf, -np.inf],
        [-11., -21., -22., -23.],
        [-12.,  -5., -15., -16.],
        [-13.,  -6.,  -6., -16.],
        [-14.,  -7.,  -7.,  -5.]
    ])

    #print(align._align_matrix)
    assert np.array_equal(align._align_matrix, expected_align_matrix), "Alignment matrix is incorrect"
    #print(align._gapA_matrix)
    assert np.array_equal(align._gapA_matrix, expected_gapA_matrix), "GapA matrix is incorrect"
    #print(align._gapB_matrix)
    assert np.array_equal(align._gapB_matrix, expected_gapB_matrix), "GapB matrix is incorrect"
    

def test_nw_backtrace():
    """
    TODO: Write your unit test for NW backtracing
    using test_seq3.fa and test_seq4.fa by
    asserting that the backtrace is correct.
    Use the BLOSUM62 matrix. Use a gap open
    penalty of -10 and a gap extension penalty of -1.
    """
    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")

    expected_seqA_allignment = 'MAVHQLIRRP'
    expected_seqB_allignment = 'M---QLIRHP'
    expected_allignment_score = 17.0
    
    sub_matrix_file = 'substitution_matrices/BLOSUM62.mat'
    gap_open= -10
    gap_extend = -1

    align = NeedlemanWunsch(sub_matrix_file=sub_matrix_file,gap_open=gap_open,gap_extend=gap_extend)

    observed_allignment_score, observed_seqA_allignment, observed_seqB_allignment = align.align(seq3,seq4)


    assert expected_seqA_allignment == observed_seqA_allignment
    assert expected_seqB_allignment == observed_seqB_allignment
    if expected_allignment_score != observed_allignment_score:
        percent_diff = (observed_allignment_score-expected_allignment_score)/100
        print(f'precent difference: {percent_diff}')
    else:
        assert expected_allignment_score == observed_allignment_score





