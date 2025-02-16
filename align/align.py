# Importing Dependencies
import numpy as np
from typing import Tuple

# Defining class for Needleman-Wunsch Algorithm for Global pairwise alignment
class NeedlemanWunsch:
    """ Class for NeedlemanWunsch Alignment

    Parameters:
        sub_matrix_file: str
            Path/filename of substitution matrix
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty

    Attributes:
        seqA_align: str
            seqA alignment
        seqB_align: str
            seqB alignment
        alignment_score: float
            Score of alignment from algorithm
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty
    """
    def __init__(self, sub_matrix_file: str, gap_open: float, gap_extend: float):
        # Init alignment and gap matrices
        self._align_matrix = None
        self._gapA_matrix = None
        self._gapB_matrix = None

        # Init matrices for backtrace procedure
        self._back = None
        self._back_A = None
        self._back_B = None

        # Init alignment_score
        self.alignment_score = 0

        # Init empty alignment attributes
        self.seqA_align = ""
        self.seqB_align = ""

        # Init empty sequences
        self._seqA = ""
        self._seqB = ""

        # Setting gap open and gap extension penalties
        self.gap_open = gap_open
        assert gap_open < 0, "Gap opening penalty must be negative."
        self.gap_extend = gap_extend
        assert gap_extend < 0, "Gap extension penalty must be negative."

        # Generating substitution matrix
        self.sub_dict = self._read_sub_matrix(sub_matrix_file) # substitution dictionary

    def _read_sub_matrix(self, sub_matrix_file):
        """
        DO NOT MODIFY THIS METHOD! IT IS ALREADY COMPLETE!

        This function reads in a scoring matrix from any matrix like file.
        Where there is a line of the residues followed by substitution matrix.
        This file also saves the alphabet list attribute.

        Parameters:
            sub_matrix_file: str
                Name (and associated path if not in current working directory)
                of the matrix file that contains the scoring matrix.

        Returns:
            dict_sub: dict
                Substitution matrix dictionary with tuple of the two residues as
                the key and score as value e.g. {('A', 'A'): 4} or {('A', 'D'): -8}
        """
        with open(sub_matrix_file, 'r') as f:
            dict_sub = {}  # Dictionary for storing scores from sub matrix
            residue_list = []  # For storing residue list
            start = False  # trigger for reading in score values
            res_2 = 0  # used for generating substitution matrix
            # reading file line by line
            for line_num, line in enumerate(f):
                # Reading in residue list
                if '#' not in line.strip() and start is False:
                    residue_list = [k for k in line.strip().upper().split(' ') if k != '']
                    start = True
                # Generating substitution scoring dictionary
                elif start is True and res_2 < len(residue_list):
                    line = [k for k in line.strip().split(' ') if k != '']
                    # reading in line by line to create substitution dictionary
                    assert len(residue_list) == len(line), "Score line should be same length as residue list"
                    for res_1 in range(len(line)):
                        dict_sub[(residue_list[res_1], residue_list[res_2])] = float(line[res_1])
                    res_2 += 1
                elif start is True and res_2 == len(residue_list):
                    break
        return dict_sub

    def align(self, seqA: str, seqB: str) -> Tuple[float, str, str]:
        
        self.seqA_align = ""
        self.seqB_align = ""
        self.alignment_score = 0
        self._seqA = seqA
        self._seqB = seqB

        # Get lengths of sequences
        n, m = len(seqA), len(seqB)

        # Initialize matrices with -infinity
        self._align_matrix = np.full((n+1, m+1), -np.inf)
        self._gapA_matrix = np.full((n+1, m+1), -np.inf)
        self._gapB_matrix = np.full((n+1, m+1), -np.inf)

        # Set starting score for main alignment matrix
        self._align_matrix[0, 0] = 0.0

        # Initialize gapA matrix (gaps in seqA)
        for j in range(1, m+1):
            self._gapA_matrix[0, j] = self.gap_open + (j) * self.gap_extend 

        # Initialize gapB matrix (gaps in seqB)
        for i in range(1, n+1):
            self._gapB_matrix[i, 0] = self.gap_open + (i) * self.gap_extend

        # Set first row and column of _align_matrix from gapA and gapB
        for j in range(1, m+1):
            self._align_matrix[0, j] = self._gapA_matrix[0, j]
        for i in range(1, n+1):
            self._align_matrix[i, 0] = self._gapB_matrix[i, 0]

        # Fill the matrices
        for i in range(1, n+1):
            for j in range(1, m+1):
                # Calculate scores for aligning residues
                match_score = self.sub_dict[(seqA[i-1], seqB[j-1])]
                align_score = self._align_matrix[i-1][j-1] + match_score

                # Calculate scores for inserting gaps in seqA
                gapA_open = self._align_matrix[i][j-1] + self.gap_open
                gapA_extend = self._gapA_matrix[i][j-1] + self.gap_extend
                gapA_score = max(gapA_open, gapA_extend)

                # Calculate scores for inserting gaps in seqB
                gapB_open = self._align_matrix[i-1][j] + self.gap_open
                gapB_extend = self._gapB_matrix[i-1][j] + self.gap_extend
                gapB_score = max(gapB_open, gapB_extend)

                # Update gap matrices
                self._gapA_matrix[i][j] = gapA_score
                self._gapB_matrix[i][j] = gapB_score

                # Update main alignment matrix
                self._align_matrix[i][j] = max(align_score, gapA_score, gapB_score)

        # Store the final alignment score
        self.alignment_score = self._align_matrix[n][m]

        return self._backtrace()


    def _backtrace(self) -> Tuple[float, str, str]:
        """
        TODO
        
        This function traces back through the back matrix created with the
        align function in order to return the final alignment score and strings.
        
        Parameters:
        	None
        
        Returns:
         	(alignment score, seqA alignment, seqB alignment) : Tuple[float, str, str]
         		the score and corresponding strings for the alignment of seqA and seqB
        """
        
        i, j = len(self._seqA), len(self._seqB) # This sets the indicies for the traceback

        alignA, alignB = [], [] # Initializze the alligned sequences

        while i > 0 or j > 0: # Continue traceback until you hit the topleft corner of the backtrace matrix
            # This should check if the current score is comming from the align matrix
            if i > 0 and j > 0 and np.isclose(
                self._align_matrix[i][j],
                self._align_matrix[i-1][j-1] + self.sub_dict[(self._seqA[i-1], self._seqB[j-1])]
        ):            

                alignA.append(self._seqA[i-1]) # Align residues from seqA and seqB
                alignB.append(self._seqB[j-1]) # Align residues from seqA and seqB
                i -= 1 # update index
                j -= 1 # update index

            # Check if the current score comes from the gapA matrix
            elif j > 0 and np.isclose(
                self._align_matrix[i][j],
                self._gapA_matrix[i][j]
        ):
                
                alignA.append('-') # Insert gap in seqA
                alignB.append(self._seqB[j-1]) # Align residues from seqB and gap
                j -= 1 # update index

            # Check if the current score comes from the gapB matrix
            elif i > 0 and np.isclose(
                self._align_matrix[i][j],
                self._gapB_matrix[i][j]
        ):
                # Insert gap in seqB
                alignA.append(self._seqA[i-1]) # Align residues from seqA and gap
                alignB.append('-') # Insert gap in seqB
                i -= 1 # update index

            else:
                # If no valid path is found, raise an error
                raise ValueError("Invalid traceback path")

        # This should reverse the alligned sequences for A and B so the order is correct (when we traceback the order is wrong)
        self.seqA_align = ''.join(reversed(alignA))
        self.seqB_align = ''.join(reversed(alignB))

        return (self.alignment_score, self.seqA_align, self.seqB_align)


def read_fasta(fasta_file: str) -> Tuple[str, str]:
    """
    DO NOT MODIFY THIS FUNCTION! IT IS ALREADY COMPLETE!

    This function reads in a FASTA file and returns the associated
    string of characters (residues or nucleotides) and the header.
    This function assumes a single protein or nucleotide sequence
    per fasta file and will only read in the first sequence in the
    file if multiple are provided.

    Parameters:
        fasta_file: str
            name (and associated path if not in current working directory)
            of the Fasta file.

    Returns:
        seq: str
            String of characters from FASTA file
        header: str
            Fasta header
    """
    assert fasta_file.endswith(".fa"), "Fasta file must be a fasta file with the suffix .fa"
    with open(fasta_file) as f:
        seq = ""  # initializing sequence
        first_header = True
        for line in f:
            is_header = line.strip().startswith(">")
            # Reading in the first header
            if is_header and first_header:
                header = line.strip()  # reading in fasta header
                first_header = False
            # Reading in the sequence line by line
            elif not is_header:
                seq += line.strip().upper()  # generating full sequence
            # Breaking if more than one header is provided in the fasta file
            elif is_header and not first_header:
                break
    return seq, header
