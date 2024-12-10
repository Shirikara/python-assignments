# # Test sequence file
# #import pytest
# import pandas as pd
# from Day05_sequence_stats import count_amino_acids, find_mode

# # Load and process sequences
# def load_sequence(file_path):
#     seq_df = pd.read_csv(file_path, header=0)
#     return "".join(seq_df.iloc[:, 0].astype(str))

# # Correct and fake sequences
# correct_sequence = load_sequence("BRCA1_protein_sequence.txt")
# fake_sequence = load_sequence("BRCA1_protein_sequence_fake.txt")

# # Expected modes
# expected_correct_mode = ['S'] 
# expected_fake_mode = ['A']

# # Test function
# def test_find_mode():
#     # Correct sequence test
#     correct_counts = count_amino_acids(correct_sequence)
#     correct_mode = find_mode(correct_counts)
#     assert correct_mode == expected_correct_mode, \
#         f"Test failed for correct sequence: Expected {expected_correct_mode}, but got {correct_mode}."
    
#     # Fake sequence test
#     fake_counts = count_amino_acids(fake_sequence)
#     fake_mode = find_mode(fake_counts)
#     assert fake_mode == expected_fake_mode, \
#         f"Test failed for fake sequence: Expected {expected_fake_mode}, but got {fake_mode}."

# # Run the tests with pytest by executing the following command in the terminal:
# # pytest -v Test_sequence_file.py


#import pytest
import pandas as pd
from Day05_sequence_stats import count_amino_acids, find_mode

# Load sequences once outside the test function
def load_sequence(file_path):
    """
    Load a protein sequence from a file.

    Parameters:
        file_path (str): Path to the file containing the protein sequence.
    
    Returns:
        str: A concatenated protein sequence string.
    """
    brca1_seq = pd.read_csv(file_path, header=0)
    return "".join(brca1_seq.iloc[:, 0].astype(str))

# Correct and fake sequences
correct_sequence = load_sequence("BRCA1_protein_sequence.txt")
fake_sequence = load_sequence("BRCA1_protein_sequence_fake.txt")

# Test functions
def test_find_mode_correct():
    correct_counts = count_amino_acids(correct_sequence)
    correct_mode = find_mode(correct_counts)
    expected_correct_mode = ['S']  # Update based on actual mode
    assert correct_mode == expected_correct_mode, f"Expected {expected_correct_mode}, but got {correct_mode}"

def test_find_mode_fake():
    fake_counts = count_amino_acids(fake_sequence)
    fake_mode = find_mode(fake_counts)
    expected_fake_mode = ['A']  # Update based on actual mode
    assert fake_mode == expected_fake_mode, f"Expected {expected_fake_mode}, but got {fake_mode}"