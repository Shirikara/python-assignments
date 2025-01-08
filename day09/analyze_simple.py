'''
This program finds the longest repeating subsequence that appears twice within the sequence
Requirments: It accepts a file path to a FASTA or GenBank format file of a gene of interest.
'''
from Bio import SeqIO

def find_longest_repeated_subsequence(sequence):
    """
    Finds the longest repeated subsequence in the given sequence.

    Parameters:
        sequence (str): DNA or protein sequence.

    Returns:
        tuple: Longest repeated subsequence and its length.
    """
    n = len(sequence)
    suffix_array = sorted((sequence[i:], i) for i in range(n))

    def lcp(s1, s2):
        """Helper function to calculate longest common prefix of two strings."""
        common_length = 0
        while common_length < len(s1) and common_length < len(s2) and s1[common_length] == s2[common_length]:
            common_length += 1
        return common_length

    longest_subseq = ""
    max_length = 0

    for i in range(n - 1):
        common_length = lcp(suffix_array[i][0], suffix_array[i + 1][0])
        if common_length > max_length:
            max_length = common_length
            longest_subseq = suffix_array[i][0][:common_length]

    return longest_subseq, max_length

def process_file(file_path):
    """
    Processes a file in Fasta or GenBank format and finds the longest repeated subsequence.

    Parameters:
        file_path (str): Path to the file.
    """
    for record in SeqIO.parse(file_path, "fasta"):
        sequence = str(record.seq)
        longest_subseq, length = find_longest_repeated_subsequence(sequence)
        print(f"Longest repeated subsequence: {longest_subseq}")
        print(f"Length: {length}")
        return  # Process only the first sequence in the file

    for record in SeqIO.parse(file_path, "genbank"):
        sequence = str(record.seq)
        longest_subseq, length = find_longest_repeated_subsequence(sequence)
        print(f"Longest repeated subsequence: {longest_subseq}")
        print(f"Length: {length}")
        return  # Process only the first sequence in the file

if __name__ == "__main__":
    file_path = input("Enter the path to the FASTA or GenBank file: ")
    process_file(file_path)

