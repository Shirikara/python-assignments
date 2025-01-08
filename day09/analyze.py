'''
This program finds the longest repeating subsequence that appears twice within the sequence.

NCBI Integration:
The program uses Biopython's Entrez module to search NCBI's nucleotide database.
User is required to enter an email address in this format to comply with NCBI's usage policy.

Fetching the Sequence:
It performs a search using the query (e.g., accession number or gene name).
Retrieves the first matching sequence in FASTA format.

Optional : TANDEM repeats feature
'''

from Bio import Entrez, SeqIO
import argparse

def fetch_gene_sequence(gene_name, email):
    """
    Fetches the sequence of a gene of interest from NCBI using Entrez.

    Parameters:
        gene_name (str): The gene name or accession number.
        email (str): The email address for NCBI Entrez query compliance.

    Returns:
        str: The sequence of the gene.
    """
    Entrez.email = email
    search_handle = Entrez.esearch(db="nucleotide", term=gene_name, retmax=1)
    search_results = Entrez.read(search_handle)
    search_handle.close()

    if search_results["Count"] == "0":
        raise ValueError(f"No gene found for {gene_name}.")

    gene_id = search_results["IdList"][0]
    fetch_handle = Entrez.efetch(db="nucleotide", id=gene_id, rettype="fasta", retmode="text")
    gene_record = SeqIO.read(fetch_handle, "fasta")
    fetch_handle.close()

    return str(gene_record.seq)

def find_tandem_repeats(sequence, min_unit=2, max_unit=10):
    """
    Finds tandem repeats (STRs) in the sequence.

    Parameters:
        sequence (str): DNA sequence.
        min_unit (int): Minimum length of the repeating sequence.
        max_unit (int): Maximum length of the repeating sequence.

    Returns:
        list: List of tuples (repeating unit, start index, repeat count).
    """
    repeats = []
    n = len(sequence)
    for unit_length in range(min_unit, max_unit + 1):
        for i in range(n - unit_length):
            unit = sequence[i:i + unit_length]
            count = 0
            while sequence[i + count * unit_length: i + (count + 1) * unit_length] == unit:
                count += 1
            if count > 1:
                repeats.append((unit, i, count))
    return repeats

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

if __name__ == "__main__":
    # Ask for user input for the email address first
    email = input("Enter your email address (required for NCBI): ")

    # Ask for user input for the gene name
    gene_name = input("Enter the gene of interest (e.g., BRCA1): ")

    try:
        # Fetch gene sequence from NCBI
        sequence = fetch_gene_sequence(gene_name, email)
        print(f"Fetched sequence (first 50 bases): {sequence[:50]}...")

        # Ask if user wants to find the longest repeated subsequence
        duplicate_option = input("Do you want to find the longest repeated subsequence? (yes/no): ").strip().lower()
        if duplicate_option == 'yes':
            longest_subseq, length = find_longest_repeated_subsequence(sequence)
            print(f"Longest repeated subsequence: {longest_subseq}")
            print(f"Length: {length}")

        # Ask if user wants to find and print tandem repeats
        tandem_option = input("Do you want to find and print tandem repeats results? (yes/no): ").strip().lower()
        if tandem_option == 'yes':
            repeats = find_tandem_repeats(sequence)
            if repeats:
                print("\nTandem repeats found:")
                for repeat in repeats:
                    print(f"Repeating unit: {repeat[0]}, Start: {repeat[1]}, Count: {repeat[2]}")
            else:
                print("No tandem repeats found.")

    except ValueError as e:
        print(e)
    except Exception as e:
        print(f"An error occurred: {e}")
