# Day05: Stats of sequences
### Here I focused on two amino acid sequences of proteins that known to be muteted in breast cancer BRCA1 and BRCA2 expecially in Ashkenazi Jews

# Import Libraries
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Define functions
# Count amino acid occurance
def count_amino_acids(sequence):
    """
    Count occurrences of each amino acid in a given sequence.

    Parameters:
        sequence (str): A string representing a protein sequence.

    Returns:
        dict: A dictionary with amino acids as keys and their counts as values.
    """
    amino_acid_counts = {}
    for char in sequence:
        if char in amino_acid_counts:
            amino_acid_counts[char] += 1
        else:
            amino_acid_counts[char] = 1
    return amino_acid_counts

# Mode
# Calculate mode
def find_mode(amino_acid_counts):
    """
    Find the amino acid(s) with the highest frequency in a dictionary.
    
    Parameters:
        amino_acid_counts (dict): Dictionary with amino acids as keys and counts as values.
    
    Returns:
        list: A list of amino acid(s) with the highest frequency. If there are ties, return the first one.
    """
    max_count = max(amino_acid_counts.values())  # Find the maximum count
    # Find the amino acid(s) with max count and return the first one in case of ties
    mode = [aa for aa, count in amino_acid_counts.items() if count == max_count]
    return mode[:1]  # Return only the first mode if there are ties

# Main execution block
if __name__ == "__main__":
    folder_path = "C:/Course_Python_2024/python-assignments/day05"
    os.chdir(folder_path)

# Read the BRCA1 prtein sequence file
brca1_seq = pd.read_csv("BRCA1_protein_sequence.txt",header= 0) # correct sequence test PASS
#IMPORTANT!!!! 
# Remove the comment in the line below if you want to fail the test: 
#brca1_seq = pd.read_csv("BRCA1_protein_sequence_fake.txt",header= 0) # wrong sequence to make the test FAIL

print(brca1_seq.head())
type(brca1_seq)

# Check DataFrame dimentions
rows, columns = brca1_seq.shape
print(f"Number of rows: {rows}, Number of columns: {columns}")

# Convert to a single concatenated string
concatenated_vector_brca1 = "".join(brca1_seq.iloc[:,0].astype(str))
print(concatenated_vector_brca1)

# Get unique characters (amino acids) from the string
amino_acids= sorted(set(concatenated_vector_brca1)) 
num_amino_acids_brca1  = len(amino_acids)
print(amino_acids,"\n the number of unique amino acids is:", num_amino_acids_brca1)

# Apply the count function
amino_acid_counts_brca1 = count_amino_acids(concatenated_vector_brca1)
amino_acid_counts_brca1

# Visualize
# Convert counts to sorted lists for plotting
amino_acids_brca1 = sorted(amino_acid_counts_brca1.keys())
counts = [amino_acid_counts_brca1[aa] for aa in amino_acids_brca1]

# Plot the counts
plt.figure(figsize=(10, 6))
sns.barplot(x=amino_acids, y=counts, palette="viridis")

# Add labels and title
plt.xlabel("Amino Acid", fontsize=12)
plt.ylabel("Count", fontsize=12)
plt.title("Amino Acid Occurrences in BRCA1 protein", fontsize=15)
plt.show()

# Get the mode (most frequent amino acid(s))
mode_amino_acids_brca1 = find_mode(amino_acid_counts_brca1)
print("The most frequent amino acid in BRCA1 is:",mode_amino_acids_brca1)
