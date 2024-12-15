# This code is a command line tool that can download data from NCBI similar to manualy searching and downloading data on the website
# https://www.ncbi.nlm.nih.gov/

#README:
'''
How to run the code in the terminal:
python NCBI.py --database DATABASE --term TERM --number NUMBER --csv CSV

List of available databases:
    ["pubmed", "protein", "nuccore", "nucleotide", "gene", "genome",
    "structure", "taxonomy", "snp", "dbvar", "biosample", "bioproject"]
'''

# Import modules
import os
import argparse
import csv
from datetime import datetime
from Bio import Entrez

# ---------------------- Non interactive example of search ----------------------#

# Entrez.email = "your_email@example.com"
# doc_id = "CAG30470.1" #SOX10 protein ID in homo sapiens
# handle = Entrez.efetch(db = "Protein", id = doc_id, rettype = "gb",retmode = "text")
# data = handle.read()
# handle.close()
# print(data)

# ---------------------- User Interactive part ----------------------#

# Define the function to insert email address
def get_email():
    """Prompt the user to input their email address for NCBI Entrez API."""
    email_address = input("Please enter your email address (required for NCBI API) \nfor example: your_email@example.com\n")
    while not email_address.strip():
        print("Email address cannot be empty.")
        email_address = input("Please enter your email address (required for NCBI API): ")
    return email_address


# Search function that retrieves the data according to the user input
def fetch_and_save_data(database, term, number):
    """Fetch data from NCBI and save each record in its own file."""
    try:
        # Search the database
        print(f"Searching {database} for term: {term}")
        handle = Entrez.esearch(db=database, term=term, retmax=number)
        record = Entrez.read(handle)
        handle.close()

        ids = record["IdList"]
        total_found = record["Count"]

        if not ids:
            print("No results found for the term.")
            return [], int(total_found)

        print(f"Found {total_found} results. Fetching up to {len(ids)} items...")
        filenames = []

        for idx, id_ in enumerate(ids):
            handle = Entrez.efetch(db=database, id=id_, rettype="fasta", retmode="text")
            data = handle.read()
            handle.close()

            # Save each record in its own file
            filename = f"{term}_{idx + 1}.fasta"
            with open(filename, "w") as f:
                f.write(data)
            filenames.append(filename)

        return filenames, int(total_found)

    except Exception as e:
        print(f"An error occurred: {e}")
        return [], 0

# Save the data in a CSV file format
def save_metadata_to_csv(csv_file, database, term, max_items, total_found):
    """Save metadata to a CSV file."""
    date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    file_exists = os.path.isfile(csv_file)

    with open(csv_file, mode="a", newline="") as f:
        writer = csv.writer(f)
        if not file_exists:
            writer.writerow(["date", "term", "database", "max", "total"])
        writer.writerow([date, term, database, max_items, total_found])

def main():
    # Ask for email at the start
    Entrez.email = get_email()

    parser = argparse.ArgumentParser(description="Download data from NCBI databases")
    parser.add_argument("--database", required=True, help="NCBI database to query (e.g., genome, nucleotide, protein)")
    parser.add_argument("--term", required=True, help="Search term for the database query")
    parser.add_argument("--number", type=int, default=10, help="Number of items to retrieve. Default is 10")
    parser.add_argument("--csv", default="metadata.csv", help="CSV file to save metadata. Default is 'metadata.csv'")

    args = parser.parse_args()

    filenames, total_found = fetch_and_save_data(args.database, args.term, args.number)

    if filenames:
        print("Files saved:")
        for filename in filenames:
            print(f"  {filename}")

    save_metadata_to_csv(args.csv, args.database, args.term, args.number, total_found)
    print(f"Metadata saved to {args.csv}")


if __name__ == "__main__":
    main()