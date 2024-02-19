import csv
import sys

def read_blacklist(blacklist_file):
    with open(blacklist_file, 'r') as file:
        blacklist = file.read().splitlines()
    return blacklist

def filter_columns(input_csv, output_csv, blacklist):
    with open(input_csv, 'r') as infile, open(output_csv, 'w', newline='') as outfile:
        reader = csv.reader(infile)
        writer = csv.writer(outfile)
        header = next(reader)  # Read the header

        # Find indices of columns to keep
        indices_to_keep = [i for i, column in enumerate(header) if not any(word in column for word in blacklist)]

        # Write updated header
        updated_header = [header[i] for i in indices_to_keep]
        writer.writerow(updated_header)

        # Write rows with filtered columns
        for row in reader:
            filtered_row = [row[i] for i in indices_to_keep]
            writer.writerow(filtered_row)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python blacklist_filter.py input_csv blacklist_file output_csv")
        sys.exit(1)

    input_csv = sys.argv[1]
    blacklist_file = sys.argv[2]
    output_csv = sys.argv[3]

    blacklist = read_blacklist(blacklist_file)
    filter_columns(input_csv, output_csv, blacklist)
    print("Filtered columns based on the blacklist successfully.")
