from Bio import SeqIO # type: ignore
import sys

def extract_subsequence(fasta_file, start, end, output_file=None):
    """
    Extracts a subsequence from a FASTA file given start and end coordinates.
    
    Args:
        fasta_file (str): Path to the FASTA file.
        start (int): Start coordinate (1-based, inclusive).
        end (int): End coordinate (1-based, inclusive).
        output_file (str, optional): File to save the subsequence. If None, prints to stdout.
    """
    # Read the FASTA
    record = SeqIO.read(fasta_file, "fasta")
    
    # Convert to Python's 0-based indexing
    subseq = record.seq[start-1:end]
    
    # FASTA header for output
    header = f"{record.id}_{start}_{end}"
    fasta_output = f">{header}\n{subseq}\n"
    
    if output_file:
        with open(output_file, "w") as f:
            f.write(fasta_output)
    else:
        print(fasta_output)

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python extract_fasta_range.py <fasta_file> <start> <end> [output_file]")
        sys.exit(1)

    fasta_file = sys.argv[1]
    start = int(sys.argv[2])
    end = int(sys.argv[3])
    output_file = sys.argv[4] if len(sys.argv) > 4 else None

    extract_subsequence(fasta_file, start, end, output_file)
