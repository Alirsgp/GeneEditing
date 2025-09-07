def fetch_full_sequence(file_path: str) -> str:
    """Reads a FASTA file and returns the continuous sequence as a string."""
    with open(file_path, "r") as f:
        return "".join(line.strip() for line in f if not line.startswith(">"))

def find_guides_in_region(seq: str, start: int, end: int) -> list:
    """
    Finds all valid 20-nt guides in [start, end] followed by NGG PAM.
    Returns list of (guide, pam, position) with 1-based positions.
    """
    region = seq[start-1:end]
    guides = []
    for i in range(len(region) - 23):  # 20 nt guide + 3 nt PAM
        guide = region[i:i+20]
        pam = region[i+20:i+23]
        if len(pam) == 3 and pam[1:] == "GG":  # NGG PAM requirement
            guides.append((guide, pam, start + i))
    return guides

def make_forward_oligo_pcas9(guide: str) -> str:
    """Cloning oligo for pCas9: forward = AAAC + guide + G."""
    return "AAAC" + guide + "G"

def make_reverse_oligo_pcas9(guide: str) -> str:
    """Cloning oligo for pCas9: reverse = reverse complement of guide + CAAA."""
    complement = str.maketrans("ATCG", "TAGC")
    revcomp = guide.translate(complement)[::-1]
    return revcomp + "CAAA"

if __name__ == "__main__":
    fasta_file = "sequence.fasta"
    full_seq = fetch_full_sequence(fasta_file)

    # Define lacZα regions (including MCS) — modify as needed
    ranges = [(238, 395), (396, 452), (455, 682)]

    all_guides = [g for r in ranges for g in find_guides_in_region(full_seq, *r)]

    print("All NGG sites in specified regions with pCas9-compatible oligos:\n")
    for guide, pam, pos in all_guides:
        print(f"Position {pos}-{pos+19}: {guide} PAM={pam}")
        print(f"  Forward (pCas9): {make_forward_oligo_pcas9(guide)}")
        print(f"  Reverse (pCas9): {make_reverse_oligo_pcas9(guide)}\n")
