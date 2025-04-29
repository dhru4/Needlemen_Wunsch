def read_fasta(filename):
    with open(filename, 'r') as f:
        sequences = []
        current_seq = ''
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_seq:
                    sequences.append(current_seq)
                    current_seq = ''
            else:
                current_seq += line
        if current_seq:
            sequences.append(current_seq)
    return sequences

def needleman_wunsch(seq1, seq2, match_score=1, mismatch_penalty=-1, gap_penalty=-2):
    m, n = len(seq1), len(seq2)

    # Initialize scoring matrix
    score_matrix = [[0 for _ in range(n + 1)] for _ in range(m + 1)]
    for i in range(m + 1):
        score_matrix[i][0] = i * gap_penalty
    for j in range(n + 1):
        score_matrix[0][j] = j * gap_penalty

    # Fill matrix
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match = score_matrix[i - 1][j - 1] + (match_score if seq1[i - 1] == seq2[j - 1] else mismatch_penalty)
            delete = score_matrix[i - 1][j] + gap_penalty
            insert = score_matrix[i][j - 1] + gap_penalty
            score_matrix[i][j] = max(match, delete, insert)

    # Traceback
    align1, align2 = '', ''
    i, j = m, n
    while i > 0 and j > 0:
        current = score_matrix[i][j]
        diagonal = score_matrix[i - 1][j - 1]
        up = score_matrix[i - 1][j]
        left = score_matrix[i][j - 1]

        if current == diagonal + (match_score if seq1[i - 1] == seq2[j - 1] else mismatch_penalty):
            align1 = seq1[i - 1] + align1
            align2 = seq2[j - 1] + align2
            i -= 1
            j -= 1
        elif current == up + gap_penalty:
            align1 = seq1[i - 1] + align1
            align2 = '-' + align2
            i -= 1
        else:
            align1 = '-' + align1
            align2 = seq2[j - 1] + align2
            j -= 1

    while i > 0:
        align1 = seq1[i - 1] + align1
        align2 = '-' + align2
        i -= 1
    while j > 0:
        align1 = '-' + align1
        align2 = seq2[j - 1] + align2
        j -= 1

    # Terminal output
    print("\nAligned Sequence 1:", align1)
    print("Aligned Sequence 2:", align2)
    print("\nSimilarity Score:", score_matrix[m][n])
    print("\nScoring Matrix:")
    for row in score_matrix:
        print(' '.join(f"{num:4d}" for num in row))

    # File output
    with open("alignment_output.txt", "w") as f:
        f.write("Aligned Sequence 1:\n" + align1 + "\n")
        f.write("Aligned Sequence 2:\n" + align2 + "\n\n")
        f.write("Similarity Score: " + str(score_matrix[m][n]) + "\n\n")
        f.write("Scoring Matrix:\n")
        for row in score_matrix:
            f.write(' '.join(f"{num:4d}" for num in row) + "\n")

    print("\n[Output has been saved to 'alignment_output.txt']")

# Main program
if __name__ == "__main__":
    print("Choose input method:")
    print("1. Manual sequence input")
    print("2. Load from FASTA file")

    method = input("Enter 1 or 2: ").strip()

    if method == '1':
        seq1 = input("Enter the first sequence: ").strip().upper()
        seq2 = input("Enter the second sequence: ").strip().upper()
        needleman_wunsch(seq1, seq2)

    elif method == '2':
        file1 = input("Enter the first FASTA filename: ").strip()
        seqs1 = read_fasta(file1)

        if len(seqs1) >= 2:
            needleman_wunsch(seqs1[0].upper(), seqs1[1].upper())
        elif len(seqs1) == 1:
            file2 = input("Only one sequence found. Enter a second FASTA filename: ").strip()
            seqs2 = read_fasta(file2)

            if len(seqs2) < 1:
                print("Second file has no sequence. Exiting.")
            else:
                needleman_wunsch(seqs1[0].upper(), seqs2[0].upper())
        else:
            print("No sequences found in the file. Exiting.")

    else:
        print("Invalid option. Please enter 1 or 2.")