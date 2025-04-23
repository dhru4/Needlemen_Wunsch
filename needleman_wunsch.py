def needleman_wunsch(seq1, seq2, match_score=1, mismatch_penalty=-1, gap_penalty=-2):
    m, n = len(seq1), len(seq2)

    # Initialize the scoring matrix
    score_matrix = [[0 for _ in range(n + 1)] for _ in range(m + 1)]

    # Initialize first row and column with gap penalties
    for i in range(m + 1):
        score_matrix[i][0] = i * gap_penalty
    for j in range(n + 1):
        score_matrix[0][j] = j * gap_penalty

    # Fill the scoring matrix
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match = score_matrix[i - 1][j - 1] + (match_score if seq1[i - 1] == seq2[j - 1] else mismatch_penalty)
            delete = score_matrix[i - 1][j] + gap_penalty
            insert = score_matrix[i][j - 1] + gap_penalty
            score_matrix[i][j] = max(match, delete, insert)

    # Traceback to get the aligned sequences
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

    # Print results
    print("Aligned Sequence 1:", align1)
    print("Aligned Sequence 2:", align2)
    print("Similarity Score:", score_matrix[m][n])
    print("\nScoring Matrix:")
    for row in score_matrix:
        print(row)

# Example usage
if __name__ == "__main__":
    seq1 = input("Enter the first sequence: ").upper()
    seq2 = input("Enter the second sequence: ").upper()
    needleman_wunsch(seq1, seq2)
