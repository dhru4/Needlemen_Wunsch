# 🧬 Pairwise Sequence Alignment - Needleman-Wunsch Algorithm

This project implements the **Needleman-Wunsch algorithm** for global sequence alignment of DNA or protein sequences.  
It compares two sequences, aligns them optimally, and outputs the aligned sequences, similarity score, and the full scoring matrix.

---

## 🚀 Features

- Manual sequence input or FASTA file input
- Full implementation of Needleman-Wunsch global alignment
- Dynamic programming matrix filling and traceback
- Outputs:
  - Aligned sequences
  - Similarity score
  - Full scoring matrix
- Results automatically saved to a structured `alignment_output.txt` file
- Clear, readable Python code structure using functions

---

## 📚 How It Works

1. **Input Sequences**  
   Users can manually type two sequences or load sequences from a FASTA file.

2. **Alignment**  
   A scoring matrix is built using:
   - Match reward (+1)
   - Mismatch penalty (-1)
   - Gap penalty (-2)
   Best alignment is found using dynamic programming.

3. **Traceback**  
   The optimal aligned sequences are reconstructed by tracing the highest scoring path through the matrix.

4. **Output**  
   Aligned sequences, similarity score, and scoring matrix are printed to the terminal and saved in `alignment_output.txt`.

---

## ⚙️ How to Run

1. Clone the repository or download the files.
2. Open a terminal and run:
   ```bash
   python needleman_wunsch.py

   or

   OPen the file in any code editor and run it.
