**DNA Sequence Alignment: README**

Authors: Varvara Esina and Allan Mukkuzhi

**Overview**

This program implements the Sequence Alignment algorithm as defined in §6 
(Dynamic Programming). It finds the minimum cost of an alignment between two DNA 
sequences given a gap penalty δ and a 4×4 similarity matrix.

**Algorithm & Implementation**

The program uses a dynamic programming approach built on memoization to avoid 
redundant computation.

**Subproblem Definition**

For substrings s₁...sᵢ and t₁...tⱼ, the minimum alignment cost is opt(i, j). We store these values in a 2D memoization table M.

**Recurrence [§6]**

  opt(i, j) is calculated as the minimum of[cite: 1052]:
    1. α_{sᵢtⱼ} + opt(i-1, j-1)  (Pairing sᵢ and tⱼ)
    2. δ + opt(i-1, j)           (Gap in sequence t)
    3. δ + opt(i, j-1)           (Gap in sequence s)

**Traceback**

  Once the table M is filled, the program reconstructs the optimal alignment 
  by backtracking from M[m, n] to M[0, 0] to identify the specific pairings 
  and gaps.

**Runtime**

  - Filling the (m+1)×(n+1) table takes O(mn) time.
  - The traceback takes O(m+n) time.
  - Total time complexity: O(mn).

**Data Structures**

  - int[][] M: The shared table used for memoizing opt(i, j).
  - ArrayList: Used to store the characters and penalties during traceback.

**Java Classes Used**

  - java.util.Scanner: For tokenized parsing of the input file.
  - java.io.FileReader/BufferedReader: For efficient file handling.
  - java.util.ArrayList/Collections: To manage and reverse alignment results.

**Steps to Run**

1. Compile the program:
   javac DNAAlign.java

2. Execute the program with an input file:
   java DNAAlign input.txt

**Input File Format**

The input file should follow this structure:
  - <delta penalty>
  - <4x4 mismatch matrix α (rows A, C, G, T)>
  - <sequence s>
  - <sequence t>
