================================================================================
DNA Sequence Alignment — README
================================================================================
Authors: [Varvara Esina] and [Allan Mukkuzhi]

--------------------------------------------------------------------------------
Overview
--------------------------------------------------------------------------------
This program computes the best (minimum edit-distance) alignment of two DNA
sequences using the classic Needleman-Wunsch O(mn) dynamic programming
algorithm.  The user supplies a gap penalty δ, a 4×4 mismatch-penalty matrix
(rows/columns ordered A, C, G, T), and the two sequences in a plain-text file.

Algorithm summary
-----------------
  1. Build an (m+1) × (n+1) DP table where dp[i][j] holds the minimum edit
     distance for aligning s[0..i-1] with t[0..j-1].
     - Base cases:  dp[i][0] = i*delta,  dp[0][j] = j*delta
     - Recurrence:  dp[i][j] = min(
           dp[i-1][j-1] + mismatchPenalty(s[i], t[j]),   // align pair
           dp[i-1][j]   + delta,                          // gap in t
           dp[i][j-1]   + delta )                         // gap in s
  2. Trace back from dp[m][n] to reconstruct the aligned strings.
  3. Print the aligned sequences, per-column penalties, and edit distance.

Time complexity: O(m*n)   Space complexity: O(m*n)

Data structures used
--------------------
  - int[][] dp       : the DP table
  - ArrayList<Character>, ArrayList<Integer>  : traceback lists (reversed at end)

Java APIs used
--------------
  java.io.BufferedReader, java.io.FileReader   — efficient file reading
  java.util.Scanner                            — tokenized input parsing
  java.util.ArrayList, java.util.Collections   — traceback list and reversal

--------------------------------------------------------------------------------
Input file format
--------------------------------------------------------------------------------
<delta>
<row A of similarity matrix>
<row C of similarity matrix>
<row G of similarity matrix>
<row T of similarity matrix>
<sequence 1>
<sequence 2>

Example (input.txt):
  2
  0 1 3 0
  1 0 1 2
  3 1 0 4
  0 2 4 0
  AAAGTCTGAC
  AACGTTTAC

--------------------------------------------------------------------------------
Compiling
--------------------------------------------------------------------------------
Make sure you have Java SE 8 (or later) installed.  From the directory
containing DNAAlign.java, run:

  javac DNAAlign.java

This produces DNAAlign.class in the same directory.

--------------------------------------------------------------------------------
Running
--------------------------------------------------------------------------------
  java DNAAlign input.txt

Replace input.txt with the path to your input file.

Expected output for the example above:
  The best alignment is

  A A A G T C T G A C
  A A C G T T T - A C
  0 0 1 0 0 2 0 2 0 0

  with the minimum edit distance of 5.

--------------------------------------------------------------------------------
Packaging as a JAR (for submission)
--------------------------------------------------------------------------------
After compiling, create the JAR (replace "xy-ab" with your initials):

  jar cvf xy-ab.jar DNAAlign.java DNAAlign.class README

Verify its contents:
  jar tf xy-ab.jar

To run directly from the JAR:
  java -cp xy-ab.jar DNAAlign input.txt

================================================================================
