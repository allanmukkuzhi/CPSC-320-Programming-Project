================================================================================
DNA Sequence Alignment — README
================================================================================
Authors: [Your Name] and [Partner Name]

--------------------------------------------------------------------------------
Overview
--------------------------------------------------------------------------------
This program finds the best (minimum edit-distance) alignment of two DNA
sequences given a gap penalty δ and a 4×4 mismatch-cost matrix α (indexed
over {A, C, G, T}).

The algorithm is the sequence-alignment dynamic programming procedure
align(s, t) from §6 (Dynamic Programming) of the course lecture notes.

--------------------------------------------------------------------------------
Algorithm — align(s, t)   [§6, Lecture Notes]
--------------------------------------------------------------------------------
Let s = s₁s₂…sₘ and t = t₁t₂…tₙ be the two input sequences.

SUBPROBLEM DEFINITION
  For i = 0, 1, …, m and j = 0, 1, …, n, define

      opt(i, j) = the minimum cost of an alignment of s₁…sᵢ and t₁…tⱼ.

  We store all values in a memoization table M: M[i][j] = opt(i, j).

BASE CASES
  opt(i, 0) = i · δ   (aligning s₁…sᵢ with the empty string needs i gaps)
  opt(0, j) = j · δ   (symmetrically)

RECURRENCE
  For i ≥ 1, j ≥ 1, any optimal alignment of s₁…sᵢ and t₁…tⱼ must end
  in one of three ways:

    1. sᵢ is paired with tⱼ (match or mismatch):
         cost = α_{sᵢtⱼ} + opt(i−1, j−1)

    2. sᵢ is unmatched — a gap is placed in t:
         cost = δ + opt(i−1, j)

    3. tⱼ is unmatched — a gap is placed in s:
         cost = δ + opt(i, j−1)

  Therefore:
    opt(i, j) = min( α_{sᵢtⱼ} + opt(i−1, j−1),
                     δ          + opt(i−1, j),
                     δ          + opt(i, j−1) )

  The answer is opt(m, n) = M[m][n].

TRACEBACK
  After filling M, we start at M[m][n] and retrace the choices that
  produced each opt value, working back to M[0][0], to reconstruct the
  full aligned strings and per-column penalty values.

RUNNING TIME
  Filling the (m+1)×(n+1) table M takes O(mn) time (each of the
  (m+1)(n+1) entries is computed in O(1) time by the recurrence).
  The traceback visits at most m+n cells: O(m+n) time.
  Total running time: O(mn).

DATA STRUCTURES
  int[][] M                       Memoization table (the DP table)
  ArrayList<Character> / <Integer> Traceback lists (reversed at the end)

JAVA APIS USED
  java.io.BufferedReader, java.io.FileReader  — efficient file reading
  java.util.Scanner                           — tokenized input parsing
  java.util.ArrayList, java.util.Collections  — traceback list + reversal

--------------------------------------------------------------------------------
Input file format
--------------------------------------------------------------------------------
  <delta>
  <row A of α>      (four space-separated integers)
  <row C of α>
  <row G of α>
  <row T of α>
  <sequence s>
  <sequence t>

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
Requires Java SE 8 or later. From the directory containing DNAAlign.java:

  javac DNAAlign.java

--------------------------------------------------------------------------------
Running
--------------------------------------------------------------------------------
  java DNAAlign input.txt

Expected output for the example above:

  The best alignment is

  A A A G T C T G A C
  A A C G T T T - A C
  0 0 1 0 0 2 0 2 0 0

  with the minimum edit distance of 5.

--------------------------------------------------------------------------------
Packaging as a JAR (for submission)
--------------------------------------------------------------------------------
  jar cvf xy-ab.jar DNAAlign.java DNAAlign.class README

Verify contents:
  jar tf xy-ab.jar

Run from the JAR:
  java -cp xy-ab.jar DNAAlign input.txt

================================================================================
