import java.io.*;
import java.util.*;

/**
 * DNAAlign.java
 *
 * Reads a gap penalty delta (δ), a 4×4 mismatch-cost matrix (α), and two DNA
 * sequences from a plain-text file, then computes and prints the best alignment
 * and its minimum edit distance.
 *
 * The algorithm is the sequence-alignment dynamic programming procedure
 * {@code align(s, t)} from §6 of the course lecture notes. We define
 * {@code opt(i, j)} as the minimum cost of an alignment of the prefixes
 * s₁...sᵢ and t₁...tⱼ, and store all values in a memoization table M
 * (so {@code M[i][j] = opt(i, j)}). After filling M bottom-up, we traceback
 * from {@code M[m][n]} to reconstruct the actual aligned strings and
 * per-column penalties.
 *
 * Usage: {@code java DNAAlign <inputFile>}
 *
 * APIs used: {@code java.io.BufferedReader}, {@code java.io.FileReader},
 * {@code java.util.Scanner}, {@code java.util.ArrayList},
 * {@code java.util.Collections}.
 *
 * @author Varvara Esina
 * @author Allan Mukkuzhi
 * @version 1.0
 */
public class DNAAlign {

    // CONSTANTS

    /**
     * Canonical ordering of nucleotides used to index the mismatch-cost
     * matrix α: A=0, C=1, G=2, T=3.
     */
    private static final String NUCLEOTIDES = "ACGT";

    // FIELDS

    /** The gap penalty δ. Every unmatched character (gap) incurs this cost. */
    private final int delta;

    /**
     * The 4×4 mismatch-cost matrix α, indexed by {@link #indexOf(char)}.
     * {@code alpha[i][j]} is the cost α_{xy} of aligning nucleotide x
     * (index i) with nucleotide y (index j). Diagonal entries are 0
     * (no cost for a match); off-diagonal entries are positive penalties.
     * The matrix is symmetric: α_{xy} = α_{yx}.
     */
    private final int[][] alpha;

    /** The first DNA sequence s = s₁s₂...sₘ. */
    private final String s;

    /** The second DNA sequence t = t₁t₂...tₙ. */
    private final String t;

    // CONSTRUCTORS

    /**
     * Constructs a DNAAlign instance with all parameters needed for alignment.
     *
     * @param delta the gap penalty δ
     * @param alpha the 4×4 mismatch-cost matrix α (rows/columns in A,C,G,T order)
     * @param s     the first DNA sequence
     * @param t     the second DNA sequence
     */
    public DNAAlign(int delta, int[][] alpha, String s, String t) {
        this.delta = delta;
        this.alpha = alpha;
        this.s = s.toUpperCase();
        this.t = t.toUpperCase();
    }

    // CORE ALGORITHM - align(s, t) from §6 of the lecture notes

    /**
     * Computes and prints the best alignment of {@code s} and {@code t}
     * using the sequence-alignment procedure {@code align(s, t)} from §6
     * of the course lecture notes.
     *
     * Subproblem definition:
     * For i = 0, 1, 2, ..., m and j = 0, 1, 2, ..., n, we define opt(i, j) to be
     * the minimum cost of an alignment of the prefixes s₁...sᵢ and t₁...tⱼ.
     * We store these values in a memoization table M: {@code M[i][j] = opt(i, j)}.
     *
     * Base cases (from the lecture notes):
     * opt(i, 0) = i · δ — aligning s₁...sᵢ with the empty string requires i gaps.
     * opt(0, j) = j · δ — symmetrically for t₁...tⱼ against empty s.
     *
     * Recurrence (from the lecture notes):
     * For i ≥ 1, j ≥ 1, the value opt(i, j) is the minimum of:
     * 1. Match/Mismatch: α_{sᵢtⱼ} + opt(i-1, j-1)
     * 2. Gap in t: δ + opt(i-1, j)
     * 3. Gap in s: δ + opt(i, j-1)
     *
     * Traceback:
     * Starting at M[m][n] and working back to M[0][0], at each cell we
     * determine which of the three cases produced the stored opt value,
     * and record the corresponding alignment column.
     *
     * Running time:
     * The table M has (m+1)(n+1) entries; each is computed in O(1) time.
     * Filling M takes O(mn) time and the traceback takes O(m+n) time.
     */
    public void align() {
        int m = s.length();
        int n = t.length();

        // Step 1: Build the memoization table M where M[i][j] = opt(i, j).
        int[][] M = new int[m + 1][n + 1];

        // Base cases: opt(i, 0) = i·δ  and  opt(0, j) = j·δ
        for (int i = 0; i <= m; i++) M[i][0] = i * delta;
        for (int j = 0; j <= n; j++) M[0][j] = j * delta;

        // Fill bottom-up using the recurrence from the lecture notes
        for (int i = 1; i <= m; i++) {
            for (int j = 1; j <= n; j++) {
                int pairCost = alphaCost(s.charAt(i - 1), t.charAt(j - 1)) + M[i - 1][j - 1];
                int gapInT   = delta + M[i - 1][j];
                int gapInS   = delta + M[i][j - 1];
                M[i][j] = Math.min(pairCost, Math.min(gapInT, gapInS));
            }
        }

        // Step 2: Traceback from M[m][n] to reconstruct the alignment.
        List<Character> alignS    = new ArrayList<>();
        List<Character> alignT    = new ArrayList<>();
        List<Integer>   penalties = new ArrayList<>();

        int i = m, j = n;
        while (i > 0 || j > 0) {
            if (i > 0 && j > 0) {
                int pen = alphaCost(s.charAt(i - 1), t.charAt(j - 1));
                if (M[i][j] == pen + M[i - 1][j - 1]) {
                    alignS.add(s.charAt(i - 1));
                    alignT.add(t.charAt(j - 1));
                    penalties.add(pen);
                    i--; j--;
                } else if (M[i][j] == delta + M[i - 1][j]) {
                    alignS.add(s.charAt(i - 1));
                    alignT.add('-');
                    penalties.add(delta);
                    i--;
                } else {
                    alignS.add('-');
                    alignT.add(t.charAt(j - 1));
                    penalties.add(delta);
                    j--;
                }
            } else if (i > 0) {
                alignS.add(s.charAt(i - 1));
                alignT.add('-');
                penalties.add(delta);
                i--;
            } else {
                alignS.add('-');
                alignT.add(t.charAt(j - 1));
                penalties.add(delta);
                j--;
            }
        }

        Collections.reverse(alignS);
        Collections.reverse(alignT);
        Collections.reverse(penalties);

        // Step 3: Print result
        printAlignment(alignS, alignT, penalties, M[m][n]);
    }

    // HELPER METHODS

    /**
     * Returns the mismatch cost α_{ab} for aligning nucleotide {@code a}
     * with nucleotide {@code b}, as given by the matrix {@link #alpha}.
     * When {@code a == b} this returns 0; otherwise it returns the penalty.
     *
     * @param a a nucleotide character (A, C, G, or T)
     * @param b a nucleotide character (A, C, G, or T)
     * @return  α_{ab}, the alignment cost for the pair (a, b)
     */
    private int alphaCost(char a, char b) {
        return alpha[indexOf(a)][indexOf(b)];
    }

    /**
     * Returns the index of a nucleotide in the canonical ordering
     * A=0, C=1, G=2, T=3, used for matrix lookup.
     *
     * @param  c a nucleotide character (A, C, G, or T)
     * @return   the index 0–3
     * @throws IllegalArgumentException if {@code c} is not a valid nucleotide
     */
    private int indexOf(char c) {
        int idx = NUCLEOTIDES.indexOf(Character.toUpperCase(c));
        if (idx == -1) {
            throw new IllegalArgumentException("Invalid nucleotide: " + c);
        }
        return idx;
    }

    /**
     * Formats and prints the aligned sequences, per-column penalties, and
     * minimum edit distance to standard output.
     *
     * @param alignS    characters from s in the alignment
     * @param alignT    characters from t in the alignment
     * @param penalties per-column cost values
     * @param minCost   the minimum edit distance opt(m, n)
     */
    private void printAlignment(List<Character> alignS,
                                List<Character> alignT,
                                List<Integer>   penalties,
                                int             minCost) {
        System.out.println("The best alignment is\n");

        StringBuilder sbS = new StringBuilder();
        StringBuilder sbT = new StringBuilder();
        StringBuilder sbP = new StringBuilder();

        for (int k = 0; k < alignS.size(); k++) {
            if (k > 0) { sbS.append(' '); sbT.append(' '); sbP.append(' '); }
            sbS.append(alignS.get(k));
            sbT.append(alignT.get(k));
            sbP.append(penalties.get(k));
        }

        System.out.println(sbS);
        System.out.println(sbT);
        System.out.println(sbP);
        System.out.println("\nwith the minimum edit distance of " + minCost + ".");
    }

    // ENTRY POINT

    /**
     * Entry point. Reads the input file, constructs a {@link DNAAlign}
     * instance, and calls {@link #align()} to compute result.
     *
     * @param args command-line arguments; {@code args[0]} is the input file
     */
    public static void main(String[] args) {
        if (args.length < 1) {
            System.err.println("Usage: java DNAAlign <inputFile>");
            System.exit(1);
        }

        try (Scanner sc = new Scanner(new BufferedReader(new FileReader(args[0])))) {
            int delta = sc.nextInt();

            int[][] alpha = new int[4][4];
            for (int r = 0; r < 4; r++) {
                for (int c = 0; c < 4; c++) {
                    alpha[r][c] = sc.nextInt();
                }
            }

            String s = sc.next();
            String t = sc.next();

            new DNAAlign(delta, alpha, s, t).align();

        } catch (FileNotFoundException e) {
            System.err.println("Error: file not found — " + args[0]);
            System.exit(1);
        }
    }
}