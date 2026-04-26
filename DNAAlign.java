import java.io.*;
import java.util.*;

/**
 * DNAAlign.java
 *
 * <p>Reads a gap penalty δ, a 4×4 similarity matrix (indexed A,C,G,T), and two
 * DNA sequences from a text file, then computes and prints the best (minimum
 * edit-distance) alignment using a standard O(mn) sequence-alignment dynamic
 * programming algorithm (Needleman–Wunsch).
 *
 * <p>Usage:  java DNAAlign &lt;inputFile&gt;
 *
 * <p>APIs used: java.io (BufferedReader, FileReader), java.util (Scanner, ArrayList,
 * Collections).
 *
 * @author  [Your Name]
 * @author  [Partner Name]
 * @version 1.0
 */
public class DNAAlign {

    // -----------------------------------------------------------------------
    // Constants
    // -----------------------------------------------------------------------

    /** The nucleotide order used to index the similarity matrix. */
    private static final String NUCLEOTIDES = "ACGT";

    // -----------------------------------------------------------------------
    // Fields
    // -----------------------------------------------------------------------

    /** Gap penalty δ. */
    private final int delta;

    /**
     * 4×4 similarity (mismatch-penalty) matrix.
     * Indexed by {@link #indexOf(char)}.
     */
    private final int[][] simMatrix;

    /** First DNA sequence. */
    private final String s;

    /** Second DNA sequence. */
    private final String t;

    // -----------------------------------------------------------------------
    // Constructor
    // -----------------------------------------------------------------------

    /**
     * Constructs a DNAAlign instance with all parameters needed for alignment.
     *
     * @param delta     the gap penalty
     * @param simMatrix the 4×4 mismatch-penalty matrix (A,C,G,T order)
     * @param s         the first DNA sequence
     * @param t         the second DNA sequence
     */
    public DNAAlign(int delta, int[][] simMatrix, String s, String t) {
        this.delta     = delta;
        this.simMatrix = simMatrix;
        this.s         = s.toUpperCase();
        this.t         = t.toUpperCase();
    }

    // -----------------------------------------------------------------------
    // Core algorithm
    // -----------------------------------------------------------------------

    /**
     * Computes the best alignment of {@code s} and {@code t} using the
     * Needleman–Wunsch O(mn) dynamic programming algorithm, then prints the
     * result to standard output.
     *
     * <p><b>Algorithm outline:</b>
     * <ol>
     *   <li>Build an (m+1) × (n+1) DP table where {@code dp[i][j]} is the
     *       minimum edit distance when aligning {@code s[0..i-1]} with
     *       {@code t[0..j-1]}.</li>
     *   <li>Base cases: {@code dp[i][0] = i*delta}, {@code dp[0][j] = j*delta}
     *       (aligning a prefix against an empty string requires all gaps).</li>
     *   <li>Recurrence:
     *       {@code dp[i][j] = min(dp[i-1][j-1] + mismatchPenalty(s[i], t[j]),
     *                             dp[i-1][j]   + delta,
     *                             dp[i][j-1]   + delta)}</li>
     *   <li>Traceback from {@code dp[m][n]} to reconstruct the aligned strings
     *       and per-column penalties.</li>
     * </ol>
     *
     * <p><b>Correctness:</b> The recurrence covers every possible case for
     * a column in the alignment: match/mismatch (diagonal move), gap in t
     * (up move), or gap in s (left move).  Because every alignment corresponds
     * to a unique path in the DP table, the globally optimal value is found
     * at {@code dp[m][n]}.
     *
     * <p><b>Running time:</b> Filling the (m+1)×(n+1) table is O(mn);
     * traceback is O(m+n).  Total: O(mn).
     */
    public void align() {
        int m = s.length();
        int n = t.length();

        // ------------------------------------------------------------------
        // 1. Build DP table
        // ------------------------------------------------------------------
        int[][] dp = new int[m + 1][n + 1];

        // Base cases
        for (int i = 0; i <= m; i++) dp[i][0] = i * delta;
        for (int j = 0; j <= n; j++) dp[0][j] = j * delta;

        // Fill table
        for (int i = 1; i <= m; i++) {
            for (int j = 1; j <= n; j++) {
                int matchCost = dp[i - 1][j - 1] + mismatchPenalty(s.charAt(i - 1), t.charAt(j - 1));
                int gapInT   = dp[i - 1][j]     + delta;  // gap inserted in t (under s[i])
                int gapInS   = dp[i][j - 1]     + delta;  // gap inserted in s (under t[j])
                dp[i][j] = Math.min(matchCost, Math.min(gapInT, gapInS));
            }
        }

        // ------------------------------------------------------------------
        // 2. Traceback to reconstruct alignment
        // ------------------------------------------------------------------
        List<Character> alignS    = new ArrayList<>();
        List<Character> alignT    = new ArrayList<>();
        List<Integer>   penalties = new ArrayList<>();

        int i = m, j = n;
        while (i > 0 || j > 0) {
            if (i > 0 && j > 0) {
                int pen = mismatchPenalty(s.charAt(i - 1), t.charAt(j - 1));
                if (dp[i][j] == dp[i - 1][j - 1] + pen) {
                    // Diagonal: match or mismatch
                    alignS.add(s.charAt(i - 1));
                    alignT.add(t.charAt(j - 1));
                    penalties.add(pen);
                    i--; j--;
                } else if (dp[i][j] == dp[i - 1][j] + delta) {
                    // Up: gap in t
                    alignS.add(s.charAt(i - 1));
                    alignT.add('-');
                    penalties.add(delta);
                    i--;
                } else {
                    // Left: gap in s
                    alignS.add('-');
                    alignT.add(t.charAt(j - 1));
                    penalties.add(delta);
                    j--;
                }
            } else if (i > 0) {
                // Remaining characters in s → gaps in t
                alignS.add(s.charAt(i - 1));
                alignT.add('-');
                penalties.add(delta);
                i--;
            } else {
                // Remaining characters in t → gaps in s
                alignS.add('-');
                alignT.add(t.charAt(j - 1));
                penalties.add(delta);
                j--;
            }
        }

        // Traceback produces the alignment in reverse; reverse all three lists.
        Collections.reverse(alignS);
        Collections.reverse(alignT);
        Collections.reverse(penalties);

        // ------------------------------------------------------------------
        // 3. Print result
        // ------------------------------------------------------------------
        printAlignment(alignS, alignT, penalties, dp[m][n]);
    }

    // -----------------------------------------------------------------------
    // Helper methods
    // -----------------------------------------------------------------------

    /**
     * Returns the mismatch penalty for aligning nucleotide {@code a} with
     * nucleotide {@code b} as given by the similarity matrix.
     *
     * @param a a nucleotide character (A, C, G, or T)
     * @param b a nucleotide character (A, C, G, or T)
     * @return  the penalty for aligning {@code a} with {@code b}
     */
    private int mismatchPenalty(char a, char b) {
        return simMatrix[indexOf(a)][indexOf(b)];
    }

    /**
     * Returns the index of a nucleotide in the canonical order A=0, C=1,
     * G=2, T=3.
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
     * @param alignS    aligned characters from the first sequence (with gaps)
     * @param alignT    aligned characters from the second sequence (with gaps)
     * @param penalties per-column penalty values
     * @param minDist   the minimum edit distance
     */
    private void printAlignment(List<Character> alignS,
                                List<Character> alignT,
                                List<Integer>   penalties,
                                int             minDist) {
        System.out.println("The best alignment is\n");

        // Build display strings (space-separated)
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
        System.out.println("\nwith the minimum edit distance of " + minDist + ".");
    }

    // -----------------------------------------------------------------------
    // Entry point
    // -----------------------------------------------------------------------

    /**
     * Entry point.  Reads the input file, constructs a {@link DNAAlign}
     * instance, and calls {@link #align()}.
     *
     * <p>Expected input-file format:
     * <pre>
     *   &lt;delta&gt;
     *   &lt;row for A: four integers&gt;
     *   &lt;row for C: four integers&gt;
     *   &lt;row for G: four integers&gt;
     *   &lt;row for T: four integers&gt;
     *   &lt;sequence s&gt;
     *   &lt;sequence t&gt;
     * </pre>
     *
     * @param args command-line arguments; {@code args[0]} must be the path to
     *             the input file
     */
    public static void main(String[] args) {
        if (args.length < 1) {
            System.err.println("Usage: java DNAAlign <inputFile>");
            System.exit(1);
        }

        try (Scanner sc = new Scanner(new BufferedReader(new FileReader(args[0])))) {

            // Read gap penalty
            int delta = sc.nextInt();

            // Read 4×4 similarity matrix
            int[][] simMatrix = new int[4][4];
            for (int r = 0; r < 4; r++) {
                for (int c = 0; c < 4; c++) {
                    simMatrix[r][c] = sc.nextInt();
                }
            }

            // Read the two sequences
            String s = sc.next();
            String t = sc.next();

            // Run alignment
            new DNAAlign(delta, simMatrix, s, t).align();

        } catch (FileNotFoundException e) {
            System.err.println("Error: file not found – " + args[0]);
            System.exit(1);
        }
    }
}