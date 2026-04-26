import java.io.*;
import java.util.*;

/**
 * DNAAlign.java
 *
 * <p>Reads a gap penalty delta (δ), a 4×4 mismatch-cost matrix (α), and two DNA
 * sequences from a plain-text file, then computes and prints the best alignment
 * and its minimum edit distance.
 *
 * <p>The algorithm is the sequence-alignment dynamic programming procedure
 * {@code align(s, t)} from §6 of the course lecture notes. We define
 * {@code opt(i, j)} as the minimum cost of an alignment of the prefixes
 * s₁…sᵢ and t₁…tⱼ, and store all values in a memoization table M
 * (so {@code M[i][j] = opt(i, j)}). After filling M bottom-up, we traceback
 * from {@code M[m][n]} to reconstruct the actual aligned strings and
 * per-column penalties.
 *
 * <p>Usage: {@code java DNAAlign <inputFile>}
 *
 * <p>APIs used: {@code java.io.BufferedReader}, {@code java.io.FileReader},
 * {@code java.util.Scanner}, {@code java.util.ArrayList},
 * {@code java.util.Collections}.
 *
 * @author [Your Name]
 * @author [Partner Name]
 * @version 1.0
 */
public class DNAAlign {

    // -----------------------------------------------------------------------
    // Constants
    // -----------------------------------------------------------------------

    /**
     * Canonical ordering of nucleotides used to index the mismatch-cost
     * matrix α: A=0, C=1, G=2, T=3.
     */
    private static final String NUCLEOTIDES = "ACGT";

    // -----------------------------------------------------------------------
    // Fields
    // -----------------------------------------------------------------------

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

    /** The first DNA sequence s = s₁s₂…sₘ. */
    private final String s;

    /** The second DNA sequence t = t₁t₂…tₙ. */
    private final String t;

    // -----------------------------------------------------------------------
    // Constructor
    // -----------------------------------------------------------------------

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

    // -----------------------------------------------------------------------
    // Core algorithm — align(s, t) from §6 of the lecture notes
    // -----------------------------------------------------------------------

    /**
     * Computes and prints the best alignment of {@code s} and {@code t}
     * using the sequence-alignment procedure {@code align(s, t)} from §6
     * of the course lecture notes.
     *
     * <p><b>Subproblem definition.</b>
     * For i = 0, 1, …, m and j = 0, 1, …, n, we define opt(i, j) to be
     * the minimum cost of an alignment of the prefixes s₁…sᵢ and t₁…tⱼ.
     * We store these values in a memoization table M: {@code M[i][j] = opt(i, j)}.
     *
     * <p><b>Base cases (from the lecture notes).</b>
     * <ul>
     *   <li>opt(i, 0) = i · δ  — aligning s₁…sᵢ with the empty string
     *       requires i gaps, each costing δ.</li>
     *   <li>opt(0, j) = j · δ  — symmetrically for t₁…tⱼ against empty s.</li>
     * </ul>
     *
     * <p><b>Recurrence (from the lecture notes).</b>
     * For i ≥ 1, j ≥ 1, consider any optimal alignment O of s₁…sᵢ and
     * t₁…tⱼ. Its last column must be one of three cases:
     * <ol>
     *   <li>sᵢ is paired with tⱼ: O contains an optimal alignment of
     *       s₁…s_{i−1} and t₁…t_{j−1}, so the cost is
     *       α_{sᵢtⱼ} + opt(i−1, j−1).</li>
     *   <li>sᵢ is unmatched (gap in t): O contains an optimal alignment
     *       of s₁…s_{i−1} and t₁…tⱼ, so the cost is δ + opt(i−1, j).</li>
     *   <li>tⱼ is unmatched (gap in s): O contains an optimal alignment
     *       of s₁…sᵢ and t₁…t_{j−1}, so the cost is δ + opt(i, j−1).</li>
     * </ol>
     * Therefore:
     * <pre>
     *   opt(i, j) = min( α_{sᵢtⱼ} + opt(i−1, j−1),
     *                    δ          + opt(i−1, j),
     *                    δ          + opt(i, j−1) )
     * </pre>
     *
     * <p><b>Traceback.</b>
     * Starting at M[m][n] and working back to M[0][0], at each cell we
     * determine which of the three cases produced the stored opt value,
     * and record the corresponding alignment column. The lists are built in
     * reverse and then reversed at the end.
     *
     * <p><b>Running time.</b>
     * The table M has (m+1)(n+1) entries; each is computed in O(1) time
     * by the recurrence. Filling M therefore takes O(mn) time. The
     * traceback from M[m][n] to M[0][0] visits at most m+n cells and
     * takes O(m+n) time. The overall running time is O(mn).
     */
    public void align() {
        int m = s.length();
        int n = t.length();

        // ------------------------------------------------------------------
        // Step 1: Build the memoization table M where M[i][j] = opt(i, j).
        // This is the align(s, t) procedure from §6 of the lecture notes.
        // ------------------------------------------------------------------
        int[][] M = new int[m + 1][n + 1];

        // Base cases: opt(i, 0) = i·δ  and  opt(0, j) = j·δ
        for (int i = 0; i <= m; i++) M[i][0] = i * delta;
        for (int j = 0; j <= n; j++) M[0][j] = j * delta;

        // Fill bottom-up using the recurrence from the lecture notes:
        //   M[i][j] = min( α_{sᵢtⱼ} + M[i-1][j-1],
        //                  δ          + M[i-1][j],
        //                  δ          + M[i][j-1]  )
        for (int i = 1; i <= m; i++) {
            for (int j = 1; j <= n; j++) {
                int pairCost = alphaCost(s.charAt(i - 1), t.charAt(j - 1)) + M[i - 1][j - 1];
                int gapInT   = delta + M[i - 1][j];   // sᵢ unmatched
                int gapInS   = delta + M[i][j - 1];   // tⱼ unmatched
                M[i][j] = Math.min(pairCost, Math.min(gapInT, gapInS));
            }
        }

        // ------------------------------------------------------------------
        // Step 2: Traceback from M[m][n] to reconstruct the alignment.
        // ------------------------------------------------------------------
        List<Character> alignS    = new ArrayList<>();
        List<Character> alignT    = new ArrayList<>();
        List<Integer>   penalties = new ArrayList<>();

        int i = m, j = n;
        while (i > 0 || j > 0) {
            if (i > 0 && j > 0) {
                int pen = alphaCost(s.charAt(i - 1), t.charAt(j - 1));
                if (M[i][j] == pen + M[i - 1][j - 1]) {
                    // Case 1: sᵢ paired with tⱼ (match or mismatch)
                    alignS.add(s.charAt(i - 1));
                    alignT.add(t.charAt(j - 1));
                    penalties.add(pen);
                    i--; j--;
                } else if (M[i][j] == delta + M[i - 1][j]) {
                    // Case 2: sᵢ unmatched — gap inserted in t
                    alignS.add(s.charAt(i - 1));
                    alignT.add('-');
                    penalties.add(delta);
                    i--;
                } else {
                    // Case 3: tⱼ unmatched — gap inserted in s
                    alignS.add('-');
                    alignT.add(t.charAt(j - 1));
                    penalties.add(delta);
                    j--;
                }
            } else if (i > 0) {
                // Consumed all of t; remaining sᵢ characters become gaps
                alignS.add(s.charAt(i - 1));
                alignT.add('-');
                penalties.add(delta);
                i--;
            } else {
                // Consumed all of s; remaining tⱼ characters become gaps
                alignS.add('-');
                alignT.add(t.charAt(j - 1));
                penalties.add(delta);
                j--;
            }
        }

        // Traceback produces the alignment in reverse; flip all three lists.
        Collections.reverse(alignS);
        Collections.reverse(alignT);
        Collections.reverse(penalties);

        // ------------------------------------------------------------------
        // Step 3: Print the result.
        // ------------------------------------------------------------------
        printAlignment(alignS, alignT, penalties, M[m][n]);
    }

    // -----------------------------------------------------------------------
    // Helper methods
    // -----------------------------------------------------------------------

    /**
     * Returns the mismatch cost α_{ab} for aligning nucleotide {@code a}
     * with nucleotide {@code b}, as given by the matrix {@link #alpha}.
     * When {@code a == b} this returns 0 (a perfect match); otherwise it
     * returns the positive penalty from the similarity matrix.
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
     * A=0, C=1, G=2, T=3, which is used to look up entries in the
     * mismatch-cost matrix α.
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
     * minimum edit distance to standard output in the required format.
     *
     * @param alignS    characters from s in the alignment (gaps shown as '-')
     * @param alignT    characters from t in the alignment (gaps shown as '-')
     * @param penalties per-column cost values corresponding to the alignment
     * @param minCost   the minimum edit distance opt(m, n) = M[m][n]
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

    // -----------------------------------------------------------------------
    // Entry point
    // -----------------------------------------------------------------------

    /**
     * Entry point. Reads the input file, constructs a {@link DNAAlign}
     * instance, and calls {@link #align()} to compute and print the result.
     *
     * @param args command-line arguments; {@code args[0]} must be the path
     *             to the input file
     */
    public static void main(String[] args) {
        if (args.length < 1) {
            System.err.println("Usage: java DNAAlign <inputFile>");
            System.exit(1);
        }

        try (Scanner sc = new Scanner(new BufferedReader(new FileReader(args[0])))) {

            // Read the gap penalty δ
            int delta = sc.nextInt();

            // Read the 4×4 mismatch-cost matrix α (rows in A, C, G, T order)
            int[][] alpha = new int[4][4];
            for (int r = 0; r < 4; r++) {
                for (int c = 0; c < 4; c++) {
                    alpha[r][c] = sc.nextInt();
                }
            }

            // Read the two DNA sequences s and t
            String s = sc.next();
            String t = sc.next();

            // Run the alignment
            new DNAAlign(delta, alpha, s, t).align();

        } catch (FileNotFoundException e) {
            System.err.println("Error: file not found — " + args[0]);
            System.exit(1);
        }
    }
}