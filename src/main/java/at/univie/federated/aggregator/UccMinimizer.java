package at.univie.federated.aggregator;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collection;
import java.util.Comparator;
import java.util.List;

/**
 * Handles minimization of UCC sets.
 * Simple responsibility: keep only subset-minimal BitSets.
 */
public class UccMinimizer {
    
    /**
     * Keeps only subset-minimal BitSets.
     * Sorts by cardinality, then keeps sets that are not dominated by any kept set.
     */
    public List<BitSet> minimalize(Collection<BitSet> sets) {
        List<BitSet> sorted = new ArrayList<>(sets.size());
        for (BitSet s : sets) {
            sorted.add((BitSet) s.clone());
        }
        
        sorted.sort(Comparator
                .comparingInt(BitSet::cardinality)
                .thenComparing(this::lexCompareBitSet));
        
        List<BitSet> kept = new ArrayList<>();
        for (BitSet candidate : sorted) {
            if (!isDominated(candidate, kept)) {
                kept.add(candidate);
            }
        }
        
        return kept;
    }
    
    private boolean isDominated(BitSet candidate, List<BitSet> kept) {
        for (BitSet k : kept) {
            if (isSubset(k, candidate)) {
                return true;
            }
        }
        return false;
    }
    
    private boolean isSubset(BitSet a, BitSet b) {
        // a ⊆ b  <=>  (a \ b) is empty
        BitSet tmp = (BitSet) a.clone();
        tmp.andNot(b);
        return tmp.isEmpty();
    }
    
    /**
     * Lexicographic comparison for deterministic ordering.
     */
    private int lexCompareBitSet(BitSet a, BitSet b) {
        int i = a.nextSetBit(0);
        int j = b.nextSetBit(0);
        
        while (i >= 0 || j >= 0) {
            if (i != j) {
                return Integer.compare(i, j);
            }
            i = a.nextSetBit(i + 1);
            j = b.nextSetBit(j + 1);
        }
        return 0;
    }
    
    /**
     * Converts BitSet to comma-separated string.
     */
    public static String toCommaList(BitSet bs) {
        StringBuilder sb = new StringBuilder();
        for (int i = bs.nextSetBit(0); i >= 0; i = bs.nextSetBit(i + 1)) {
            if (sb.length() > 0) {
                sb.append(',');
            }
            sb.append(i);
        }
        return sb.toString();
    }
}

