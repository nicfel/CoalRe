package coalre.util;

import java.util.BitSet;
import java.util.stream.IntStream;

public class BitSetFun {

    public static BitSet or(BitSet a, BitSet b) {
        BitSet result = (BitSet) a.clone();
        result.or(b);
        return result;
    }

    public static BitSet and(BitSet a, BitSet b) {
        BitSet result = (BitSet) a.clone();
        result.and(b);
        return result;
    }

    public static BitSet xor(BitSet a, BitSet b) {
        BitSet result = (BitSet) a.clone();
        result.xor(b);
        return result;
    }

    public static BitSet not(BitSet a) {
        BitSet result = (BitSet) a.clone();
        result.flip(0, a.length());
        return result;
    }

    public static BitSet full(int n) {
        BitSet result = new BitSet(n);
        result.set(0, n);
        return result;
    }
    
    public static BitSet empty(int n) {
        return new BitSet(n);
    }

    public static BitSet copy(BitSet a) {
        return (BitSet) a.clone();
    }

    public static BitSet andNot(BitSet a, BitSet b) {
        BitSet result = (BitSet) a.clone();
        result.andNot(b);
        return result;
    }
    
    public static boolean disjoint(BitSet a, BitSet b) {
        return and(a, b).isEmpty();
    }

    public static BitSet singleton(int i) {
        BitSet result = new BitSet();
        result.set(i);
        return result;
    }

    public static Iterable<Integer> iterate(BitSet bitset) {
        return bitset.stream()::iterator;
    }

    public static String toBitString(BitSet bitset, int nbits) {
        final StringBuilder buffer = new StringBuilder(nbits);
        IntStream.range(0, nbits)
                .mapToObj(i -> bitset.get(i) ? '1' : '0')
                .forEach(buffer::append);
        return buffer.toString();
    }


}
