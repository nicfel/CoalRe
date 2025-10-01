package coalre.util;

import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class Util {

    public static <T> boolean allElementsUnique(List<T> list) {
        Set<T> set = new HashSet<>(list);
        return set.size() == list.size();
    }

    public static Iterable<Integer> range(int endExclusive) {
        return IntStream.range(0, endExclusive)::iterator;
    }

    public static Iterable<Integer> range(int startInclusive, int endExclusive) {
        return IntStream.range(startInclusive, endExclusive)::iterator;
    }

    public static Double sum(Iterable<Double> xs) {
        Double s = 0.0;
        for (Double x : xs) {
            s += x;
        }
        return s;
    }

    public static List<Double> normalize(List<Double> xs) {
        double sum = sum(xs);
        return xs.stream()
            .map(x -> x / sum)
            .collect(Collectors.toList());
    }

    public static <T> void writeToColumn(T[][] matrix, int col, List<T> values) {
        for (int i : range(values.size())) {
            matrix[i][col] = values.get(i);
        }
    }

    public static <T> void writeToColumn(T[][] matrix, int col, T[] values) {
        for (int i : range(values.length)) {
            matrix[i][col] = values[i];
        }
    }
}
