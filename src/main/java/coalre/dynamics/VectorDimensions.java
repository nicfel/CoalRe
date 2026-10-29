package coalre.dynamics;

import beast.base.core.BEASTInterface;
import beast.base.spec.inference.parameter.RealVectorParam;
import beast.base.spec.type.RealVector;

final class VectorDimensions {

    private VectorDimensions() {
    }

    /**
     * Gives a vector input the dimension its CalcNode needs, typically one entry per rate-shift
     * interval. Only a parameter can be resized; any other vector must already have that size.
     */
    static void ensureDimension(RealVector<?> vector, int dimension) {
        if (vector.size() == dimension)
            return;

        if (vector instanceof RealVectorParam<?> param) {
            param.setDimension(dimension);
        } else {
            String id = vector instanceof BEASTInterface o ? o.getID() : vector.toString();
            throw new IllegalArgumentException(id + " has dimension " + vector.size()
                    + " but expected " + dimension);
        }
    }
}
