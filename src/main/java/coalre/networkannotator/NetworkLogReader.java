package coalre.networkannotator;

import coalre.network.Network;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public interface NetworkLogReader extends Iterable<Network> {

    int getNetworkCount();
    int getCorrectedNetworkCount();
}
