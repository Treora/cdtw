"""This file shows a simple usage example of (constrained) Dynamic Time Warping. Given a query sequence, the goal is to find the sequence in the data set that provides the closest match."""

import sys
from operator import itemgetter

from numpy import *
try:
    from matplotlib import pyplot as plt
except ImportError:
    enable_plots = False
else:
    enable_plots = True

from cdtw import cdtw_sakoe_chiba


def main(argv):
    # Generate some arbitrary sequences
    t = linspace(0, 6*pi, 40)
    noise = lambda: random.normal(0,0.1,len(t))
    dataset = array([
        cos(t) + noise(),
        abs(cos(t)) + noise(),
        sign(cos(t+1.5)),
        noise()
    ])

    # The query sequence that will be matched to the dataset's sequences
    query_sequence = cos(t+2) + noise()

    # Determine value of r
    if len(argv) >= 2:
        r_param = argv[1]
    else:
        r_param = "10%"
        print("Using default: r={0}".format(r_param))
    sequence_length = len(query_sequence)
    if r_param.endswith('%'):
        r = int(float(r_param[:-1])/100 * sequence_length)
        print("Sequence length is {0}, so r={1}".format(sequence_length, r))
    else:
        r = int(r_param)
    params = {'r': r}
    
    # Find the closest match
    best_index, distances = nearest_neighbour(query_sequence, dataset, params)
    
    set_printoptions(precision=2)
    print("\nData set:")
    for seq_id, sequence in enumerate(dataset):
        print("<{0}>: {1}".format(seq_id, sequence))
    print("\nQuery senquence:")
    print("{0}".format(query_sequence))
    
    print("\nDistances:")
    for seq_id, sequence in enumerate(dataset):
        print("<{0}>: {1}".format(seq_id, distances[seq_id]))
    
    print("Best match was sequence <{0}>".format(best_index))

    if enable_plots:
        for seq_id, sequence in enumerate(dataset):
            plt.subplot(len(dataset), 1, seq_id+1)
            plt.ylim(-1.5, 1.5)
            plt.plot(sequence)
        plt.subplot(len(dataset), 1, best_index+1)
        plt.plot(query_sequence)
        plt.show()


def nearest_neighbour(query_sequence, training_set, params):
    distances = [cdtw_sakoe_chiba(query_sequence, training_sequence, params['r']) for training_sequence in training_set]
    best_seq_index = argmin(distances)
    return best_seq_index, distances

def argmin(iterator):
    return min(enumerate(iterator), key=itemgetter(1))[0]


if __name__ == "__main__":
    main(sys.argv)


def test_cdtw():
    """Simple test to run cDTW with different values of r"""
    print("\n\n== Testing cDTW ==\n")
    a = array([1,2,2,2,3,4,5,6,7], dtype="float64")
    b = array([1,2,3,4,5,6,6,6,7], dtype="float64")
    
    results = [(r, cdtw_sakoe_chiba(a, b, r)) for r in [0,1,2]]
    print("a: {0}".format(a))
    print("b: {0}".format(b))
    print("cDTW results: (should be 8, 4, and 0)")
    for result in results:
        print("r={0[0]}: {0[1]}".format(result))
    print("\n")
