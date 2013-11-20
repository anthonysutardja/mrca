#!/usr/bin/env python
"""
CS176 - Problem Set 5
HMM implementation for finding most recent common ancestor.

Anthony Sutardja
Kevin Tee
"""
from collections import namedtuple

# Where m: P(Q=k), a: transition prob, e: emission prob
Theta = namedtuple('Theta', ['m', 'a', 'e'])

INITIAL_MU_PARAMS = Theta(
    # Marginal Probabilities
    {
        1: 0.603154,
        2: 0.357419,
        3: 0.0388879,
        4: 0.000539295
    },
    # Transition Probabilities
    {
        1: {1: 0.999916, 2: 0.0000760902, 3: 8.27877e-6, 4: 1.14809e-7},
        2: {1: 0.000128404, 2: 0.999786, 3: 0.0000848889, 4: 1.17723e-6},
        3: {1: 0.000128404, 2: 0.000780214, 3: 0.999068, 4: 0.0000235507},
        4: {1: 0.000128404, 2: 0.000780214, 3: 0.00169821, 4: 0.997393},
    },
    # Emission Probabilities
    {
        1: {'I': 0.999608, 'D': 0.000391695},
        2: {'I': 0.998334, 'D': 0.00166636},
        3: {'I': 0.995844, 'D': 0.00415567},
        4: {'I': 0.991548, 'D': 0.008452},
    }
)


def read_fasta_sequences_to_str(filename):
    """
    Reads a FASTA formatted file to a list of sequences.

    Returns a list of sequences. (should be just 2)
    """
    with open(filename) as f:
        lines = [line.strip() for line in f.readlines()]
        sequences = []
        text = ''

        for line in lines:
            if line[0] == '>':
                if len(text) > 0:
                    sequences.append(text)
                text = ''
            else:
                if len(line):
                    text += line
        if len(text) > 0:
            sequences.append(text)

        return sequences


def observe_differences(seq1, seq2):
    """
    Checks if each index of seq1 and seq2 are the same. If it is the same, then
    we indicate it with a 1, otherwise 0.

    Returns a list of difference observations.
    >>> observe_differences('hello', 'jello')
    ['D', 'I', 'I', 'I', 'I']
    >>> observe_differences('abc', 'abcdefg')
    ['I', 'I', 'I', 'D', 'D', 'D', 'D']
    """
    i, j = 0, 0
    observations = []
    while i < len(seq1) and j < len(seq2):
        if seq1[i] == seq2[j]:
            observations.append('I')
        else:
            observations.append('D')
        i += 1
        j += 1

    # process tails
    while i < len(seq1):
        observations.append('D')
        i += 1
    while j < len(seq2):
        observations.append('D')
        j += 1

    return observations


def forward_algorithm(theta, observations):
    """Performs forward algorithm for computing table values of P( x , q )."""
    forward_table = {}
    return forward_table


def backward_algorithm(theta, observations):
    """Performs backward algorithm for computing table values of P( x | q )."""
    backward_table = {}
    return backward_table


class EM(object):

    def __init__(self, observations, initial_params=None, thresh=0.001):
        self.x = observations
        self.thresh = thresh
        # let self.theta be a dictionary of params
        self.theta = initial_params

    def estimate_params(self):
        # perform the iteration until some value
        while True:
            self.iteration()

        return self.theta

    def iteration(self):
        self.process_forward_algorithm()
        self.process_backward_algorithm()

        # do all the actual crap here...

        pass

    def process_forward_algorithm(self):
        self.forward = forward_algorithm(self.theta, self.x)

    def process_backward_algorithm(self):
        self.backward = backward_algorithm(self.theta, self.x)

    def get_forward_value(self, idx, state):
        return self.forward[idx][state]

    def get_backward_value(self, idx, state):
        return self.backward[idx][state]


class Decoding(object):
    
    def __init__(self, observations, params):
        self.x = observations
        self.theta = params

    def viterbi(self):
        return []

    def marginal(self):
        return []


def main():
    pass


if __name__ == '__main__':
    main()
