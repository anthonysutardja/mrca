#!/usr/bin/env python
"""
CS176 - Problem Set 5
HMM implementation for finding most recent common ancestor.

Anthony Sutardja
Kevin Tee
"""


def read_fasta_sequence_to_str(filename):
    """
    Reads a FASTA formatted file to a string.

    Returns a tuple containing the sequence's name and the sequence.
    """
    with open(filename) as f:
        lines = [line.strip() for line in f.readlines()]
        sequence_name = lines.pop(0).replace('>', '')
        return (sequence_name, ''.join(lines))


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
        self.thresh = 0.001
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
