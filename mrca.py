#!/usr/bin/env python
"""
CS176 - Problem Set 5
HMM implementation for finding most recent common ancestor.

Anthony Sutardja
Kevin Tee
"""
from params import Theta
from params import INITIAL_MU_PARAMS, INITIAL_2MU_PARAMS, INITIAL_5MU_PARAMS

from math import e
from math import log

"""
Access to the parameter of Theta (marginal, transition, emission) can be
done with the following:
>>> theta = INITIAL_MU_PARAMS
>>> theta.m[1]  # marginal probability of state k=1, Pr(Q_i = k)
0.603154
>>> theta.a[2][1]  # transition probability from 2 to 1, Pr(Q_j = 1 | Q_i = 2)
0.000128404
>>> theta.e[3]['D']  # emission probability of different char while in state 3
0.00415567
"""

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
    """
    Performs forward alg for computing table values of log P( x , q ).

    i.e. f_t(k) can be accessed by forward_table[k][t]
    """
    forward_table = {}

    states = theta.m.keys()
    for t in range(len(observations)):
        for k in states:
            if t == 0:
                forward_table[k] = [log(theta.e[k][observations[t]]) + \
                    log(theta.m[k])]
            else:
                offset = max([forward_table[i][t - 1] for i in states])
                forward_table[k].append(log(theta.e[k][observations[t]]) + \
                    offset + log(sum([e**(forward_table[i][t - 1] - offset) * \
                    theta.a[i][k] for i in states])))

    return forward_table


def backward_algorithm(theta, observations):
    """Performs backward alg for computing table values of log P( x | q )."""
    backward_table = {}
    
    states = theta.m.keys()
    for k in states:
        backward_table[k] = [0] * len(observations)

    for t in reversed(range(len(observations))):
        for k in states:
            if t == len(observations) - 1:
                backward_table[k][t] = 0
            else:
                offset = max([backward_table[i][t + 1] for i in states])
                backward_table[k][t] = offset + log(sum([theta.a[k][j] * \
                    theta.e[j][observations[t + 1]] * \
                    e**(backward_table[j][t + 1] - offset) for j in states]))

    return backward_table


class EM(object):

    def __init__(self, observ, init_params, thresh=1e-5, max_iter=100):
        self.x = observ
        self.thresh = thresh
        # let self.theta be a dictionary of params
        self.theta = init_params
        self.max_iter = max_iter

    def estimate_params(self):
        # perform the iteration until some value
        i = 0

        # perform first iteration manually
        old_lhood = float('-inf')
        self.process_forward_algorithm()
        self.process_backward_algorithm()
        self.lhood = self.calculate_likelihood()

        # iterate until improvement meets threshold
        while self.check_improvement(self.lhood, old_lhood) and i < self.max_iter:
            print self.lhood
            old_lhood = self.lhood

            # self.iteration() will always have access to the correct
            # forward and backward tables based on the current self.theta
            self.theta = self.iteration()

            # Update forward and backward tables
            self.process_forward_algorithm()
            self.process_backward_algorithm()
            self.lhood = self.calculate_likelihood()
            i += 1

        return self.theta

    def check_improvement(self, new_lhood, old_lhood):
        return new_lhood - old_lhood > self.thresh

    def iteration(self):
        """Generate a new theta."""
        states = self.theta.m.keys()
        
        # do all the actual crap here...
        marginal = {}
        for k in states:
            marginal[k] = self.calculate_marginal_for_state(k)
        marginal = self.theta.m

        transition = {}
        for i in states:
            for j in states:
                pass
        transition = self.theta.a

        emission = {}
        symbols = ['I', 'D']
        for k in states:
            for symbol in symbols:
                pass
        # need to normalize
        emission = self.theta.e

        return Theta(marginal, transition, emission)

    def calculate_marginal_for_state(self, k):
        return e**(self.f(k)[0] + self.b(k)[0] - self.lhood)

    def calculate_likelihood(self, theta=None):
        """Return the log-likelihood of the data given theta."""
        if theta:
            forward = forward_algorithm(theta, self.x)
        else:
            theta = self.theta
            forward = self.forward

        states = theta.m.keys()
        D = max([forward[k][-1] for k in states])

        return log(sum([e**(forward[k][-1] - D) for k in states])) + D

    def process_forward_algorithm(self):
        self.forward = forward_algorithm(self.theta, self.x)

    def process_backward_algorithm(self):
        self.backward = backward_algorithm(self.theta, self.x)

    def f(self, state):
        return self.forward[state]

    def b(self, state):
        return self.backward[state]


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
