#!/usr/bin/env python
"""
CS176 - Problem Set 5
HMM implementation for finding most recent common ancestor.

Anthony Sutardja
Kevin Tee
"""
from collections import namedtuple
from math import e
from math import log

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
    """Performs forward algorithm for computing table values of P( x , q )."""
    forward_table = {}
    """
    TODO(kevintee):

    See example on line 39 to see how to access `theta`.

    forward_table should be a dict with the keys being the state of Q, and the
    value being a list of forward values with length up to the observation seq.

    i.e. f_t(k) can be accessed by forward_table[k][t]
    """
    states = theta.m.keys()
    for t in len(range(observations)):
        for k in states:
            if t == 0:
                forward_table[k] = [log(theta.e[k][observations[t]] + \
                    log(theta.m[k]))]
            else:
                offset = max([forward_table[i][t - 1] for i in states])
                forward_table[k].append(log(theta.e[k][observations[t]]) + \
                    offset + log(sum([e**(forward_table[i][t - 1] - offset) * \
                    theta.a[i][k] for i in states])))

    return forward_table


def backward_algorithm(theta, observations):
    """Performs backward algorithm for computing table values of P( x | q )."""
    backward_table = {}
    """
    TODO(kevintee):
    See forward_algorithm comment. Same biznis.
    i.e. b_t(k) can be accessed by backward_table[k][t]

    Would recommend initializing the entire length of the array first since
    you're working backwards.
    """
    states = theta.m.keys()
    for k in states:
        backward_table = [0] * len(range(observations))

    for t in reversed(len(range(observations))):
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

    def __init__(self, observ, init_params=None, thresh=1e-5, max_iter=100):
        self.x = observ
        self.thresh = thresh
        # let self.theta be a dictionary of params
        self.theta = init_params
        self.max_iter = max_iter

    def estimate_params(self):
        # perform the iteration until some value
        i = 0

        # perform first iteration manually
        old_theta = self.theta
        self.theta = self.iteration()
        i += 1

        # iterate until improvement meets threshold
        while check_improvement(self.theta, old_theta) and i < self.max_iter:
            old_theta = self.theta
            self.theta = self.iteration()
            i += 1

        return self.theta

    def check_improvement(self, new_theta, old_theta):
        return (self.calculate_likelihood(new_theta, self.x) \
                - self.calculate_likelihood(old_theta, self.x)) > self.thresh

    def iteration(self):
        """Generate a new theta."""
        self.process_forward_algorithm()
        self.process_backward_algorithm()

        theta_prime = {}
        # do all the actual crap here...

        return theta_prime

    def calculate_likelihood(self, theta):
        """Return the log-likelihood of the data given theta."""
        forward = forward_algorithm(theta, self.x)

        # do something with logs....
        return 0.0

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
