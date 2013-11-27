#!/usr/bin/env python
"""
CS176 - Problem Set 5
HMM implementation for finding most recent common ancestor.

Anthony Sutardja
Kevin Tee
"""
from params import Theta
from params import INITIAL_MU_PARAMS, INITIAL_2MU_PARAMS, INITIAL_5MU_PARAMS

from params import INITIAL_4MU_PARAMS

from hmm import Decoding, EM

from util import time_it


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

"""
STATE_TO_TIME, SEQUENCE_TYPE_TO_THETA, SEQUENCE_TYPE_TO_THETA are all constants
that are specific to the tmrca problem.
"""

"""
STATE_TO_TIME is translates states k in {1,2,3,4} to it's respective times.
"""
STATE_TO_TIME = {
    1: 0.32,
    2: 1.75,
    3: 4.54,
    4: 9.40
}

"""
SEQUENCE_TYPE_TO_FILE provides the correct paths to the data files.
"""
SEQUENCE_TYPE_TO_FILE = {
    'mu': 'data/sequences_mu.fasta',
    '2mu': 'data/sequences_2mu.fasta',
    '5mu': 'data/sequences_5mu.fasta',
    '4mu': 'example/sequences_4mu.fasta',
}

"""
SEQUENCE_TYPE_TO_THETA provides the correct variable names for the init params.
"""
SEQUENCE_TYPE_TO_THETA = {
    'mu': INITIAL_MU_PARAMS,
    '2mu': INITIAL_2MU_PARAMS,
    '4mu': INITIAL_4MU_PARAMS,
    '5mu': INITIAL_5MU_PARAMS,
}


def convert_state_seq_to_time_seq(q_sequence):
    """
    Returns a list of times converted from state numbers.
    >>> convert_state_seq_to_time_seq([1, 2, 2])
    [0.32, 1.75, 1.75]
    """
    return [STATE_TO_TIME[q] for q in q_sequence]


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


@time_it
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


class TSequence(object):

    def __init__(self, s):
        """
        Initialize TSequence to the dataset you are interested in investigating
        Can only be called with the following strings `s`:
        'mu', '2mu', or '5mu'

        e.g. TSequence('mu')
        """
        if s not in SEQUENCE_TYPE_TO_THETA:
            raise Error('Invalid mu type')

        self.sequences = read_fasta_sequences_to_str(SEQUENCE_TYPE_TO_FILE[s])
        self.obs = observe_differences(self.sequences[0], self.sequences[1])
        self.theta = SEQUENCE_TYPE_TO_THETA[s]
        self.estimate = None
        self.likelihood = None

    def estimate_params(self, thresh=1e-3, max_iter=15):
        """Performs EM on this dataset and initial parameters."""
        em = EM(self.obs, self.theta, thresh=thresh, max_iter=max_iter)
        self.estimate = em.estimate_params()
        self.likelihood = em.lhood

    def decode_initial(self):
        """Performs marginal posterior decoding and viterbi decoding."""
        decode = Decoding(self.obs, self.theta)
        # need to also return Viterbi
        return decode.posterior(), decode.viterbi(), decode.expectation()

    def decode_estimate(self):
        """Must run estimate_params() on object first."""
        if self.estimate == None:
            raise Error('Must run estimate_params first!')
        decode = Decoding(self.obs, self.estimate)
        # need to also return viterbi
        return decode.posterior(), decode.viterbi(), decode.expectation()

    def theta_to_str(self):
        return ''


if __name__ == '__main__':
    pass
