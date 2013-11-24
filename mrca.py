#!/usr/bin/env python
"""
CS176 - Problem Set 5
HMM implementation for finding most recent common ancestor.

Anthony Sutardja
Kevin Tee
"""
from params import Theta
from params import INITIAL_MU_PARAMS, INITIAL_2MU_PARAMS, INITIAL_5MU_PARAMS
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

STATE_TO_TIME = {
    1: 0.32,
    2: 1.75,
    3: 4.54,
    4: 9.40
}


def convert_state_seq_to_time_seq(q_sequence):
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



SEQUENCE_TYPE_TO_FILE = {
    'mu': 'data/sequences_mu.fasta',
    '2mu': 'data/sequences_2mu.fasta',
    '5mu': 'data/sequences_5mu.fasta',
}


SEQUENCE_TYPE_TO_THETA = {
    'mu': INITIAL_MU_PARAMS,
    '2mu': INITIAL_2MU_PARAMS,
    '5mu': INITIAL_5MU_PARAMS,
}


class TSequence(object):

    def __init__(self, s):
        if s not in SEQUENCE_TYPE_TO_THETA:
            raise Error('Invalid mu type')

        self.sequences = read_fasta_sequences_to_str(SEQUENCE_TYPE_TO_FILE[s])
        self.obs = observe_differences(self.sequences[0], self.sequences[1])
        self.theta = SEQUENCE_TYPE_TO_THETA[s]
        self.estimate = None

    def initial_decoding(self):
        decode = Decoding(self.obs, self.theta)
        # need to also return Viterbi
        return decode.posterior()

    def estimate_params(self):
        em = EM(self.obs, self.theta, max_iter=10)
        self.estimate = em.estimate_params()

    def estimate_decoding(self):
        if self.estimate == None:
            raise Error('Must run estimate_params first!')
        decode = Decoding(self.obs, self.estimate)
        # need to also return viterbi
        return decode.posterior()



if __name__ == '__main__':
    pass
