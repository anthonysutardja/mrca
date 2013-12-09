#!/usr/bin/env python
"""
CS176 - Problem Set 5
HMM implementation for finding most recent common ancestor.

Anthony Sutardja
Kevin Tee
"""
from argparse import ArgumentParser
from params import Theta
from params import INITIAL_MU_PARAMS, INITIAL_2MU_PARAMS, INITIAL_5MU_PARAMS

from params import INITIAL_4MU_PARAMS
from params import parse_params

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
    'mu': 'data/initial_parameters_mu.txt',
    '2mu': 'data/initial_parameters_2mu.txt',
    '4mu': 'example/initial_parameters_4mu.txt',
    '5mu': 'data/initial_parameters_5mu.txt',
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

    def __init__(self, initial_param_file, fasta_file):
        """
        Initialize TSequence to the dataset you are interested in investigating
        Can only be called with the following strings `s`:
        'mu', '2mu', or '5mu'

        e.g. TSequence('mu')
        """
        self.sequences = read_fasta_sequences_to_str(fasta_file)
        self.obs = observe_differences(self.sequences[0], self.sequences[1])
        self.theta = parse_params(initial_param_file)
        self.estimate = None
        self.likelihood = None
        self.initial_likelihood = None

    def estimate_params(self, thresh=1e-3, max_iter=15):
        """Performs EM on this dataset and initial parameters."""
        em = EM(self.obs, self.theta, thresh=thresh, max_iter=max_iter)
        self.estimate = em.estimate_params()
        self.likelihood = em.lhood
        self.initial_likelihood = em.calculate_likelihood(theta=self.theta)

    def decode_initial(self):
        """Performs marginal posterior decoding and viterbi decoding."""
        decode = Decoding(self.obs, self.theta)
        # need to also return Viterbi
        return (
            decode.posterior(),
            decode.viterbi(),
            decode.expectation(STATE_TO_TIME)
        )

    def decode_estimate(self):
        """Must run estimate_params() on object first."""
        if self.estimate == None:
            raise Error('Must run estimate_params first!')
        decode = Decoding(self.obs, self.estimate)
        # need to also return viterbi
        return (
            decode.posterior(),
            decode.viterbi(),
            decode.expectation(STATE_TO_TIME)
        )

    def theta_to_str(self):
        if self.estimate == None:
            raise Error('Must run estimate_params first!')

        output = []
        output.extend([
            'Marginal Probabilities\n',
            '========================================\n',
        ])

        for k in range(1, 5):
            output.append('{:e}'.format(self.estimate.m[k]))
            output.append('\n')

        output.append('\n')

        output.extend([
            'Transition Probabilities\n',
            '========================================\n',
        ])

        for i in range(1, 5):
            probs = []
            for j in range(1, 5):
                probs.append('{:e}'.format(self.estimate.a[i][j]))
            output.append(' '.join(probs) + '\n')

        output.append('\n')

        output.extend([
            'Emission Probabilities\n',
            '========================================\n',
        ])
        for k in range(1,5):
            probs = []
            for symbol in ['I', 'D']:
                probs.append('{:e}'.format(self.estimate.e[k][symbol]))
            output.append(' '.join(probs) + '\n')

        return ''.join(output) + '\n'

    def decoding_to_str(self, flag):
        if flag not in ('initial', 'estimate'):
            raise Error('Must call on initial or estimate')

        if flag == 'initial':
            posterior, viterbi, expectation = self.decode_initial()
        else:
            posterior, viterbi, expectation = self.decode_estimate()

        lines = []
        lines.append(
            '# Viterbi_decoding posterior_decoding posterior_mean'
        )
        for i in range(len(posterior)):
            lines.append(
                '%.2f %.2f %.6f' % (
                    STATE_TO_TIME[viterbi[i]],
                    STATE_TO_TIME[posterior[i]],
                    expectation[i]
                )
            )
        return '\n'.join(lines) + '\n'

    def likelihood_to_str(self):
        if self.likelihood == None or self.initial_likelihood == None:
            raise Error('Need to run estimate_params first')
        lines = [
            '# Likelihood of initial (first) and estimated (last)',
            '%.6f' % self.initial_likelihood,
            '%.6f' % self.likelihood
        ]
        return '\n'.join(lines) + '\n'

def main():
    parser = ArgumentParser(description='Parser for Arguments')

    parser.add_argument('initial_param_file', type=str, metavar='input1',
                        help='initial input file')
    parser.add_argument('fasta_file', type=str, metavar='input2',
                        help='input fasta sequence file')
    parser.add_argument('output_estimate', type=str, metavar='output1',
                        help='output file for estimated paramters')
    parser.add_argument('output_initial_decoding', type=str, metavar='output2',
                        help='output file for initial decodings')
    parser.add_argument('output_estimate_decoding', type=str, metavar='output3',
                        help='output file for estimated decodings')
    parser.add_argument('output_likelihood', type=str, metavar='output4',
                        help='output file for log likelihoods')
    args = parser.parse_args()

    tmrca = TSequence(args.initial_param_file, args.fasta_file)
    with open(args.output_initial_decoding, 'w') as f:
        f.write(tmrca.decoding_to_str('initial'))
        f.close()

    tmrca.estimate_params(max_iter=15)

    with open(args.output_estimate, 'w') as f:
        f.write(tmrca.theta_to_str())
        f.close()

    with open(args.output_estimate_decoding, 'w') as f:
        f.write(tmrca.decoding_to_str('estimate'))
        f.close()

    with open(args.output_likelihood, 'w') as f:
        f.write(tmrca.likelihood_to_str())
        f.close()

if __name__ == '__main__':
    main()
