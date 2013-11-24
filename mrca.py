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
from multiprocessing import Pool, Array

from util import time_it

# Necessary for allowing pickling of instance methods in concurrent processes
import copy_reg
import types

def reduce_method(m):
    return (getattr, (m.__self__, m.__func__.__name__))

copy_reg.pickle(types.MethodType, reduce_method)

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


@time_it
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


@time_it
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

    def __init__(self, observ, init_params, thresh=1e-3, max_iter=50):
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
        self.join_forward_backward()
        self.lhood = self.calculate_likelihood()

        # iterate until improvement meets threshold
        while self.check_improvement(self.lhood, old_lhood) and i < self.max_iter:
            print ' ' + str(i) + '    ' + str(self.lhood)
            old_lhood = self.lhood

            # self.iteration() will always have access to the correct
            # forward and backward tables based on the current self.theta
            self.theta = self.iteration()

            # Update forward and backward tables
            self.process_forward_algorithm()
            self.process_backward_algorithm()
            self.join_forward_backward()
            self.lhood = self.calculate_likelihood()
            i += 1

        return self.theta

    def check_improvement(self, new_lhood, old_lhood):
        return new_lhood - old_lhood > self.thresh

    @time_it
    def iteration(self):
        """Generate a new theta."""
        print self.theta.m
        states = self.theta.m.keys()
        
        # do all the actual crap here...
        marginal = {}
        marginal_norm = 0
        for k in states:
            marginal[k] = self.calculate_marginal_for_state(k)
            marginal_norm += marginal[k]
        # need to normalize
        for state, val in marginal.items():
            marginal[state] = val / marginal_norm

        transition = {}

        # Run these calculations asynchronously
        async_results = []
        for i in states:
            for j in states:
                async_results.append(
                    pool.apply_async(self.calculate_transition, (i,j))
                )
        async_e_results = []
        symbols = self.theta.e[states[0]].keys()

        for k in states:
            for symbol in symbols:
                async_e_results.append(
                    pool.apply_async(self.calculate_emission, (k, symbol))
                )

        # Gather the asynchronous results
        for i in states:
            transition[i] = {}
            trans_norm = 0
            for j in states:
                transition[i][j] = async_results.pop(0).get()
                #transition[i][j] = self.calculate_transition(i,j)
                trans_norm += transition[i][j]
            # need to normalize
            for j, val in transition[i].items():
                transition[i][j] = val / trans_norm

        emission = {}

        for k in states:
            emission[k] = {}
            emission_norm = 0
            for symbol in symbols:
                emission[k][symbol] = async_e_results.pop(0).get()
                emission_norm += emission[k][symbol]
            # need to normalize
            for symbol, val in emission[k].items():
                emission[k][symbol] = val / emission_norm

        return Theta(marginal, transition, emission)

    @time_it
    def calculate_marginal_for_state(self, k):
        """Return the marginal expectation for being in state k given theta."""
        return e**(self.f(k)[0] + self.b(k)[0] - self.lhood)

    @time_it
    def calculate_transition(self, i, j):
        """Return the transition expectation for A_{ij} given theta."""
        num_samp = len(self.x) - 1

        D = max([self.f(i)[t] + self.b(j)[t+1] for t in range(num_samp)])
        inside = lambda s: (e**(self.f(i)[s] + self.b(j)[s + 1] - D)) * \
                    self.theta.e[j][self.x[s + 1]]

        term = log(sum([inside(t) for t in range(num_samp)])) + D - self.lhood
        return self.theta.a[i][j] * (e**term)

    @time_it
    def calculate_emission(self, k, sigma):
        D = max([self.f(k)[t] + self.b(k)[t] for t in range(len(self.x)) \
                if self.x[t] == sigma])
        term = log(sum([e**(self.f(k)[t] + self.b(k)[t] - D) \
                    for t in range(len(self.x)) if self.x[t] == sigma]))
        return e**(term + D - self.lhood)

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
        #self.async_forward = pool.apply_async(
        #    forward_algorithm, (self.theta, self.x)
        #)
        self.forward = forward_algorithm(self.theta, self.x)

    def process_backward_algorithm(self):
        self.backward = backward_algorithm(self.theta, self.x)

    def join_forward_backward(self):
        #self.forward = self.async_forward.get()
        pass

    def f(self, state):
        return self.forward[state]

    def b(self, state):
        return self.backward[state]


class Decoding(object):

    def __init__(self, observations, params):
        self.x = observations
        self.theta = params
        self.forward = forward_algorithm(self.theta, self.x)
        self.backward = backward_algorithm(self.theta, self.x)
        self.lhood = self.calculate_likelihood()

    def calculate_likelihood(self, theta=None):
        """Return the log-likelihood of the data given theta."""
        states = self.theta.m.keys()
        D = max([self.forward[k][-1] for k in states])

        return log(sum([e**(self.forward[k][-1] - D) for k in states])) + D

    def posterior(self):
        states = self.theta.m.keys()
        ck = lambda x: x[1]
        results = []
        for t in range(len(self.x)):
            q, prob = max([(k, self.calc_marg(k, t)) for k in states], key=ck)
            results.append(q)
        return results

    def calc_marg(self, k, t):
        return e**(self.f(k)[t] + self.b(k)[t] - self.lhood)

    def viterbi(self):
        return []

    def f(self, k):
        return self.forward[k]

    def b(self, k):
        return self.backward[k]


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
        self.obs = observe_differences(sequences[0], sequences[1])
        self.theta = SEQUENCE_TYPE_TO_THETA[s]
        self.estimate = None

    def initial_decoding(self):
        decode = Decoding(obs, self.theta)
        # need to also return Viterbi
        return decode.posterior()

    def estimate_params(self):
        em = EM(obs, self.theta, max_iter=10)
        self.estimate = em.estimate_params()

    def estimate_decoding(self):
        if self.estimate == None:
            raise Error('Must run estimate_params first!')
        decode = Decoding(obs, self.estimate)
        # need to also return viterbi
        return decode.posterior()


pool = Pool(processes=3)

if __name__ == '__main__':
    pass
