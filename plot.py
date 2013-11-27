import numpy as np
import matplotlib.pyplot as plt

from mrca import convert_state_seq_to_time_seq
from mrca import TSequence

def read_true_tmrca():
    with open('data/true_tmrca.txt') as f:
        lines = f.readlines()
        return [float(line) for line in lines]


def read_example_true_tmrca():
    with open('example/true_tmrca_test.txt') as f:
        lines = f.readlines()
        return [float(line) for line in lines]


TRUE_TMRCA_SEQ = read_true_tmrca()
EX_TRUE_TMRCA_SEQ = read_example_true_tmrca()


def plot_initial(variant, sequences=None):
    if variant not in ('mu', '2mu', '4mu', '5mu'):
        raise Error('plot_initial not ran with correct variant')
    if sequences:
        posterior = sequences[0]
        viterbi = sequences[1]
        mean = sequences[2]
    else:
        t = TSequence(variant)
        posterior, viterbi, mean = t.decode_initial()
    posterior_decoding = convert_state_seq_to_time_seq(
        posterior
    )

    # perform plotting
    line1 = plt.plot(posterior_decoding)
    line3 = plt.plot(mean)
    if variant == '4mu':
        true_line = plt.plot(EX_TRUE_TMRCA_SEQ)
    else:
        true_line = plt.plot(TRUE_TMRCA_SEQ)
    plt.legend(
        (line1[0], line3[0], true_line[0]),
        ('Initial Posterior', 'Mean', 'True Tmrca')
    )
    plt.show()


def plot_estimate(variant, sequences=None):
    if variant not in ('mu', '2mu', '4mu', '5mu'):
        raise Error('plot_estimate not ran with correct variant')
    if sequences:
        posterior = sequences[0]
        viterbi = sequences[1]
        mean = sequences[2]
    else:
        t = TSequence(variant)
        t.estimate_params(max_iter=15)
        posterior, viterbi, mean = t.decode_estimate()
    posterior_decoding = convert_state_seq_to_time_seq(
        posterior
    )

    # perform plotting
    line1 = plt.plot(posterior_decoding)
    line3 = plt.plot(mean)
    if variant == '4mu':
        true_line = plt.plot(EX_TRUE_TMRCA_SEQ)
    else:
        true_line = plt.plot(TRUE_TMRCA_SEQ)
    plt.legend(
        (line1[0], line3[0], true_line[0]),
        ('Estimate Posterior', 'Mean', 'True Tmrca')
    )
    plt.show()
