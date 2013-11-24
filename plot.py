import numpy as np
import matplotlib.pyplot as plt

from mrca import convert_state_seq_to_time_seq
from mrca import TSequence

def read_true_tmrca():
    with open('data/true_tmrca.txt') as f:
        lines = f.readlines()
        return [float(line) for line in lines]


TRUE_TMRCA_SEQ = read_true_tmrca()


def plot_initial(variant, sequence=None):
    if variant not in ('mu', '2mu', '5mu'):
        raise Error('plot_initial not ran with correct variant')
    if sequence:
        posterior = sequences[0]
        viterbi = sequences[1]
    else:
        t = TSequence(variant)
        posterior, viterbi = t.decode_initial()
    posterior_decoding = convert_state_seq_to_time_seq(
        posterior
    )
    line1 = plt.plot(posterior_decoding)
    true_line = plt.plot(TRUE_TMRCA_SEQ)
    plt.legend(
        (line1[0], true_line[0]),
        ('Initial Posterior', 'True Tmrca')
    )
    plt.show()


def plot_estimate(variant, sequences=None):
    if variant not in ('mu', '2mu', '5mu'):
        raise Error('plot_estimate not ran with correct variant')
    if sequences:
        posterior = sequences[0]
        viterbi = sequences[1]
    else:
        t = TSequence(variant)
        t.estimate_params(max_iter=10)
        posterior, viterbi = t.decode_estimate()
    posterior_decoding = convert_state_seq_to_time_seq(
        posterior
    )
    line1 = plt.plot(posterior_decoding)
    true_line = plt.plot(TRUE_TMRCA_SEQ)
    plt.legend(
        (line1[0], true_line[0]),
        ('Estimate Posterior', 'True Tmrca')
    )
    plt.show()
