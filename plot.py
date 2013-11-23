import numpy as np
import matplotlib.pyplot as plt

from mrca import run_mu_initial_decoding
from mrca import run_mu_estimate_decoding
from mrca import convert_state_seq_to_time_seq

def read_true_tmrca():
    with open('data/true_tmrca.txt') as f:
        lines = f.readlines()
        return [float(line) for line in lines]


TRUE_TMRCA_SEQ = read_true_tmrca()


def plot_initial_mu():
    posterior_decoding = convert_state_seq_to_time_seq(
        run_mu_initial_decoding()
    )
    line1 = plt.plot(posterior_decoding)
    true_line = plt.plot(TRUE_TMRCA_SEQ)
    plt.legend(
        (line1[0], true_line[0]),
        ('Initial Posterior Decoding', 'True Tmrca')
    )
    plt.show()

def plot_estimate_mu():
    init_posterior_decoding = convert_state_seq_to_time_seq(
        run_mu_initial_decoding()
    )
    posterior_decoding = convert_state_seq_to_time_seq(
        run_mu_estimate_decoding()
    )
    line1 = plt.plot(init_posterior_decoding)
    line2 = plt.plot(posterior_decoding)
    true_line = plt.plot(TRUE_TMRCA_SEQ)
    plt.legend(
        (line1[0], line2[0], true_line[0]),
        ('Initial Posterior', 'Estimate Posterior', 'True Tmrca')
    )
    plt.show()
