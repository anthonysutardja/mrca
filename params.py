"""
Module containing all the initial parameters, as well as Theta class.
"""
from collections import namedtuple

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


INITIAL_2MU_PARAMS = Theta(
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
        1: {'I': 0.999217, 'D': 0.000782947},
        2: {'I': 0.996674, 'D': 0.00332648},
        3: {'I': 0.991725, 'D': 0.00827535},
        4: {'I': 0.983241, 'D': 0.0167592},
    }
)


INITIAL_5MU_PARAMS = Theta(
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
        1: {'I': 0.998046, 'D': 0.00195405},
        2: {'I': 0.99173, 'D': 0.00826963},
        3: {'I': 0.979578, 'D': 0.0204217},
        4: {'I': 0.959163, 'D': 0.040837},
    }
)


# Example for testing

INITIAL_4MU_PARAMS = Theta(
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


def parse_params(file_name):
    with open(file_name) as f:
        lines = [line.strip() for line in f.readlines()]
        lines = [line for line in lines \
            if len(line) != 0 and line[0] != '#']

        # Parse marginal probabilities
        m = {}
        for i in range(4):
            state, prob = lines.pop(0).split()
            m[int(state)] = float(prob)

        # Parse transitional probabilities
        t = {}
        for i in range(1, 5):
            a, b, c, d = lines.pop(0).split()
            t[i] = {}

            t[i][1] = float(a)
            t[i][2] = float(b)
            t[i][3] = float(c)
            t[i][4] = float(d)

        # Parse emission probabilities
        e = {}
        for i in range(1,5):
            k, I, D = lines.pop(0).split()
            e[int(k)] = {}
            e[int(k)]['I'] = float(I)
            e[int(k)]['D'] = float(D)
        
        return Theta(m, t, e)
