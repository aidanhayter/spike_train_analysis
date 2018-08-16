import numpy as np

"""spike_train creates random data in a Poisson model or MMPP model"""

def create_spikes(model):

    num_samples = 100000
    s_train = np.zeros(num_samples)

    if model == 'Poisson':
        p = 0.01
        r = np.random.rand(len(s_train))
        s_train = r > (1-p)
        spike_train = s_train.astype(int)

    elif model == 'MMPP':   # Creates array of bursts and non bursts based on Markov process
        nb2b = .003  # probability of going from non burst to burst
        b2nb = .05  # probability of going from burst to non burst
        nbfr = .001  # non burst fire rate
        bfr = .1    # burst firing rate
        state_train = np.zeros(num_samples) # Array of zeroes, of size num_samples
        state1_p = nb2b / b2nb  # 1/5

        if np.random.rand() >= state1_p: # If random number is greater than or equal to 1/5
            state1 = 0
        else:
            state1 = 1

        state_train[0] = state1    # Assigns 1 or 0 to first index of array


        for i_state in range (1, num_samples):  # Assigns 1 or 0 to next index, based on previous value
            if state_train[i_state - 1]:    # Burst state
                p_state = b2nb      # Probability is for burst to non burst
                if np.random.rand() >= p_state:
                    state_train[i_state] = 1    # Stays as burst
                else:
                    state_train[i_state] = 0    # Moves to non burst

            else:
                p_state = nb2b      # Probability is for non burst to burst
                if np.random.rand() >= p_state:
                    state_train[i_state] = 0    # Stays as non burst
                else:
                    state_train[i_state] = 1    # Moves to burst

        # Creates array opposite to state_train
        nb_state_train = state_train + 1
        nb_state_train = nb_state_train - (2 * state_train)

        # Non burst probability
        nb_p = nbfr * nb_state_train
        temp_rand_array = np.random.rand(num_samples)
        nb_spike_train = temp_rand_array < nb_p

        # Burst probability
        b_p = bfr * state_train
        b_spike_train = temp_rand_array < b_p

        spike_train = nb_spike_train + b_spike_train

    return spike_train