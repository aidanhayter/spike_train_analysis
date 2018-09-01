import numpy as np
import matplotlib.pyplot as plt
from spike_train import create_spikes

"""histogram_data calls the function from spike_train to create random data,
as a Poisson model or MMPP model.
This code then creates statistical histograms from the data."""

plt.switch_backend('Qt5Agg')
model = 'Poisson'
spike_train = create_spikes(model)

# Spike Times
spike_times = np.where(spike_train)

# Spike Intervals
spike_intervals = np.diff(spike_times[0])

# Plots the IDC for the random data
np.var(spike_train) / np.mean(spike_train)
tau_min = 1
tau_max = 1000      # CORRECT: This is in time bins, need in ms
idc = [0] * tau_max
num_samples = 100000

for tau in range(tau_min, tau_max):
    sptr_tau, t_tau = np.histogram(spike_times, bins=int(num_samples/tau))
    idc[tau] = np.var(sptr_tau) / np.mean(sptr_tau)

plt.subplot(3, 1, 1)
plt.title('IDC (FF)')
plt.xlabel('Bin Width (ms/10)')
plt.ylabel('IDC')
plt.plot(idc[1:])

# Plots the analytic data
lambda1 = 30    # Number must tbe adjusted
lambda2 = 30    # Number must tbe adjusted
r1 = 101        # Number must tbe adjusted
r2 = 50         # Number must tbe adjusted
T_min = int(input('T_min (in ms):  '))  # Need input to be int
T_length = int(input('T_length:  '))
# T_length = 0.5
T_array = T_min * np.arange(1, T_length)

# Next section is for pgf
Pi = 1/(r1+r2) * np.array([r2, r1])
e = np.array([1, 1])
R = np.array([[-r1, r1], [r2, -r2]])
V = np.array([[11, 0],[0, 12]])

IDC_array = np.zeros(T_length-1)
ExpArray = np.zeros(T_length-1)
pgf = np.zeros(T_length-1)

for i_win_size in range(1, T_length-1):
    t = (i_win_size*T_min)/1000     # Previously we called this winsize, but it's put in seconds
    A_num = 2*(11-12)**2*r1*r2
    A_denom = (r1 + r2)**2*(11*r2 + 12*r1)
    B_num = A_num
    B_denom = (r1 + r2)**3*(11*r2 + 12*r1)*t

    A = A_num/A_denom
    B = B_num/B_denom
    C = 1 - np.exp(-(r1+r2)*t)

    IDC_array[i_win_size] = 1 + A - (B * C)
    ExpArray[i_win_size] = (11*r2 + 12*r1) / (r1 + r2)*t

plt.subplot(3, 1, 2)
plt.plot(T_array, IDC_array)
plt.ylabel('IDC')
plt.title('Index for the dispersion of counts')

plt.subplot(3, 1, 3)
plt.plot(T_array, ExpArray)
plt.xlabel('Window size (ms)')
plt.ylabel('Expected number of spikes')
plt.title('Expected number of spikes')