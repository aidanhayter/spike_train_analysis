import numpy as np
import matplotlib.pyplot as plt
from spike_train import create_spikes

"""histogram_data calls the function from spike_train to create random data,
as a Poisson model or MMPP model.
This code then creates statistical histograms from the data."""

plt.switch_backend('Qt5Agg')
model = 'MMPP'
spike_train = create_spikes(model)

# Spike Times
spike_times = np.where(spike_train)

# Spike Intervals
spike_intervals = np.diff(spike_times[0])

# Plotting PSTH
psth_array = [1] * len(spike_times[0])
plt.title('Peristimulus Time Histogram')
plt.ylabel('Spikes')
plt.xlabel('Time')
y, bin_edges = np.histogram(spike_times, bins=100)
bin_centers = (bin_edges - np.diff(bin_edges)[0]/2)[1:]
plt.plot(bin_centers, y)

# Train of spikes
plt.title('Spike Train')
plt.ylabel('Spikes')
plt.xlabel('Time')
plt.yticks([])
plt.scatter(spike_times[0], psth_array, marker='|')

# Plotting IIH
iih_array = [1] * len(spike_intervals)
plt.title('Interspike Interval Histogram')
plt.ylabel('Spikes')
plt.xlabel('Time')
y, bin_edges = np.histogram(spike_intervals, bins=100)
bin_centres = (bin_edges - np.diff(bin_edges)[0]/2)[1:]
plt.plot(bin_centres, y)

# Plotting IDC (FF)
np.var(spike_train) / np.mean(spike_train)
tau_min = 1
tau_max = 1000      # CORRECT: This is in time bins, need in ms
idc = [0] * tau_max
num_samples = 100000

for tau in range(tau_min, tau_max):
    sptr_tau, t_tau = np.histogram(spike_times, bins=int(num_samples/tau))
    idc[tau] = np.var(sptr_tau) / np.mean(sptr_tau)

plt.title('IDC (FF)')
plt.xlabel('Bin Width (ms/10)')
plt.ylabel('IDC')
plt.plot(idc[1:])
