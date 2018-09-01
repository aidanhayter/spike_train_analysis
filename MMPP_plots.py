import numpy as np
import matplotlib.pyplot as plt
from spike_train import create_spikes
import math
plt.switch_backend('Qt5Agg')

"""histogram_data calls the function from spike_train to create random data,
as a Poisson model or MMPP model.
This code then creates statistical histograms from the data."""

model = 'MMPP'
spike_train = create_spikes(model)

# Spike Times
spike_times = np.where(spike_train)

# Spike Intervals
spike_intervals = np.diff(spike_times[0])

# Plotting IIH for MMPP process random data
plt.subplot(2, 1, 1)
iih_array = [1] * len(spike_intervals)
plt.title('Interspike Interval Histogram')
plt.ylabel('Spikes')
plt.xlabel('Time')
y, bin_edges = np.histogram(spike_intervals, bins=100)
bin_centres = (bin_edges - np.diff(bin_edges)[0]/2)[1:]
plt.plot(bin_centres, y)

# Plots analytic data of MMPP as IIH
trans1 = 101
trans2 = 50
lambda1 = 30
lambda2 = 30
t_length = int(input('t_length:  '))
n_bins = int(input('n_bins:  '))
int_length = t_length / n_bins
a = trans1
b = trans2
c1 = lambda1
c2 = lambda2
t = np.arange(t_length/n_bins, t_length, int_length)    # Begins at 0.0005 and goes to 0.5 with step 0.0005


A = a**2 + 2*a*c1 + 2*a*b - 2*c2*a + c1**2 - 2*c1*b - 2*c1*c2 + b**2 + 2*b*c2 + c2**2
B = a + c1 + b + c2 + math.sqrt(A)
C = a + c1 + b + c2 - math.sqrt(A)

L1 = c1**2*b**2*math.sqrt(A) + 2*c1*b*c2*math.sqrt(A)*a + c2**3*a*c1 - c2**2*a**2*b
L2 = -2*c2**2*a**2*c1 - c1**2*b**2*a + c1**3*b**2*a + c1**3+b**2 + c2**2*a**3 - c1**2*b**3
L3 = c1**2*b*c2*math.sqrt(A) - 2*c1*b*a**2*c2 + c1**3*b*c2 - c1**2*b*c2**2 - 2*c1*b**2*c2*a
L4 = -c2*a*c1**2*b - c2**2*a*c1**2 - 2*c1**2*b**2*c2 - c1*b*c2**2*a + c2**2*a**2*math.sqrt(A)
L5 = c2**2*a*math.sqrt(A)*c1
L6 = c1**2*b**2*math.sqrt(A) + c1**2*b*c2*math.sqrt(A) - c1**3*b**2 + c1**2*b**3 - c2**3*a**2 + c2**2*a**3 + 2*c2**2*a**2*c1
L7 = 2*c1*b*c2*math.sqrt(A)*a - c1**3*b*c2 + 2*c1*b*a**2*c2 + c1*b*c2*a + c2**2*a*math.sqrt(A)*c1
L8 = c1**2*b*c2**2 + c2**2*a**2*b + c2**2*a*c1**2 + c2*a*c1**2*b - c2**3*a*c1 + c1**2*b**2*a
L9 = 2*c1**2*b**2*c2 + 2*c1*b**2*c2*a + c2**2*a**2*math.sqrt(A)

part1 = L1 + L2 + L3 + L4 + L5
part2 = L6 + L7 + L8 + L9
denom = (c1*b + c2*a)*C*math.sqrt(A)*B

prob = -2*(part1*np.exp(-0.5*B*t) + part2*part2*np.exp(-0.5*C*t)*(-0.5*C)) / denom
plt.subplot(2, 1, 2)
plt.ylabel('Spikes')
plt.xlabel('Time')
plt.plot(t, prob)
