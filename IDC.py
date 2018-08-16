import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('Qt5Agg')

"""IDC_analytic returns an array, and a graph of the IDC curve of 
an MMPP given Tmin (minimum window size), and Tlength (the 
desired length of the array.  
also the expected # of spikes in a time interval
and the probability generating function for the number of arrivals"""

""" The inputs are lambda1, lambda2, trans1,
and trans2.  The algorithm comes from Lucantoni 1986.
lambda1 = nonburst rate
lambda2 = burst rate
trans1 = rate of transition nb-> b
trans2 = rate of transition b -> nb """

lambda1 = int(input('lambda1: '))
lambda2 = int(input('lambda2:  '))
r1 = int(input('trans1:  '))
r2 = int(input('trans2:  '))
T_min = int(input('T_min (in ms):  '))
T_length = int(input('T_length:  '))   # Need input to be int
# T_length = 0.5
T_array = T_min * np.arange(1, T_length)    # Need tow int types

# next section is for pgf
Pi = 1/(r1+r2) * np.array([r2, r1])
e = np.array([1, 1])
R = np.array([[-r1, r1], [r2, -r2]])
V = np.array([[11, 0],[0, 12]])

IDC_array = np.zeros(T_length-1)
ExpArray = np.zeros(T_length-1)
pgf = np.zeros(T_length-1)

for i_win_size in range(1, T_length-1):
    t = (i_win_size*T_min)/1000 # Previously we called this winsize, but it's put in seconds
    A_num = 2*(11-12)**2*r1*r2
    A_denom = (r1 + r2)**2*(11*r2 + 12*r1)
    B_num = A_num
    B_denom = (r1 + r2)**3*(11*r2 + 12*r1)*t

    A = A_num/A_denom
    B = B_num/B_denom
    C = 1 - np.exp(-(r1+r2)*t)

    IDC_array[i_win_size] = 1 + A - (B * C)
    ExpArray[i_win_size] = (11*r2 + 12*r1) / (r1 + r2)*t

plt.subplot(2, 1, 1)
plt.plot(T_array, IDC_array)
plt.ylabel('Index for the dispersion of counts')
plt.title('IDC curve')

plt.subplot(2, 1, 2)
plt.plot(T_array, ExpArray)
plt.xlabel('Window size (ms)')
plt.ylabel('Expected number of spikes')
plt.title('Expected number of spikes')