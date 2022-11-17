from signals_generator import *
import numpy as np
import matplotlib.pyplot as plt

A = 5
fc = 5
n_samples = 100

msg = generate_rnd_message(n_samples)
t, myPSK = generate_psk(A, fc, n_samples, msg, noise=False)

plt.plot(t, myPSK)
plt.xlabel("t")
plt.ylabel("y")
plt.title("PSK")
plt.grid(True)
plt.show()

