from compute_threshold import *
import matplotlib.pyplot as plt

N = 15000
sigma_n = 1
SNRp = range(-30, -6, 2)
Pf_fixed = 0.10
Pd_fixed = 0.90
Pm_fixed = 1 - Pd_fixed

op_threshold = []
lf, lm, le = [], [], []
i = 0

for gamma_p in SNRp:
    gamma_p = dB_to_times(gamma_p)
    lf.append(get_lf(N, Pf_fixed, sigma_n))
    lm.append(get_lm(N, Pm_fixed, sigma_n, gamma_p))
    le.append(get_le(N, sigma_n, gamma_p))
    op_threshold.append(op_threshold_sel(lf[i], lm[i], le[i], Pf_fixed, Pd_fixed, sigma_n, gamma_p))
    i += 1

print(op_threshold)

plt.plot(SNRp, lf)
plt.plot(SNRp, lm)
plt.plot(SNRp, le)
plt.plot(SNRp, op_threshold)
plt.xlabel("SNRp [dB]")
plt.ylabel("threshold")
plt.title("test")
plt.grid(True)
plt.show()

