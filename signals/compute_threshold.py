import numpy as np
from math import sqrt, log
from scipy import special


def q_func(val):
    aux = 0.5 - 0.5 * special.erf(val / sqrt(2))
    return aux


def q_func_inv(val):
    aux = sqrt(2) * special.erfinv(1 - 2 * val)
    return aux


def dB_to_times(db):
    return 10 ** (db / 10)


def get_lf(N: int, pf_fix, sigma_n):  # false alarm probability, var = varianza
    lf = (sqrt(2 / N) * special.erfcinv(2 * pf_fix) + 1) * N * sigma_n
    return int(lf)


def get_lm(N: int, pm_fix, sigma_n, gamma_p):  # miss detection probability
    lm = N * sigma_n * (sqrt(2 / N) * (1 + gamma_p) * special.erfcinv(2 * (1 - pm_fix)) + (1 + gamma_p))
    return int(lm)


def get_le(N: int, sigma_n, gamma_p):  # error probability
    le = (N * sigma_n / 2) * (1 + sqrt(1 + (2 * (2 + gamma_p) * log(1 + gamma_p)) / (N * gamma_p))) * (
                (1 + gamma_p) / (1 + gamma_p / 2))
    return int(le)


def op_threshold_sel(lf: int,
                     lm: int,
                     le: int,
                     pf_fix: float,
                     pd_fix: float,
                     sigma_n: float,
                     gamma_p: float):
    if lf <= lm:
        if lf <= le <= lm:
            lopt = le
        else:
            if le < lf:
                lopt = lf
            else:
                lopt = lm
    else:
        # no hay valor optimo de lf para predefinir el valor de N
        # computar el valor N*
        N_aux = (1 / gamma_p ** 2)*(q_func_inv(pf_fix) - q_func_inv(pd_fix) * sqrt(2 * gamma_p + 1)) ** 2

        # encontrar el valor de lf*, lm*, le* en N* para el valor correspondiente a SNRp
        lf_aux = get_lf(N_aux, pf_fix, sigma_n)
        lm_aux = get_lm(N_aux, 1 - pd_fix, sigma_n, gamma_p)
        le_aux = get_le(N_aux, sigma_n, gamma_p)

        # lm_aux = le_aux
        # lf_aux = lm_aux
        lopt = op_threshold_sel(lf_aux, lm_aux, le_aux, pf_fix, pd_fix, sigma_n, gamma_p)

    return lopt

