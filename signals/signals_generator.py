import numpy as np
import matplotlib.pyplot as plt


def generate_rnd_message(N: int) -> np.array:
    msg = np.array([])

    for i in range(N):
        aleatorio = np.random.random()
        msg = np.append(msg, 1) if aleatorio >= 0.5 else np.append(msg, 0)

    return msg


def add_noise(msg: np.array,
              N: int,
              sigma_w: float = 0.2) -> np.array:
    """
        Add noise to signal
    :param msg:     message to add noise
    :param sigma_w: standard deviation for noise
    :param N:       number of samples
    :return:        np. Array of length N
    """
    j = 0
    rdm_noise = np.random.normal(0, sigma_w, N)
    for i in rdm_noise:
        np.put(msg, j, msg[j] + i)
        j += 1

    return msg


def generate_psk(fc: int,
                 A: int,
                 N: int,
                 msg=np.array([]),
                 noise: bool = True,
                 ):
    """
        Generate PSK modulated signal
    :param noise: determine if PSK signal has noise or not
    :param msg: message to be transmitted
    :param N: number of samples
    :param fc: carrier sine wave frequency
    :param A:  Amplitude
    :return: simple modulated PSK signal in time
    """
    steps = 1/N
    t = np.arange(0, 1, steps)

    if len(msg) < 1:
        msg = np.zeros(0, dtype=int)
        b = np.arange(0.2, 1.2, 0.2)
        s = [1]

        for i in t:
            msg = np.append(msg, s)
            if i == b[0]:
                b = np.delete(b, 0)
                s[0] = 1 if s[0] == 0 else 0

    psk = np.zeros(0, dtype=int)
    for i in range(N):
        if msg[i] == 1:
            psk = np.append(psk, A * np.sin(2 * np.pi * fc * t[i]))
        else:
            psk = np.append(psk, A * np.sin(2 * np.pi * fc * t[i]) * -1)

    if noise:
        psk = add_noise(msg, N)

    return t, psk


def generate_qpsk(sigma_s: float,
                  fc: int,
                  fs: int,
                  Rs: int,
                  N: int) -> np.array:
    """
        Generate QPSK modulated signal
    :param sigma_s: standard deviation
    :param fc:  carry frequency
    :param fs:  sample frequency
    :param Rs:  simbol rate bits/s
    :param N:   number of samples
    :return:    modulated QPSK signal
    """
    Ns = int(fs / Rs)       # Muestras por simbolo
    Nbits = int(N / Ns)     # Numero total de muestras
    t = np.r_[0.0: N] / fs  # Tiempos de muestreo

    inputBits = np.random.randn(Nbits, 1) > 0

    # Señales de portadora. carrier2 desfasada 90° respecto carrier1 -> modulacion IQ
    carrier1 = np.cos(2 * np.pi * fc * t)
    carrier2 = np.sin(2 * np.pi * fc * t)

    # Serial-to-Parallel Converter (1/2 of data rate)
    I_bits = inputBits[0::2]
    Q_bits = inputBits[1::2]

    # DAC
    I_signal = (np.tile(I_bits * 2 - 1, (1, 2 * Ns))).ravel()
    Q_signal = (np.tile(Q_bits * 2 - 1, (1, 2 * Ns))).ravel()

    # Multiplicador
    I_signal_modulated = I_signal * carrier1
    Q_signal_modulated = Q_signal * carrier2

    # Señal modulada
    signal = I_signal_modulated + Q_signal_modulated

    # Normalizamos la señal
    signal = (sigma_s / np.std(signal)) * signal

    # Señal QPSK
    return signal




