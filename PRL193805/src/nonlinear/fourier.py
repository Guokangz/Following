"""Fourier utilities with the Supplemental Material's positive sign.

The convention is

``f(omega) = (1/sqrt(2pi)) integral dt exp(+i omega t) f(t)``,

matching SM Eq. (S.20).
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray


def positive_frequency_fft(
    time: ArrayLike,
    signal: ArrayLike,
    *,
    zero_padding_factor: int = 1,
) -> tuple[NDArray[np.float64], NDArray[np.complex128]]:
    """Return the uniformly sampled positive-sign Fourier transform."""

    time_array = np.asarray(time, dtype=float)
    signal_array = np.asarray(signal, dtype=complex)
    if time_array.ndim != 1 or signal_array.ndim != 1 or time_array.shape != signal_array.shape:
        raise ValueError("time and signal must be matching 1D arrays")
    if time_array.size < 2:
        raise ValueError("at least two time samples are required")
    spacing = np.diff(time_array)
    if not np.allclose(spacing, spacing[0], rtol=1e-10, atol=1e-13):
        raise ValueError("time grid must be uniform")
    if zero_padding_factor < 1 or int(zero_padding_factor) != zero_padding_factor:
        raise ValueError("zero_padding_factor must be a positive integer")

    dt = float(spacing[0])
    padded_size = int(zero_padding_factor) * time_array.size
    padded = np.zeros(padded_size, dtype=complex)
    padded[: signal_array.size] = signal_array
    frequency = 2.0 * np.pi * np.fft.fftfreq(padded_size, d=dt)
    transform = (
        padded_size
        * dt
        / np.sqrt(2.0 * np.pi)
        * np.fft.ifft(padded)
        * np.exp(1j * frequency * time_array[0])
    )
    return np.fft.fftshift(frequency), np.fft.fftshift(transform)


def direct_positive_frequency_transform(
    time: ArrayLike,
    signal: ArrayLike,
    frequency: ArrayLike,
) -> NDArray[np.complex128]:
    """Direct trapezoidal transform used as a low-dimensional benchmark."""

    time_array = np.asarray(time, dtype=float)
    signal_array = np.asarray(signal, dtype=complex)
    frequency_array = np.asarray(frequency, dtype=float)
    phase = np.exp(1j * frequency_array[..., None] * time_array)
    return np.trapezoid(phase * signal_array, x=time_array, axis=-1) / np.sqrt(2.0 * np.pi)

