from __future__ import annotations

import unittest

import numpy as np

from src.nonlinear.fourier import (
    direct_positive_frequency_transform,
    positive_frequency_fft,
)


class TestPositiveFrequencyFourierTransform(unittest.TestCase):
    def test_gaussian_matches_sm_eq_s20(self) -> None:
        dt = 0.002
        time = np.arange(-4.0, 8.0, dt)
        center = 1.3
        width = 0.2
        signal = np.exp(-0.5 * ((time - center) / width) ** 2)
        frequency, transform = positive_frequency_fft(
            time,
            signal,
            zero_padding_factor=4,
        )
        analytic = width * np.exp(-0.5 * (frequency * width) ** 2) * np.exp(
            1j * frequency * center
        )
        mask = np.abs(frequency) < 15.0
        self.assertLess(np.max(np.abs(transform[mask] - analytic[mask])), 2e-13)

    def test_fft_matches_direct_transform(self) -> None:
        dt = 0.01
        time = np.arange(0.0, 20.0, dt)
        signal = np.exp(-0.4 * time) * np.exp(-1.7j * time)
        frequency, transform = positive_frequency_fft(time, signal)
        indices = np.array([900, 950, 1000, 1050, 1100])
        direct = direct_positive_frequency_transform(time, signal, frequency[indices])
        # FFT uses a rectangle sum, while the direct benchmark uses trapezoids;
        # the difference is first-order in dt at the nonzero initial endpoint.
        self.assertLess(np.max(np.abs(transform[indices] - direct)), 0.005)


if __name__ == "__main__":
    unittest.main()

