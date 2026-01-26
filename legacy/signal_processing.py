"""
Signal Processing Module (Legacy)

NOTE: This is legacy code preserved for reference.
Adds signal smoothing functionality to hydrological data.
"""

import numpy as np
from scipy.ndimage import uniform_filter1d, gaussian_filter1d
from scipy import signal
import matplotlib.pyplot as plt


class SignalProcessor:
    """
    Signal processing utilities for hydrological time series.
    Can be mixed into existing classes or used standalone.
    """

    def apply_moving_average_filter(self, data, window_size):
        """
        Apply moving average filter to smooth signals.

        Parameters
        ----------
        data : array-like
            Input signal data
        window_size : int
            Size of the moving window

        Returns
        -------
        numpy.ndarray
            Filtered signal
        """
        if window_size < 1:
            raise ValueError("Window size must be >= 1")
        if window_size > len(data):
            raise ValueError("Window size cannot be larger than data length")

        filtered_signal = uniform_filter1d(data, size=window_size, mode='nearest')
        return filtered_signal

    def apply_gaussian_filter(self, data, sigma):
        """
        Apply Gaussian filter for smooth signal filtering.

        Parameters
        ----------
        data : array-like
            Input signal data
        sigma : float
            Standard deviation of Gaussian kernel

        Returns
        -------
        numpy.ndarray
            Gaussian-filtered signal
        """
        if sigma <= 0:
            raise ValueError("Sigma must be positive")

        filtered_signal = gaussian_filter1d(data, sigma=sigma)
        return filtered_signal

    def apply_butterworth_filter(self, data, cutoff_freq, filter_order=4):
        """
        Apply Butterworth low-pass filter.

        Parameters
        ----------
        data : array-like
            Input signal data
        cutoff_freq : float
            Cutoff frequency (normalized: 0 < cutoff_freq < 1)
        filter_order : int, optional
            Filter order (default: 4)

        Returns
        -------
        numpy.ndarray
            Butterworth-filtered signal
        """
        if not (0 < cutoff_freq < 1):
            raise ValueError("Cutoff frequency must be between 0 and 1")
        if filter_order < 1:
            raise ValueError("Filter order must be >= 1")

        b, a = signal.butter(filter_order, cutoff_freq, btype='low')
        filtered_signal = signal.filtfilt(b, a, data)
        return filtered_signal

    def compare_signals(self, original, filtered, time_array=None, title="Signal Comparison"):
        """
        Compare original and filtered signals.

        Parameters
        ----------
        original : array-like
            Original signal
        filtered : array-like
            Filtered signal
        time_array : array-like, optional
            Time array for x-axis
        title : str, optional
            Plot title
        """
        if time_array is None:
            time_array = np.arange(len(original))

        fig, axes = plt.subplots(2, 1, figsize=(12, 8))

        axes[0].plot(time_array, original, 'b-', linewidth=2, label='Original')
        axes[0].plot(time_array, filtered, 'r-', linewidth=2, label='Filtered')
        axes[0].set_ylabel('Signal Value')
        axes[0].set_title(title)
        axes[0].legend()
        axes[0].grid(True, alpha=0.3)

        difference = original - filtered
        axes[1].plot(time_array, difference, 'g-', linewidth=1)
        axes[1].set_ylabel('Difference')
        axes[1].set_xlabel('Time')
        axes[1].set_title('Signal Difference (Original - Filtered)')
        axes[1].grid(True, alpha=0.3)

        plt.tight_layout()
        plt.show()

        rmse = np.sqrt(np.mean(difference**2))
        max_diff = np.max(np.abs(difference))
        print(f"Filter performance metrics:")
        print(f"  RMSE: {rmse:.4f}")
        print(f"  Max difference: {max_diff:.4f}")
        print(f"  Signal reduction: {(1 - np.std(filtered)/np.std(original))*100:.1f}%")
