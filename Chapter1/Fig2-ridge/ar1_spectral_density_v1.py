"""Generate the AR(1) spectral-density illustration used in Chapter 3.

The plot compares the theoretical spectral density of an AR(1) process with
the raw and smoothed periodogram from one simulated realization.

Run from the book root. The figure is saved under ./figures/.
"""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt


def main() -> None:
    output_dir = Path("figures")
    output_dir.mkdir(parents=True, exist_ok=True)

    rng = np.random.default_rng(20260510)
    T = 512
    phi = 0.8
    sigma_eps = 1.0

    x = np.empty(T)
    x[0] = rng.normal(scale=sigma_eps / np.sqrt(1 - phi**2))
    eps = rng.normal(scale=sigma_eps, size=T)

    for t in range(1, T):
        x[t] = phi * x[t - 1] + eps[t]

    x_centered = x - x.mean()
    fft = np.fft.fft(x_centered)

    omega = 2 * np.pi * np.arange(T // 2 + 1) / T
    periodogram = (np.abs(fft[: T // 2 + 1]) ** 2) / (2 * np.pi * T)
    f_theory = sigma_eps**2 / (2 * np.pi * (1 + phi**2 - 2 * phi * np.cos(omega)))

    m = 7
    kernel = np.ones(2 * m + 1) / (2 * m + 1)
    periodogram_smooth = np.convolve(periodogram, kernel, mode="same")

    plt.figure(figsize=(7.2, 4.6))
    plt.plot(omega, f_theory, linewidth=2, label="Theoretical spectral density")
    plt.plot(omega, periodogram, linewidth=0.8, alpha=0.45, label="Raw periodogram")
    plt.plot(omega, periodogram_smooth, linewidth=1.5, label="Smoothed periodogram")
    plt.xlabel(r"Frequency $\omega$")
    plt.ylabel(r"Spectral density")
    plt.title(r"AR(1) spectral density and periodogram, $\phi=0.8$")
    plt.legend()
    plt.tight_layout()

    plt.savefig(output_dir / "esttest_ar1_spectral_density.pdf")
    plt.savefig(output_dir / "esttest_ar1_spectral_density.png", dpi=200)
    plt.close()


if __name__ == "__main__":
    main()
