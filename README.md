# Statistical and Numerical Precision Analysis of Power Spectral Density Estimators

## Overview

This project provides a C-based implementation for a deep analysis of Power Spectral Density (PSD) estimation techniques. It focuses on a signal composed of two sinusoids embedded in additive white Gaussian noise (AWGN) and explores the trade-offs between different estimators, with a primary focus on the impact of numerical precision on the results.

The project is divided into four key analyses:

1.  **Comparative Analysis of PSD Estimators**: It generates and plots a comparison between four key methods: the Periodogram (with rectangular and Hamming windows), the Multitaper method, and Welch's method. This plot highlights the fundamental trade-offs between spectral resolution, variance, and spectral leakage.

2.  **Statistical Distribution Analysis**: It investigates the probability distribution of the Welch PSD estimate at a specific frequency. It calculates the **skewness** and **kurtosis** of this distribution across many realizations for different noise levels (`sigma`). Furthermore, it uses the **Kolmogorov-Smirnov (K-S) test** to numerically quantify the deviation of the resulting distribution from a Gaussian model.

3.  **Numerical Precision Analysis**: This is the core of the project. The entire simulation is run at three distinct levels of numerical precision to quantify the impact on accuracy and performance:
    *   **Standard Precision**: Using the standard `double` type (64-bit IEEE 754).
    *   **High Precision**: Using the `long double` type (80-bit x86 extended precision).
    *   **Arbitrary Precision**: Using the **GNU MPFR library** with a high precision of **4096 bits** as a "gold standard" reference.

4.  **Performance and Resource Analysis**: The project meticulously measures and reports the **CPU execution time** and **maximum memory usage** for each of the three precision levels, demonstrating the significant computational cost associated with higher precision.

All simulations are performed by the compiled C code, which leverages the FFTW3 library for Fourier transforms. The project then generates comprehensive visualizations and reports using a Python/Matplotlib script.

**System Information:** This project was developed and tested on a machine running a Linux-based OS (Linux Mint).

## Theoretical Background

### 1. The Periodogram

The Periodogram is derived directly from the Discrete Fourier Transform (DFT). For a signal $x[n]$ of length $N$ with sampling frequency $F_S$, the DFT is:

$$
X[k] = \sum_{n=0}^{N-1} x[n] \cdot e^{-j \frac{2\pi k n}{N}}
$$

The standard Periodogram is the squared magnitude of the DFT, normalized by the signal length and sampling frequency:

$$
P_p(f_k) = \frac{1}{F_S \cdot N} |X[k]|^2
$$

-   **Windowed Periodogram**: The standard periodogram implicitly uses a rectangular window, which causes significant spectral leakage. To mitigate this, the signal is first multiplied by a window function $w[n]$ (e.g., Hamming). The resulting estimate must be normalized by the window's power, given by $U$:

    $$
    U = \frac{1}{N} \sum_{n=0}^{N-1} w[n]^2
    $$

    The windowed periodogram is then:

    $$
    P_{wp}(f_k) = \frac{1}{F_S \cdot U} \left|\text{DFT}(w[n] \cdot x[n])\right|^2
    $$

    Applying a window reduces leakage at the cost of widening the main lobe, thereby decreasing frequency resolution.

### 2. Welch's Method

Welch's method reduces the high variance of the Periodogram by averaging. The original signal is divided into $L$ possibly overlapping segments. A windowed periodogram is calculated for each segment, and the results are averaged:

$$
P_w(f_k) = \frac{1}{L} \sum_{i=0}^{L-1} P_i(f_k)
$$

where $P_i(f_k)$ is the windowed periodogram of the i-th segment. This averaging significantly smooths the spectrum (reduces variance), but because each segment is shorter than the original signal, the frequency resolution is lower.

### 3. The Multitaper Method

The Multitaper method seeks to reduce variance while minimizing the loss of frequency resolution. It computes $K$ independent PSD estimates by applying a set of $K$ specially designed orthogonal window functions (tapers), $w_k[n]$, to the signal. The optimal tapers are Discrete Prolate Spheroidal Sequences (DPSS), but this project uses **sine tapers**, which are a computationally efficient and effective alternative:

$$
w_k[n] = \sqrt{\frac{2}{N+1}} \sin\left(\frac{\pi(k+1)(n+1)}{N+1}\right) \quad \text{for } n = 0, \dots, N-1
$$

The final PSD estimate is the average of the $K$ individual periodograms (eigenspectra):

$$
P_{mt}(f) = \frac{1}{K} \sum_{k=0}^{K-1} P_k(f)
$$

This method provides excellent leakage suppression while preserving better resolution than Welch's method for a given variance reduction.

### 4. Statistical Analysis of the Welch Estimator

The PSD estimate at any given frequency is a random variable. For a signal dominated by Gaussian noise, the Welch estimate $P_w(f)$ at a frequency bin $f$ follows a **scaled chi-squared distribution**:

$$
P_w(f) \approx \frac{\sigma_n^2}{2 \cdot DF} \chi^2(2 \cdot DF)
$$

where $\sigma_n^2$ is the true noise variance and $DF$ is the degrees of freedom. This project investigates the shape of this distribution by calculating key statistical moments and performing a goodness-of-fit test.

*   **Skewness ($\gamma_1$)**: Measures the asymmetry of the distribution. A positive value indicates a tail extending to the right, which is characteristic of the $\chi^2$ distribution.
    $$
    \gamma_1 = E\left[\left(\frac{X - \mu}{\sigma}\right)^3\right]
    $$

*   **Excess Kurtosis ($\gamma_2$)**: Measures the "tailedness" or "peakedness" relative to a Gaussian distribution (which has an excess kurtosis of 0).
    $$
    \gamma_2 = E\left[\left(\frac{X - \mu}{\sigma}\right)^4\right] - 3
    $$

*   **Kolmogorov-Smirnov (K-S) Test**: This is a non-parametric test used to determine if a sample comes from a specific distribution. We use it to numerically evaluate the similarity of our PSD distribution to a Gaussian distribution. The test statistic, $D_n$, is defined as the largest absolute difference between the empirical distribution function (EDF) of the sample, $F_n(x)$, and the cumulative distribution function (CDF) of the reference distribution, $F(x)$.

    $$
    D_n = \sup_x |F_n(x) - F(x)|
    $$

    The test yields a **p-value**, which is the probability of observing a test statistic as extreme as $D_n$ under the null hypothesis ($H_0$) that the data *does* follow the reference distribution.
    *   A **high p-value** (e.g., > 0.05) suggests that the data is consistent with a Gaussian distribution.
    *   A **low p-value** (e.g., < 0.05) provides evidence to reject the null hypothesis, indicating the distribution is statistically different from a Gaussian.

### 5. Numerical Precision Analysis

This project compares three levels of floating-point precision to analyze their impact on both results and performance. The precision of a floating-point type is often characterized by its **machine epsilon ($\epsilon$)**, which is the smallest number that, when added to 1, gives a result greater than 1.

*   **Standard Precision (`double`)**: 64-bit IEEE 754 format.
    *   $\epsilon \approx 2.22 \times 10^{-16}$

*   **High Precision (`long double`)**: Typically implemented as 80-bit x86 extended precision.
    *   $\epsilon \approx 1.08 \times 10^{-19}$

*   **Arbitrary Precision (`mpfr_t`)**: A software-based implementation via the GNU MPFR library, not limited by hardware. For this project, we use a high precision of **4096 bits**.
    *   $\epsilon \approx 2.25 \times 10^{-1233}$

While higher precision offers a more accurate representation of real numbers, it comes at a significant computational cost. A key goal of this project is to determine if this extra cost yields a meaningful improvement in the final PSD estimates.

#### Practical Limits of Arbitrary Precision

While MPFR theoretically supports immense precision levels, practical limits are imposed by system hardware. The primary constraint after CPU time is **system RAM**. Each 4096-bit number requires 512 bytes of storage. The simulation simultaneously allocates memory for the input signal, multiple windowing arrays, and accumulator arrays, totaling thousands of high-precision variables.

For a typical system with 8 GB of RAM (with ~5 GB usable by a single process), the practical maximum precision for this specific program is estimated to be around **1 million bits**. Exceeding this limit would likely cause the system to rely on disk swap space, leading to a catastrophic drop in performance or a program crash.

Therefore, **4096 bits** was chosen as a value that is exceptionally high compared to standard hardware types but remains feasible to execute in a reasonable timeframe, effectively demonstrating the project's conclusions without requiring specialized hardware.

## Project Structure

```
psd_analysis/
├── Makefile
├── README.md
├── build.sh
|
├── src/
│   ├── main.c                # High Precision (long double)
│   ├── main_double.c         # Standard Precision (double)
│   └── main_mpfr.c           # Arbitrary Precision (MPFR)
|
├── scripts/
│   └── plotter.py
|
├── bin/                  # (Generated) Executables
├── data/                 # (Generated) Data files (.csv, .log)
└── plots/                # (Generated) Plots (.png)
```

## Prerequisites

### C / Gnuplot Dependencies
**On Debian / Ubuntu / Linux Mint:**
```bash
sudo apt-get update
sudo apt-get install build-essential libfftw3-dev libfftw3-bin gnuplot libgmp-dev libmpfr-dev
```
**On Fedora / RHEL / CentOS:**
```bash
sudo dnf groupinstall "Development Tools"
sudo dnf install fftw-devel gnuplot gmp-devel mpfr-devel
```

### Python Dependencies
```bash
pip3 install pandas matplotlib scipy
```

## Compilation and Execution

The `build.sh` script automates the entire process.

**1. Make the Script Executable (run once):**
```bash
chmod +x build.sh
```

**2. Compile, Run, and Plot:**
```bash
./build.sh run
```
> **WARNING:** The arbitrary-precision simulation (`psd_analyzer_mpfr`) is **computationally intensive** due to the 4096-bit precision. Expect it to take significantly longer to complete than the standard and high-precision versions.

**3. Clean the Project:**
```bash
./build.sh clean
```

## Interpreting the Results

### The Performance and Precision Report (Text Window)

This window provides a quantitative summary of the entire project.
*   **Mean Squared Error (MSE) vs. Ideal Signal**: This table shows the accuracy of each method at each precision level against a theoretical ideal. **Crucially, you will observe that the MSE does not improve with higher precision.** This is the project's key finding: the error is dominated by the statistical noise and the inherent limitations of the PSD methods, not by the floating-point precision.
*   **Execution Time Comparison**: This table clearly shows the performance cost of precision. The `Arbitrary Precision (4096-bit)` column will have drastically higher execution times, demonstrating the immense computational burden.
*   **Resource Usage**: This summarizes the maximum memory and total CPU time for each of the three main executables, further reinforcing the performance cost.

### Figure 1: Statistical Moments and Gaussian Similarity vs. Noise

This plot provides a deep dive into the statistical nature of the PSD estimate.
*   **Left Y-Axis (Moments)**: Shows how skewness and kurtosis evolve with noise (`sigma`). For low noise, they are near zero (Gaussian-like). As noise increases, they converge towards the positive theoretical values for a Chi-squared distribution, proving the distribution is non-Gaussian.
*   **Right Y-Axis (P-value)**: This shows the result of the K-S test. The p-value starts high (cannot reject that it's Gaussian) and **plummets below the red significance line** as `sigma` increases. This numerically pinpoints where the distribution becomes statistically different from a Gaussian.

### Figure 2: Statistical Moments Precision Comparison

This new plot directly compares the calculated skewness and kurtosis from the Standard, High, and Arbitrary precision simulations.
*   **Observation**: The lines for all three precision levels will be almost perfectly superimposed.
*   **Conclusion**: This graphically proves that even for sensitive statistical calculations, increasing numerical precision beyond `double` has a negligible impact on the final result when statistical noise is present.

### Figure 3: PSD Histograms (from Arbitrary Precision Data)

These histograms visualize the probability distribution of the PSD estimate at a single frequency, generated from the most accurate **4096-bit simulation**. They visually confirm the trend seen in Figure 1: the distribution is asymmetric and non-Gaussian, especially at higher noise levels.

### Figure 4: Full Precision Comparison Grid

This is the main visual result of the project. Each of the four subplots shows a PSD estimation method and compares the results from all three precision levels.
*   **Observation**: In every plot, the lines for Standard, High, and Arbitrary Precision will be virtually indistinguishable from one another.
*   **Conclusion**: This is the ultimate demonstration that the choice of numerical precision has no visible impact on the final PSD estimate for this signal. The `Arbitrary Precision (4096-bit)` line serves as a definitive reference, confirming that the `double` and `long double` results were already as accurate as they could be.

## Project Conclusions

1.  **Statistical Error Dominates Numerical Error**: The primary finding of this project is that for PSD estimation of noisy signals, the errors introduced by statistical variance (the randomness of noise) and the methodological limitations of the estimators (e.g., spectral leakage, resolution loss) are orders of magnitude larger than the errors introduced by floating-point precision.

2.  **Point of Diminishing Returns**: Increasing precision from `double` (64-bit) to `long double` (80-bit) and further to an extreme `mpfr_t` (4096-bit) results in a **massive increase in computational cost** (CPU time and memory) but yields **no discernible improvement in the accuracy** of the final PSD estimate, as measured by the MSE and visual comparison.

3.  **Confirmation of Theory**: The statistical analysis confirms theoretical expectations. The distribution of the Welch PSD estimator is shown to be non-Gaussian under noisy conditions, converging towards a Chi-squared distribution as predicted by theory. The K-S test provides a clear numerical threshold for when this deviation becomes statistically significant.