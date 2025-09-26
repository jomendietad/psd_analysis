# Statistical Analysis and Comparison of Power Spectral Density Estimators

## Overview

This project provides a C-based implementation for analyzing and comparing several Power Spectral Density (PSD) estimation techniques. It focuses on a signal composed of two sinusoids embedded in additive white Gaussian noise (AWGN).

The project is divided into two main parts:

1.  **Comparative Analysis of PSD Estimators**: It generates and plots a comparison between four key methods: the Periodogram (with a rectangular window), the Periodogram with a Hamming window, the Multitaper method, and Welch's method. This plot highlights the fundamental trade-offs between spectral resolution, variance, and spectral leakage.
2.  **Statistical Distribution Analysis**: It investigates the probability distribution of the Welch PSD estimate at a specific frequency. It calculates the skewness and kurtosis of this distribution across many realizations and for different noise levels (`sigma`). The goal is to demonstrate how the estimator's distribution behaves as the signal-to-noise ratio (SNR) changes.

All simulations and data generation are performed by the compiled C code, which is optimized for numerical precision using `long double` data types and the `FFTW3` library. The project then generates visualizations using both a native C/Gnuplot plotter and a more feature-rich Python/Matplotlib script, which saves the results as `.png` files for easy access.

**System Information:** This project was developed and tested on a machine running `Linux Mint`. Instructions for other operating systems are provided below.

## Theoretical Background

### 1. The Periodogram

The Periodogram is derived directly from the Discrete Fourier Transform (DFT). For a signal $x[n]$ of length $N$ with sampling frequency $F_S$, the DFT is:

$X[k] = \sum_{n=0}^{N-1} x[n] \cdot e^{-j \frac{2\pi k n}{N}}$

The standard Periodogram is the squared magnitude of the DFT, normalized by the signal length and sampling frequency:

$P_p(f_k) = \frac{1}{F_S \cdot N} |X[k]|^2$

-   **Windowed Periodogram**: The standard periodogram implicitly uses a rectangular window, which causes significant spectral leakage. To mitigate this, the signal is first multiplied by a window function $w[n]$ (e.g., Hamming). The resulting estimate must be normalized by the window's power, given by $U$:

    $U = \frac{1}{N} \sum_{n=0}^{N-1} w[n]^2$

    $P_{wp}(f_k) = \frac{1}{F_S \cdot U} \left|\text{DFT}(w[n] \cdot x[n])\right|^2$

    Applying a window reduces leakage at the cost of widening the main lobe, thereby decreasing frequency resolution.

### 2. Welch's Method

Welch's method reduces the high variance of the Periodogram by averaging. The original signal is divided into $L$ possibly overlapping segments. A windowed periodogram is calculated for each segment, and the results are averaged:

$P_w(f_k) = \frac{1}{L} \sum_{i=0}^{L-1} P_i(f_k)$

where $P_i(f_k)$ is the windowed periodogram of the i-th segment. This averaging significantly smooths the spectrum (reduces variance), but because each segment is shorter than the original signal, the frequency resolution is lower.

### 3. The Multitaper Method

The Multitaper method seeks to reduce variance while minimizing the loss of frequency resolution. It computes $K$ independent PSD estimates by applying a set of $K$ specially designed orthogonal window functions (tapers), $w_k[n]$, to the signal. The optimal tapers are Discrete Prolate Spheroidal Sequences (DPSS), but this project uses **sine tapers**, which are a computationally efficient and effective alternative:

$w_k[n] = \sqrt{\frac{2}{N}} \sin\left(\frac{\pi(k+1)n}{N}\right) \quad \text{for } n = 0, \dots, N-1$

The final PSD estimate is the average of the $K$ individual periodograms (eigenspectra):

$P_{mt}(f) = \frac{1}{K} \sum_{k=1}^{K} P_k(f)$

This method provides excellent leakage suppression while preserving better resolution than Welch's method for a given variance reduction.

### 4. Statistical Analysis of the Welch Estimator

The PSD estimate at any given frequency is a random variable. For a signal dominated by Gaussian noise, the Welch estimate $P_w(f)$ at a frequency bin $f$ follows a **scaled chi-squared distribution**:

$P_w(f) \approx \frac{\sigma_n^2}{2 \cdot DF} \chi^2(2 \cdot DF)$

where $\sigma_n^2$ is the true noise variance and $DF$ is the degrees of freedom (approximately the number of averaged segments). This project investigates the shape of this distribution by calculating two key statistical moments:

*   **Skewness ($\gamma_1$)**: Measures the asymmetry of the distribution. A positive value indicates a tail extending to the right, which is characteristic of the $\chi^2$ distribution.
    $\gamma_1 = E\left[\left(\frac{X - \mu}{\sigma}\right)^3\right]$

*   **Excess Kurtosis ($\gamma_2$)**: Measures the "tailedness" or "peakedness" relative to a Gaussian distribution (which has an excess kurtosis of 0).
    $\gamma_2 = E\left[\left(\frac{X - \mu}{\sigma}\right)^4\right] - 3$

## Project Structure

```
psd_analysis/
├── Makefile              # Compiles the C source files.
├── README.md             # Main project documentation.
├── build.sh              # Master script to compile, run, and plot.
|
├── src/                  # Folder for all C source code.
│   ├── main.c
│   └── plotter.c
|
├── scripts/              # Folder for the Python plotting script.
│   └── plotter.py
|
├── bin/                  # (Generated) Output directory for executables.
├── data/                 # (Generated) Output directory for data files (.csv, .log).
└── plots/                # (Generated) Output directory for saved plots (.png).
```

## Prerequisites

### C / Gnuplot Dependencies
**On Debian / Ubuntu / Linux Mint:**
```bash
sudo apt-get update
sudo apt-get install build-essential libfftw3-dev libfftw3-bin gnuplot
```
**On Fedora / RHEL / CentOS:**
```bash
sudo dnf groupinstall "Development Tools"
sudo dnf install fftw-devel gnuplot
```

### Python Dependencies
Install Python 3 and the following libraries using `pip`:
```bash
pip3 install pandas matplotlib scipy
```

### Instructions for Windows Users

#### Option 1 (Recommended): Using WSL
The Windows Subsystem for Linux (WSL) provides a full Linux environment directly within Windows, ensuring complete compatibility.

1.  **Install WSL**: Open **PowerShell as an Administrator** and run `wsl --install`. Restart your computer when prompted. An Ubuntu terminal will open to complete the setup (create a username and password).
2.  **Install Dependencies**: Open the Ubuntu terminal and run the installation commands for both C/Gnuplot and Python as shown for Debian/Ubuntu.
3.  **Access Project Files**: Your Windows drives are mounted under `/mnt/`. For example, to access a project in `C:\Users\YourUser\Documents\psd_analysis`, you would run:
    ```bash
    cd /mnt/c/Users/YourUser/Documents/psd_analysis
    ```
4.  **Run the Project**: Follow the `Compilation and Execution` steps below from within the WSL terminal.

#### Option 2 (Alternative): Using MSYS2 and MinGW
This method uses a toolchain to build native Windows executables without a full Linux VM. It requires using a specific Bash-like terminal.

1.  **Install MSYS2**: Download and run the installer from the [official MSYS2 website](https://www.msys2.org/). Follow the on-screen instructions.
2.  **Install Build Tools**: Open the **MSYS2 MINGW64** shell from the Start Menu. First, update the package database by running `pacman -Syu` (you may need to close the terminal and run it a second time). Then, install the C compiler, Make, and FFTW library:
    ```bash
    pacman -S mingw-w64-x86_64-toolchain mingw-w64-x86_64-make mingw-w64-x86_64-fftw
    ```
3.  **Install Python**: Download and install Python for Windows from [python.org](https://www.python.org/). **Crucially, check the "Add Python to PATH" option** during installation.
4.  **Install Gnuplot**: Download and install Gnuplot for Windows from the [official Gnuplot website](http://www.gnuplot.info/).
5.  **Run the Project**: Navigate to your project directory inside the **MSYS2 MINGW64 shell** and follow the execution steps below.

## Compilation and Execution

The included `build.sh` script automates the entire process.

### 1. Make the Script Executable
This command only needs to be run once from the project's root directory:
```bash
chmod +x build.sh
```

### 2. Compile, Run, and Plot
To compile the C code, generate the data files, and generate all plots, run the following from the root directory:
```bash
./build.sh run```
This will:
1.  Compile the C source code into the `bin/` directory.
2.  Run the `psd_analyzer` executable to generate data files in the `data/` directory.
3.  Execute the C-based plotter, which will attempt to open three **Gnuplot** windows.
4.  Execute the Python plotter, which will **save three `.png` images** to the `plots/` directory and then display them in interactive windows.

**Note for Windows Users:** On Windows 11 with WSL, graphical applications usually work out-of-the-box. On Windows 10 or with MSYS2, you may need an X server (like VcXsrv) for the Gnuplot windows to appear. However, the Python script will always save the plots to the `plots/` folder, which you can view from the Windows File Explorer.

### 3. Clean the Project
To remove all generated files and directories (`bin/`, `data/`, `plots/`), use the `clean` argument:
```bash
./build.sh clean
```

## Interpreting the Results

### Window 0: Statistical Analysis (Skewness & Kurtosis vs. Sigma)
This plot shows the evolution of skewness and excess kurtosis as a function of the noise standard deviation (`sigma`).
*   For **low sigma** (e.g., 0.1), the signal is dominant, and the estimator's distribution is nearly symmetric, with skewness and kurtosis values close to 0.
*   As **sigma increases**, the noise component begins to dominate, and the distribution transitions towards the theoretical **chi-squared distribution**.
*   The plot correctly shows that both skewness and kurtosis converge to stable, **positive** values for high sigma, demonstrating the non-Gaussian nature of the estimator in noisy conditions.

### Window 1: PSD Histograms
These subplots visually confirm the trend from the first plot.
*   At **low sigma (0.5 and 3.5)**, the distribution is clearly asymmetric and non-Gaussian.
*   At **high sigma (50.0)**, the distribution becomes much more defined and approaches a symmetric, bell-like shape, where a Gaussian fit is more appropriate.

### Window 2: Comparison of PSD Estimators
This plot compares the four PSD estimates for the two-tone signal.
1.  **Periodogram (Rect & Hamming Windows)**: Show sharp peaks (high resolution) but have a very noisy floor (high variance).
2.  **Multitaper Method**: Provides a great balance. The variance is reduced, leakage is very low, and frequency resolution is well-preserved.
3.  **Welch's Method**: Produces the smoothest estimate (lowest variance) at the cost of the lowest frequency resolution, smearing the two peaks.

This graph effectively illustrates the classic trade-offs in spectral estimation.