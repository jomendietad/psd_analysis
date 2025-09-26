import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os

# Definir directorios de entrada y salida
DATA_DIR = 'data'
PLOTS_DIR = 'plots'

# Asegurar que el directorio de plots exista
os.makedirs(PLOTS_DIR, exist_ok=True)

def plot_statistics():
    print("--- Python: Generating and saving Statistics Plot ---")
    try:
        stats_file = os.path.join(DATA_DIR, 'statistics.csv')
        stats_df = pd.read_csv(stats_file)
        
        plt.figure(figsize=(10, 6))
        plt.plot(stats_df['sigma'], stats_df['skewness'], marker='o', linestyle='-', label='Skewness')
        plt.plot(stats_df['sigma'], stats_df['kurtosis'], marker='s', linestyle='-', label='Excess Kurtosis')
        plt.axhline(0, color='black', linestyle='--', linewidth=1, label='Ideal Gaussian')
        plt.xscale('log')
        plt.xlabel('Sigma (Noise Standard Deviation)')
        plt.ylabel('Statistical Moment Value')
        plt.title('Python Plot: Skewness and Kurtosis vs. Noise Deviation')
        plt.grid(True, which="both", ls="--")
        plt.legend()
        
        # Guardar la figura en un archivo .png
        output_path = os.path.join(PLOTS_DIR, 'stats_py.png')
        plt.savefig(output_path)
        print(f"Plot saved to '{output_path}'")
        
    except FileNotFoundError:
        print(f"Error: '{stats_file}' not found.")

def plot_histograms():
    print("--- Python: Generating and saving Histograms Plot ---")
    def gaussian(x, mu, sigma, A):
        return A * np.exp(- (x - mu)**2 / (2 * sigma**2))

    try:
        hist_file = os.path.join(DATA_DIR, 'histograms.csv')
        hist_df = pd.read_csv(hist_file)
        sigmas_to_plot = [0.5, 3.5, 50.0]
        fig, axes = plt.subplots(3, 1, figsize=(10, 15))
        fig.suptitle('Python Plot: Distribution of PSD at f0', fontsize=16)

        for i, sigma_val in enumerate(sigmas_to_plot):
            ax = axes[i]
            subset = hist_df[hist_df['sigma'] == sigma_val].copy()
            if subset.empty:
                print(f"Warning: No data found for sigma={sigma_val}")
                continue
                
            subset['width'] = subset['bin_center'].diff().fillna(method='bfill')
            ax.bar(subset['bin_center'], subset['count'], width=subset['width']*0.9, 
                   alpha=0.7, label='Histogram', color='darkviolet')
            
            # Bloque de ajuste robusto
            params = None
            try:
                # Estimaciones iniciales
                x_data = subset['bin_center']
                y_data = subset['count']
                mu0 = np.average(x_data, weights=y_data)
                sigma0 = np.sqrt(np.average((x_data - mu0)**2, weights=y_data))
                A0 = y_data.max()
                p0 = [mu0, sigma0, A0]
                
                params, _ = curve_fit(gaussian, x_data, y_data, p0=p0, maxfev=5000)
            
            except RuntimeError as e:
                print(f"Gaussian fit failed for sigma={sigma_val}: {e}")

            if params is not None:
                x_fit = np.linspace(x_data.min(), x_data.max(), 200)
                y_fit = gaussian(x_fit, *params)
                ax.plot(x_fit, y_fit, color='red', linewidth=2, label='Gaussian Fit')

            ax.set_title(f'Sigma = {sigma_val:.2f}')
            ax.set_xlabel('PSD [Units/Hz]')
            ax.set_ylabel('Occurrence Count')
            ax.legend()
            ax.grid(alpha=0.3)
        
        plt.tight_layout(rect=[0, 0, 1, 0.96])

        # Guardar la figura en un archivo .png
        output_path = os.path.join(PLOTS_DIR, 'histograms_py.png')
        plt.savefig(output_path)
        print(f"Plot saved to '{output_path}'")

    except FileNotFoundError:
        print(f"Error: '{hist_file}' not found.")

def plot_psd_comparison():
    print("--- Python: Generating and saving PSD Estimators Comparison Plot ---")
    try:
        psd_file = os.path.join(DATA_DIR, 'results.csv')
        psd_df = pd.read_csv(psd_file)
        
        FS = 1000.0
        T = 1024
        NSEG = 256
        
        NFFT_RES = T
        NFFT_WE = NSEG
        
        freq_res = np.fft.rfftfreq(NFFT_RES, 1/FS)
        freq_we = np.fft.rfftfreq(NFFT_WE, 1/FS)

        n_pts_res = NFFT_RES // 2 + 1
        n_pts_we = NFFT_WE // 2 + 1
        
        col_periodogram = psd_df.columns[1]
        col_periodogram_hamm = psd_df.columns[2]
        col_multitaper = psd_df.columns[3]
        col_welch = psd_df.columns[5]
        
        periodogram_db = 10 * np.log10(psd_df[col_periodogram].dropna().values[:n_pts_res])
        periodogram_hamm_db = 10 * np.log10(psd_df[col_periodogram_hamm].dropna().values[:n_pts_res])
        multitaper_db = 10 * np.log10(psd_df[col_multitaper].dropna().values[:n_pts_res])
        welch_db = 10 * np.log10(psd_df[col_welch].dropna().values[:n_pts_we])

        plt.figure(figsize=(14, 8))
        plt.plot(freq_res, periodogram_db, label='Periodogram (Rect Window)', lw=1.0, alpha=0.8, color='skyblue')
        plt.plot(freq_res, periodogram_hamm_db, label='Periodogram (Hamming)', lw=1.5, alpha=0.9, color='chocolate')
        plt.plot(freq_res, multitaper_db, label='Multitaper Method', lw=2.0, color='darkviolet')
        plt.plot(freq_we, welch_db, label='Welch Method (Hamming)', lw=2.5, color='red')

        plt.ylim(-60, 20)
        plt.xlim(0, FS / 2)
        plt.xlabel('Frequency [Hz]')
        plt.ylabel('PSD [dB/Hz]')
        plt.title('Python Plot: Comparison of PSD Estimators')
        plt.legend(loc='upper right')
        plt.grid(True, linestyle='--', alpha=0.6)
        
        # Guardar la figura en un archivo .png
        output_path = os.path.join(PLOTS_DIR, 'psd_comp_py.png')
        plt.savefig(output_path)
        print(f"Plot saved to '{output_path}'")
        
    except FileNotFoundError:
        print(f"Error: '{psd_file}' not found.")
    except Exception as e:
        print(f"An unexpected error occurred while processing 'results.csv': {e}")

if __name__ == "__main__":
    plot_statistics()
    plot_histograms()
    plot_psd_comparison()
    print("\n--- Python plots generated and saved. Showing all plots now. ---")
    plt.show()