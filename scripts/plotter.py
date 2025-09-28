import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import stats
import os
import tkinter as tk
from tkinter import scrolledtext

# --- CONFIGURATION ---
DATA_DIR = 'data'
PLOTS_DIR = 'plots'
FS = 1000.0
T = 1024
M_COMP = 100
SIGMA_COMP = 3.5
NSEG = 256
OVERLAP = 128
HOP = NSEG - OVERLAP
SEGMENTS_PER_REALIZATION = (T - NSEG) // HOP + 1
DF_THEORETICAL = 2 * SEGMENTS_PER_REALIZATION
MPFR_BITS = 4096 # Define the bit precision for labels

os.makedirs(PLOTS_DIR, exist_ok=True)

# --- REPORTING FUNCTIONS ---
def calculate_and_write_reports():
    """
    Calculates Mean Squared Error (MSE) for each precision level against a 
    theoretical ideal signal and writes the results to a log file.
    """
    print("--- Python: Calculating and writing performance and precision reports ---")
    try:
        df_d = pd.read_csv(os.path.join(DATA_DIR, 'results_double.csv'))
        df_ld = pd.read_csv(os.path.join(DATA_DIR, 'results_long_double.csv'))
        try:
            df_mpfr = pd.read_csv(os.path.join(DATA_DIR, 'results_mpfr.csv'))
        except FileNotFoundError:
            df_mpfr = None
            print("Warning: MPFR results file not found. MSE report will be incomplete.")

        A1, A2, F1, F2 = 1.0, 0.7, 120.0, 250.0
        len_res, len_we = len(np.fft.rfftfreq(T, 1/FS)), len(np.fft.rfftfreq(NSEG, 1/FS))
        
        freq_res = np.fft.rfftfreq(T, 1/FS)
        ideal_psd_res = np.full(len(freq_res), SIGMA_COMP**2 / FS)
        k1_res, k2_res = np.argmin(np.abs(freq_res - F1)), np.argmin(np.abs(freq_res - F2))
        ideal_psd_res[k1_res] += (A1**2 / 2)
        ideal_psd_res[k2_res] += (A2**2 / 2)

        freq_we = np.fft.rfftfreq(NSEG, 1/FS)
        ideal_psd_we = np.full(len(freq_we), SIGMA_COMP**2 / FS)
        k1_we, k2_we = np.argmin(np.abs(freq_we - F1)), np.argmin(np.abs(freq_we - F2))
        ideal_psd_we[k1_we] += (A1**2 / 2)
        ideal_psd_we[k2_we] += (A2**2 / 2)

        mse_vals = {}
        for name, col, ideal, length in [
            ('Periodogram (Rect)', 'PSD_periodogram', ideal_psd_res, len_res),
            ('Periodogram (Hamming)', 'PSD_periodogram_hamm', ideal_psd_res, len_res),
            ('Multitaper', 'PSD_mt', ideal_psd_res, len_res),
            ('Welch', 'PSD_welch', ideal_psd_we, len_we)
        ]:
            mse_vals[name] = {
                'd': np.mean((df_d[col].iloc[:length] - ideal)**2),
                'ld': np.mean((df_ld[col].iloc[:length] - ideal)**2),
                'mpfr': np.mean((df_mpfr[col].iloc[:length].astype(float) - ideal)**2) if df_mpfr is not None else 'N/A'
            }

        with open(os.path.join(DATA_DIR, 'precision_metrics.log'), 'w') as f:
            f.write("--- Precision Error Metrics ---\n\n")
            f.write("--- Mean Squared Error (MSE) vs. Ideal Signal ---\n")
            f.write(f"{'Method':<25} {'MSE (Standard Precision)':>28} {'MSE (High Precision)':>25} {'MSE (Arbitrary - ' + str(MPFR_BITS) + '-bit)':>30}\n")
            f.write("-" * 112 + "\n")
            for name, mses in mse_vals.items():
                mpfr_str = f"{mses['mpfr']:30.4e}" if isinstance(mses['mpfr'], float) else mses['mpfr'].rjust(30)
                f.write(f"{name:<25} {mses['d']:28.4e} {mses['ld']:25.4e} {mpfr_str}\n")

    except Exception as e:
        print(f"CRITICAL: Could not calculate precision metrics: {e}")

def show_text_reports():
    """
    Reads all timing and resource logs and displays them in a formatted text window.
    """
    try:
        with open(os.path.join(DATA_DIR, 'precision_metrics.log'), 'r') as f:
            precision_report = f.read()
        
        df_d_time = pd.read_csv(os.path.join(DATA_DIR, 'timing_log_d.csv')).set_index('method')
        df_ld_time = pd.read_csv(os.path.join(DATA_DIR, 'timing_log_ld.csv')).set_index('method')
        try:
            df_mpfr_time = pd.read_csv(os.path.join(DATA_DIR, 'timing_log_mpfr.csv')).set_index('method')
        except FileNotFoundError:
            df_mpfr_time = pd.DataFrame()

        timing_report = "--- Execution Time Comparison (avg ms per realization) ---\n"
        timing_report += f"{'Method':<25} {'Standard Precision':>22} {'High Precision':>20} {'Arbitrary Precision (' + str(MPFR_BITS) + '-bit)':>28}\n"
        timing_report += "-"*100 + "\n"
        
        all_methods = sorted(list(set(df_d_time.index) | set(df_ld_time.index) | set(df_mpfr_time.index)))

        for method in all_methods:
            time_d = df_d_time.loc[method, 'avg_cpu_time_per_realization_ms'] if method in df_d_time.index else 'N/A'
            time_ld = df_ld_time.loc[method, 'avg_cpu_time_per_realization_ms'] if method in df_ld_time.index else 'N/A'
            time_mpfr = df_mpfr_time.loc[method, 'avg_cpu_time_per_realization_ms'] if method in df_mpfr_time.index else 'N/A'
            
            time_d_str = f"{time_d:.4f}" if isinstance(time_d, float) else time_d
            time_ld_str = f"{time_ld:.4f}" if isinstance(time_ld, float) else time_ld
            time_mpfr_str = f"{time_mpfr:.4f}" if isinstance(time_mpfr, float) else time_mpfr

            timing_report += f"{method:<25} {time_d_str:>22} {time_ld_str:>20} {time_mpfr_str:>28}\n"

        resource_report = "\n--- Resource Usage (Summary) ---\n"
        for precision, label in [('d', 'Standard Precision'), ('ld', 'High Precision'), ('mpfr', f'Arbitrary Precision ({MPFR_BITS}-bit)')]:
            try:
                with open(os.path.join(DATA_DIR, f'resource_usage_{precision}.txt'), 'r') as f:
                    lines = f.readlines()
                    mem = [l for l in lines if 'Maximum resident set size' in l][0].split(':')[-1].strip()
                    cpu = [l for l in lines if 'User time' in l][0].split(':')[-1].strip()
                    resource_report += f"-> {label}:\n   Max Memory: {mem} KB\n   User CPU Time: {cpu} s\n"
            except (FileNotFoundError, IndexError):
                 resource_report += f"-> {label}: Resource file not found or incomplete.\n"

        instruction = "\n\n--- NOTE: Close this window to show interactive plots. ---"
        root = tk.Tk()
        root.title("Performance and Precision Report")
        st = scrolledtext.ScrolledText(root, wrap=tk.WORD, font=("Monospace", 10), width=115, height=35)
        st.pack(padx=10, pady=10, fill=tk.BOTH, expand=True)
        full_report = f"{precision_report}\n{timing_report}\n{resource_report}\n{instruction}"
        st.insert(tk.INSERT, full_report)
        st.configure(state='disabled')
        root.mainloop()
    except Exception as e:
        print(f"CRITICAL: Could not generate report window: {e}")

# --- PLOTTING FUNCTIONS ---
def plot_statistics_gaussian_similarity():
    """
    Plots skewness, kurtosis, and the K-S test p-value against noise sigma
    using the high-precision data as the reference.
    """
    print("--- Python: Generating Gaussian Similarity Statistics Plot ---")
    try:
        stats_df = pd.read_csv(os.path.join(DATA_DIR, 'statistics.csv'))
        realizations_df = pd.read_csv(os.path.join(DATA_DIR, 'psd_realizations.csv'))
        
        p_values = []
        for sigma_val in stats_df['sigma']:
            data = realizations_df[realizations_df['sigma'] == sigma_val]['psd_value']
            if len(data) > 20:
                z_scores = (data - data.mean()) / data.std()
                _, p_value = stats.kstest(z_scores, 'norm')
                p_values.append(p_value)
            else:
                p_values.append(np.nan)
        
        stats_df['p_value'] = p_values
        th_skew, th_kurt = np.sqrt(8 / DF_THEORETICAL), 12 / DF_THEORETICAL
        
        fig, ax1 = plt.subplots(figsize=(12, 7))
        fig.canvas.manager.set_window_title('Figure 1: Gaussian Similarity Analysis')

        ax1.set_xlabel('Sigma (Noise Standard Deviation)')
        ax1.set_ylabel('Statistical Moment Value', color='black')
        ax1.plot(stats_df['sigma'], stats_df['skewness'], 'o-', color='C0', label='Simulated Skewness')
        ax1.plot(stats_df['sigma'], stats_df['kurtosis'], 's-', color='C1', label='Simulated Excess Kurtosis')
        ax1.axhline(0, color='k', ls='--', lw=1, label='Ideal Gaussian (Skew=0, Kurt=0)')
        ax1.axhline(th_skew, color='blue', ls=':', lw=2, label=f'Theoretical Chi-Sq Skewness ({th_skew:.2f})')
        ax1.axhline(th_kurt, color='orangered', ls=':', lw=2, label=f'Theoretical Chi-Sq Kurtosis ({th_kurt:.2f})')
        ax1.set_xscale('log'); ax1.tick_params(axis='y', labelcolor='black'); ax1.grid(True, which="both", ls="--")

        ax2 = ax1.twinx()
        ax2.set_ylabel('K-S Test P-value (vs. Normal)', color='darkgreen')
        ax2.plot(stats_df['sigma'], stats_df['p_value'], 'd-', color='darkgreen', label='P-value (Gaussian Similarity)')
        ax2.axhline(0.05, color='red', ls=':', lw=2, label='Significance Level (α=0.05)')
        ax2.set_yscale('log'); ax2.tick_params(axis='y', labelcolor='darkgreen')
        
        fig.suptitle('Python Plot: Statistical Moments and Gaussian Similarity vs. Noise')
        fig.legend(loc='upper left', bbox_to_anchor=(0.1, 0.9)); fig.tight_layout(rect=[0, 0, 1, 0.96])
        plt.savefig(os.path.join(PLOTS_DIR, 'stats_gaussian_similarity.png'))
        print(f"Plot saved to 'plots/stats_gaussian_similarity.png'")
    except Exception as e: print(f"Error plotting statistics: {e}")

def plot_statistics_precision_comparison():
    """
    Generates a new plot comparing the calculated skewness and kurtosis
    across all three precision levels.
    """
    print("--- Python: Generating Statistical Moments Precision Comparison Plot ---")
    try:
        df_d = pd.read_csv(os.path.join(DATA_DIR, 'statistics_double.csv'))
        df_ld = pd.read_csv(os.path.join(DATA_DIR, 'statistics.csv'))
        df_mpfr = pd.read_csv(os.path.join(DATA_DIR, 'statistics_mpfr.csv'))

        fig, axs = plt.subplots(2, 1, figsize=(12, 10), sharex=True, constrained_layout=True)
        fig.canvas.manager.set_window_title('Figure 2: Statistical Moments Precision Comparison')
        fig.suptitle('Comparison of Calculated Statistical Moments Across Precisions')

        # Skewness Subplot
        axs[0].plot(df_d['sigma'], df_d['skewness'], 'b--o', label='Standard Precision')
        axs[0].plot(df_ld['sigma'], df_ld['skewness'], 'r-s', label='High Precision')
        axs[0].plot(df_mpfr['sigma'], df_mpfr['skewness'], 'g:d', label=f'Arbitrary Precision ({MPFR_BITS}-bit)')
        axs[0].set_ylabel('Skewness')
        axs[0].set_title('Skewness vs. Noise')
        axs[0].grid(True, which="both", ls="--"); axs[0].legend()

        # Kurtosis Subplot
        axs[1].plot(df_d['sigma'], df_d['kurtosis'], 'b--o', label='Standard Precision')
        axs[1].plot(df_ld['sigma'], df_ld['kurtosis'], 'r-s', label='High Precision')
        axs[1].plot(df_mpfr['sigma'], df_mpfr['kurtosis'], 'g:d', label=f'Arbitrary Precision ({MPFR_BITS}-bit)')
        axs[1].set_xlabel('Sigma (Noise Standard Deviation)')
        axs[1].set_ylabel('Excess Kurtosis')
        axs[1].set_title('Excess Kurtosis vs. Noise')
        axs[1].grid(True, which="both", ls="--"); axs[1].legend()
        
        axs[1].set_xscale('log')
        plt.savefig(os.path.join(PLOTS_DIR, 'stats_precision_comparison.png'))
        print(f"Plot saved to 'plots/stats_precision_comparison.png'")
    except Exception as e:
        print(f"Error plotting statistical precision comparison: {e}")

def plot_histograms():
    """
    Plots the distribution of the PSD estimate for different noise levels,
    using data from the arbitrary-precision simulation.
    """
    print("--- Python: Generating and saving Histograms Plot (from Arbitrary Precision data) ---")
    def gaussian(x, mu, sigma, A): return A*np.exp(-(x-mu)**2/(2*sigma**2))
    try:
        hist_df = pd.read_csv(os.path.join(DATA_DIR, 'histograms_mpfr.csv'))
        fig = plt.figure("Figure 3: PSD Histograms", figsize=(10, 15))
        subfigs = fig.subfigures(3, 1)
        fig.suptitle(f'Distribution of PSD at f0 (Arbitrary Precision, {MPFR_BITS}-bit)', fontsize=16)
        for i, sigma_val in enumerate([0.5, 3.5, 50.0]):
            ax = subfigs[i].subplots()
            subset = hist_df[hist_df['sigma'] == sigma_val].copy()
            if subset.empty: continue
            ax.bar(subset['bin_center'], subset['count'], width=subset['bin_center'].diff().median(), 
                   alpha=0.7, label='Histogram', color='darkviolet')
            try:
                x_data, y_data = subset['bin_center'], subset['count']
                p0 = [np.average(x_data,weights=y_data), np.sqrt(np.average((x_data-np.average(x_data,weights=y_data))**2,weights=y_data)), y_data.max()]
                params, _ = curve_fit(gaussian, x_data, y_data, p0=p0, maxfev=5000)
                x_fit = np.linspace(x_data.min(), x_data.max(), 200)
                ax.plot(x_fit, gaussian(x_fit, *params), 'r-', lw=2, label='Gaussian Fit')
            except Exception as e: print(f"Fit failed for sigma={sigma_val}: {e}")
            ax.set_title(f'Sigma = {sigma_val:.2f}'); ax.set_xlabel('PSD [Units/Hz]'); ax.set_ylabel('Occurrence Count'); ax.legend(); ax.grid(alpha=0.3)
        plt.savefig(os.path.join(PLOTS_DIR, 'histograms_py.png'))
        print(f"Plot saved to 'plots/histograms_py.png'")
    except Exception as e: print(f"Error plotting histograms: {e}")

def plot_precision_comparison_grid():
    """
    Generates a 2x2 grid comparing all precision levels for each PSD method.
    This is the main comparison plot.
    """
    print("--- Python: Generating Full Precision Comparison Grid Plot ---")
    try:
        df_d = pd.read_csv(os.path.join(DATA_DIR, 'results_double.csv'))
        df_ld = pd.read_csv(os.path.join(DATA_DIR, 'results_long_double.csv'))
        try:
            df_mpfr = pd.read_csv(os.path.join(DATA_DIR, 'results_mpfr.csv'))
        except FileNotFoundError:
            df_mpfr = None
            print("Warning: MPFR results not found, will be skipped in comparison grid.")

        freq_res, freq_we = np.fft.rfftfreq(T, 1/FS), np.fft.rfftfreq(NSEG, 1/FS)
        noise_floor_db = 10 * np.log10((SIGMA_COMP**2) / FS)

        fig, axs = plt.subplots(2, 2, figsize=(20, 14), sharex=True, sharey=True, constrained_layout=True)
        fig.canvas.manager.set_window_title('Figure 4: Full Precision Comparison')
        fig.suptitle('Precision Comparison of PSD Estimators vs. Ideal Reference', fontsize=18)
        
        methods = {
            'Periodogram (Rect)': ('PSD_periodogram', freq_res),
            'Periodogram (Hamming)': ('PSD_periodogram_hamm', freq_res),
            'Multitaper': ('PSD_mt', freq_res),
            'Welch': ('PSD_welch', freq_we)
        }
        
        for i, (name, (col, freq)) in enumerate(methods.items()):
            ax = axs.flatten()[i]
            len_f = len(freq)
            
            psd_d = 10*np.log10(df_d[col].iloc[:len_f])
            psd_ld = 10*np.log10(df_ld[col].iloc[:len_f])
            
            ax.plot(freq, psd_d, 'b--', lw=1, label='Standard Precision')
            ax.plot(freq, psd_ld, 'r-', lw=1.5, label='High Precision')
            
            if df_mpfr is not None and col in df_mpfr.columns:
                psd_mpfr = 10*np.log10(df_mpfr[col].iloc[:len_f].astype(float))
                ax.plot(freq, psd_mpfr, 'g:', lw=2.5, label=f'Arbitrary Precision ({MPFR_BITS}-bit)')

            ax.axhline(noise_floor_db, color='grey', ls=':', lw=2, label=f'Ideal Noise Floor ({noise_floor_db:.2f} dB/Hz)')
            ax.axvline(120.0, color='k', ls='-.', lw=1.5, label='Tone Frequencies')
            ax.axvline(250.0, color='k', ls='-.', lw=1.5)
            ax.set(title=name, xlabel='Frequency [Hz]', ylabel='PSD [dB/Hz]', xlim=(0, FS/2), ylim=(-60, 20))
            ax.grid(True, ls='--', alpha=0.6)
            ax.legend()

        plt.savefig(os.path.join(PLOTS_DIR, 'precision_comparison_grid.png'))
        print(f"Plot saved to 'plots/precision_comparison_grid.png'")
    except Exception as e:
        print(f"Error generating precision comparison grid plot: {e}")

# --- Main Execution ---
if __name__ == "__main__":
    calculate_and_write_reports()
    show_text_reports()
    
    plt.close('all')

    # Generate all plots
    plot_statistics_gaussian_similarity()
    plot_statistics_precision_comparison()
    plot_histograms()
    plot_precision_comparison_grid()
    
    print("\n--- All Python plots generated. Now showing interactive plots. ---")
    plt.show()