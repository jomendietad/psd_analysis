// Se define _GNU_SOURCE para asegurar la disponibilidad de M_PIL y hypotl
#define _GNU_SOURCE
#include <math.h>

#include <fftw3.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

// ============================================================================
// ESTRUCTURAS Y FUNCIONES DE ANÁLISIS ESTADÍSTICO
// ============================================================================

typedef struct {
    long double mean;
    long double variance;
    long double skewness;
    long double kurtosis;
} Stats;

Stats calculate_statistics(const long double* data, int n) {
    Stats result = {0.0L, 0.0L, 0.0L, 0.0L};
    if (n < 4) return result;
    long double sum = 0.0L;
    for (int i = 0; i < n; ++i) sum += data[i];
    result.mean = sum / n;
    long double sum_sq_diff = 0.0L, sum_cub_diff = 0.0L, sum_quar_diff = 0.0L;
    for (int i = 0; i < n; ++i) {
        long double diff = data[i] - result.mean;
        sum_sq_diff += diff * diff;
        sum_cub_diff += diff * diff * diff;
        sum_quar_diff += diff * diff * diff * diff;
    }
    result.variance = sum_sq_diff / n;
    if (result.variance > 1e-30L) {
        long double std_dev = sqrtl(result.variance);
        result.skewness = (sum_cub_diff / n) / (std_dev * std_dev * std_dev);
        result.kurtosis = (sum_quar_diff / n) / (result.variance * result.variance) - 3.0L;
    }
    return result;
}

void generate_histogram(FILE* file, const long double* data, int n, long double sigma_val, int num_bins) {
    if (n == 0) return;
    long double min_val = data[0], max_val = data[0];
    for (int i = 1; i < n; ++i) {
        if (data[i] < min_val) min_val = data[i];
        if (data[i] > max_val) max_val = data[i];
    }
    if (fabsl(max_val - min_val) < 1e-30L) {
        fprintf(file, "%.8Lg,%.17Lg,%d\n", sigma_val, min_val, n);
        return;
    }
    long double bin_width = (max_val - min_val) / num_bins;
    int* bins = (int*) calloc(num_bins, sizeof(int));
    for (int i = 0; i < n; ++i) {
        int bin_index = (int)((data[i] - min_val) / bin_width);
        if (bin_index >= num_bins) bin_index = num_bins - 1;
        if (bin_index < 0) bin_index = 0;
        bins[bin_index]++;
    }
    for (int i = 0; i < num_bins; ++i) {
        long double bin_center = min_val + (i + 0.5L) * bin_width;
        fprintf(file, "%.8Lg,%.17Lg,%d\n", sigma_val, bin_center, bins[i]);
    }
    free(bins);
}

// ============================================================================
// FUNCIONES AUXILIARES DE PROCESAMIENTO DE SEÑALES
// ============================================================================

static inline long double randu() { return (rand() + 1.0L) / (RAND_MAX + 2.0L); }
static long double randn_boxmuller() {
    long double u1 = randu(), u2 = randu();
    long double r = sqrtl(-2.0L * logl(u1 + 1e-18L));
    long double th = 2.0L * M_PIl * u2;
    return r * cosl(th);
}
static void make_hamming(long double* w, int N) { for (int n=0; n<N; ++n) w[n] = 0.54L - 0.46L * cosl(2.0L * M_PIl * n / (N - 1)); }
static long double window_U(const long double* w, int N) {
    long double s2 = 0.0L;
    for (int n=0; n<N; ++n) s2 += w[n] * w[n];
    return s2 / (long double)N;
}

// ============================================================================
// IMPLEMENTACIÓN DEL MÉTODO MULTITAPER
// ============================================================================

void generate_sine_tapers(int N, int K, long double** tapers) {
    for (int k = 0; k < K; k++) {
        long double norm_factor = sqrtl(2.0L / (long double)(N));
        for (int i = 0; i < N; i++) {
            tapers[k][i] = norm_factor * sinl(M_PIl * (long double)(k + 1) * (long double)(i) / (long double)(N));
        }
    }
}

// ============================================================================
// FUNCIÓN PRINCIPAL DE COMPARACIÓN DE PSD (CORREGIDA Y OPTIMIZADA)
// ============================================================================
void generate_psd_comparison_data(
    long double fs, int T, int M, long double sigma, 
    long double A1, long double f1, long double A2, long double f2, 
    int Nseg, int overlap, long double NW_mt) {
    
    fprintf(stderr, "Generating 'data/results.csv' for method comparison...\n");
    fprintf(stderr, "Signal: f1=%.1Lf Hz (A=%.1Lf), f2=%.1Lf Hz (A=%.1Lf), sigma=%.2Lf\n", f1, A1, f2, A2, sigma);
    
    const int Nfft = T;
    const int Nfft_we = Nseg;
    const int K_mt = (int)(2 * NW_mt) - 1;

    long double* wT_hamm = (long double*) malloc(sizeof(long double) * T);
    make_hamming(wT_hamm, T);
    long double U_hamm = window_U(wT_hamm, T);
    
    long double** dpss_tapers = (long double**) malloc(K_mt * sizeof(long double*));
    for(int i = 0; i < K_mt; ++i) dpss_tapers[i] = (long double*) malloc(T * sizeof(long double));
    generate_sine_tapers(T, K_mt, dpss_tapers);
    fprintf(stderr, "Multitaper method configured with NW=%.1Lf, K=%d tapers.\n", NW_mt, K_mt);

    fftwl_complex* in_fft = (fftwl_complex*) fftwl_malloc(sizeof(fftwl_complex) * Nfft);
    fftwl_complex* out_fft = (fftwl_complex*) fftwl_malloc(sizeof(fftwl_complex) * Nfft);
    fftwl_plan plan_fft = fftwl_plan_dft_1d(Nfft, in_fft, out_fft, FFTW_FORWARD, FFTW_PATIENT);

    long double* Periodogram_acc = (long double*) calloc(Nfft, sizeof(long double));
    long double* Periodogram_hamm_acc = (long double*) calloc(Nfft, sizeof(long double));
    long double* MT_acc = (long double*) calloc(Nfft, sizeof(long double));
    long double* Welch_acc = (long double*) calloc(Nfft_we, sizeof(long double));
    long double* x = (long double*) malloc(sizeof(long double) * T);

    long double* wS_welch = (long double*)malloc(sizeof(long double)*Nseg);
    make_hamming(wS_welch, Nseg);
    long double Useg_welch = window_U(wS_welch, Nseg);
    fftwl_complex* in_we = (fftwl_complex*)fftwl_malloc(sizeof(fftwl_complex)*Nfft_we);
    fftwl_complex* out_we = (fftwl_complex*)fftwl_malloc(sizeof(fftwl_complex)*Nfft_we);
    fftwl_plan plan_we = fftwl_plan_dft_1d(Nfft_we, in_we, out_we, FFTW_FORWARD, FFTW_ESTIMATE);
    int hop = Nseg - overlap;
    int S_per_real = (T - Nseg) / hop + 1;

    for (int m = 0; m < M; ++m) {
        for (int n = 0; n < T; ++n) {
            x[n] = A1 * cosl(2.0L * M_PIl * f1 * n / fs) + 
                   A2 * cosl(2.0L * M_PIl * f2 * n / fs) + 
                   sigma * randn_boxmuller();
        }

        for(int i=0; i<T; ++i) { in_fft[i][0] = x[i]; in_fft[i][1] = 0.0L; }
        fftwl_execute(plan_fft);
        for(int i=0; i<Nfft; ++i) Periodogram_acc[i] += (out_fft[i][0]*out_fft[i][0] + out_fft[i][1]*out_fft[i][1]);

        for(int i=0; i<T; ++i) { in_fft[i][0] = x[i] * wT_hamm[i]; in_fft[i][1] = 0.0L; }
        fftwl_execute(plan_fft);
        for(int i=0; i<Nfft; ++i) Periodogram_hamm_acc[i] += (out_fft[i][0]*out_fft[i][0] + out_fft[i][1]*out_fft[i][1]);

        long double* k_spectrum = calloc(Nfft, sizeof(long double));
        for (int k = 0; k < K_mt; ++k) {
            for(int i=0; i<T; ++i) { in_fft[i][0] = x[i] * dpss_tapers[k][i]; in_fft[i][1] = 0.0L; }
            fftwl_execute(plan_fft);
            for(int i=0; i<Nfft; ++i) k_spectrum[i] += (out_fft[i][0]*out_fft[i][0] + out_fft[i][1]*out_fft[i][1]);
        }
        for(int i=0; i<Nfft; ++i) MT_acc[i] += k_spectrum[i] / K_mt;
        free(k_spectrum);
        
        for (int sidx=0; sidx<S_per_real; ++sidx){
            int start=sidx*hop;
            for (int n=0; n<Nseg; ++n){ in_we[n][0]=x[start+n]*wS_welch[n]; in_we[n][1]=0.0L; }
            fftwl_execute(plan_we);
            for (int k=0; k<Nfft_we; ++k) Welch_acc[k] += (out_we[k][0]*out_we[k][0]+out_we[k][1]*out_we[k][1]);
        }
    }
    
    for (int k=0; k<Nfft; ++k) {
        Periodogram_acc[k] /= (fs * T * M);
        Periodogram_hamm_acc[k] /= (fs * T * U_hamm * M);
        MT_acc[k] /= (fs * T * M);
    }
    for (int k=0; k<Nfft_we; ++k) {
        Welch_acc[k] /= (fs * Nseg * Useg_welch * S_per_real * M);
    }
    
    FILE* res_file = fopen("data/results.csv", "w");
    fprintf(res_file, "k_res,PSD_periodogram,PSD_periodogram_hamm,PSD_mt,k_we,PSD_welch\n");
    for (int i=0; i < Nfft; ++i){
        fprintf(res_file, "%d,%.21Lg,%.21Lg,%.21Lg,%d,%.21Lg\n",
                i, Periodogram_acc[i], Periodogram_hamm_acc[i], MT_acc[i],
                (i<Nfft_we)?i:-1, (i<Nfft_we)?Welch_acc[i]:NAN);
    }
    fclose(res_file);
    fprintf(stderr, "'data/results.csv' generated successfully.\n");
    
    free(wT_hamm); free(x); free(wS_welch);
    free(Periodogram_acc); free(Periodogram_hamm_acc); free(MT_acc); free(Welch_acc);
    for(int i = 0; i < K_mt; ++i) free(dpss_tapers[i]); free(dpss_tapers);
    fftwl_destroy_plan(plan_fft); fftwl_destroy_plan(plan_we);
    fftwl_free(in_fft); fftwl_free(out_fft); fftwl_free(in_we); fftwl_free(out_we);
}

int main() {
    const long double fs = 1000.0L;
    const int T = 1024;        
    const int M_comp = 100;
    const int M_stats = 1000;
    const int Nseg = 256;      
    const int overlap = 128;   

    const long double sigma_comp = 3.5L;
    const long double f1 = 120.0L, A1 = 1.0L;    
    const long double f2 = 250.0L, A2 = 0.7L;    
    const long double NW_mt = 4.0L; 

    const long double f0_stats = 250.0L, A_stats = 1.0L;    

    srand(time(NULL));
    fftwl_init_threads();
    fftwl_plan_with_nthreads(2);

    generate_psd_comparison_data(fs, T, M_comp, sigma_comp, A1, f1, A2, f2, Nseg, overlap, NW_mt);

    fprintf(stderr, "\nStarting statistical analysis over multiple sigmas...\n");
    const long double sigmas[] = {0.1L, 0.5L, 1.0L, 3.5L, 10.0L, 25.0L, 50.0L};
    const int num_sigmas = sizeof(sigmas) / sizeof(sigmas[0]);
    const int k_f0_we = (int)(f0_stats * Nseg / fs);
    const int hop = Nseg - overlap;
    const int S_per_real = (T - Nseg) / hop + 1;

    long double* wS = (long double*)malloc(sizeof(long double) * Nseg);
    make_hamming(wS, Nseg);
    const long double Useg = window_U(wS, Nseg);
    
    const long double norm_welch_stats = fs * Nseg * Useg * S_per_real;

    FILE* stats_file = fopen("data/statistics.csv", "w");
    FILE* hist_file = fopen("data/histograms.csv", "w");
    fprintf(stats_file, "sigma,mean_psd,variance_psd,skewness,kurtosis\n");
    fprintf(hist_file, "sigma,bin_center,count\n");

    clock_t t_total_start = clock();

    fftwl_complex* in_we = (fftwl_complex*)fftwl_malloc(sizeof(fftwl_complex)*Nseg);
    fftwl_complex* out_we = (fftwl_complex*)fftwl_malloc(sizeof(fftwl_complex)*Nseg);
    fftwl_plan plan_we = fftwl_plan_dft_1d(Nseg, in_we, out_we, FFTW_FORWARD, FFTW_ESTIMATE);

    for (int s_idx = 0; s_idx < num_sigmas; ++s_idx) {
        const long double current_sigma = sigmas[s_idx];
        fprintf(stderr, "Processing sigma = %.4Lf ...\n", current_sigma);

        long double* realization_psd_values = (long double*)malloc(sizeof(long double) * M_stats);
        long double* x = (long double*)malloc(sizeof(long double) * T);

        for (int m = 0; m < M_stats; ++m) {
            for (int n = 0; n < T; ++n) x[n] = A_stats * cosl(2.0L * M_PIl * f0_stats * n / fs) + current_sigma * randn_boxmuller();
            
            long double Pxx_realization_sum = 0.0L;
            for (int sidx = 0; sidx < S_per_real; ++sidx) {
                int start = sidx * hop;
                for (int n = 0; n < Nseg; ++n){ in_we[n][0] = x[start+n]*wS[n]; in_we[n][1] = 0.0L; }
                fftwl_execute(plan_we);
                long double re = out_we[k_f0_we][0];
                long double im = out_we[k_f0_we][1];
                Pxx_realization_sum += (re*re + im*im);
            }
            realization_psd_values[m] = (2.0L * Pxx_realization_sum) / norm_welch_stats;
        }
        
        Stats stats = calculate_statistics(realization_psd_values, M_stats);
        fprintf(stats_file, "%.8Lg,%.17Lg,%.17Lg,%.17Lg,%.17Lg\n", current_sigma, stats.mean, stats.variance, stats.skewness, stats.kurtosis);
        generate_histogram(hist_file, realization_psd_values, M_stats, current_sigma, 40);
        
        free(realization_psd_values); free(x);
    }
    
    fftwl_destroy_plan(plan_we); fftwl_free(in_we); fftwl_free(out_we);
    
    clock_t t_total_end = clock();
    double elapsed_total = (double)(t_total_end-t_total_start)/CLOCKS_PER_SEC;
    fprintf(stderr, "\nStatistical analysis completed.\nTotal time: %.4f s\n", elapsed_total);

    fclose(stats_file); fclose(hist_file);
    free(wS);
    fftwl_cleanup_threads();
    return 0;
}