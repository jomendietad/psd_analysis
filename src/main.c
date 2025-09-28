#define _GNU_SOURCE
#include <math.h>
#include <locale.h>
#include <fftw3.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <pthread.h>

// A structure to hold key statistical moments
typedef struct { long double mean, variance, skewness, kurtosis; } Stats;

// Calculates the first four statistical moments for a given dataset.
Stats calculate_statistics(const long double* data, int n) {
    Stats result = {0.0L, 0.0L, 0.0L, 0.0L};
    if (n < 4) return result;

    long double sum = 0.0L;
    for (int i = 0; i < n; ++i) sum += data[i];
    result.mean = sum / n;

    long double s2 = 0.0L, s3 = 0.0L, s4 = 0.0L;
    for (int i = 0; i < n; ++i) {
        long double delta = data[i] - result.mean;
        s2 += delta * delta;
        s3 += delta * delta * delta;
        s4 += delta * delta * delta * delta;
    }
    result.variance = s2 / n;

    if (result.variance > 1e-30L) {
        long double std_dev = sqrtl(result.variance);
        result.skewness = (s3 / n) / (std_dev * std_dev * std_dev);
        result.kurtosis = (s4 / n) / (result.variance * result.variance) - 3.0L;
    }
    return result;
}

void generate_histogram(FILE* file, const long double* data, int n, long double sigma, int nbins) {
    if (n == 0) return;
    long double min_v = data[0], max_v = data[0];
    for (int i = 1; i < n; ++i) {
        if (data[i] < min_v) min_v = data[i];
        if (data[i] > max_v) max_v = data[i];
    }
    if (fabsl(max_v - min_v) < 1e-30L) {
        fprintf(file, "%.8Lg,%.17Lg,%d\n", sigma, min_v, n);
        return;
    }
    long double width = (max_v - min_v) / nbins;
    int *bins = calloc(nbins, sizeof(int));
    for (int i = 0; i < n; ++i) {
        int idx = (int)((data[i] - min_v) / width);
        if (idx >= nbins) idx = nbins - 1;
        bins[idx]++;
    }
    for (int i = 0; i < nbins; ++i) {
        fprintf(file, "%.8Lg,%.17Lg,%d\n", sigma, min_v + (i + 0.5L) * width, bins[i]);
    }
    free(bins);
}

static inline long double rand_uniform() { return (rand() + 1.0L) / (RAND_MAX + 2.0L); }
static long double randn_box_muller() { return sqrtl(-2.0L * logl(rand_uniform() + 1e-18L)) * cosl(2.0L * M_PIl * rand_uniform()); }
static void make_hamming(long double *w, int N) { for (int n = 0; n < N; ++n) w[n] = 0.54L - 0.46L * cosl(2 * M_PIl * n / (N - 1)); }
static long double window_power(const long double *w, int N) { long double s2 = 0; for (int n = 0; n < N; ++n) s2 += w[n] * w[n]; return s2 / N; }
void generate_sine_tapers(int N, int K, long double **tapers) { for (int k = 0; k < K; k++) { long double norm = sqrtl(2.0L / (N + 1)); for (int i = 0; i < N; i++) tapers[k][i] = norm * sinl(M_PIl * (k + 1) * (i + 1) / (N + 1)); } }

void generate_psd_data(long double fs, int T, int M, long double sigma, long double A1, long double f1, long double A2, long double f2, int Nseg, int overlap, long double NW_mt) {
    struct timespec start, end;
    double time_p = 0, time_ph = 0, time_m = 0, time_w = 0;
    fprintf(stderr, "Generating 'data/results_long_double.csv'...\n");

    const int Nfft = T, Nfft_we = Nseg, K_mt = (int)(2 * NW_mt) - 1, hop = Nseg - overlap, num_segs = (T - Nseg) / hop + 1;
    long double *x = malloc(T * sizeof(*x)), *wT_h = malloc(T * sizeof(*wT_h)), *wS_w = malloc(Nseg * sizeof(*wS_w));
    long double **sine_tapers = malloc(K_mt * sizeof(long double*));
    for (int i = 0; i < K_mt; ++i) { sine_tapers[i] = malloc(T * sizeof(long double)); }

    fftwl_complex *in_fft = fftwl_malloc(Nfft * sizeof(*in_fft)), *out_fft = fftwl_malloc(Nfft * sizeof(*out_fft));
    fftwl_complex *in_welch = fftwl_malloc(Nfft_we * sizeof(*in_welch)), *out_welch = fftwl_malloc(Nfft_we * sizeof(*out_welch));
    fftwl_plan plan_fft = fftwl_plan_dft_1d(Nfft, in_fft, out_fft, FFTW_FORWARD, FFTW_ESTIMATE);
    fftwl_plan plan_welch = fftwl_plan_dft_1d(Nfft_we, in_welch, out_welch, FFTW_FORWARD, FFTW_ESTIMATE);

    make_hamming(wT_h, T);
    generate_sine_tapers(T, K_mt, sine_tapers);
    make_hamming(wS_w, Nseg);
    long double U_h = window_power(wT_h, T), U_w = window_power(wS_w, Nseg);
    long double *P_rect = calloc(Nfft, sizeof(*P_rect)), *P_hamm = calloc(Nfft, sizeof(*P_hamm)), *P_mt = calloc(Nfft, sizeof(*P_mt)), *P_welch = calloc(Nfft_we, sizeof(*P_welch));

    for (int m = 0; m < M; ++m) {
        for (int n = 0; n < T; ++n) x[n] = A1 * cosl(2 * M_PIl * f1 * n / fs) + A2 * cosl(2 * M_PIl * f2 * n / fs) + sigma * randn_box_muller();
        
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
        for (int i = 0; i < T; ++i) { in_fft[i][0] = x[i]; in_fft[i][1] = 0; }
        fftwl_execute(plan_fft);
        for (int i = 0; i < Nfft; ++i) P_rect[i] += out_fft[i][0] * out_fft[i][0] + out_fft[i][1] * out_fft[i][1];
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
        time_p += (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
        for (int i = 0; i < T; ++i) { in_fft[i][0] = x[i] * wT_h[i]; in_fft[i][1] = 0; }
        fftwl_execute(plan_fft);
        for (int i = 0; i < Nfft; ++i) P_hamm[i] += out_fft[i][0] * out_fft[i][0] + out_fft[i][1] * out_fft[i][1];
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
        time_ph += (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
        long double *eigenspectra_sum = calloc(Nfft, sizeof(*eigenspectra_sum));
        for (int k = 0; k < K_mt; ++k) {
            for (int i = 0; i < T; ++i) { in_fft[i][0] = x[i] * sine_tapers[k][i]; in_fft[i][1] = 0; }
            fftwl_execute(plan_fft);
            for (int i = 0; i < Nfft; ++i) eigenspectra_sum[i] += out_fft[i][0] * out_fft[i][0] + out_fft[i][1] * out_fft[i][1];
        }
        for (int i = 0; i < Nfft; ++i) P_mt[i] += eigenspectra_sum[i] / K_mt;
        free(eigenspectra_sum);
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
        time_m += (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
        for (int i = 0; i < num_segs; ++i) {
            int seg_start = i * hop;
            for (int n = 0; n < Nseg; ++n) { in_welch[n][0] = x[seg_start + n] * wS_w[n]; in_welch[n][1] = 0; }
            fftwl_execute(plan_welch);
            for (int k = 0; k < Nfft_we; ++k) P_welch[k] += out_welch[k][0] * out_welch[k][0] + out_welch[k][1] * out_welch[k][1];
        }
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
        time_w += (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    }

    for (int k = 0; k < Nfft; ++k) {
        P_rect[k] /= (fs * T * M);
        P_hamm[k] /= (fs * T * U_h * M);
        P_mt[k] /= (fs * M);
    }
    for (int k = 0; k < Nfft_we; ++k) P_welch[k] /= (fs * Nseg * U_w * num_segs * M);

    FILE *results_file = fopen("data/results_long_double.csv", "w"), *timing_file = fopen("data/timing_log_ld.csv", "w");
    fprintf(results_file, "k_res,PSD_periodogram,PSD_periodogram_hamm,PSD_mt,k_we,PSD_welch\n");
    for (int i = 0; i < Nfft; ++i) fprintf(results_file, "%d,%.21Lg,%.21Lg,%.21Lg,%d,%.21Lg\n", i, P_rect[i], P_hamm[i], P_mt[i], (i < Nfft_we) ? i : -1, (i < Nfft_we) ? P_welch[i] : NAN);
    fprintf(timing_file, "method,total_cpu_time_s,avg_cpu_time_per_realization_ms\nPeriodogram,%.6f,%.6f\nPeriodogram (Hamming),%.6f,%.6f\nMultitaper,%.6f,%.6f\nWelch,%.6f,%.6f\n", time_p, time_p / M * 1e3, time_ph, time_ph / M * 1e3, time_m, time_m / M * 1e3, time_w, time_w / M * 1e3);
    fclose(results_file); fclose(timing_file);

    free(x); free(wT_h); free(wS_w); for (int i = 0; i < K_mt; ++i) { free(sine_tapers[i]); } free(sine_tapers);
    free(P_rect); free(P_hamm); free(P_mt); free(P_welch);
    fftwl_destroy_plan(plan_fft); fftwl_destroy_plan(plan_welch);
    fftwl_free(in_fft); fftwl_free(out_fft); fftwl_free(in_welch); fftwl_free(out_welch);
}

int main() {
    setlocale(LC_NUMERIC, "C");
    const long double fs = 1000.0L, sigma_comp = 3.5L, f1 = 120.0L, A1 = 1.0L, f2 = 250.0L, A2 = 0.7L, NW_mt = 4.0L, f0_stats = 250.0L, A_stats = 1.0L;
    const int T = 1024, M_comp = 100, M_stats = 1000, Nseg = 256, overlap = 128;
    srand(time(NULL));
    fftwl_init_threads();
    fftwl_plan_with_nthreads(4);

    generate_psd_data(fs, T, M_comp, sigma_comp, A1, f1, A2, f2, Nseg, overlap, NW_mt);

    fprintf(stderr, "\nStarting statistical analysis...\n");
    const long double sigmas[] = {0.1L, 0.5L, 1.0L, 3.5L, 10.0L, 25.0L, 50.0L};
    const int num_sigmas = sizeof(sigmas) / sizeof(sigmas[0]);
    const int kf0_welch = (int)(f0_stats * Nseg / fs), hop = Nseg - overlap, num_segs = (T - Nseg) / hop + 1;
    long double *wS = malloc(Nseg * sizeof(*wS));
    make_hamming(wS, Nseg);
    const long double Useg = window_power(wS, Nseg), norm_welch = fs * Nseg * Useg * num_segs;

    FILE *stats_file = fopen("data/statistics.csv", "w"), *hist_file = fopen("data/histograms.csv", "w"), *realizations_file = fopen("data/psd_realizations.csv", "w");
    fprintf(stats_file, "sigma,mean_psd,variance_psd,skewness,kurtosis\n");
    fprintf(hist_file, "sigma,bin_center,count\n");
    fprintf(realizations_file, "sigma,realization_index,psd_value\n");

    fftwl_complex *in_w = fftwl_malloc(Nseg * sizeof(*in_w)), *out_w = fftwl_malloc(Nseg * sizeof(*out_w));
    fftwl_plan p_w = fftwl_plan_dft_1d(Nseg, in_w, out_w, FFTW_FORWARD, FFTW_ESTIMATE);

    for (int i = 0; i < num_sigmas; ++i) {
        fprintf(stderr, "Processing sigma = %.4Lf ...\n", sigmas[i]);
        long double *psd_vals = malloc(M_stats * sizeof(*psd_vals)), *x = malloc(T * sizeof(*x));
        for (int m = 0; m < M_stats; ++m) {
            for (int n = 0; n < T; ++n) x[n] = A_stats * cosl(2 * M_PIl * f0_stats * n / fs) + sigmas[i] * randn_box_muller();
            long double Pxx_sum = 0;
            for (int s = 0; s < num_segs; ++s) {
                int start = s * hop;
                for (int n = 0; n < Nseg; ++n) { in_w[n][0] = x[start + n] * wS[n]; in_w[n][1] = 0; }
                fftwl_execute(p_w);
                Pxx_sum += out_w[kf0_welch][0] * out_w[kf0_welch][0] + out_w[kf0_welch][1] * out_w[kf0_welch][1];
            }
            psd_vals[m] = (2.0L * Pxx_sum) / norm_welch;
            fprintf(realizations_file, "%.8Lg,%d,%.17Lg\n", sigmas[i], m, psd_vals[m]);
        }
        Stats stats = calculate_statistics(psd_vals, M_stats);
        fprintf(stats_file, "%.8Lg,%.17Lg,%.17Lg,%.17Lg,%.17Lg\n", sigmas[i], stats.mean, stats.variance, stats.skewness, stats.kurtosis);
        generate_histogram(hist_file, psd_vals, M_stats, sigmas[i], 40);
        free(psd_vals);
        free(x);
    }
    fftwl_destroy_plan(p_w);
    fftwl_free(in_w);
    fftwl_free(out_w);
    fclose(stats_file);
    fclose(hist_file);
    fclose(realizations_file);
    free(wS);
    fftwl_cleanup_threads();
    return 0;
}