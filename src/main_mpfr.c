#define _GNU_SOURCE
#include <math.h>
#include <locale.h>
#include <fftw3.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <pthread.h>
#include <mpfr.h>

// --- CONFIGURATION ---
const mpfr_prec_t PREC = 4096; // Set the bits for arbitrary precision

// --- MPFR DATA STRUCTURES & UTILITIES ---
typedef struct { mpfr_t mean, variance, skewness, kurtosis; } Stats_mpfr;

void randn_box_muller_mpfr(mpfr_t result) {
    mpfr_t u1, u2, log_u1, two_pi, cos_u2, sqrt_term, temp;
    mpfr_init2(u1, PREC); mpfr_init2(u2, PREC); mpfr_init2(log_u1, PREC);
    mpfr_init2(two_pi, PREC); mpfr_init2(cos_u2, PREC); mpfr_init2(sqrt_term, PREC);
    mpfr_init2(temp, PREC);
    mpfr_set_d(u1, (double)rand() / (double)RAND_MAX, MPFR_RNDN);
    mpfr_set_d(u2, (double)rand() / (double)RAND_MAX, MPFR_RNDN);
    mpfr_log(log_u1, u1, MPFR_RNDN);
    mpfr_mul_si(sqrt_term, log_u1, -2, MPFR_RNDN);
    mpfr_sqrt(sqrt_term, sqrt_term, MPFR_RNDN);
    mpfr_const_pi(two_pi, PREC);
    mpfr_mul_ui(two_pi, two_pi, 2, MPFR_RNDN);
    mpfr_mul(temp, two_pi, u2, MPFR_RNDN);
    mpfr_cos(cos_u2, temp, MPFR_RNDN);
    mpfr_mul(result, sqrt_term, cos_u2, MPFR_RNDN);
    mpfr_clears(u1, u2, log_u1, two_pi, cos_u2, sqrt_term, temp, (mpfr_ptr)0);
}

Stats_mpfr calculate_statistics_mpfr(mpfr_t* data, int n) {
    Stats_mpfr result;
    mpfr_init2(result.mean, PREC); mpfr_init2(result.variance, PREC);
    mpfr_init2(result.skewness, PREC); mpfr_init2(result.kurtosis, PREC);
    mpfr_set_zero(result.mean, 1); mpfr_set_zero(result.variance, 1);
    mpfr_set_zero(result.skewness, 1); mpfr_set_zero(result.kurtosis, 1);
    if (n < 4) return result;

    mpfr_t sum, s2, s3, s4, delta, std_dev, temp;
    mpfr_init2(sum, PREC); mpfr_init2(s2, PREC); mpfr_init2(s3, PREC); mpfr_init2(s4, PREC);
    mpfr_init2(delta, PREC); mpfr_init2(std_dev, PREC); mpfr_init2(temp, PREC);
    mpfr_set_zero(sum, 1);

    for (int i = 0; i < n; ++i) mpfr_add(sum, sum, data[i], MPFR_RNDN);
    mpfr_div_ui(result.mean, sum, n, MPFR_RNDN);

    mpfr_set_zero(s2, 1); mpfr_set_zero(s3, 1); mpfr_set_zero(s4, 1);
    for (int i = 0; i < n; ++i) {
        mpfr_sub(delta, data[i], result.mean, MPFR_RNDN);
        mpfr_pow_ui(temp, delta, 2, MPFR_RNDN); mpfr_add(s2, s2, temp, MPFR_RNDN);
        mpfr_pow_ui(temp, delta, 3, MPFR_RNDN); mpfr_add(s3, s3, temp, MPFR_RNDN);
        mpfr_pow_ui(temp, delta, 4, MPFR_RNDN); mpfr_add(s4, s4, temp, MPFR_RNDN);
    }
    mpfr_div_ui(result.variance, s2, n, MPFR_RNDN);

    if (mpfr_cmp_ui(result.variance, 0) > 0) {
        mpfr_sqrt(std_dev, result.variance, MPFR_RNDN);
        mpfr_pow_ui(temp, std_dev, 3, MPFR_RNDN);
        mpfr_div_ui(result.skewness, s3, n, MPFR_RNDN);
        mpfr_div(result.skewness, result.skewness, temp, MPFR_RNDN);
        
        mpfr_pow_ui(temp, result.variance, 2, MPFR_RNDN);
        mpfr_div_ui(result.kurtosis, s4, n, MPFR_RNDN);
        mpfr_div(result.kurtosis, result.kurtosis, temp, MPFR_RNDN);
        mpfr_sub_ui(result.kurtosis, result.kurtosis, 3, MPFR_RNDN);
    }
    mpfr_clears(sum, s2, s3, s4, delta, std_dev, temp, (mpfr_ptr)0);
    return result;
}

void generate_histogram_mpfr(FILE* file, mpfr_t* data, int n, mpfr_t sigma, int nbins) {
    if (n == 0) return;
    mpfr_t min_v, max_v, width, bin_center, temp_idx;
    mpfr_init2(min_v, PREC); mpfr_init2(max_v, PREC); mpfr_init2(width, PREC);
    mpfr_init2(bin_center, PREC); mpfr_init2(temp_idx, PREC);
    mpfr_set(min_v, data[0], MPFR_RNDN);
    mpfr_set(max_v, data[0], MPFR_RNDN);

    for (int i = 1; i < n; ++i) {
        if (mpfr_cmp(data[i], min_v) < 0) mpfr_set(min_v, data[i], MPFR_RNDN);
        if (mpfr_cmp(data[i], max_v) > 0) mpfr_set(max_v, data[i], MPFR_RNDN);
    }
    
    mpfr_sub(width, max_v, min_v, MPFR_RNDN);
    if (mpfr_zero_p(width)) {
        mpfr_fprintf(file, "%.8Rg,%.50Rg,%d\n", sigma, min_v, n);
        mpfr_clears(min_v, max_v, width, bin_center, temp_idx, (mpfr_ptr)0);
        return;
    }
    mpfr_div_ui(width, width, nbins, MPFR_RNDN);
    int *bins = calloc(nbins, sizeof(int));
    for (int i = 0; i < n; ++i) {
        mpfr_sub(temp_idx, data[i], min_v, MPFR_RNDN);
        mpfr_div(temp_idx, temp_idx, width, MPFR_RNDN);
        int idx = mpfr_get_si(temp_idx, MPFR_RNDD);
        if (idx >= nbins) idx = nbins - 1;
        if (idx < 0) idx = 0;
        bins[idx]++;
    }
    for (int i = 0; i < nbins; ++i) {
        mpfr_mul_ui(bin_center, width, i, MPFR_RNDN);
        mpfr_add(bin_center, bin_center, min_v, MPFR_RNDN);
        mpfr_mul_d(temp_idx, width, 0.5, MPFR_RNDN);
        mpfr_add(bin_center, bin_center, temp_idx, MPFR_RNDN);
        mpfr_fprintf(file, "%.8Rg,%.50Rg,%d\n", sigma, bin_center, bins[i]);
    }
    free(bins);
    mpfr_clears(min_v, max_v, width, bin_center, temp_idx, (mpfr_ptr)0);
}

void make_hamming_mpfr(mpfr_t *w, int N) {
    mpfr_t n_mpfr, N_minus_1, two_pi, term, cos_term;
    mpfr_init2(n_mpfr, PREC); mpfr_init2(N_minus_1, PREC); mpfr_init2(two_pi, PREC);
    mpfr_init2(term, PREC); mpfr_init2(cos_term, PREC);
    mpfr_set_ui(N_minus_1, N - 1, MPFR_RNDN);
    mpfr_const_pi(two_pi, PREC);
    mpfr_mul_ui(two_pi, two_pi, 2, MPFR_RNDN);
    for (int n = 0; n < N; ++n) {
        mpfr_set_ui(n_mpfr, n, MPFR_RNDN);
        mpfr_mul(term, two_pi, n_mpfr, MPFR_RNDN);
        mpfr_div(term, term, N_minus_1, MPFR_RNDN);
        mpfr_cos(cos_term, term, MPFR_RNDN);
        mpfr_mul_d(cos_term, cos_term, -0.46, MPFR_RNDN);
        mpfr_add_d(w[n], cos_term, 0.54, MPFR_RNDN);
    }
    mpfr_clears(n_mpfr, N_minus_1, two_pi, term, cos_term, (mpfr_ptr)0);
}

void window_power_mpfr(mpfr_t result, mpfr_t *w, int N) {
    mpfr_t sum_sq, temp;
    mpfr_init2(sum_sq, PREC); mpfr_init2(temp, PREC);
    mpfr_set_zero(sum_sq, 1);
    for (int n = 0; n < N; ++n) {
        mpfr_mul(temp, w[n], w[n], MPFR_RNDN);
        mpfr_add(sum_sq, sum_sq, temp, MPFR_RNDN);
    }
    mpfr_div_ui(result, sum_sq, N, MPFR_RNDN);
    mpfr_clears(sum_sq, temp, (mpfr_ptr)0);
}

void generate_sine_tapers_mpfr(int N, int K, mpfr_t **tapers) {
    mpfr_t norm, pi, N_plus_1, k_plus_1, i_plus_1, arg, sin_arg;
    mpfr_init2(norm, PREC); mpfr_init2(pi, PREC); mpfr_init2(N_plus_1, PREC);
    mpfr_init2(k_plus_1, PREC); mpfr_init2(i_plus_1, PREC); mpfr_init2(arg, PREC);
    mpfr_init2(sin_arg, PREC);
    mpfr_const_pi(pi, PREC);
    mpfr_set_ui(N_plus_1, N + 1, MPFR_RNDN);
    mpfr_set_d(norm, 2.0, MPFR_RNDN);
    mpfr_div(norm, norm, N_plus_1, MPFR_RNDN);
    mpfr_sqrt(norm, norm, MPFR_RNDN);
    for (int k = 0; k < K; k++) {
        mpfr_set_ui(k_plus_1, k + 1, MPFR_RNDN);
        for (int i = 0; i < N; i++) {
            mpfr_set_ui(i_plus_1, i + 1, MPFR_RNDN);
            mpfr_mul(arg, pi, k_plus_1, MPFR_RNDN);
            mpfr_mul(arg, arg, i_plus_1, MPFR_RNDN);
            mpfr_div(arg, arg, N_plus_1, MPFR_RNDN);
            mpfr_sin(sin_arg, arg, MPFR_RNDN);
            mpfr_mul(tapers[k][i], norm, sin_arg, MPFR_RNDN);
        }
    }
    mpfr_clears(norm, pi, N_plus_1, k_plus_1, i_plus_1, arg, sin_arg, (mpfr_ptr)0);
}

void generate_psd_data_mpfr(double fs_d, int T, int M, double sigma_d, double A1_d, double f1_d, double A2_d, double f2_d, int Nseg, int overlap, double NW_mt_d) {
    struct timespec start, end;
    double time_p = 0, time_ph = 0, time_m = 0, time_w = 0;
    fprintf(stderr, "Generating 'data/results_mpfr.csv' with %ld-bit precision...\n", (long)PREC);

    mpfr_t fs, sigma, A1, f1, A2, f2;
    mpfr_init2(fs, PREC); mpfr_set_d(fs, fs_d, MPFR_RNDN);
    mpfr_init2(sigma, PREC); mpfr_set_d(sigma, sigma_d, MPFR_RNDN);
    mpfr_init2(A1, PREC); mpfr_set_d(A1, A1_d, MPFR_RNDN);
    mpfr_init2(f1, PREC); mpfr_set_d(f1, f1_d, MPFR_RNDN);
    mpfr_init2(A2, PREC); mpfr_set_d(A2, A2_d, MPFR_RNDN);
    mpfr_init2(f2, PREC); mpfr_set_d(f2, f2_d, MPFR_RNDN);

    const int Nfft = T, Nfft_we = Nseg, K_mt = (int)(2 * NW_mt_d) - 1, hop = Nseg - overlap, num_segs = (T - Nseg) / hop + 1;
    
    mpfr_t *x = malloc(T * sizeof(mpfr_t));
    mpfr_t *wT_h = malloc(T * sizeof(mpfr_t));
    mpfr_t *wS_w = malloc(Nseg * sizeof(mpfr_t));
    mpfr_t **sine_tapers = malloc(K_mt * sizeof(mpfr_t*));
    for(int i=0; i<T; ++i) { mpfr_init2(x[i], PREC); mpfr_init2(wT_h[i], PREC); }
    for(int i=0; i<Nseg; ++i) mpfr_init2(wS_w[i], PREC);
    for(int i=0; i<K_mt; ++i) {
        sine_tapers[i] = malloc(T * sizeof(mpfr_t));
        for(int j=0; j<T; ++j) mpfr_init2(sine_tapers[i][j], PREC);
    }

    fftwl_complex *in_fft = fftwl_malloc(Nfft * sizeof(*in_fft)), *out_fft = fftwl_malloc(Nfft * sizeof(*out_fft));
    fftwl_complex *in_welch = fftwl_malloc(Nfft_we * sizeof(*in_welch)), *out_welch = fftwl_malloc(Nfft_we * sizeof(*out_welch));
    fftwl_plan plan_fft = fftwl_plan_dft_1d(Nfft, in_fft, out_fft, FFTW_FORWARD, FFTW_ESTIMATE);
    fftwl_plan plan_welch = fftwl_plan_dft_1d(Nfft_we, in_welch, out_welch, FFTW_FORWARD, FFTW_ESTIMATE);

    make_hamming_mpfr(wT_h, T);
    make_hamming_mpfr(wS_w, Nseg);
    generate_sine_tapers_mpfr(T, K_mt, sine_tapers);
    
    mpfr_t U_h, U_w;
    mpfr_init2(U_h, PREC); mpfr_init2(U_w, PREC);
    window_power_mpfr(U_h, wT_h, T);
    window_power_mpfr(U_w, wS_w, Nseg);

    mpfr_t *P_rect = malloc(Nfft * sizeof(mpfr_t)), *P_hamm = malloc(Nfft * sizeof(mpfr_t)), *P_mt = malloc(Nfft * sizeof(mpfr_t)), *P_welch = malloc(Nfft_we * sizeof(mpfr_t));
    for(int i=0; i<Nfft; ++i) { mpfr_init2(P_rect[i], PREC); mpfr_set_zero(P_rect[i], 1); mpfr_init2(P_hamm[i], PREC); mpfr_set_zero(P_hamm[i], 1); mpfr_init2(P_mt[i], PREC); mpfr_set_zero(P_mt[i], 1); }
    for(int i=0; i<Nfft_we; ++i) { mpfr_init2(P_welch[i], PREC); mpfr_set_zero(P_welch[i], 1); }

    mpfr_t term1, term2, noise, temp, two_pi_fs, n_mpfr;
    mpfr_init2(term1, PREC); mpfr_init2(term2, PREC); mpfr_init2(noise, PREC);
    mpfr_init2(temp, PREC); mpfr_init2(two_pi_fs, PREC); mpfr_init2(n_mpfr, PREC);
    mpfr_const_pi(two_pi_fs, PREC);
    mpfr_mul_ui(two_pi_fs, two_pi_fs, 2, MPFR_RNDN);
    mpfr_div(two_pi_fs, two_pi_fs, fs, MPFR_RNDN);

    for (int m = 0; m < M; ++m) {
        for (int n = 0; n < T; ++n) {
            mpfr_set_ui(n_mpfr, n, MPFR_RNDN);
            mpfr_mul(term1, f1, n_mpfr, MPFR_RNDN); mpfr_mul(term1, term1, two_pi_fs, MPFR_RNDN); mpfr_cos(term1, term1, MPFR_RNDN); mpfr_mul(term1, term1, A1, MPFR_RNDN);
            mpfr_mul(term2, f2, n_mpfr, MPFR_RNDN); mpfr_mul(term2, term2, two_pi_fs, MPFR_RNDN); mpfr_cos(term2, term2, MPFR_RNDN); mpfr_mul(term2, term2, A2, MPFR_RNDN);
            randn_box_muller_mpfr(noise); mpfr_mul(noise, noise, sigma, MPFR_RNDN);
            mpfr_add(x[n], term1, term2, MPFR_RNDN); mpfr_add(x[n], x[n], noise, MPFR_RNDN);
        }

        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
        for (int i = 0; i < T; ++i) { in_fft[i][0] = mpfr_get_ld(x[i], MPFR_RNDN); in_fft[i][1] = 0; }
        fftwl_execute(plan_fft);
        for (int i = 0; i < Nfft; ++i) { mpfr_set_ld(temp, out_fft[i][0]*out_fft[i][0] + out_fft[i][1]*out_fft[i][1], MPFR_RNDN); mpfr_add(P_rect[i], P_rect[i], temp, MPFR_RNDN); }
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end); time_p += (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
        for (int i = 0; i < T; ++i) { mpfr_mul(temp, x[i], wT_h[i], MPFR_RNDN); in_fft[i][0] = mpfr_get_ld(temp, MPFR_RNDN); in_fft[i][1] = 0; }
        fftwl_execute(plan_fft);
        for (int i = 0; i < Nfft; ++i) { mpfr_set_ld(temp, out_fft[i][0]*out_fft[i][0] + out_fft[i][1]*out_fft[i][1], MPFR_RNDN); mpfr_add(P_hamm[i], P_hamm[i], temp, MPFR_RNDN); }
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end); time_ph += (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
        mpfr_t *e_sum = malloc(Nfft * sizeof(mpfr_t)); for(int i=0; i<Nfft; ++i) { mpfr_init2(e_sum[i], PREC); mpfr_set_zero(e_sum[i], 1); }
        for (int k = 0; k < K_mt; ++k) {
            for (int i = 0; i < T; ++i) { mpfr_mul(temp, x[i], sine_tapers[k][i], MPFR_RNDN); in_fft[i][0] = mpfr_get_ld(temp, MPFR_RNDN); in_fft[i][1] = 0; }
            fftwl_execute(plan_fft);
            for (int i = 0; i < Nfft; ++i) { mpfr_set_ld(temp, out_fft[i][0]*out_fft[i][0] + out_fft[i][1]*out_fft[i][1], MPFR_RNDN); mpfr_add(e_sum[i], e_sum[i], temp, MPFR_RNDN); }
        }
        for (int i = 0; i < Nfft; ++i) { mpfr_div_ui(temp, e_sum[i], K_mt, MPFR_RNDN); mpfr_add(P_mt[i], P_mt[i], temp, MPFR_RNDN); }
        for(int i=0; i<Nfft; ++i) { mpfr_clear(e_sum[i]); } free(e_sum);
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end); time_m += (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
        for (int i = 0; i < num_segs; ++i) {
            for (int n = 0; n < Nseg; ++n) { mpfr_mul(temp, x[i*hop + n], wS_w[n], MPFR_RNDN); in_welch[n][0] = mpfr_get_ld(temp, MPFR_RNDN); in_welch[n][1] = 0; }
            fftwl_execute(plan_welch);
            for (int k = 0; k < Nfft_we; ++k) { mpfr_set_ld(temp, out_welch[k][0]*out_welch[k][0] + out_welch[k][1]*out_welch[k][1], MPFR_RNDN); mpfr_add(P_welch[k], P_welch[k], temp, MPFR_RNDN); }
        }
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end); time_w += (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    }

    mpfr_t norm_p, norm_ph, norm_m, norm_w;
    mpfr_init2(norm_p, PREC); mpfr_init2(norm_ph, PREC); mpfr_init2(norm_m, PREC); mpfr_init2(norm_w, PREC);
    mpfr_mul_ui(norm_p, fs, T * M, MPFR_RNDN);
    mpfr_mul(norm_ph, norm_p, U_h, MPFR_RNDN);
    mpfr_mul_ui(norm_m, fs, M, MPFR_RNDN);
    mpfr_mul_ui(norm_w, fs, Nseg * num_segs * M, MPFR_RNDN);
    mpfr_mul(norm_w, norm_w, U_w, MPFR_RNDN);

    for (int k = 0; k < Nfft; ++k) { mpfr_div(P_rect[k], P_rect[k], norm_p, MPFR_RNDN); mpfr_div(P_hamm[k], P_hamm[k], norm_ph, MPFR_RNDN); mpfr_div(P_mt[k], P_mt[k], norm_m, MPFR_RNDN); }
    for (int k = 0; k < Nfft_we; ++k) mpfr_div(P_welch[k], P_welch[k], norm_w, MPFR_RNDN);

    FILE *res_f = fopen("data/results_mpfr.csv", "w"), *time_f = fopen("data/timing_log_mpfr.csv", "w");
    fprintf(res_f, "k_res,PSD_periodogram,PSD_periodogram_hamm,PSD_mt,k_we,PSD_welch\n");
    for (int i = 0; i < Nfft; ++i) {
        fprintf(res_f, "%d,", i);
        mpfr_fprintf(res_f, "%.50Rg,%.50Rg,%.50Rg,", P_rect[i], P_hamm[i], P_mt[i]);
        if (i < Nfft_we) mpfr_fprintf(res_f, "%d,%.50Rg\n", i, P_welch[i]);
        else fprintf(res_f, "-1,nan\n");
    }
    fprintf(time_f, "method,total_cpu_time_s,avg_cpu_time_per_realization_ms\n");
    fprintf(time_f, "Periodogram,%.6f,%.6f\n", time_p, time_p / M * 1e3);
    fprintf(time_f, "Periodogram (Hamming),%.6f,%.6f\n", time_ph, time_ph / M * 1e3);
    fprintf(time_f, "Multitaper,%.6f,%.6f\n", time_m, time_m / M * 1e3);
    fprintf(time_f, "Welch,%.6f,%.6f\n", time_w, time_w / M * 1e3);
    fclose(res_f); fclose(time_f);

    mpfr_clears(fs, sigma, A1, f1, A2, f2, U_h, U_w, norm_p, norm_ph, norm_m, norm_w, (mpfr_ptr)0);
    mpfr_clears(term1, term2, noise, temp, two_pi_fs, n_mpfr, (mpfr_ptr)0);
    for(int i=0; i<T; ++i) { mpfr_clear(x[i]); mpfr_clear(wT_h[i]); } free(x); free(wT_h);
    for(int i=0; i<Nseg; ++i) { mpfr_clear(wS_w[i]); } free(wS_w);
    for(int i=0; i<K_mt; ++i) { for(int j=0; j<T; ++j) mpfr_clear(sine_tapers[i][j]); free(sine_tapers[i]); } free(sine_tapers);
    for(int i=0; i<Nfft; ++i) { mpfr_clear(P_rect[i]); mpfr_clear(P_hamm[i]); mpfr_clear(P_mt[i]); } free(P_rect); free(P_hamm); free(P_mt);
    for(int i=0; i<Nfft_we; ++i) { mpfr_clear(P_welch[i]); } free(P_welch);
    fftwl_destroy_plan(plan_fft); fftwl_destroy_plan(plan_welch);
    fftwl_free(in_fft); fftwl_free(out_fft); fftwl_free(in_welch); fftwl_free(out_welch);
}

int main() {
    setlocale(LC_NUMERIC, "C");
    const double fs_d = 1000.0, sigma_comp_d = 3.5, f1_d = 120.0, A1_d = 1.0, f2_d = 250.0, A2_d = 0.7, NW_mt_d = 4.0, f0_stats_d = 250.0, A_stats_d = 1.0;
    const int T = 1024, M_comp = 100, M_stats = 1000, Nseg = 256, overlap = 128;
    srand(time(NULL));
    fftwl_init_threads();
    fftwl_plan_with_nthreads(4);

    generate_psd_data_mpfr(fs_d, T, M_comp, sigma_comp_d, A1_d, f1_d, A2_d, f2_d, Nseg, overlap, NW_mt_d);

    fprintf(stderr, "\nStarting statistical analysis (MPFR)...\n");
    const double sigmas_d[] = {0.1, 0.5, 1.0, 3.5, 10.0, 25.0, 50.0};
    const int num_s = sizeof(sigmas_d) / sizeof(sigmas_d[0]);
    const int kf0_w = (int)(f0_stats_d * Nseg / fs_d), hop = Nseg - overlap, segs = (T - Nseg) / hop + 1;

    mpfr_t fs, f0_stats, A_stats, wS[Nseg], Useg, norm_w;
    mpfr_init2(fs, PREC); mpfr_set_d(fs, fs_d, MPFR_RNDN);
    mpfr_init2(f0_stats, PREC); mpfr_set_d(f0_stats, f0_stats_d, MPFR_RNDN);
    mpfr_init2(A_stats, PREC); mpfr_set_d(A_stats, A_stats_d, MPFR_RNDN);
    for(int i=0; i<Nseg; ++i) mpfr_init2(wS[i], PREC);
    mpfr_init2(Useg, PREC); mpfr_init2(norm_w, PREC);
    
    make_hamming_mpfr(wS, Nseg);
    window_power_mpfr(Useg, wS, Nseg);
    mpfr_mul_ui(norm_w, fs, Nseg * segs, MPFR_RNDN);
    mpfr_mul(norm_w, norm_w, Useg, MPFR_RNDN);

    FILE *sf = fopen("data/statistics_mpfr.csv", "w"), *hf = fopen("data/histograms_mpfr.csv", "w");
    fprintf(sf, "sigma,mean_psd,variance_psd,skewness,kurtosis\n");
    fprintf(hf, "sigma,bin_center,count\n");

    fftwl_complex *in_w = fftwl_malloc(Nseg * sizeof(*in_w)), *out_w = fftwl_malloc(Nseg * sizeof(*out_w));
    fftwl_plan p_w = fftwl_plan_dft_1d(Nseg, in_w, out_w, FFTW_FORWARD, FFTW_ESTIMATE);
    
    mpfr_t term1, noise, temp, two_pi_fs, n_mpfr;
    mpfr_init2(term1, PREC); mpfr_init2(noise, PREC); mpfr_init2(temp, PREC);
    mpfr_init2(two_pi_fs, PREC); mpfr_init2(n_mpfr, PREC);
    mpfr_const_pi(two_pi_fs, PREC);
    mpfr_mul_ui(two_pi_fs, two_pi_fs, 2, MPFR_RNDN);
    mpfr_div(two_pi_fs, two_pi_fs, fs, MPFR_RNDN);

    for (int i = 0; i < num_s; ++i) {
        mpfr_t current_sigma;
        mpfr_init2(current_sigma, PREC);
        mpfr_set_d(current_sigma, sigmas_d[i], MPFR_RNDN);
        fprintf(stderr, "Processing sigma = %.4f ...\n", sigmas_d[i]);

        mpfr_t *psd_vals = malloc(M_stats * sizeof(mpfr_t));
        mpfr_t *x = malloc(T * sizeof(mpfr_t));
        for(int j=0; j<M_stats; ++j) mpfr_init2(psd_vals[j], PREC);
        for(int j=0; j<T; ++j) mpfr_init2(x[j], PREC);

        for(int m=0; m<M_stats; ++m){
            for(int n=0; n<T; ++n) {
                mpfr_set_ui(n_mpfr, n, MPFR_RNDN);
                mpfr_mul(term1, f0_stats, n_mpfr, MPFR_RNDN);
                mpfr_mul(term1, term1, two_pi_fs, MPFR_RNDN);
                mpfr_cos(term1, term1, MPFR_RNDN);
                mpfr_mul(term1, term1, A_stats, MPFR_RNDN);
                randn_box_muller_mpfr(noise);
                mpfr_mul(noise, noise, current_sigma, MPFR_RNDN);
                mpfr_add(x[n], term1, noise, MPFR_RNDN);
            }
            mpfr_t Pxx_sum;
            mpfr_init2(Pxx_sum, PREC);
            mpfr_set_zero(Pxx_sum, 1);
            for(int s=0; s<segs; ++s){
                int start = s*hop;
                for(int n=0; n<Nseg; ++n){
                    mpfr_mul(temp, x[start+n], wS[n], MPFR_RNDN);
                    in_w[n][0] = mpfr_get_ld(temp, MPFR_RNDN);
                    in_w[n][1] = 0;
                }
                fftwl_execute(p_w);
                long double mag_sq = out_w[kf0_w][0]*out_w[kf0_w][0] + out_w[kf0_w][1]*out_w[kf0_w][1];
                mpfr_set_ld(temp, mag_sq, MPFR_RNDN);
                mpfr_add(Pxx_sum, Pxx_sum, temp, MPFR_RNDN);
            }
            mpfr_mul_ui(psd_vals[m], Pxx_sum, 2, MPFR_RNDN);
            mpfr_div(psd_vals[m], psd_vals[m], norm_w, MPFR_RNDN);
            mpfr_clear(Pxx_sum);
        }
        
        Stats_mpfr stats = calculate_statistics_mpfr(psd_vals, M_stats);
        mpfr_fprintf(sf, "%.8Rg,%.50Rg,%.50Rg,%.50Rg,%.50Rg\n", current_sigma, stats.mean, stats.variance, stats.skewness, stats.kurtosis);
        generate_histogram_mpfr(hf, psd_vals, M_stats, current_sigma, 40);
        
        mpfr_clear(current_sigma);
        for(int j=0; j<M_stats; ++j) { mpfr_clear(psd_vals[j]); } free(psd_vals);
        for(int j=0; j<T; ++j) { mpfr_clear(x[j]); } free(x);
        mpfr_clears(stats.mean, stats.variance, stats.skewness, stats.kurtosis, (mpfr_ptr)0);
    }

    fclose(sf); fclose(hf);
    fftwl_destroy_plan(p_w); fftwl_free(in_w); fftwl_free(out_w);
    mpfr_clears(fs, f0_stats, A_stats, Useg, norm_w, (mpfr_ptr)0);
    mpfr_clears(term1, noise, temp, two_pi_fs, n_mpfr, (mpfr_ptr)0);
    for(int i=0; i<Nseg; ++i) mpfr_clear(wS[i]);

    fftwl_cleanup_threads();
    mpfr_free_cache();
    return 0;
}