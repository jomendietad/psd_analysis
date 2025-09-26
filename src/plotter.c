#include <stdio.h>
#include <stdlib.h>

const char* STATS_FILENAME = "data/statistics.csv";
const char* HIST_FILENAME = "data/histograms.csv";
const char* RESULTS_FILENAME = "data/results.csv"; 
const char* FIT_LOG_FILENAME = "data/fit.log";

const double sigmas_to_plot[] = {0.5, 3.5, 50.0};
const int num_sigmas_to_plot = sizeof(sigmas_to_plot) / sizeof(sigmas_to_plot[0]);

void plot_statistics() {
    FILE *gnuplotPipe = popen("gnuplot -persistent", "w");
    if (!gnuplotPipe) { fprintf(stderr, "Error Gnuplot.\n"); return; }
    fprintf(gnuplotPipe, "set terminal wxt 0 title 'Statistical Analysis'\n");
    fprintf(gnuplotPipe, "set title 'Skewness and Kurtosis vs. Noise Deviation (sigma)'\n");
    fprintf(gnuplotPipe, "set xlabel 'Sigma (Noise Standard Deviation)'\n");
    fprintf(gnuplotPipe, "set ylabel 'Statistical Moment Value'\n");
    fprintf(gnuplotPipe, "set grid\n");
    fprintf(gnuplotPipe, "set key top right\n");
    fprintf(gnuplotPipe, "set datafile separator ','\n");
    fprintf(gnuplotPipe, "set logscale x\n");
    fprintf(gnuplotPipe, "set arrow from graph 0, first 0 to graph 1, first 0 nohead dashtype 2 lc 'black' lw 1\n");
    fprintf(gnuplotPipe, "set label 'Ideal Gaussian' at graph 0.02, first 0.1 tc 'black'\n");
    fprintf(gnuplotPipe, "plot '%s' u 1:4 w lp pt 7 ps 1 lw 2 t 'Skewness', '' u 1:5 w lp pt 5 ps 1 lw 2 t 'Excess Kurtosis'\n", STATS_FILENAME);
    pclose(gnuplotPipe);
    printf("C/Gnuplot: Statistics plot generated (Window 0).\n");
}

void plot_histograms() {
    FILE *gnuplotPipe = popen("gnuplot -persistent", "w");
    if (!gnuplotPipe) { fprintf(stderr, "Error Gnuplot.\n"); return; }

    fprintf(gnuplotPipe, "set fit logfile '%s'\n", FIT_LOG_FILENAME);

    fprintf(gnuplotPipe, "set terminal wxt 1 size 1200, 900 title 'PSD Histograms'\n");
    fprintf(gnuplotPipe, "set multiplot layout %d,1 title 'Distribution of PSD at f0 for Different Noise Levels' font ',14'\n", num_sigmas_to_plot);
    fprintf(gnuplotPipe, "set datafile separator ','\n");
    fprintf(gnuplotPipe, "set xlabel 'PSD [Units/Hz]' font ',12'\n");
    fprintf(gnuplotPipe, "set ylabel 'Occurrence Count' font ',12'\n");
    for (int i = 0; i < num_sigmas_to_plot; ++i) {
        double current_sigma = sigmas_to_plot[i];
        fprintf(gnuplotPipe, "set title 'Sigma = %.2f' font ',12'\n", current_sigma);
        fprintf(gnuplotPipe, "g(x, mu, s, A) = A * exp(-(x-mu)**2 / (2*s**2))\n");
        fprintf(gnuplotPipe, "fit g(x, mu, s, A) '%s' using (column(1)==%f ? $2 : 1/0):3 via mu, s, A\n", HIST_FILENAME, current_sigma);
        fprintf(gnuplotPipe, "plot '%s' using (column(1)==%f ? $2 : 1/0):3 with boxes fill solid 0.5 title 'Histogram', g(x, mu, s, A) with lines lw 2 lc 'red' title 'Gaussian Fit'\n", HIST_FILENAME, current_sigma);
    }
    fprintf(gnuplotPipe, "unset multiplot\n");
    pclose(gnuplotPipe);
    printf("C/Gnuplot: Histograms plot generated (Window 1).\n");
}

void plot_psd_comparison() {
    const double FS = 1000.0;
    const int T = 1024;
    const int NFFT_RES = T; 
    const int NFFT_WE = 256;
    
    FILE *gnuplotPipe = popen("gnuplot -persistent", "w");
    if (!gnuplotPipe) { fprintf(stderr, "Error Gnuplot.\n"); return; }

    fprintf(gnuplotPipe, "set terminal wxt 2 size 1024, 768 title 'PSD Estimator Comparison'\n");
    fprintf(gnuplotPipe, "set title 'Comparison of Power Spectral Density (PSD) Estimators for a Two-Tone Signal'\n");
    fprintf(gnuplotPipe, "set xlabel 'Frequency [Hz]'\n");
    fprintf(gnuplotPipe, "set ylabel 'PSD [dB/Hz]'\n");
    fprintf(gnuplotPipe, "set grid xtics ytics dashtype 2 lc 'gray'\n");
    fprintf(gnuplotPipe, "set key top right\n");
    fprintf(gnuplotPipe, "set datafile separator ','\n");
    fprintf(gnuplotPipe, "set xrange [0:%f]\n", FS / 2.0);
    fprintf(gnuplotPipe, "set yrange [-60:20]\n");

    fprintf(gnuplotPipe,
        "plot '%s' using ($1 * %f / %d):(10*log10($2)) skip 1 with lines lw 1.0 lc 'light-blue' title 'Periodogram (Rect Window)', \\\n"
        "     ''    using ($1 * %f / %d):(10*log10($3)) skip 1 with lines lw 1.5 lc 'dark-orange' title 'Periodogram (Hamming)', \\\n"
        "     ''    using ($1 * %f / %d):(10*log10($4)) skip 1 with lines lw 2.0 lc 'dark-violet' title 'Multitaper Method', \\\n"
        "     ''    using ($5 * %f / %d):(10*log10($6)) skip 1 with lines lw 2.5 lc 'red' title 'Welch Method (Hamming)'\n",
        RESULTS_FILENAME, 
        FS, NFFT_RES,
        FS, NFFT_RES,
        FS, NFFT_RES,
        FS, NFFT_WE
    );
    
    pclose(gnuplotPipe);
    printf("C/Gnuplot: PSD comparison plot generated (Window 2).\n");
}

int main() {
    plot_statistics();
    plot_histograms();
    plot_psd_comparison();
    return 0;
}