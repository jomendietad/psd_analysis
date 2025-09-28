#define _GNU_SOURCE
#include <math.h>
#include <locale.h>
#include <fftw3.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <pthread.h>

typedef struct { double mean, variance, skewness, kurtosis; } Stats;

Stats calculate_statistics(const double* data, int n) {
    Stats result = {0.0, 0.0, 0.0, 0.0}; if (n < 4) return result;
    double sum = 0.0; for (int i = 0; i < n; ++i) sum += data[i]; result.mean = sum / n;
    double s2=0, s3=0, s4=0;
    for (int i = 0; i < n; ++i) { double d=data[i]-result.mean; s2+=d*d; s3+=d*d*d; s4+=d*d*d*d; }
    result.variance = s2/n;
    if (result.variance > 1e-30) {
        double std_dev = sqrt(result.variance);
        result.skewness = (s3/n)/(std_dev*std_dev*std_dev);
        result.kurtosis = (s4/n)/(result.variance*result.variance) - 3.0;
    }
    return result;
}

void generate_histogram(FILE* file, const double* data, int n, double sigma, int nbins) {
    if (n==0) return;
    double min_v=data[0], max_v=data[0];
    for (int i=1; i<n; ++i) { if(data[i]<min_v) min_v=data[i]; if(data[i]>max_v) max_v=data[i]; }
    if (fabs(max_v-min_v) < 1e-30) { fprintf(file, "%.8g,%.17g,%d\n", sigma, min_v, n); return; }
    double width=(max_v-min_v)/nbins; int *bins=calloc(nbins,sizeof(int));
    for (int i=0; i<n; ++i) { int idx=(int)((data[i]-min_v)/width); if(idx>=nbins)idx=nbins-1; bins[idx]++; }
    for (int i=0; i<nbins; ++i) fprintf(file,"%.8g,%.17g,%d\n",sigma,min_v+(i+0.5)*width,bins[i]);
    free(bins);
}

static inline double rand_u(){return(rand()+1.0)/(RAND_MAX+2.0);}
static double randn_bm(){return sqrt(-2.0*log(rand_u()+1e-18))*cos(2.0*M_PI*rand_u());}
static void make_hamming(double *w, int N){for(int n=0;n<N;++n)w[n]=0.54-0.46*cos(2.0*M_PI*n/(N-1));}
static double win_pow(const double *w,int N){double s2=0;for(int n=0;n<N;++n)s2+=w[n]*w[n];return s2/N;}
void gen_sine_tapers(int N, int K, double **t){for(int k=0;k<K;k++){double norm=sqrt(2.0/(N+1)); for(int i=0;i<N;i++)t[k][i]=norm*sin(M_PI*(k+1)*(i+1)/(N+1));}}

void generate_psd_data_double(double fs, int T, int M, double sigma, double A1, double f1, double A2, double f2, int Nseg, int overlap, double NW_mt) {
    struct timespec start, end;
    double time_p = 0, time_ph = 0, time_m = 0, time_w = 0;
    fprintf(stderr, "Generating 'data/results_double.csv'...\n");
    const int Nfft=T, Nfft_we=Nseg, K_mt=(int)(2*NW_mt)-1, hop=Nseg-overlap, segs=(T-Nseg)/hop+1;
    double *x=malloc(T*sizeof(*x)), *wT_h=malloc(T*sizeof(*wT_h)), *wS_w=malloc(Nseg*sizeof(*wS_w));
    double **s_t=malloc(K_mt*sizeof(double*)); for(int i=0;i<K_mt;++i){s_t[i]=malloc(T*sizeof(double));}
    fftw_complex *in_f=fftw_malloc(Nfft*sizeof(*in_f)), *out_f=fftw_malloc(Nfft*sizeof(*out_f)), *in_w=fftw_malloc(Nfft_we*sizeof(*in_w)), *out_w=fftw_malloc(Nfft_we*sizeof(*out_w));
    fftw_plan p_f=fftw_plan_dft_1d(Nfft,in_f,out_f,FFTW_FORWARD,FFTW_ESTIMATE), p_w=fftw_plan_dft_1d(Nfft_we,in_w,out_w,FFTW_FORWARD,FFTW_ESTIMATE);
    make_hamming(wT_h, T); gen_sine_tapers(T, K_mt, s_t); make_hamming(wS_w, Nseg);
    double U_h=win_pow(wT_h,T), U_w=win_pow(wS_w,Nseg);
    double *P=calloc(Nfft,sizeof(*P)), *Ph=calloc(Nfft,sizeof(*Ph)), *MT=calloc(Nfft,sizeof(*MT)), *W=calloc(Nfft_we,sizeof(*W));

    for(int m=0;m<M;++m){
        for(int n=0;n<T;++n) x[n]=A1*cos(2*M_PI*f1*n/fs)+A2*cos(2*M_PI*f2*n/fs)+sigma*randn_bm();
        
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
        for(int i=0;i<T;++i){in_f[i][0]=x[i];in_f[i][1]=0;} fftw_execute(p_f); for(int i=0;i<Nfft;++i)P[i]+=out_f[i][0]*out_f[i][0]+out_f[i][1]*out_f[i][1];
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end); time_p += (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
        for(int i=0;i<T;++i){in_f[i][0]=x[i]*wT_h[i];in_f[i][1]=0;} fftw_execute(p_f); for(int i=0;i<Nfft;++i)Ph[i]+=out_f[i][0]*out_f[i][0]+out_f[i][1]*out_f[i][1];
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end); time_ph += (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
        double *k_s=calloc(Nfft,sizeof(*k_s)); for(int k=0;k<K_mt;++k){for(int i=0;i<T;++i){in_f[i][0]=x[i]*s_t[k][i];in_f[i][1]=0;} fftw_execute(p_f); for(int i=0;i<Nfft;++i)k_s[i]+=out_f[i][0]*out_f[i][0]+out_f[i][1]*out_f[i][1];} for(int i=0;i<Nfft;++i)MT[i]+=k_s[i]/K_mt; free(k_s);
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end); time_m += (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
        for(int i=0;i<segs;++i){int start=i*hop; for(int n=0;n<Nseg;++n){in_w[n][0]=x[start+n]*wS_w[n];in_w[n][1]=0;} fftw_execute(p_w); for(int k=0;k<Nfft_we;++k)W[k]+=out_w[k][0]*out_w[k][0]+out_w[k][1]*out_w[k][1];}
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end); time_w += (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    }
    
    for(int k=0;k<Nfft;++k){P[k]/=(fs*T*M); Ph[k]/=(fs*T*U_h*M); MT[k]/=(fs*M);}
    for(int k=0;k<Nfft_we;++k) W[k]/=(fs*Nseg*U_w*segs*M);

    FILE *res=fopen("data/results_double.csv","w"), *time_log=fopen("data/timing_log_d.csv","w");
    fprintf(res,"k_res,PSD_periodogram,PSD_periodogram_hamm,PSD_mt,k_we,PSD_welch\n");
    for(int i=0;i<Nfft;++i) fprintf(res,"%d,%.15g,%.15g,%.15g,%d,%.15g\n",i,P[i],Ph[i],MT[i],(i<Nfft_we)?i:-1,(i<Nfft_we)?W[i]:NAN);
    fprintf(time_log,"method,total_cpu_time_s,avg_cpu_time_per_realization_ms\nPeriodogram,%.6f,%.6f\nPeriodogram (Hamming),%.6f,%.6f\nMultitaper,%.6f,%.6f\nWelch,%.6f,%.6f\n",time_p,time_p/M*1e3,time_ph,time_ph/M*1e3,time_m,time_m/M*1e3,time_w,time_w/M*1e3);
    fclose(res); fclose(time_log);

    free(x);free(wT_h);free(wS_w);for(int i=0;i<K_mt;++i){free(s_t[i]);} free(s_t); free(P);free(Ph);free(MT);free(W);
    fftw_destroy_plan(p_f);fftw_destroy_plan(p_w);fftw_free(in_f);fftw_free(out_f);fftw_free(in_w);fftw_free(out_w);
}

int main() {
    setlocale(LC_NUMERIC,"C");
    const double fs=1000, sigma_comp=3.5, f1=120, A1=1, f2=250, A2=0.7, NW_mt=4, f0_stats=250, A_stats=1;
    const int T=1024, M_comp=100, M_stats=1000, Nseg=256, overlap=128;
    srand(time(NULL)); fftw_init_threads(); fftw_plan_with_nthreads(4);
    
    generate_psd_data_double(fs, T, M_comp, sigma_comp, A1, f1, A2, f2, Nseg, overlap, NW_mt);
    
    fprintf(stderr, "\nStarting statistical analysis (double)...\n");
    const double sigmas[]={0.1,0.5,1.0,3.5,10.0,25.0,50.0}; const int num_s=sizeof(sigmas)/sizeof(sigmas[0]);
    const int kf0_w=(int)(f0_stats*Nseg/fs), hop=Nseg-overlap, segs=(T-Nseg)/hop+1;
    double *wS=malloc(Nseg*sizeof(*wS)); make_hamming(wS,Nseg); const double Useg=win_pow(wS,Nseg), norm_w=fs*Nseg*Useg*segs;
    FILE *sf=fopen("data/statistics_double.csv","w"),*hf=fopen("data/histograms_double.csv","w");
    fprintf(sf,"sigma,mean_psd,variance_psd,skewness,kurtosis\n"); fprintf(hf,"sigma,bin_center,count\n");
    fftw_complex *in_w=fftw_malloc(Nseg*sizeof(*in_w)), *out_w=fftw_malloc(Nseg*sizeof(*out_w));
    fftw_plan p_w=fftw_plan_dft_1d(Nseg,in_w,out_w,FFTW_FORWARD,FFTW_ESTIMATE);

    for(int i=0;i<num_s;++i){
        fprintf(stderr,"Processing sigma = %.4f ...\n",sigmas[i]);
        double *psd_vals=malloc(M_stats*sizeof(*psd_vals)), *x=malloc(T*sizeof(*x));
        for(int m=0;m<M_stats;++m){
            for(int n=0;n<T;++n) x[n]=A_stats*cos(2*M_PI*f0_stats*n/fs)+sigmas[i]*randn_bm();
            double Pxx_sum=0;
            for(int s=0;s<segs;++s){int start=s*hop; for(int n=0;n<Nseg;++n){in_w[n][0]=x[start+n]*wS[n];in_w[n][1]=0;} fftw_execute(p_w); Pxx_sum+=out_w[kf0_w][0]*out_w[kf0_w][0]+out_w[kf0_w][1]*out_w[kf0_w][1];}
            psd_vals[m]=(2.0*Pxx_sum)/norm_w;
        }
        Stats stats=calculate_statistics(psd_vals,M_stats);
        fprintf(sf,"%.8g,%.17g,%.17g,%.17g,%.17g\n",sigmas[i],stats.mean,stats.variance,stats.skewness,stats.kurtosis);
        generate_histogram(hf,psd_vals,M_stats,sigmas[i],40);
        free(psd_vals); free(x);
    }
    fftw_destroy_plan(p_w);fftw_free(in_w);fftw_free(out_w);
    fclose(sf); fclose(hf); free(wS);
    
    fftw_cleanup_threads();
    return 0;
}