#include "bayesm.h"
 
//The functions below are used to print the output from MCMC draws for many of the bayesm functions

time_t itime;
char buf[100];

void startMcmcTimer() {
    itime = time(NULL);
    Rcout << " MCMC Iteration (est time to end - min) \n";
}

void infoMcmcTimer(int rep, int R) {
    time_t ctime = time(NULL);    
    char buf[32];
    
    double timetoend = difftime(ctime, itime) / 60.0 * (R - rep - 1) / (rep+1);
//  sprintf(buf, " %d (%.1f)\n", rep+1, timetoend);  changed 11/30/22 P. Rossi
    snprintf(buf,32, " %d (%.1f)\n", rep+1, timetoend);
    Rcout <<  buf;
}

void endMcmcTimer() {
    time_t ctime = time(NULL);
    char buf[32];

 //   sprintf(buf, " Total Time Elapsed: %.2f \n", difftime(ctime, itime) / 60.0); 
 //     changed P. Rossi 11/30/2022
    snprintf(buf,32, " Total Time Elapsed: %.2f \n", difftime(ctime, itime) / 60.0);  
    Rcout << buf;

    itime = 0;
}
