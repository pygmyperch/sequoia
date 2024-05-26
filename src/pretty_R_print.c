#include <R.h>


/* 
* print nicely formatted messages & data from Fortran to R console
*/


/* 
* print timestamp, iteration number, step, no. assigned dams + sires, total likelihood 
* in a table-like format. Print time + iter + step at start of round/step, 
* and no. assigned dams & sires + total likelihood at end of round/step.
*/
void F77_SUB(rprint_status_tbl_header)() {
    Rprintf("%8s | %2s | %10s | %10s | %5s | %5s | %10s \n",
            "Time", "R", "Step", "progress", "dams", "sires", "Total LL"); 
    Rprintf("-------- | -- | ---------- | ---------- | ----- | ----- | ----------\n");
}



void F77_SUB(rprint_status_tbl_entry)(int *time, int *iter, char* step, int *n_parents, double *total_LL) {
    Rprintf("%02d:%02d:%02d | %2d | %.10s |            | %5d | %5d | %10.1f \n", 
            time[0], time[1], time[2], *iter, step, n_parents[0], n_parents[1], *total_LL);
}


void F77_SUB(rprint_status_tbl_a)(int *time, int *iter, char* step) { 
    Rprintf("%02d:%02d:%02d | %2d | %.10s | ", time[0], time[1], time[2], *iter, step);
}

void F77_SUB(rprint_status_tbl_b)(int *n_parents, double *total_LL) {
    Rprintf(" | %5d | %5d | %10.1f \n", n_parents[0], n_parents[1], *total_LL);
}

/* in the for loop of each step, print progress dots at every 10% */
void F77_SUB(rprint_status_tbl_dot)() {
    Rprintf(".");
}
 
void F77_SUB(rprint_status_tbl_no_dots)() {
    Rprintf("          ");
}

void F77_SUB(rprint_status_tbl_eol)() {
    Rprintf(" | \n");
}
 


/* progress bar for GetMaybeRel() etc */
void F77_SUB(rprint_progbar_header)() {
    Rprintf("\n");
    Rprintf(" 0   10  20  30  40  50  60  70  80  90  100%% \n");
    Rprintf(" |   |   |   |   |   |   |   |   |   |   |\n  ");
}

void F77_SUB(rprint_progbar_dot)() {
    Rprintf("*");
}

void F77_SUB(rprint_eol)() {
    Rprintf("\n");
}