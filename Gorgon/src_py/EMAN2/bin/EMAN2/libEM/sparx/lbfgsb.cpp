#include <cstdio>
#include <cstring>
#include <cmath>
#include <cstdlib>

#define TRUE_ (1)
#define FALSE_ (0)
#define min_(a,b) ((a) <= (b) ? (a) : (b))
#define max_(a,b) ((a) >= (b) ? (a) : (b))
#define abs_(x) ((x) >= 0 ? (x) : -(x))

typedef struct
{	long int cierr;
	long int ciunit;
	long int ciend;
	char *cifmt;
	long int cirec;
} cilist;

typedef struct
{	long int oerr;
	long int ounit;
	char *ofnm;
	long int ofnmlen;
	char *osta;
	char *oacc;
	char *ofm;
	long int orl;
	char *oblnk;
} olist;

/* Table of constant values */

static double c_b9 = 0.;
static long int c__1 = 1;
//static long int c__9 = 9;
static long int c__11 = 11;
//static long int c__3 = 3;
static double c_b275 = .001;
static double c_b276 = .9;
static double c_b277 = .1;
//static long int c__5 = 5;

const int SIXTY=60;

// function prototype for L-BFGS-B routine
// NOTES: All arguments must be passed as pointers. Must add an underscore
//  to the function name. 'extern "C"' keeps c++ from mangling the function
//  name. Be aware that fortran indexes 1d arrays [1,n] not [0,n-1], and
//  expects 2d arrays to be in col-major order not row-major order.
/*int setulb_(long int *n, long int *m, double *x, double *l, double *u, long int *nbd, double *f, double *g, 
		   double *factr, double *pgtol, double *wa, long int *iwa, char *task, long int *iprint, char *csave, long int *lsave,
		   long int *isave, double *dsave, long int aaa1, long int aaa2);	   */
	   
//int MAIN__() { abort(); return 1; }

/*int setulb_(long int *, long int *, double *, double *, double *, long int *, double *, double *, 
		   double *, double *, double *, long int *, char *, long int *, char *, long int *,
		   long int *, double *, long int , long int);*/



/*int main() {
  
  //c     This simple driver demonstrates how to call the L-BFGS-B code to
  //c       solve a sample problem (the extended Rosenbrock function 
  //c       subject to bounds on the variables). The dimension n of this
  //c       problem is variable.


  int nmax=1024, mmax=17;
  //c        nmax is the dimension of the largest problem to be solved.
  //c        mmax is the maximum number of limited memory corrections.
 
  //c     Declare the variables needed by the code.
  //c       A description of all these variables is given at the end of 
  //c       the driver.
  char task[SIXTY], csave[SIXTY];
  long int lsave[4];
  long int n, m, iprint, nbd[nmax], iwa[3*nmax], isave[44];
  double f, factr, pgtol, x[nmax], l[nmax], u[nmax], g[nmax], dsave[29], 
    wa[2*mmax*nmax+4*nmax+12*mmax*mmax+12*mmax];

  //c     Declare a few additional variables for this sample problem.
  double t1, t2;
  int i;
 
  //c     We wish to have output at every iteration.
  iprint = -1;

  //c     We specify the tolerances in the stopping criteria.
  factr=1.0e1;
  pgtol=1.0e-9;

  //c     We specify the dimension n of the sample problem and the number
  //c        m of limited memory corrections stored.  (n and m should not
  //c        exceed the limits nmax and mmax respectively.)
  n=25;
  m=5;
 
  //c     We now provide nbd which defines the bounds on the variables:
  //c                    l   specifies the lower bounds,
  //c                    u   specifies the upper bounds. 
 
  //c     First set bounds on the odd-numbered variables. (*even-numbered in c*)
  for (i=0; i<n; i+=2){
    nbd[i]=2;
    l[i]=1.0;
    u[i]=100.0;
  }

  //c     Next set bounds on the even-numbered variables. (*odd-numbered in c*)
  for (i=1; i<n; i+=2){
    nbd[i]=2;
    l[i]=-100.0;
    u[i]=100.0;
  }

  //c     We now define the starting point.
  for (i=0; i<n; i++)
    x[i]=3.0;
 
  printf("Solving sample problem.\n");
  printf(" (f=0.0 at the optimal solution.)\n");


  //c     We start the iteration by initializing task.
  // (**MUST clear remaining chars in task with spaces (else crash)!**)
  strcpy(task,"START");
  for (i=5;i<SIXTY;i++)
    task[i]=' ';

  //c     This is the call to the L-BFGS-B code.
  // (* call the L-BFGS-B routine with task='START' once before loop *)
  setulb_(&n,&m,x,l,u,nbd,&f,g,&factr,&pgtol,wa,iwa,task,&iprint,
	 csave,lsave,isave,dsave,(long int)60,(long int)60);

  // (* while routine returns "FG" or "NEW_X" in task, keep calling it *)
  while (strncmp(task,"FG",2)==0 || strncmp(task,"NEW_X",5)==0) {

    if (strncmp(task,"FG",2)==0) {
      //c   the minimization routine has returned to request the
      //c   function f and gradient g values at the current x

      //c        Compute function value f for the sample problem.

      f = 0.25 * (x[0]-1.0)*(x[0]-1.0);
      for (i=1; i<n; i++)
	f += pow(x[i] - x[i-1]*x[i-1], 2);
      f *= 4.0;

      //c        Compute gradient g for the sample problem.
      //t1 = (x[1]-x[0])*(x[1]-x[0]);
      t1 = x[1]-x[0]*x[0];
      g[0] = 2.0*(x[0]-1.0) - 16.0*x[0]*t1;
      for (i=1; i<n-1; i++) {
	t2=t1;
	//t1=(x[i+1]-x[i])*(x[i+1]-x[i]);
	t1 = x[i+1]-x[i]*x[i];
	g[i]=8.0*t2-16.0*x[i]*t1;
      }
      g[n-1]=8.0*t1;

    }
    else {
      // the minimization routine has returned with a new iterate,
      // and we have opted to continue the iteration (do nothing here)
    }

    //c          go back to the minimization routine.
    setulb_(&n,&m,x,l,u,nbd,&f,g,&factr,&pgtol,wa,iwa,task,&iprint,
	   csave,lsave,isave,dsave,(long int)60,(long int)60);
  }

  for (i=0; i<25; i++) { printf("i=%2d     x[i]=%f\n", i, x[i]); }
  //c     If task is neither FG nor NEW_X we terminate execution.
  // (* the minimization routine has returned with one of
  //  {CONV, ABNO, ERROR} *)

  return 0;
}*/

long int s_cmp(char *str1, const char *const str2, long int l1, long int l2) {
	long int l;
	int i = 0;
	long int flag = 0;
	
	l = min_(l1,l2);
	
	while (i<l && flag==0) {
		if (str1[i] != str2[i]) { flag = 1; }
		i++;
	}
	return flag;
}

int s_copy(char *str1, const char *const str2, long int l1, long int l2) {
	long int l;
	int i = 0;
	l = min_(l1,l2);
	
	while (i<l) {
		str1[i]=str2[i];
		i++;
	}	
	return 0;
}

/*******************************************************************
c     --------------------------------------------------------------
c             DESCRIPTION OF THE VARIABLES IN L-BFGS-B
c     --------------------------------------------------------------
c
c     n is an long int variable that must be set by the user to the
c       number of variables.  It is not altered by the routine.
c
c     m is an long int variable that must be set by the user to the
c       number of corrections used in the limited memory matrix.
c       It is not altered by the routine.  Values of m < 3  are
c       not recommended, and large values of m can result in excessive
c       computing time. The range  3 <= m <= 20 is recommended. 
c
c     x is a DOUBLE PRECISION array of length n.  On initial entry
c       it must be set by the user to the values of the initial
c       estimate of the solution vector.  Upon successful exit, it
c       contains the values of the variables at the best point
c       found (usually an approximate solution).
c
c     l is a DOUBLE PRECISION array of length n that must be set by
c       the user to the values of the lower bounds on the variables. If
c       the i-th variable has no lower bound, l(i) need not be defined.
c
c     u is a DOUBLE PRECISION array of length n that must be set by
c       the user to the values of the upper bounds on the variables. If
c       the i-th variable has no upper bound, u(i) need not be defined.
c
c     nbd is an long int array of dimension n that must be set by the
c       user to the type of bounds imposed on the variables:
c       nbd(i)=0 if x(i) is unbounded,
c              1 if x(i) has only a lower bound,
c              2 if x(i) has both lower and upper bounds, 
c              3 if x(i) has only an upper bound.
c
c     f is a DOUBLE PRECISION variable.  If the routine setulb returns
c       with task(1:2)= 'FG', then f must be set by the user to
c       contain the value of the function at the point x.
c
c     g is a DOUBLE PRECISION array of length n.  If the routine setulb
c       returns with taskb(1:2)= 'FG', then g must be set by the user to
c       contain the components of the gradient at the point x.
c
c     factr is a DOUBLE PRECISION variable that must be set by the user.
c       It is a tolerance in the termination test for the algorithm.
c       The iteration will stop when
c
c        (f^k - f^{k+1})/max{|f^k|,|f^{k+1}|,1} <= factr*epsmch
c
c       where epsmch is the machine precision which is automatically
c       generated by the code. Typical values for factr on a computer
c       with 15 digits of accuracy in double precision are:
c       factr=1.d+12 for low accuracy;
c             1.d+7  for moderate accuracy; 
c             1.d+1  for extremely high accuracy.
c       The user can suppress this termination test by setting factr=0.
c
c     pgtol is a double precision variable.
c       On entry pgtol >= 0 is specified by the user.  The iteration
c         will stop when
c
c                 max{|proj g_i | i = 1, ..., n} <= pgtol
c
c         where pg_i is the ith component of the projected gradient.
c       The user can suppress this termination test by setting pgtol=0.
c
c     wa is a DOUBLE PRECISION  array of length 
c       (2mmax + 4)nmax + 12mmax^2 + 12mmax used as workspace.
c       This array must not be altered by the user.
c
c     iwa is an long int  array of length 3nmax used as
c       workspace. This array must not be altered by the user.
c
c     task is a CHARACTER string of length 60.
c       On first entry, it must be set to 'START'.
c       On a return with task(1:2)='FG', the user must evaluate the
c         function f and gradient g at the returned value of x.
c       On a return with task(1:5)='NEW_X', an iteration of the
c         algorithm has concluded, and f and g contain f(x) and g(x)
c         respectively.  The user can decide whether to continue or stop
c         the iteration. 
c       When
c         task(1:4)='CONV', the termination test in L-BFGS-B has been 
c           satisfied;
c         task(1:4)='ABNO', the routine has terminated abnormally
c           without being able to satisfy the termination conditions,
c           x contains the best approximation found,
c           f and g contain f(x) and g(x) respectively;
c         task(1:5)='ERROR', the routine has detected an error in the
c           input parameters;
c       On exit with task = 'CONV', 'ABNO' or 'ERROR', the variable task
c         contains additional information that the user can print.
c       This array should not be altered unless the user wants to
c          stop the run for some reason.  See driver2 or driver3
c          for a detailed explanation on how to stop the run 
c          by assigning task(1:4)='STOP' in the driver.
c
c     iprint is an long int variable that must be set by the user.
c       It controls the frequency and type of output generated:
c        iprint<0    no output is generated;
c        iprint=0    print only one line at the last iteration;
c        0<iprint<99 print also f and |proj g| every iprint iterations;
c        iprint=99   print details of every iteration except n-vectors;
c        iprint=100  print also the changes of active set and final x;
c        iprint>100  print details of every iteration including x and g;
c       When iprint > 0, the file iterate.dat will be created to
c                        summarize the iteration.
c
c     csave  is a CHARACTER working array of length 60.
c
c     lsave is a long int working array of dimension 4.
c       On exit with task = 'NEW_X', the following information is
c         available:
c       lsave(1) = .true.  the initial x did not satisfy the bounds;
c       lsave(2) = .true.  the problem contains bounds;
c       lsave(3) = .true.  each variable has upper and lower bounds.
c
c     isave is an long int working array of dimension 44.
c       On exit with task = 'NEW_X', it contains information that
c       the user may want to access:
c         isave(30) = the current iteration number;
c         isave(34) = the total number of function and gradient
c                         evaluations;
c         isave(36) = the number of function value or gradient
c                                  evaluations in the current iteration;
c         isave(38) = the number of free variables in the current
c                         iteration;
c         isave(39) = the number of active constraints at the current
c                         iteration;
c
c         see the subroutine setulb.f for a description of other 
c         information contained in isave
c
c     dsave is a DOUBLE PRECISION working array of dimension 29.
c       On exit with task = 'NEW_X', it contains information that
c         the user may want to access:
c         dsave(2) = the value of f at the previous iteration;
c         dsave(5) = the machine precision epsmch generated by the code;
c         dsave(13) = the infinity norm of the projected gradient;
c
c         see the subroutine setulb.f for a description of other 
c         information contained in dsave
c
c     --------------------------------------------------------------
c           END OF THE DESCRIPTION OF THE VARIABLES IN L-BFGS-B
c     --------------------------------------------------------------
c
c     << An example of subroutine 'timer' for AIX Version 3.2 >>
c
c     subroutine timer(ttime)
c     double precision ttime
c     long int itemp, long int mclock
c     itemp = mclock()
c     ttime = dble(itemp)*1.0d-2
c     return
c     end

********************************************************************/

/* Subroutine */ int setulb_(long int *n, long int *m, double *x, double *l, double *u, long int *nbd, double *f, double *g, 
			 double *factr, double *pgtol, double *wa, long int *iwa, char *task, long int *iprint, char *csave, long int *lsave,
			 long int *isave, double *dsave, long int task_len, long int csave_len)
/*long int *n, *m;
double *x, *l, *u;
long int *nbd;
double *f, *g, *factr, *pgtol, *wa;
long int *iwa;
char *task;
long int *iprint;
char *csave;
long int *lsave;
long int *isave;
double *dsave;
long int task_len;
long int csave_len;*/
{
    /* System generated locals */
    long int i__1;

    /* Builtin functions */
    long int s_cmp(char *, const char *const, long int, long int);

    /* Local variables */
    static long int lsnd, lsgo, lygo, l1, l2, l3, ld, lr, lt;
    extern /* Subroutine */ int mainlb_(long int *n, long int *m, double *x, double *l, double *u, long int *nbd, double *f, 
			 double *g, double *factr, double *pgtol, double *ws, double *wy,
			 double *sy, double *ss, double *yy, double *wt, double *wn, double *snd, 
			 double *z__, double *r__, double *d__, double *t, double *wa, double *sg, 
			 double *sgo, double *yg, double *ygo, long int *index, long int *iwhere, long int *indx2, 
			 char *task, long int *iprint, char *csave, long int *lsave, long int *isave, double *dsave, 
			 long int task_len, long int csave_len);
    static long int lz, lwa, lsg, lyg, lwn, lss, lws, lwt, lsy, lwy, lyy;

/*     ************ */

/*     Subroutine setulb */

/*     This subroutine partitions the working arrays wa and iwa, and */
/*       then uses the limited memory BFGS method to solve the bound */
/*       constrained optimization problem by calling mainlb. */
/*       (The direct method will be used in the subspace minimization.) */

/*     n is an long int variable. */
/*       On entry n is the dimension of the problem. */
/*       On exit n is unchanged. */

/*     m is an long int variable. */
/*       On entry m is the maximum number of variable metric corrections */
/*         used to define the limited memory matrix. */
/*       On exit m is unchanged. */

/*     x is a double precision array of dimension n. */
/*       On entry x is an approximation to the solution. */
/*       On exit x is the current approximation. */

/*     l is a double precision array of dimension n. */
/*       On entry l is the lower bound on x. */
/*       On exit l is unchanged. */

/*     u is a double precision array of dimension n. */
/*       On entry u is the upper bound on x. */
/*       On exit u is unchanged. */

/*     nbd is an long int array of dimension n. */
/*       On entry nbd represents the type of bounds imposed on the */
/*         variables, and must be specified as follows: */
/*         nbd(i)=0 if x(i) is unbounded, */
/*                1 if x(i) has only a lower bound, */
/*                2 if x(i) has both lower and upper bounds, and */
/*                3 if x(i) has only an upper bound. */
/*       On exit nbd is unchanged. */

/*     f is a double precision variable. */
/*       On first entry f is unspecified. */
/*       On final exit f is the value of the function at x. */

/*     g is a double precision array of dimension n. */
/*       On first entry g is unspecified. */
/*       On final exit g is the value of the gradient at x. */

/*     factr is a double precision variable. */
/*       On entry factr >= 0 is specified by the user.  The iteration */
/*         will stop when */

/*         (f^k - f^{k+1})/max{|f^k|,|f^{k+1}|,1} <= factr*epsmch */

/*         where epsmch is the machine precision, which is automatically */
/*         generated by the code. Typical values for factr: 1.d+12 for */
/*         low accuracy; 1.d+7 for moderate accuracy; 1.d+1 for extremely */
/*         high accuracy. */
/*       On exit factr is unchanged. */

/*     pgtol is a double precision variable. */
/*       On entry pgtol >= 0 is specified by the user.  The iteration */
/*         will stop when */

/*                 max{|proj g_i | i = 1, ..., n} <= pgtol */

/*         where pg_i is the ith component of the projected gradient. */
/*       On exit pgtol is unchanged. */

/*     wa is a double precision working array of length */
/*       (2mmax + 4)nmax + 12mmax^2 + 12mmax. */

/*     iwa is an long int working array of length 3nmax. */

/*     task is a working string of characters of length 60 indicating */
/*       the current job when entering and quitting this subroutine. */

/*     iprint is an long int variable that must be set by the user. */
/*       It controls the frequency and type of output generated: */
/*        iprint<0    no output is generated; */
/*        iprint=0    print only one line at the last iteration; */
/*        0<iprint<99 print also f and |proj g| every iprint iterations; */
/*        iprint=99   print details of every iteration except n-vectors; */
/*        iprint=100  print also the changes of active set and final x; */
/*        iprint>100  print details of every iteration including x and g; */
/*       When iprint > 0, the file iterate.dat will be created to */
/*                        summarize the iteration. */

/*     csave is a working string of characters of length 60. */

/*     lsave is a long int working array of dimension 4. */
/*       On exit with 'task' = NEW_X, the following information is */
/*                                                             available: */
/*         If lsave(1) = .true.  then  the initial X has been replaced by */
/*                                     its projection in the feasible set; */
/*         If lsave(2) = .true.  then  the problem is constrained; */
/*         If lsave(3) = .true.  then  each variable has upper and lower */
/*                                     bounds; */

/*     isave is an long int working array of dimension 44. */
/*       On exit with 'task' = NEW_X, the following information is */
/*                                                             available: */
/*         isave(22) = the total number of intervals explored in the */
/*                         search of Cauchy points; */
/*         isave(26) = the total number of skipped BFGS updates before */
/*                         the current iteration; */
/*         isave(30) = the number of current iteration; */
/*         isave(31) = the total number of BFGS updates prior the current */
/*                         iteration; */
/*         isave(33) = the number of intervals explored in the search of */
/*                         Cauchy point in the current iteration; */
/*         isave(34) = the total number of function and gradient */
/*                         evaluations; */
/*         isave(36) = the number of function value or gradient */
/*                                  evaluations in the current iteration; */
/*         if isave(37) = 0  then the subspace argmin is within the box; */
/*         if isave(37) = 1  then the subspace argmin is beyond the box; */
/*         isave(38) = the number of free variables in the current */
/*                         iteration; */
/*         isave(39) = the number of active constraints in the current */
/*                         iteration; */
/*         n + 1 - isave(40) = the number of variables leaving the set of */
/*                           active constraints in the current iteration; */
/*         isave(41) = the number of variables entering the set of active */
/*                         constraints in the current iteration. */

/*     dsave is a double precision working array of dimension 29. */
/*       On exit with 'task' = NEW_X, the following information is */
/*                                                             available: */
/*         dsave(1) = current 'theta' in the BFGS matrix; */
/*         dsave(2) = f(x) in the previous iteration; */
/*         dsave(3) = factr*epsmch; */
/*         dsave(4) = 2-norm of the line search direction vector; */
/*         dsave(5) = the machine precision epsmch generated by the code; */
/*         dsave(7) = the accumulated time spent on searching for */
/*                                                         Cauchy points; */
/*         dsave(8) = the accumulated time spent on */
/*                                                 subspace minimization; */
/*         dsave(9) = the accumulated time spent on line search; */
/*         dsave(11) = the slope of the line search function at */
/*                                  the current point of line search; */
/*         dsave(12) = the maximum relative step length imposed in */
/*                                                           line search; */
/*         dsave(13) = the infinity norm of the projected gradient; */
/*         dsave(14) = the relative step length in the line search; */
/*         dsave(15) = the slope of the line search function at */
/*                                 the starting point of the line search; */
/*         dsave(16) = the square of the 2-norm of the line search */
/*                                                      direction vector. */

/*     Subprograms called: */

/*       L-BFGS-B Library ... mainlb. */


/*     References: */

/*       [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited */
/*       memory algorithm for bound constrained optimization'', */
/*       SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208. */

/*       [2] C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, ``L-BFGS-B: a */
/*       limited memory FORTRAN code for solving bound constrained */
/*       optimization problems'', Tech. Report, NAM-11, EECS Department, */
/*       Northwestern University, 1994. */

/*       (Postscript files of these papers are available via anonymous */
/*        ftp to eecs.nwu.edu in the directory pub/lbfgs/lbfgs_bcm.) */

/*                           *  *  * */

/*     NEOS, November 1994. (Latest revision June 1996.) */
/*     Optimization Technology Center. */
/*     Argonne National Laboratory and Northwestern University. */
/*     Written by */
/*                        Ciyou Zhu */
/*     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal. */


/*     ************ */
    /* Parameter adjustments */
    --iwa;
    --g;
    --nbd;
    --u;
    --l;
    --x;
    --wa;
    --lsave;
    --isave;
    --dsave;

    /* Function Body */
    if (s_cmp(task, "START", (long int)60, (long int)5) == 0) {
	isave[1] = *m * *n;
/* Computing 2nd power */
	i__1 = *m;
	isave[2] = i__1 * i__1;
/* Computing 2nd power */
	i__1 = *m;
	isave[3] = i__1 * i__1 << 2;
	isave[4] = 1;
	isave[5] = isave[4] + isave[1];
	isave[6] = isave[5] + isave[1];
	isave[7] = isave[6] + isave[2];
	isave[8] = isave[7] + isave[2];
	isave[9] = isave[8] + isave[2];
	isave[10] = isave[9] + isave[2];
	isave[11] = isave[10] + isave[3];
	isave[12] = isave[11] + isave[3];
	isave[13] = isave[12] + *n;
	isave[14] = isave[13] + *n;
	isave[15] = isave[14] + *n;
	isave[16] = isave[15] + *n;
	isave[17] = isave[16] + (*m << 3);
	isave[18] = isave[17] + *m;
	isave[19] = isave[18] + *m;
	isave[20] = isave[19] + *m;
    }
    l1 = isave[1];
    l2 = isave[2];
    l3 = isave[3];
    lws = isave[4];
    lwy = isave[5];
    lsy = isave[6];
    lss = isave[7];
    lyy = isave[8];
    lwt = isave[9];
    lwn = isave[10];
    lsnd = isave[11];
    lz = isave[12];
    lr = isave[13];
    ld = isave[14];
    lt = isave[15];
    lwa = isave[16];
    lsg = isave[17];
    lsgo = isave[18];
    lyg = isave[19];
    lygo = isave[20];
    mainlb_(n, m, &x[1], &l[1], &u[1], &nbd[1], f, &g[1], factr, pgtol, &wa[
	    lws], &wa[lwy], &wa[lsy], &wa[lss], &wa[lyy], &wa[lwt], &wa[lwn], 
	    &wa[lsnd], &wa[lz], &wa[lr], &wa[ld], &wa[lt], &wa[lwa], &wa[lsg],
	     &wa[lsgo], &wa[lyg], &wa[lygo], &iwa[1], &iwa[*n + 1], &iwa[(*n 
	    << 1) + 1], task, iprint, csave, &lsave[1], &isave[22], &dsave[1],
	     (long int)60, (long int)60);
    return 0;
} /* setulb_ */

/* ======================= The end of setulb ============================= */
/* Subroutine */ int mainlb_(long int *n, long int *m, double *x, double *l, double *u, long int *nbd, double *f, 
			 double *g, double *factr, double *pgtol, double *ws, double *wy,
			 double *sy, double *ss, double *yy, double *wt, double *wn, double *snd, 
			 double *z__, double *r__, double *d__, double *t, double *wa, double *sg, 
			 double *sgo, double *yg, double *ygo, long int *index, long int *iwhere, long int *indx2, 
			 char *task, long int *iprint, char *csave, long int *lsave, long int *isave, double *dsave, 
			 long int task_len, long int csave_len)
/*long int *n, *m;
double *x, *l, *u;
long int *nbd;
double *f, *g, *factr, *pgtol, *ws, *wy, *sy, *ss, *yy, *wt, *wn, *snd, *
	z__, *r__, *d__, *t, *wa, *sg, *sgo, *yg, *ygo;
long int *index, *iwhere, *indx2;
char *task;
long int *iprint;
char *csave;
long int *lsave;
long int *isave;
double *dsave;
long int task_len;
long int csave_len;*/
{
    /* Format strings *//*
    static char fmt_1002[] = "(/,\002At iterate\002,i5,4x,\002f= \002,1p,d12\
.5,4x,\002|proj g|= \002,1p,d12.5)";
    static char fmt_1003[] = "(2(1x,i4),5x,\002-\002,5x,\002-\002,3x,\002\
-\002,5x,\002-\002,5x,\002-\002,8x,\002-\002,3x,1p,2(1x,d10.3))";
    static char fmt_1001[] = "(//,\002ITERATION \002,i5)";
    static char fmt_1005[] = "(/,\002 Singular triangular system detected\
;\002,/,\002   refresh the lbfgs memory and restart the iteration.\002)";
    static char fmt_1006[] = "(/,\002 Nonpositive definiteness in Cholesky f\
actorization in formk;\002,/,\002   refresh the lbfgs memory and restart the\
 iteration.\002)";
    static char fmt_1008[] = "(/,\002 Bad direction in the line search;\002,\
/,\002   refresh the lbfgs memory and restart the iteration.\002)";
    static char fmt_1004[] = "(\002  ys=\002,1p,e10.3,\002  -gs=\002,1p,e10.\
3,\002 BFGS update SKIPPED\002)";
    static char fmt_1007[] = "(/,\002 Nonpositive definiteness in Cholesky f\
actorization in formt;\002,/,\002   refresh the lbfgs memory and restart the\
 iteration.\002)";*/

    /* System generated locals */
    long int ws_dim1, ws_offset, wy_dim1, wy_offset, sy_dim1, sy_offset, 
	    ss_dim1, ss_offset, yy_dim1, yy_offset, wt_dim1, wt_offset, 
	    wn_dim1, wn_offset, snd_dim1, snd_offset, i__1;
    double d__1, d__2;
//    olist o__1;

    /* Builtin functions */
    long int s_cmp(char *, const char *const, long int, long int);
    /* Subroutine */ int s_copy(char *, const char *const, long int, long int);
    //long int f_open(olist *), s_wsfe(cilist *), do_fio(long int*, char*, long int), e_wsfe();

    /* Local variables */
    static long int head;
    static double fold;
    static long int nact;
    static double ddum;
    extern double ddot_(long int *n, double *dx, long int *incx, double *dy, long int *incy);
    static long int info;
    static double time;
    static long int nfgv, ifun, iter, nint;
    static char word[3];
    static double time1, time2;
    static long int i__, iback, k;
    extern /* Subroutine */ int dscal_(long int *n, double *da, double *dx, long int *incx);
    static double gdold;
    static long int nfree;
    static long int boxed;
    static long int itail;
    static double theta;
    extern /* Subroutine */ int freev_(long int *n, long int *nfree, long int *index, long int *nenter, long int *ileave, long int *indx2, long int *iwhere, 
	    		         long int *wrk, long int *updatd, long int *cnstnd, long int *iprint, long int *iter),
				 dcopy_(long int *n, double *dx, long int *incx, double *dy, long int *incy);
    static double dnorm;
    extern /* Subroutine */ int timer_(double *ttime), 
    			  formk_(long int *n, long int *nsub, long int *ind, long int *nenter, long int *ileave, long int *indx2, long int *iupdat, 
			         long int *updatd, double *wn, double *wn1, long int *m, double *ws, double *wy, double *sy, 
			 	 double *theta, long int *col, long int *head, long int *info);
    static long int nskip, iword;
    extern /* Subroutine */ int formt_(long int *m, double *wt, double *sy, double *ss, long int *col, double *theta, long int *info), 
    			   subsm_(long int *n, long int *m, long int *nsub, long int *ind, double *l, double *u, long int *nbd,
			          double *x, double *d__, double *ws, double *wy, double *theta, long int *col,
				  long int *head, long int *iword, double *wv, double *wn, long int *iprint, long int *info);
    static double xstep, stpmx;
    extern /* Subroutine */ int prn1lb_(long int*, long int*, double*, double*, double*, long int*, long int*, double*),
    			   prn2lb_(long int *n, double *x, double *f, double *g, long int *iprint, long int *itfile,
			 long int *iter, long int *nfgv, long int *nact, double *sbgnrm, long int *nint, char *word,
			 long int *iword, long int *iback, double *stp, double *xstep, long int word_len), 
			 prn3lb_(long int *n, double *x, double *f, char *task, long int *iprint, long int *info, 
			 long int *itfile, long int *iter, long int *nfgv, long int *nintol, long int *nskip, long int *nact,
			 double *sbgnrm, double *time, long int *nint, char *word, long int *iback, double *stp,
			 double *xstep, long int *k, double *cachyt, double *sbtime, double *lnscht,
			 long int task_len, long int word_len);
    static double gd, dr, rr;
    static long int ileave;
    extern /* Subroutine */ int errclb_(long int *n, long int *m, double *factr, double *l, double *u, long int *nbd, char *task,
			 long int *info, long int *k, long int task_len);
    static long int itfile;
    static double cachyt, epsmch;
    static long int updatd;
    static double sbtime;
    extern /* Subroutine */ int active_(long int *n, double *l, double *u, long int *nbd, double *x, long int *iwhere, 
			 long int *iprint, long int *prjctd, long int *cnstnd, long int *boxed);
    static long int prjctd;
    static long int iupdat;
    extern double dpmeps_();
    static long int cnstnd;
    static double sbgnrm;
    static long int nenter;
    static double lnscht;
    extern /* Subroutine */ int cauchy_(long int *n, double *x, double *l, double *u, long int *nbd, double *g, long int *iorder,
			 long int *iwhere, double *t, double *d__, double *xcp, long int *m, double *wy, double *ws,
			 double *sy, double *wt, double *theta, long int *col, long int *head, double *p, double *c__,
			 double *wbp, double *v, long int *nint, double *sg, double *yg, long int *iprint, double *sbgnrm, long int *info),
			  cmprlb_(long int *n, long int *m, double *x, double *g, double *ws, double *wy, double *sy,
			 double *wt, double *z__, double *r__, double *wa, long int *index, double *theta,
			 long int *col, long int *head, long int *nfree, long int *cnstnd, long int *info), 
			  lnsrlb_(long int *n, double *l, double *u, long int *nbd, double *x, double *f, double *fold,
			 double *gd, double *gdold, double *g, double *d__, double *r__, double *t, 
			 double *z__, double *stp, double *dnorm, double *dtd, double *xstep, double *stpmx,
			 long int *iter, long int *ifun, long int *iback, long int *nfgv, long int *info, char *task, long int *boxed, 
			 long int *cnstnd, char *csave, long int *isave, double *dsave, long int task_len, long int csave_len), 
			  matupd_(long int *n, long int *m, double *ws, double *wy, double *sy, double *ss,
			 double *d__, double *r__, long int *itail, long int *iupdat, long int *col, long int *head,
			 double *theta, double *rr, double *dr, double *stp, double *dtd);
    static long int nintol;
    extern /* Subroutine */ int projgr_(long int *n, double *l, double *u, long int *nbd, double *x, double *g, double *sbgnrm);
    static double dtd;
    static long int col;
    static double tol;
    static long int wrk;
    static double stp, cpu1, cpu2;

    /* Fortran I/O blocks */
/*    static cilist io___62 = { 0, 6, 0, fmt_1002, 0 };
    static cilist io___63 = { 0, 0, 0, fmt_1003, 0 };
    static cilist io___64 = { 0, 6, 0, fmt_1001, 0 };
    static cilist io___66 = { 0, 6, 0, fmt_1005, 0 };
    static cilist io___68 = { 0, 6, 0, fmt_1006, 0 };
    static cilist io___69 = { 0, 6, 0, fmt_1005, 0 };
    static cilist io___71 = { 0, 6, 0, fmt_1008, 0 };
    static cilist io___75 = { 0, 6, 0, fmt_1004, 0 };
    static cilist io___76 = { 0, 6, 0, fmt_1007, 0 };*/


/*     ************ */

/*     Subroutine mainlb */

/*     This subroutine solves bound constrained optimization problems by */
/*       using the compact formula of the limited memory BFGS updates. */

/*     n is an long int variable. */
/*       On entry n is the number of variables. */
/*       On exit n is unchanged. */

/*     m is an long int variable. */
/*       On entry m is the maximum number of variable metric */
/*          corrections allowed in the limited memory matrix. */
/*       On exit m is unchanged. */

/*     x is a double precision array of dimension n. */
/*       On entry x is an approximation to the solution. */
/*       On exit x is the current approximation. */

/*     l is a double precision array of dimension n. */
/*       On entry l is the lower bound of x. */
/*       On exit l is unchanged. */

/*     u is a double precision array of dimension n. */
/*       On entry u is the upper bound of x. */
/*       On exit u is unchanged. */

/*     nbd is an long int array of dimension n. */
/*       On entry nbd represents the type of bounds imposed on the */
/*         variables, and must be specified as follows: */
/*         nbd(i)=0 if x(i) is unbounded, */
/*                1 if x(i) has only a lower bound, */
/*                2 if x(i) has both lower and upper bounds, */
/*                3 if x(i) has only an upper bound. */
/*       On exit nbd is unchanged. */

/*     f is a double precision variable. */
/*       On first entry f is unspecified. */
/*       On final exit f is the value of the function at x. */

/*     g is a double precision array of dimension n. */
/*       On first entry g is unspecified. */
/*       On final exit g is the value of the gradient at x. */

/*     factr is a double precision variable. */
/*       On entry factr >= 0 is specified by the user.  The iteration */
/*         will stop when */

/*         (f^k - f^{k+1})/max{|f^k|,|f^{k+1}|,1} <= factr*epsmch */

/*         where epsmch is the machine precision, which is automatically */
/*         generated by the code. */
/*       On exit factr is unchanged. */

/*     pgtol is a double precision variable. */
/*       On entry pgtol >= 0 is specified by the user.  The iteration */
/*         will stop when */

/*                 max{|proj g_i | i = 1, ..., n} <= pgtol */

/*         where pg_i is the ith component of the projected gradient. */
/*       On exit pgtol is unchanged. */

/*     ws, wy, sy, and wt are double precision working arrays used to */
/*       store the following information defining the limited memory */
/*          BFGS matrix: */
/*          ws, of dimension n x m, stores S, the matrix of s-vectors; */
/*          wy, of dimension n x m, stores Y, the matrix of y-vectors; */
/*          sy, of dimension m x m, stores S'Y; */
/*          ss, of dimension m x m, stores S'S; */
/* 	   yy, of dimension m x m, stores Y'Y; */
/*          wt, of dimension m x m, stores the Cholesky factorization */
/*                                  of (theta*S'S+LD^(-1)L'); see eq. */
/*                                  (2.26) in [3]. */

/*     wn is a double precision working array of dimension 2m x 2m */
/*       used to store the LEL^T factorization of the indefinite matrix */
/*                 K = [-D -Y'ZZ'Y/theta     L_a'-R_z'  ] */
/*                     [L_a -R_z           theta*S'AA'S ] */

/*       where     E = [-I  0] */
/*                     [ 0  I] */

/*     snd is a double precision working array of dimension 2m x 2m */
/*       used to store the lower triangular part of */
/*                 N = [Y' ZZ'Y   L_a'+R_z'] */
/*                     [L_a +R_z  S'AA'S   ] */

/*     z(n),r(n),d(n),t(n),wa(8*m) are double precision working arrays. */
/*       z is used at different times to store the Cauchy point and */
/*       the Newton point. */

/*     sg(m),sgo(m),yg(m),ygo(m) are double precision working arrays. */

/*     index is an long int working array of dimension n. */
/*       In subroutine freev, index is used to store the free and fixed */
/*          variables at the Generalized Cauchy Point (GCP). */

/*     iwhere is an long int working array of dimension n used to record */
/*       the status of the vector x for GCP computation. */
/*       iwhere(i)=0 or -3 if x(i) is free and has bounds, */
/*                 1       if x(i) is fixed at l(i), and l(i) .ne. u(i) */
/*                 2       if x(i) is fixed at u(i), and u(i) .ne. l(i) */
/*                 3       if x(i) is always fixed, i.e.,  u(i)=x(i)=l(i) */
/*                -1       if x(i) is always free, i.e., no bounds on it. */

/*     indx2 is an long int working array of dimension n. */
/*       Within subroutine cauchy, indx2 corresponds to the array iorder. */
/*       In subroutine freev, a list of variables entering and leaving */
/*       the free set is stored in indx2, and it is passed on to */
/*       subroutine formk with this information. */

/*     task is a working string of characters of length 60 indicating */
/*       the current job when entering and leaving this subroutine. */

/*     iprint is an long int variable that must be set by the user. */
/*       It controls the frequency and type of output generated: */
/*        iprint<0    no output is generated; */
/*        iprint=0    print only one line at the last iteration; */
/*        0<iprint<99 print also f and |proj g| every iprint iterations; */
/*        iprint=99   print details of every iteration except n-vectors; */
/*        iprint=100  print also the changes of active set and final x; */
/*        iprint>100  print details of every iteration including x and g; */
/*       When iprint > 0, the file iterate.dat will be created to */
/*                        summarize the iteration. */

/*     csave is a working string of characters of length 60. */

/*     lsave is a long int working array of dimension 4. */

/*     isave is an long int working array of dimension 23. */

/*     dsave is a double precision working array of dimension 29. */


/*     Subprograms called */

/*       L-BFGS-B Library ... cauchy, subsm, lnsrlb, formk, */

/*        errclb, prn1lb, prn2lb, prn3lb, active, projgr, */

/*        freev, cmprlb, matupd, formt. */

/*       Minpack2 Library ... timer, dpmeps. */

/*       Linpack Library ... dcopy, ddot. */


/*     References: */

/*       [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited */
/*       memory algorithm for bound constrained optimization'', */
/*       SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208. */

/*       [2] C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, ``L-BFGS-B: FORTRAN */
/*       Subroutines for Large Scale Bound Constrained Optimization'' */
/*       Tech. Report, NAM-11, EECS Department, Northwestern University, */
/*       1994. */

/*       [3] R. Byrd, J. Nocedal and R. Schnabel "Representations of */
/*       Quasi-Newton Matrices and their use in Limited Memory Methods'', */
/*       Mathematical Programming 63 (1994), no. 4, pp. 129-156. */

/*       (Postscript files of these papers are available via anonymous */
/*        ftp to eecs.nwu.edu in the directory pub/lbfgs/lbfgs_bcm.) */

/*                           *  *  * */

/*     NEOS, November 1994. (Latest revision June 1996.) */
/*     Optimization Technology Center. */
/*     Argonne National Laboratory and Northwestern University. */
/*     Written by */
/*                        Ciyou Zhu */
/*     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal. */


/*     ************ */
    /* Parameter adjustments */
    --indx2;
    --iwhere;
    --index;
    --t;
    --d__;
    --r__;
    --z__;
    --g;
    --nbd;
    --u;
    --l;
    --x;
    --ygo;
    --yg;
    --sgo;
    --sg;
    --wa;
    snd_dim1 = 2 * *m;
    snd_offset = 1 + snd_dim1 * 1;
    snd -= snd_offset;
    wn_dim1 = 2 * *m;
    wn_offset = 1 + wn_dim1 * 1;
    wn -= wn_offset;
    wt_dim1 = *m;
    wt_offset = 1 + wt_dim1 * 1;
    wt -= wt_offset;
    yy_dim1 = *m;
    yy_offset = 1 + yy_dim1 * 1;
    yy -= yy_offset;
    ss_dim1 = *m;
    ss_offset = 1 + ss_dim1 * 1;
    ss -= ss_offset;
    sy_dim1 = *m;
    sy_offset = 1 + sy_dim1 * 1;
    sy -= sy_offset;
    wy_dim1 = *n;
    wy_offset = 1 + wy_dim1 * 1;
    wy -= wy_offset;
    ws_dim1 = *n;
    ws_offset = 1 + ws_dim1 * 1;
    ws -= ws_offset;
    --lsave;
    --isave;
    --dsave;


    /* Function Body */
    if (s_cmp(task, "START", (long int)60, (long int)5) == 0) {
	timer_(&time1);
/*        Generate the current machine precision. */
	epsmch = dpmeps_();
/*        Initialize counters and scalars when task='START'. */
/*           for the limited memory BFGS matrices: */
	col = 0;
	head = 1;
	theta = 1.;
	iupdat = 0;
	updatd = FALSE_;
/*           for operation counts: */
	iter = 0;
	nfgv = 0;
	nint = 0;
	nintol = 0;
	nskip = 0;
	nfree = *n;
/*           for stopping tolerance: */
	tol = *factr * epsmch;
/*           for measuring running time: */
	cachyt = 0.;
	sbtime = 0.;
	lnscht = 0.;
/*           'word' records the status of subspace solutions. */
	s_copy(word, "---", (long int)3, (long int)3);
/*           'info' records the termination information. */
	info = 0;
	if (*iprint >= 1) {
/*                                open a summary file 'iterate.dat' */
/*	    o__1.oerr = 0;
	    o__1.ounit = 8;
	    o__1.ofnmlen = 11;
	    o__1.ofnm = "iterate.dat";
	    o__1.orl = 0;
	    o__1.osta = "unknown";
	    o__1.oacc = 0;
	    o__1.ofm = 0;
	    o__1.oblnk = 0;
	    f_open(&o__1); */
	    itfile = 8;
	}
/*        Check the input arguments for errors. */
	errclb_(n, m, factr, &l[1], &u[1], &nbd[1], task, &info, &k, (long int)60);
	if (s_cmp(task, "ERROR", (long int)5, (long int)5) == 0) {
	    prn3lb_(n, &x[1], f, task, iprint, &info, &itfile, &iter, &nfgv, &
		    nintol, &nskip, &nact, &sbgnrm, &c_b9, &nint, word, &
		    iback, &stp, &xstep, &k, &cachyt, &sbtime, &lnscht, (
		    long int)60, (long int)3);
	    return 0;
	}
	prn1lb_(n, m, &l[1], &u[1], &x[1], iprint, &itfile, &epsmch);
/*        Initialize iwhere & project x onto the feasible set. */
	active_(n, &l[1], &u[1], &nbd[1], &x[1], &iwhere[1], iprint, &prjctd, 
		&cnstnd, &boxed);
/*        The end of the initialization. */
    } else {
/*          restore local variables. */
	prjctd = lsave[1];
	cnstnd = lsave[2];
	boxed = lsave[3];
	updatd = lsave[4];
	nintol = isave[1];
	itfile = isave[3];
	iback = isave[4];
	nskip = isave[5];
	head = isave[6];
	col = isave[7];
	itail = isave[8];
	iter = isave[9];
	iupdat = isave[10];
	nint = isave[12];
	nfgv = isave[13];
	info = isave[14];
	ifun = isave[15];
	iword = isave[16];
	nfree = isave[17];
	nact = isave[18];
	ileave = isave[19];
	nenter = isave[20];
	theta = dsave[1];
	fold = dsave[2];
	tol = dsave[3];
	dnorm = dsave[4];
	epsmch = dsave[5];
	cpu1 = dsave[6];
	cachyt = dsave[7];
	sbtime = dsave[8];
	lnscht = dsave[9];
	time1 = dsave[10];
	gd = dsave[11];
	stpmx = dsave[12];
	sbgnrm = dsave[13];
	stp = dsave[14];
	gdold = dsave[15];
	dtd = dsave[16];
/*        After returning from the driver go to the point where execution */
/*        is to resume. */
	if (s_cmp(task, "FG_LN", (long int)5, (long int)5) == 0) {
	    goto L666;
	}
	if (s_cmp(task, "NEW_X", (long int)5, (long int)5) == 0) {
	    goto L777;
	}
	if (s_cmp(task, "FG_ST", (long int)5, (long int)5) == 0) {
	    goto L111;
	}
	if (s_cmp(task, "STOP", (long int)4, (long int)4) == 0) {
	    if (s_cmp(task + 6, "CPU", (long int)3, (long int)3) == 0) {
/*                                          restore the previous iterate. */
		dcopy_(n, &t[1], &c__1, &x[1], &c__1);
		dcopy_(n, &r__[1], &c__1, &g[1], &c__1);
		*f = fold;
	    }
	    goto L999;
	}
    }
/*     Compute f0 and g0. */
    s_copy(task, "FG_START", (long int)60, (long int)8);
/*          return to the driver to calculate f and g; reenter at 111. */
    goto L1000;
L111:
    nfgv = 1;
/*     Compute the infinity norm of the (-) projected gradient. */
    projgr_(n, &l[1], &u[1], &nbd[1], &x[1], &g[1], &sbgnrm);
    if (*iprint >= 1) {
/*	s_wsfe(&io___62);
	do_fio(&c__1, (char *)&iter, (long int)sizeof(long int));
	do_fio(&c__1, (char *)&(*f), (long int)sizeof(double));
	do_fio(&c__1, (char *)&sbgnrm, (long int)sizeof(double));
	e_wsfe();
	io___63.ciunit = itfile;
	s_wsfe(&io___63);
	do_fio(&c__1, (char *)&iter, (long int)sizeof(long int));
	do_fio(&c__1, (char *)&nfgv, (long int)sizeof(long int));
	do_fio(&c__1, (char *)&sbgnrm, (long int)sizeof(double));
	do_fio(&c__1, (char *)&(*f), (long int)sizeof(double));
	e_wsfe();*/
    }
    if (sbgnrm <= *pgtol) {
/*                                terminate the algorithm. */
	s_copy(task, "CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL", (
		long int)60, (long int)48);
	goto L999;
    }
/* ----------------- the beginning of the loop -------------------------- */
L222:
    if (*iprint >= 99) {
/*	s_wsfe(&io___64);
	i__1 = iter + 1;
	do_fio(&c__1, (char *)&i__1, (long int)sizeof(long int));
	e_wsfe();*/
    }
    iword = -1;

    if (! cnstnd && col > 0) {
/*                                            skip the search for GCP. */
	dcopy_(n, &x[1], &c__1, &z__[1], &c__1);
	wrk = updatd;
	nint = 0;
	goto L333;
    }
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */

/*     Compute the Generalized Cauchy Point (GCP). */

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
    timer_(&cpu1);
    cauchy_(n, &x[1], &l[1], &u[1], &nbd[1], &g[1], &indx2[1], &iwhere[1], &t[
	    1], &d__[1], &z__[1], m, &wy[wy_offset], &ws[ws_offset], &sy[
	    sy_offset], &wt[wt_offset], &theta, &col, &head, &wa[1], &wa[(*m 
	    << 1) + 1], &wa[(*m << 2) + 1], &wa[*m * 6 + 1], &nint, &sg[1], &
	    yg[1], iprint, &sbgnrm, &info);
    if (info != 0) {
/*         singular triangular system detected; refresh the lbfgs memory. */
	if (*iprint >= 1) {
/*	    s_wsfe(&io___66);
	    e_wsfe();*/
	}
	info = 0;
	col = 0;
	head = 1;
	theta = 1.;
	iupdat = 0;
	updatd = FALSE_;
	timer_(&cpu2);
	cachyt = cachyt + cpu2 - cpu1;
	goto L222;
    }
    timer_(&cpu2);
    cachyt = cachyt + cpu2 - cpu1;
    nintol += nint;
/*     Count the entering and leaving variables for iter > 0; */
/*     find the index set of free and active variables at the GCP. */
    freev_(n, &nfree, &index[1], &nenter, &ileave, &indx2[1], &iwhere[1], &
	    wrk, &updatd, &cnstnd, iprint, &iter);
    nact = *n - nfree;
L333:
/*     If there are no free variables or B=I, then */
/*                                        skip the subspace minimization. */
    if (nfree == 0 || col == 0) {
	goto L555;
    }
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */

/*     Subspace minimization. */

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
    timer_(&cpu1);
/*     Form  the LEL^T factorization of the indefinite */
/*       matrix    K = [-D -Y'ZZ'Y/theta     L_a'-R_z'  ] */
/*                     [L_a -R_z           theta*S'AA'S ] */
/*       where     E = [-I  0] */
/*                     [ 0  I] */
    if (wrk) {
	formk_(n, &nfree, &index[1], &nenter, &ileave, &indx2[1], &iupdat, &
		updatd, &wn[wn_offset], &snd[snd_offset], m, &ws[ws_offset], &
		wy[wy_offset], &sy[sy_offset], &theta, &col, &head, &info);
    }
    if (info != 0) {
/*          nonpositive definiteness in Cholesky factorization; */
/*          refresh the lbfgs memory and restart the iteration. */
	if (*iprint >= 1) {
	   /* s_wsfe(&io___68);
	    e_wsfe();*/
	}
	info = 0;
	col = 0;
	head = 1;
	theta = 1.;
	iupdat = 0;
	updatd = FALSE_;
	timer_(&cpu2);
	sbtime = sbtime + cpu2 - cpu1;
	goto L222;
    }
/*        compute r=-Z'B(xcp-xk)-Z'g (using wa(2m+1)=W'(xcp-x) */
/*                                                   from 'cauchy'). */
    cmprlb_(n, m, &x[1], &g[1], &ws[ws_offset], &wy[wy_offset], &sy[sy_offset]
	    , &wt[wt_offset], &z__[1], &r__[1], &wa[1], &index[1], &theta, &
	    col, &head, &nfree, &cnstnd, &info);
    if (info != 0) {
	goto L444;
    }
/*       call the direct method. */
    subsm_(n, m, &nfree, &index[1], &l[1], &u[1], &nbd[1], &z__[1], &r__[1], &
	    ws[ws_offset], &wy[wy_offset], &theta, &col, &head, &iword, &wa[1]
	    , &wn[wn_offset], iprint, &info);
L444:
    if (info != 0) {
/*          singular triangular system detected; */
/*          refresh the lbfgs memory and restart the iteration. */
	if (*iprint >= 1) {
/*	    s_wsfe(&io___69);
	    e_wsfe();*/
	}
	info = 0;
	col = 0;
	head = 1;
	theta = 1.;
	iupdat = 0;
	updatd = FALSE_;
	timer_(&cpu2);
	sbtime = sbtime + cpu2 - cpu1;
	goto L222;
    }
    timer_(&cpu2);
    sbtime = sbtime + cpu2 - cpu1;
L555:
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */

/*     Line search and optimality tests. */

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*     Generate the search direction d:=z-x. */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	d__[i__] = z__[i__] - x[i__];
/* L40: */
    }
    timer_(&cpu1);
L666:
    lnsrlb_(n, &l[1], &u[1], &nbd[1], &x[1], f, &fold, &gd, &gdold, &g[1], &
	    d__[1], &r__[1], &t[1], &z__[1], &stp, &dnorm, &dtd, &xstep, &
	    stpmx, &iter, &ifun, &iback, &nfgv, &info, task, &boxed, &cnstnd, 
	    csave, &isave[22], &dsave[17], (long int)60, (long int)60);
    if (info != 0 || iback >= 20) {
/*          restore the previous iterate. */
	dcopy_(n, &t[1], &c__1, &x[1], &c__1);
	dcopy_(n, &r__[1], &c__1, &g[1], &c__1);
	*f = fold;
	if (col == 0) {
/*             abnormal termination. */
	    if (info == 0) {
		info = -9;
/*                restore the actual number of f and g evaluations etc. */
		--nfgv;
		--ifun;
		--iback;
	    }
	    s_copy(task, "ABNORMAL_TERMINATION_IN_LNSRCH", (long int)60, (
		    long int)30);
	    ++iter;
	    goto L999;
	} else {
/*             refresh the lbfgs memory and restart the iteration. */
	    if (*iprint >= 1) {
		/*s_wsfe(&io___71);
		e_wsfe();*/
	    }
	    if (info == 0) {
		--nfgv;
	    }
	    info = 0;
	    col = 0;
	    head = 1;
	    theta = 1.;
	    iupdat = 0;
	    updatd = FALSE_;
	    s_copy(task, "RESTART_FROM_LNSRCH", (long int)60, (long int)19);
	    timer_(&cpu2);
	    lnscht = lnscht + cpu2 - cpu1;
	    goto L222;
	}
    } else if (s_cmp(task, "FG_LN", (long int)5, (long int)5) == 0) {
/*          return to the driver for calculating f and g; reenter at 666. */
	goto L1000;
    } else {
/*          calculate and print out the quantities related to the new X. */
	timer_(&cpu2);
	lnscht = lnscht + cpu2 - cpu1;
	++iter;
/*        Compute the infinity norm of the projected (-)gradient. */
	projgr_(n, &l[1], &u[1], &nbd[1], &x[1], &g[1], &sbgnrm);
/*        Print iteration information. */
	prn2lb_(n, &x[1], f, &g[1], iprint, &itfile, &iter, &nfgv, &nact, &
		sbgnrm, &nint, word, &iword, &iback, &stp, &xstep, (long int)3);
	goto L1000;
    }
L777:
/*     Test for termination. */
    if (sbgnrm <= *pgtol) {
/*                                terminate the algorithm. */
	s_copy(task, "CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL", (
		long int)60, (long int)48);
	goto L999;
    }
/* Computing MAX */
    d__1 = abs_(fold), d__2 = abs_(*f), d__1 = max_(d__1,d__2);
    ddum = max_(d__1,1.);
    if (fold - *f <= tol * ddum) {
/*                                        terminate the algorithm. */
	s_copy(task, "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH", (
		long int)60, (long int)47);
	if (iback >= 10) {
	    info = -5;
	}
/*           i.e., to issue a warning if iback>10 in the line search. */
	goto L999;
    }
/*     Compute d=newx-oldx, r=newg-oldg, rr=y'y and dr=y's. */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	r__[i__] = g[i__] - r__[i__];
/* L42: */
    }
    rr = ddot_(n, &r__[1], &c__1, &r__[1], &c__1);
    if (stp == 1.) {
	dr = gd - gdold;
	ddum = -gdold;
    } else {
	dr = (gd - gdold) * stp;
	dscal_(n, &stp, &d__[1], &c__1);
	ddum = -gdold * stp;
    }
    if (dr <= epsmch * ddum) {
/*                            skip the L-BFGS update. */
	++nskip;
	updatd = FALSE_;
	if (*iprint >= 1) {
	   /* s_wsfe(&io___75);
	    do_fio(&c__1, (char *)&dr, (long int)sizeof(double));
	    do_fio(&c__1, (char *)&ddum, (long int)sizeof(double));
	    e_wsfe();*/
	}
	goto L888;
    }
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */

/*     Update the L-BFGS matrix. */

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
    updatd = TRUE_;
    ++iupdat;
/*     Update matrices WS and WY and form the middle matrix in B. */
    matupd_(n, m, &ws[ws_offset], &wy[wy_offset], &sy[sy_offset], &ss[
	    ss_offset], &d__[1], &r__[1], &itail, &iupdat, &col, &head, &
	    theta, &rr, &dr, &stp, &dtd);
/*     Form the upper half of the pds T = theta*SS + L*D^(-1)*L'; */
/*        Store T in the upper triangular of the array wt; */
/*        Cholesky factorize T to J*J' with */
/*           J' stored in the upper triangular of wt. */
    formt_(m, &wt[wt_offset], &sy[sy_offset], &ss[ss_offset], &col, &theta, &
	    info);
    if (info != 0) {
/*          nonpositive definiteness in Cholesky factorization; */
/*          refresh the lbfgs memory and restart the iteration. */
	if (*iprint >= 1) {
	    /*s_wsfe(&io___76);
	    e_wsfe();*/
	}
	info = 0;
	col = 0;
	head = 1;
	theta = 1.;
	iupdat = 0;
	updatd = FALSE_;
	goto L222;
    }
/*     Now the inverse of the middle matrix in B is */
/*       [  D^(1/2)      O ] [ -D^(1/2)  D^(-1/2)*L' ] */
/*       [ -L*D^(-1/2)   J ] [  0        J'          ] */
L888:
/* -------------------- the end of the loop ----------------------------- */
    goto L222;
L999:
    timer_(&time2);
    time = time2 - time1;
    prn3lb_(n, &x[1], f, task, iprint, &info, &itfile, &iter, &nfgv, &nintol, 
	    &nskip, &nact, &sbgnrm, &time, &nint, word, &iback, &stp, &xstep, 
	    &k, &cachyt, &sbtime, &lnscht, (long int)60, (long int)3);
L1000:
/*     Save local variables. */
    lsave[1] = prjctd;
    lsave[2] = cnstnd;
    lsave[3] = boxed;
    lsave[4] = updatd;
    isave[1] = nintol;
    isave[3] = itfile;
    isave[4] = iback;
    isave[5] = nskip;
    isave[6] = head;
    isave[7] = col;
    isave[8] = itail;
    isave[9] = iter;
    isave[10] = iupdat;
    isave[12] = nint;
    isave[13] = nfgv;
    isave[14] = info;
    isave[15] = ifun;
    isave[16] = iword;
    isave[17] = nfree;
    isave[18] = nact;
    isave[19] = ileave;
    isave[20] = nenter;
    dsave[1] = theta;
    dsave[2] = fold;
    dsave[3] = tol;
    dsave[4] = dnorm;
    dsave[5] = epsmch;
    dsave[6] = cpu1;
    dsave[7] = cachyt;
    dsave[8] = sbtime;
    dsave[9] = lnscht;
    dsave[10] = time1;
    dsave[11] = gd;
    dsave[12] = stpmx;
    dsave[13] = sbgnrm;
    dsave[14] = stp;
    dsave[15] = gdold;
    dsave[16] = dtd;
    return 0;
} /* mainlb_ */

/* ======================= The end of mainlb ============================= */
/* Subroutine */ int active_(long int *n, double *l, double *u, long int *nbd, double *x, long int *iwhere, 
			 long int *iprint, long int *prjctd, long int *cnstnd, long int *boxed)
/*long int *n;
double *l, *u;
long int *nbd;
double *x;
long int *iwhere, *iprint;
long int *prjctd, *cnstnd, *boxed;*/
{
    /* Format strings *//*
    static char fmt_1001[] = "(/,\002At X0 \002,i9,\002 variables are exactl\
y at the bounds\002)";*/

    /* System generated locals */
    long int i__1;

    /* Builtin functions */
    //long int s_wsle(cilist *), do_lio(long int *, long int *, char *, long int), e_wsle(), s_wsfe(cilist *), do_fio(long int*, char*, long int), e_wsfe();

    /* Local variables */
    static long int nbdd, i__;

    /* Fortran I/O blocks */
/*    static cilist io___81 = { 0, 6, 0, 0, 0 };
    static cilist io___82 = { 0, 6, 0, 0, 0 };
    static cilist io___83 = { 0, 6, 0, fmt_1001, 0 };*/


/*     ************ */

/*     Subroutine active */

/*     This subroutine initializes iwhere and projects the initial x to */
/*       the feasible set if necessary. */

/*     iwhere is an long int array of dimension n. */
/*       On entry iwhere is unspecified. */
/*       On exit iwhere(i)=-1  if x(i) has no bounds */
/*                         3   if l(i)=u(i) */
/*                         0   otherwise. */
/*       In cauchy, iwhere is given finer gradations. */


/*                           *  *  * */

/*     NEOS, November 1994. (Latest revision June 1996.) */
/*     Optimization Technology Center. */
/*     Argonne National Laboratory and Northwestern University. */
/*     Written by */
/*                        Ciyou Zhu */
/*     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal. */


/*     ************ */
/*     Initialize nbdd, prjctd, cnstnd and boxed. */
    /* Parameter adjustments */
    --iwhere;
    --x;
    --nbd;
    --u;
    --l;

    /* Function Body */
    nbdd = 0;
    *prjctd = FALSE_;
    *cnstnd = FALSE_;
    *boxed = TRUE_;
/*     Project the initial x to the easible set if necessary. */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (nbd[i__] > 0) {
	    if (nbd[i__] <= 2 && x[i__] <= l[i__]) {
		if (x[i__] < l[i__]) {
		    *prjctd = TRUE_;
		    x[i__] = l[i__];
		}
		++nbdd;
	    } else if (nbd[i__] >= 2 && x[i__] >= u[i__]) {
		if (x[i__] > u[i__]) {
		    *prjctd = TRUE_;
		    x[i__] = u[i__];
		}
		++nbdd;
	    }
	}
/* L10: */
    }
/*     Initialize iwhere and assign values to cnstnd and boxed. */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (nbd[i__] != 2) {
	    *boxed = FALSE_;
	}
	if (nbd[i__] == 0) {
/*                                this variable is always free */
	    iwhere[i__] = -1;
/*           otherwise set x(i)=mid(x(i), u(i), l(i)). */
	} else {
	    *cnstnd = TRUE_;
	    if (nbd[i__] == 2 && u[i__] - l[i__] <= 0.) {
/*                   this variable is always fixed */
		iwhere[i__] = 3;
	    } else {
		iwhere[i__] = 0;
	    }
	}
/* L20: */
    }
    if (*iprint >= 0) {
/*	if (*prjctd) {
	    s_wsle(&io___81);
	    do_lio(&c__9, &c__1, "The initial X is infeasible.  Restart with\
 its projection.", (long int)58);
	    e_wsle();
	}
	if (! (*cnstnd)) {
	    s_wsle(&io___82);
	    do_lio(&c__9, &c__1, "This problem is unconstrained.", (long int)30)
		    ;
	    e_wsle();
	}*/
    }
    if (*iprint > 0) {
/*	s_wsfe(&io___83);
	do_fio(&c__1, (char *)&nbdd, (long int)sizeof(long int));
	e_wsfe();*/
    }
    return 0;
} /* active_ */

/* ======================= The end of active ============================= */
/* Subroutine */ int bmv_(long int *m, double *sy, double *wt, long int *col, double *v, double *p, long int *info)
/*long int *m;
double *sy, *wt;
long int *col;
double *v, *p;
long int *info;*/
{
    /* System generated locals */
    long int sy_dim1, sy_offset, wt_dim1, wt_offset, i__1, i__2;

    /* Builtin functions */
    //double sqrt();

    /* Local variables */
    static long int i__, k;
    extern /* Subroutine */ int dtrsl_(double *t, long int *ldt, long int *n, double *b, long int *job, long int *info);
    static long int i2;
    static double sum;

/*     ************ */

/*     Subroutine bmv */

/*     This subroutine computes the product of the 2m x 2m middle matrix */
/* 	in the compact L-BFGS formula of B and a 2m vector v; */
/* 	it returns the product in p. */

/*     m is an long int variable. */
/*       On entry m is the maximum number of variable metric corrections */
/*         used to define the limited memory matrix. */
/*       On exit m is unchanged. */

/*     sy is a double precision array of dimension m x m. */
/*       On entry sy specifies the matrix S'Y. */
/*       On exit sy is unchanged. */

/*     wt is a double precision array of dimension m x m. */
/*       On entry wt specifies the upper triangular matrix J' which is */
/*         the Cholesky factor of (thetaS'S+LD^(-1)L'). */
/*       On exit wt is unchanged. */

/*     col is an long int variable. */
/*       On entry col specifies the number of s-vectors (or y-vectors) */
/*         stored in the compact L-BFGS formula. */
/*       On exit col is unchanged. */

/*     v is a double precision array of dimension 2col. */
/*       On entry v specifies vector v. */
/*       On exit v is unchanged. */

/*     p is a double precision array of dimension 2col. */
/*       On entry p is unspecified. */
/*       On exit p is the product Mv. */

/*     info is an long int variable. */
/*       On entry info is unspecified. */
/*       On exit info = 0       for normal return, */
/*                    = nonzero for abnormal return when the system */
/*                                to be solved by dtrsl is singular. */

/*     Subprograms called: */

/*       Linpack ... dtrsl. */


/*                           *  *  * */

/*     NEOS, November 1994. (Latest revision June 1996.) */
/*     Optimization Technology Center. */
/*     Argonne National Laboratory and Northwestern University. */
/*     Written by */
/*                        Ciyou Zhu */
/*     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal. */


/*     ************ */
    /* Parameter adjustments */
    wt_dim1 = *m;
    wt_offset = 1 + wt_dim1 * 1;
    wt -= wt_offset;
    sy_dim1 = *m;
    sy_offset = 1 + sy_dim1 * 1;
    sy -= sy_offset;
    --p;
    --v;

    /* Function Body */
    if (*col == 0) {
	return 0;
    }
/*     PART I: solve [  D^(1/2)      O ] [ p1 ] = [ v1 ] */
/*                   [ -L*D^(-1/2)   J ] [ p2 ]   [ v2 ]. */
/* 	solve Jp2=v2+LD^(-1)v1. */
    p[*col + 1] = v[*col + 1];
    i__1 = *col;
    for (i__ = 2; i__ <= i__1; ++i__) {
	i2 = *col + i__;
	sum = 0.;
	i__2 = i__ - 1;
	for (k = 1; k <= i__2; ++k) {
	    sum += sy[i__ + k * sy_dim1] * v[k] / sy[k + k * sy_dim1];
/* L10: */
	}
	p[i2] = v[i2] + sum;
/* L20: */
    }
/*     Solve the triangular system */
    dtrsl_(&wt[wt_offset], m, col, &p[*col + 1], &c__11, info);
    if (*info != 0) {
	return 0;
    }
/*     	solve D^(1/2)p1=v1. */
    i__1 = *col;
    for (i__ = 1; i__ <= i__1; ++i__) {
	p[i__] = v[i__] / sqrt(sy[i__ + i__ * sy_dim1]);
/* L30: */
    }
/*     PART II: solve [ -D^(1/2)   D^(-1/2)*L'  ] [ p1 ] = [ p1 ] */
/*                    [  0         J'           ] [ p2 ]   [ p2 ]. */
/*       solve J^Tp2=p2. */
    dtrsl_(&wt[wt_offset], m, col, &p[*col + 1], &c__1, info);
    if (*info != 0) {
	return 0;
    }
/*       compute p1=-D^(-1/2)(p1-D^(-1/2)L'p2) */
/*                 =-D^(-1/2)p1+D^(-1)L'p2. */
    i__1 = *col;
    for (i__ = 1; i__ <= i__1; ++i__) {
	p[i__] = -p[i__] / sqrt(sy[i__ + i__ * sy_dim1]);
/* L40: */
    }
    i__1 = *col;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sum = 0.;
	i__2 = *col;
	for (k = i__ + 1; k <= i__2; ++k) {
	    sum += sy[k + i__ * sy_dim1] * p[*col + k] / sy[i__ + i__ * 
		    sy_dim1];
/* L50: */
	}
	p[i__] += sum;
/* L60: */
    }
    return 0;
} /* bmv_ */

/* ======================== The end of bmv =============================== */
/* Subroutine */ int cauchy_(long int *n, double *x, double *l, double *u, long int *nbd, double *g, long int *iorder,
			 long int *iwhere, double *t, double *d__, double *xcp, long int *m, double *wy, double *ws,
			 double *sy, double *wt, double *theta, long int *col, long int *head, double *p, double *c__,
			 double *wbp, double *v, long int *nint, double *sg, double *yg, long int *iprint, double *sbgnrm, long int *info)
/*long int *n;
double *x, *l, *u;
long int *nbd;
double *g;
long int *iorder, *iwhere;
double *t, *d__, *xcp;
long int *m;
double *wy, *ws, *sy, *wt, *theta;
long int *col, *head;
double *p, *c__, *wbp, *v;
long int *nint;
double *sg, *yg;
long int *iprint;
double *sbgnrm;
long int *info;*/
{
    /* Format strings *//*
    static char fmt_3010[] = "(/,\002---------------- CAUCHY entered--------\
-----------\002)";
    static char fmt_1010[] = "(\002Cauchy X =  \002,/,(4x,1p,6(1x,d11.4)))";
    static char fmt_4011[] = "(/,\002Piece    \002,i3,\002 --f1, f2 at start\
 point \002,1p,2(1x,d11.4))";
    static char fmt_5010[] = "(\002Distance to the next break point =  \002,\
1p,d11.4)";
    static char fmt_6010[] = "(\002Distance to the stationary point =  \002,\
1p,d11.4)";
    static char fmt_4010[] = "(\002Piece    \002,i3,\002 --f1, f2 at start p\
oint \002,1p,2(1x,d11.4))";
    static char fmt_2010[] = "(/,\002---------------- exit CAUCHY-----------\
-----------\002,/)";*/

    /* System generated locals */
    long int wy_dim1, wy_offset, ws_dim1, ws_offset, sy_dim1, sy_offset, 
	    wt_dim1, wt_offset, i__1, i__2;
    double d__1;

    /* Builtin functions */
//    long int s_wsle(cilist *), do_lio(long int *, long int *, char *, long int), e_wsle(), s_wsfe(cilist *), e_wsfe(), do_fio(long int*, char*, long int);

    /* Local variables */
    static double dibp;
    extern double ddot_(long int *n, double *dx, long int *incx, double *dy, long int *incy);
    static long int iter;
    static double zibp, tsum, dibp2;
    static long int i__, j;
    static long int bnded;
    extern /* Subroutine */ int dscal_(long int *n, double *da, double *dx, long int *incx);
    static double neggi;
    static long int nfree;
    static double bkmin;
    static long int nleft;
    extern /* Subroutine */ int dcopy_(long int *n, double *dx, long int *incx, double *dy, long int *incy),
    			  daxpy_(long int *n, double *da, double *dx, long int *incx, double *dy, long int *incy);
    static double f1, f2, dt, tj, tl;
    static long int nbreak, ibkmin;
    static double tu;
    extern /* Subroutine */ int hpsolb_(long int *n, double *t, long int *iorder, long int *iheap);
    static long int pointr;
    static double tj0;
    static long int xlower, xupper;
    static long int ibp;
    static double dtm;
    extern /* Subroutine */ int bmv_(long int *m, double *sy, double *wt, long int *col, double *v, double *p, long int *info);
    static double wmc, wmp, wmw;
    static long int col2;

    /* Fortran I/O blocks */
/*    static cilist io___88 = { 0, 6, 0, 0, 0 };
    static cilist io___96 = { 0, 6, 0, fmt_3010, 0 };
    static cilist io___105 = { 0, 6, 0, fmt_1010, 0 };
    static cilist io___109 = { 0, 6, 0, 0, 0 };
    static cilist io___116 = { 0, 6, 0, fmt_4011, 0 };
    static cilist io___117 = { 0, 6, 0, fmt_5010, 0 };
    static cilist io___118 = { 0, 6, 0, fmt_6010, 0 };
    static cilist io___121 = { 0, 6, 0, 0, 0 };
    static cilist io___126 = { 0, 6, 0, 0, 0 };
    static cilist io___127 = { 0, 6, 0, 0, 0 };
    static cilist io___128 = { 0, 6, 0, fmt_4010, 0 };
    static cilist io___129 = { 0, 6, 0, fmt_6010, 0 };
    static cilist io___130 = { 0, 6, 0, fmt_1010, 0 };
    static cilist io___131 = { 0, 6, 0, fmt_2010, 0 };*/


/*     ************ */

/*     Subroutine cauchy */

/*     For given x, l, u, g (with sbgnrm > 0), and a limited memory */
/*       BFGS matrix B defined in terms of matrices WY, WS, WT, and */
/*       scalars head, col, and theta, this subroutine computes the */
/*       generalized Cauchy point (GCP), defined as the first local */
/*       minimizer of the quadratic */

/*                  Q(x + s) = g's + 1/2 s'Bs */

/*       along the projected gradient direction P(x-tg,l,u). */
/*       The routine returns the GCP in xcp. */

/*     n is an long int variable. */
/*       On entry n is the dimension of the problem. */
/*       On exit n is unchanged. */

/*     x is a double precision array of dimension n. */
/*       On entry x is the starting point for the GCP computation. */
/*       On exit x is unchanged. */

/*     l is a double precision array of dimension n. */
/*       On entry l is the lower bound of x. */
/*       On exit l is unchanged. */

/*     u is a double precision array of dimension n. */
/*       On entry u is the upper bound of x. */
/*       On exit u is unchanged. */

/*     nbd is an long int array of dimension n. */
/*       On entry nbd represents the type of bounds imposed on the */
/*         variables, and must be specified as follows: */
/*         nbd(i)=0 if x(i) is unbounded, */
/*                1 if x(i) has only a lower bound, */
/*                2 if x(i) has both lower and upper bounds, and */
/*                3 if x(i) has only an upper bound. */
/*       On exit nbd is unchanged. */

/*     g is a double precision array of dimension n. */
/*       On entry g is the gradient of f(x).  g must be a nonzero vector. */
/*       On exit g is unchanged. */

/*     iorder is an long int working array of dimension n. */
/*       iorder will be used to store the breakpoints in the piecewise */
/*       linear path and free variables encountered. On exit, */
/*         iorder(1),...,iorder(nleft) are indices of breakpoints */
/*                                which have not been encountered; */
/*         iorder(nleft+1),...,iorder(nbreak) are indices of */
/*                                     encountered breakpoints; and */
/*         iorder(nfree),...,iorder(n) are indices of variables which */
/*                 have no bound constraits along the search direction. */

/*     iwhere is an long int array of dimension n. */
/*       On entry iwhere indicates only the permanently fixed (iwhere=3) */
/*       or free (iwhere= -1) components of x. */
/*       On exit iwhere records the status of the current x variables. */
/*       iwhere(i)=-3  if x(i) is free and has bounds, but is not moved */
/*                 0   if x(i) is free and has bounds, and is moved */
/*                 1   if x(i) is fixed at l(i), and l(i) .ne. u(i) */
/*                 2   if x(i) is fixed at u(i), and u(i) .ne. l(i) */
/*                 3   if x(i) is always fixed, i.e.,  u(i)=x(i)=l(i) */
/*                 -1  if x(i) is always free, i.e., it has no bounds. */

/*     t is a double precision working array of dimension n. */
/*       t will be used to store the break points. */

/*     d is a double precision array of dimension n used to store */
/*       the Cauchy direction P(x-tg)-x. */

/*     xcp is a double precision array of dimension n used to return the */
/*       GCP on exit. */

/*     m is an long int variable. */
/*       On entry m is the maximum number of variable metric corrections */
/*         used to define the limited memory matrix. */
/*       On exit m is unchanged. */

/*     ws, wy, sy, and wt are double precision arrays. */
/*       On entry they store information that defines the */
/*                             limited memory BFGS matrix: */
/*         ws(n,m) stores S, a set of s-vectors; */
/*         wy(n,m) stores Y, a set of y-vectors; */
/*         sy(m,m) stores S'Y; */
/*         wt(m,m) stores the */
/*                 Cholesky factorization of (theta*S'S+LD^(-1)L'). */
/*       On exit these arrays are unchanged. */

/*     theta is a double precision variable. */
/*       On entry theta is the scaling factor specifying B_0 = theta I. */
/*       On exit theta is unchanged. */

/*     col is an long int variable. */
/*       On entry col is the actual number of variable metric */
/*         corrections stored so far. */
/*       On exit col is unchanged. */

/*     head is an long int variable. */
/*       On entry head is the location of the first s-vector (or y-vector) */
/*         in S (or Y). */
/*       On exit col is unchanged. */

/*     p is a double precision working array of dimension 2m. */
/*       p will be used to store the vector p = W^(T)d. */

/*     c is a double precision working array of dimension 2m. */
/*       c will be used to store the vector c = W^(T)(xcp-x). */

/*     wbp is a double precision working array of dimension 2m. */
/*       wbp will be used to store the row of W corresponding */
/*         to a breakpoint. */

/*     v is a double precision working array of dimension 2m. */

/*     nint is an long int variable. */
/*       On exit nint records the number of quadratic segments explored */
/*         in searching for the GCP. */

/*     sg and yg are double precision arrays of dimension m. */
/*       On entry sg  and yg store S'g and Y'g correspondingly. */
/*       On exit they are unchanged. */

/*     iprint is an long int variable that must be set by the user. */
/*       It controls the frequency and type of output generated: */
/*        iprint<0    no output is generated; */
/*        iprint=0    print only one line at the last iteration; */
/*        0<iprint<99 print also f and |proj g| every iprint iterations; */
/*        iprint=99   print details of every iteration except n-vectors; */
/*        iprint=100  print also the changes of active set and final x; */
/*        iprint>100  print details of every iteration including x and g; */
/*       When iprint > 0, the file iterate.dat will be created to */
/*                        summarize the iteration. */

/*     sbgnrm is a double precision variable. */
/*       On entry sbgnrm is the norm of the projected gradient at x. */
/*       On exit sbgnrm is unchanged. */

/*     info is an long int variable. */
/*       On entry info is 0. */
/*       On exit info = 0       for normal return, */
/*                    = nonzero for abnormal return when the the system */
/*                              used in routine bmv is singular. */

/*     Subprograms called: */

/*       L-BFGS-B Library ... hpsolb, bmv. */

/*       Linpack ... dscal dcopy, daxpy. */


/*     References: */

/*       [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited */
/*       memory algorithm for bound constrained optimization'', */
/*       SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208. */

/*       [2] C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, ``L-BFGS-B: FORTRAN */
/*       Subroutines for Large Scale Bound Constrained Optimization'' */
/*       Tech. Report, NAM-11, EECS Department, Northwestern University, */
/*       1994. */

/*       (Postscript files of these papers are available via anonymous */
/*        ftp to eecs.nwu.edu in the directory pub/lbfgs/lbfgs_bcm.) */

/*                           *  *  * */

/*     NEOS, November 1994. (Latest revision June 1996.) */
/*     Optimization Technology Center. */
/*     Argonne National Laboratory and Northwestern University. */
/*     Written by */
/*                        Ciyou Zhu */
/*     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal. */


/*     ************ */
/*     Check the status of the variables, reset iwhere(i) if necessary; */
/*       compute the Cauchy direction d and the breakpoints t; initialize */
/*       the derivative f1 and the vector p = W'd (for theta = 1). */
    /* Parameter adjustments */
    --xcp;
    --d__;
    --t;
    --iwhere;
    --iorder;
    --g;
    --nbd;
    --u;
    --l;
    --x;
    --yg;
    --sg;
    --v;
    --wbp;
    --c__;
    --p;
    wt_dim1 = *m;
    wt_offset = 1 + wt_dim1 * 1;
    wt -= wt_offset;
    sy_dim1 = *m;
    sy_offset = 1 + sy_dim1 * 1;
    sy -= sy_offset;
    ws_dim1 = *n;
    ws_offset = 1 + ws_dim1 * 1;
    ws -= ws_offset;
    wy_dim1 = *n;
    wy_offset = 1 + wy_dim1 * 1;
    wy -= wy_offset;

    /* Function Body */
    if (*sbgnrm <= 0.) {
	if (*iprint >= 0) {
/*	    s_wsle(&io___88);
	    do_lio(&c__9, &c__1, "Subgnorm = 0.  GCP = X.", (long int)23);
	    e_wsle();*/
	}
	dcopy_(n, &x[1], &c__1, &xcp[1], &c__1);
	return 0;
    }
    bnded = TRUE_;
    nfree = *n + 1;
    nbreak = 0;
    ibkmin = 0;
    bkmin = 0.;
    col2 = *col << 1;
    f1 = 0.;
    if (*iprint >= 99) {
/*	s_wsfe(&io___96);
	e_wsfe();*/
    }
/*     We set p to zero and build it up as we determine d. */
    i__1 = col2;
    for (i__ = 1; i__ <= i__1; ++i__) {
	p[i__] = 0.;
/* L20: */
    }
/*     In the following loop we determine for each variable its bound */
/*        status and its breakpoint, and update p accordingly. */
/*        Smallest breakpoint is identified. */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	neggi = -g[i__];
	if (iwhere[i__] != 3 && iwhere[i__] != -1) {
/*             if x(i) is not a constant and has bounds, */
/*             compute the difference between x(i) and its bounds. */
	    if (nbd[i__] <= 2) {
		tl = x[i__] - l[i__];
	    }
	    if (nbd[i__] >= 2) {
		tu = u[i__] - x[i__];
	    }
/*           If a variable is close enough to a bound */
/*             we treat it as at bound. */
	    xlower = nbd[i__] <= 2 && tl <= 0.;
	    xupper = nbd[i__] >= 2 && tu <= 0.;
/*              reset iwhere(i). */
	    iwhere[i__] = 0;
	    if (xlower) {
		if (neggi <= 0.) {
		    iwhere[i__] = 1;
		}
	    } else if (xupper) {
		if (neggi >= 0.) {
		    iwhere[i__] = 2;
		}
	    } else {
		if (abs_(neggi) <= 0.) {
		    iwhere[i__] = -3;
		}
	    }
	}
	pointr = *head;
	if (iwhere[i__] != 0 && iwhere[i__] != -1) {
	    d__[i__] = 0.;
	} else {
	    d__[i__] = neggi;
	    f1 -= neggi * neggi;
/*             calculate p := p - W'e_i* (g_i). */
	    i__2 = *col;
	    for (j = 1; j <= i__2; ++j) {
		p[j] += wy[i__ + pointr * wy_dim1] * neggi;
		p[*col + j] += ws[i__ + pointr * ws_dim1] * neggi;
		pointr = pointr % *m + 1;
/* L40: */
	    }
	    if (nbd[i__] <= 2 && nbd[i__] != 0 && neggi < 0.) {
/*                                 x(i) + d(i) is bounded; compute t(i). */
		++nbreak;
		iorder[nbreak] = i__;
		t[nbreak] = tl / (-neggi);
		if (nbreak == 1 || t[nbreak] < bkmin) {
		    bkmin = t[nbreak];
		    ibkmin = nbreak;
		}
	    } else if (nbd[i__] >= 2 && neggi > 0.) {
/*                                 x(i) + d(i) is bounded; compute t(i). */
		++nbreak;
		iorder[nbreak] = i__;
		t[nbreak] = tu / neggi;
		if (nbreak == 1 || t[nbreak] < bkmin) {
		    bkmin = t[nbreak];
		    ibkmin = nbreak;
		}
	    } else {
/*                x(i) + d(i) is not bounded. */
		--nfree;
		iorder[nfree] = i__;
		if (abs_(neggi) > 0.) {
		    bnded = FALSE_;
		}
	    }
	}
/* L50: */
    }
/*     The indices of the nonzero components of d are now stored */
/*       in iorder(1),...,iorder(nbreak) and iorder(nfree),...,iorder(n). */
/*       The smallest of the nbreak breakpoints is in t(ibkmin)=bkmin. */
    if (*theta != 1.) {
/*                   complete the initialization of p for theta not= one. */
	dscal_(col, theta, &p[*col + 1], &c__1);
    }
/*     Initialize GCP xcp = x. */
    dcopy_(n, &x[1], &c__1, &xcp[1], &c__1);
    if (nbreak == 0 && nfree == *n + 1) {
/*                  is a zero vector, return with the initial xcp as GCP. */
	if (*iprint > 100) {
/*	    s_wsfe(&io___105);
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		do_fio(&c__1, (char *)&xcp[i__], (long int)sizeof(double));
	    }
	    e_wsfe();*/
	}
	return 0;
    }
/*     Initialize c = W'(xcp - x) = 0. */
    i__1 = col2;
    for (j = 1; j <= i__1; ++j) {
	c__[j] = 0.;
/* L60: */
    }
/*     Initialize derivative f2. */
    f2 = -(*theta) * f1;
    if (*col > 0) {
	bmv_(m, &sy[sy_offset], &wt[wt_offset], col, &p[1], &v[1], info);
	if (*info != 0) {
	    return 0;
	}
	f2 -= ddot_(&col2, &v[1], &c__1, &p[1], &c__1);
    }
    dtm = -f1 / f2;
    tsum = 0.;
    *nint = 1;
    if (*iprint >= 99) {
/*	s_wsle(&io___109);
	do_lio(&c__9, &c__1, "There are ", (long int)10);
	do_lio(&c__3, &c__1, (char *)&nbreak, (long int)sizeof(long int));
	do_lio(&c__9, &c__1, "  breakpoints ", (long int)14);
	e_wsle();*/
    }
/*     If there are no breakpoints, locate the GCP and return. */
    if (nbreak == 0) {
	goto L888;
    }
    nleft = nbreak;
    iter = 1;
    tj = 0.;
/* ------------------- the beginning of the loop ------------------------- */
L777:
/*     Find the next smallest breakpoint; */
/*       compute dt = t(nleft) - t(nleft + 1). */
    tj0 = tj;
    if (iter == 1) {
/*         Since we already have the smallest breakpoint we need not do */
/*         heapsort yet. Often only one breakpoint is used and the */
/*         cost of heapsort is avoided. */
	tj = bkmin;
	ibp = iorder[ibkmin];
    } else {
	if (iter == 2) {
/*             Replace the already used smallest breakpoint with the */
/*             breakpoint numbered nbreak > nlast, before heapsort call. */
	    if (ibkmin != nbreak) {
		t[ibkmin] = t[nbreak];
		iorder[ibkmin] = iorder[nbreak];
	    }
/*        Update heap structure of breakpoints */
/*           (if iter=2, initialize heap). */
	}
	i__1 = iter - 2;
	hpsolb_(&nleft, &t[1], &iorder[1], &i__1);
	tj = t[nleft];
	ibp = iorder[nleft];
    }
    dt = tj - tj0;
    if (dt != 0. && *iprint >= 100) {
/*	s_wsfe(&io___116);
	do_fio(&c__1, (char *)&(*nint), (long int)sizeof(long int));
	do_fio(&c__1, (char *)&f1, (long int)sizeof(double));
	do_fio(&c__1, (char *)&f2, (long int)sizeof(double));
	e_wsfe();
	s_wsfe(&io___117);
	do_fio(&c__1, (char *)&dt, (long int)sizeof(double));
	e_wsfe();
	s_wsfe(&io___118);
	do_fio(&c__1, (char *)&dtm, (long int)sizeof(double));
	e_wsfe();*/
    }
/*     If a minimizer is within this interval, locate the GCP and return. */
    if (dtm < dt) {
	goto L888;
    }
/*     Otherwise fix one variable and */
/*       reset the corresponding component of d to zero. */
    tsum += dt;
    --nleft;
    ++iter;
    dibp = d__[ibp];
    d__[ibp] = 0.;
    if (dibp > 0.) {
	zibp = u[ibp] - x[ibp];
	xcp[ibp] = u[ibp];
	iwhere[ibp] = 2;
    } else {
	zibp = l[ibp] - x[ibp];
	xcp[ibp] = l[ibp];
	iwhere[ibp] = 1;
    }
    if (*iprint >= 100) {
/*	s_wsle(&io___121);
	do_lio(&c__9, &c__1, "Variable  ", (long int)10);
	do_lio(&c__3, &c__1, (char *)&ibp, (long int)sizeof(long int));
	do_lio(&c__9, &c__1, "  is fixed.", (long int)11);
	e_wsle();*/
    }
    if (nleft == 0 && nbreak == *n) {
/*                                             all n variables are fixed, */
/*                                                return with xcp as GCP. */
	dtm = dt;
	goto L999;
    }
/*     Update the derivative information. */
    ++(*nint);
/* Computing 2nd power */
    d__1 = dibp;
    dibp2 = d__1 * d__1;
/*     Update f1 and f2. */
/*        temporarily set f1 and f2 for col=0. */
    f1 = f1 + dt * f2 + dibp2 - *theta * dibp * zibp;
    f2 -= *theta * dibp2;
    if (*col > 0) {
/*                          update c = c + dt*p. */
	daxpy_(&col2, &dt, &p[1], &c__1, &c__[1], &c__1);
/*           choose wbp, */
/*           the row of W corresponding to the breakpoint encountered. */
	pointr = *head;
	i__1 = *col;
	for (j = 1; j <= i__1; ++j) {
	    wbp[j] = wy[ibp + pointr * wy_dim1];
	    wbp[*col + j] = *theta * ws[ibp + pointr * ws_dim1];
	    pointr = pointr % *m + 1;
/* L70: */
	}
/*           compute (wbp)Mc, (wbp)Mp, and (wbp)M(wbp)'. */
	bmv_(m, &sy[sy_offset], &wt[wt_offset], col, &wbp[1], &v[1], info);
	if (*info != 0) {
	    return 0;
	}
	wmc = ddot_(&col2, &c__[1], &c__1, &v[1], &c__1);
	wmp = ddot_(&col2, &p[1], &c__1, &v[1], &c__1);
	wmw = ddot_(&col2, &wbp[1], &c__1, &v[1], &c__1);
/*           update p = p - dibp*wbp. */
	d__1 = -dibp;
	daxpy_(&col2, &d__1, &wbp[1], &c__1, &p[1], &c__1);
/*           complete updating f1 and f2 while col > 0. */
	f1 += dibp * wmc;
	f2 = f2 + dibp * 2. * wmp - dibp2 * wmw;
    }
    if (nleft > 0) {
	dtm = -f1 / f2;
	goto L777;
/*                 to repeat the loop for unsearched intervals. */
    } else if (bnded) {
	f1 = 0.;
	f2 = 0.;
	dtm = 0.;
    } else {
	dtm = -f1 / f2;
    }
/* ------------------- the end of the loop ------------------------------- */
L888:
    if (*iprint >= 99) {
/*	s_wsle(&io___126);
	e_wsle();
	s_wsle(&io___127);
	do_lio(&c__9, &c__1, "GCP found in this segment", (long int)25);
	e_wsle();
	s_wsfe(&io___128);
	do_fio(&c__1, (char *)&(*nint), (long int)sizeof(long int));
	do_fio(&c__1, (char *)&f1, (long int)sizeof(double));
	do_fio(&c__1, (char *)&f2, (long int)sizeof(double));
	e_wsfe();
	s_wsfe(&io___129);
	do_fio(&c__1, (char *)&dtm, (long int)sizeof(double));
	e_wsfe();*/
    }
    if (dtm <= 0.) {
	dtm = 0.;
    }
    tsum += dtm;
/*     Move free variables (i.e., the ones w/o breakpoints) and */
/*       the variables whose breakpoints haven't been reached. */
    daxpy_(n, &tsum, &d__[1], &c__1, &xcp[1], &c__1);
L999:
/*     Update c = c + dtm*p = W'(x^c - x) */
/*       which will be used in computing r = Z'(B(x^c - x) + g). */
    if (*col > 0) {
	daxpy_(&col2, &dtm, &p[1], &c__1, &c__[1], &c__1);
    }
    if (*iprint > 100) {
/*	s_wsfe(&io___130);
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&xcp[i__], (long int)sizeof(double));
	}
	e_wsfe();*/
    }
    if (*iprint >= 99) {
/*	s_wsfe(&io___131);
	e_wsfe();*/
    }
    return 0;
} /* cauchy_ */

/* ====================== The end of cauchy ============================== */
/* Subroutine */ int cmprlb_(long int *n, long int *m, double *x, double *g, double *ws, double *wy, double *sy,
			 double *wt, double *z__, double *r__, double *wa, long int *index, double *theta,
			 long int *col, long int *head, long int *nfree, long int *cnstnd, long int *info)
/*long int *n, *m;
double *x, *g, *ws, *wy, *sy, *wt, *z__, *r__, *wa;
long int *index;
double *theta;
long int *col, *head, *nfree;
long int *cnstnd;
long int *info;*/
{
    /* System generated locals */
    long int ws_dim1, ws_offset, wy_dim1, wy_offset, sy_dim1, sy_offset, 
	    wt_dim1, wt_offset, i__1, i__2;

    /* Local variables */
    static long int i__, j, k;
    static double a1, a2;
    static long int pointr;
    extern /* Subroutine */ int bmv_(long int *m, double *sy, double *wt, long int *col, double *v, double *p, long int *info);

/*     ************ */

/*     Subroutine cmprlb */

/*       This subroutine computes r=-Z'B(xcp-xk)-Z'g by using */
/*         wa(2m+1)=W'(xcp-x) from subroutine cauchy. */

/*     Subprograms called: */

/*       L-BFGS-B Library ... bmv. */


/*                           *  *  * */

/*     NEOS, November 1994. (Latest revision June 1996.) */
/*     Optimization Technology Center. */
/*     Argonne National Laboratory and Northwestern University. */
/*     Written by */
/*                        Ciyou Zhu */
/*     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal. */


/*     ************ */
    /* Parameter adjustments */
    --index;
    --r__;
    --z__;
    --g;
    --x;
    --wa;
    wt_dim1 = *m;
    wt_offset = 1 + wt_dim1 * 1;
    wt -= wt_offset;
    sy_dim1 = *m;
    sy_offset = 1 + sy_dim1 * 1;
    sy -= sy_offset;
    wy_dim1 = *n;
    wy_offset = 1 + wy_dim1 * 1;
    wy -= wy_offset;
    ws_dim1 = *n;
    ws_offset = 1 + ws_dim1 * 1;
    ws -= ws_offset;

    /* Function Body */
    if (! (*cnstnd) && *col > 0) {
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    r__[i__] = -g[i__];
/* L26: */
	}
    } else {
	i__1 = *nfree;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    k = index[i__];
	    r__[i__] = -(*theta) * (z__[k] - x[k]) - g[k];
/* L30: */
	}
	bmv_(m, &sy[sy_offset], &wt[wt_offset], col, &wa[(*m << 1) + 1], &wa[
		1], info);
	if (*info != 0) {
	    *info = -8;
	    return 0;
	}
	pointr = *head;
	i__1 = *col;
	for (j = 1; j <= i__1; ++j) {
	    a1 = wa[j];
	    a2 = *theta * wa[*col + j];
	    i__2 = *nfree;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		k = index[i__];
		r__[i__] = r__[i__] + wy[k + pointr * wy_dim1] * a1 + ws[k + 
			pointr * ws_dim1] * a2;
/* L32: */
	    }
	    pointr = pointr % *m + 1;
/* L34: */
	}
    }
    return 0;
} /* cmprlb_ */

/* ======================= The end of cmprlb ============================= */
/* Subroutine */ int errclb_(long int *n, long int *m, double *factr, double *l, double *u, long int *nbd, char *task,
			 long int *info, long int *k, long int task_len)
/*long int *n, *m;
double *factr, *l, *u;
long int *nbd;
char *task;
long int *info, *k;
long int task_len;*/
{
    /* System generated locals */
    long int i__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, const char *const, long int, long int);

    /* Local variables */
    static long int i__;

/*     ************ */

/*     Subroutine errclb */

/*     This subroutine checks the validity of the input data. */


/*                           *  *  * */

/*     NEOS, November 1994. (Latest revision June 1996.) */
/*     Optimization Technology Center. */
/*     Argonne National Laboratory and Northwestern University. */
/*     Written by */
/*                        Ciyou Zhu */
/*     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal. */


/*     ************ */
/*     Check the input arguments for errors. */
    /* Parameter adjustments */
    --nbd;
    --u;
    --l;

    /* Function Body */
    if (*n <= 0) {
	s_copy(task, "ERROR: N .LE. 0", (long int)60, (long int)15);
    }
    if (*m <= 0) {
	s_copy(task, "ERROR: M .LE. 0", (long int)60, (long int)15);
    }
    if (*factr < 0.) {
	s_copy(task, "ERROR: FACTR .LT. 0", (long int)60, (long int)19);
    }
/*     Check the validity of the arrays nbd(i), u(i), and l(i). */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (nbd[i__] < 0 || nbd[i__] > 3) {
/*                                                   return */
	    s_copy(task, "ERROR: INVALID NBD", (long int)60, (long int)18);
	    *info = -6;
	    *k = i__;
	}
	if (nbd[i__] == 2) {
	    if (l[i__] > u[i__]) {
/*                                    return */
		s_copy(task, "ERROR: NO FEASIBLE SOLUTION", (long int)60, (
			long int)27);
		*info = -7;
		*k = i__;
	    }
	}
/* L10: */
    }
    return 0;
} /* errclb_ */

/* ======================= The end of errclb ============================= */
/* Subroutine */ int formk_(long int *n, long int *nsub, long int *ind, long int *nenter, long int *ileave, long int *indx2, long int *iupdat, 
			long int *updatd, double *wn, double *wn1, long int *m, double *ws, double *wy, double *sy, 
			double *theta, long int *col, long int *head, long int *info)
/*long int *n, *nsub, *ind, *nenter, *ileave, *indx2, *iupdat;
long int *updatd;
double *wn, *wn1;
long int *m;
double *ws, *wy, *sy, *theta;
long int *col, *head, *info;*/
{
    /* System generated locals */
    long int wn_dim1, wn_offset, wn1_dim1, wn1_offset, ws_dim1, ws_offset, 
	    wy_dim1, wy_offset, sy_dim1, sy_offset, i__1, i__2, i__3;

    /* Local variables */
    static long int dend, pend;
    extern double ddot_(long int *n, double *dx, long int *incx, double *dy, long int *incy);
    static long int upcl;
    static double temp1, temp2, temp3, temp4;
    static long int i__, k;
    extern /* Subroutine */ int dpofa_(double *a, long int *lda, long int *n, long int *info), 
    			dcopy_(long int *n, double *dx, long int *incx, double *dy, long int *incy),
			dtrsl_(double *t, long int *ldt, long int *n, double *b, long int *job, long int *info);
    static long int ipntr, jpntr, k1, m2, dbegin, is, js, iy, jy, pbegin, is1, 
	    js1, col2;

/*     ************ */

/*     Subroutine formk */

/*     This subroutine forms  the LEL^T factorization of the indefinite */

/*       matrix    K = [-D -Y'ZZ'Y/theta     L_a'-R_z'  ] */
/*                     [L_a -R_z           theta*S'AA'S ] */
/*                                                    where E = [-I  0] */
/*                                                              [ 0  I] */
/*     The matrix K can be shown to be equal to the matrix M^[-1]N */
/*       occurring in section 5.1 of [1], as well as to the matrix */
/*       Mbar^[-1] Nbar in section 5.3. */

/*     n is an long int variable. */
/*       On entry n is the dimension of the problem. */
/*       On exit n is unchanged. */

/*     nsub is an long int variable */
/*       On entry nsub is the number of subspace variables in free set. */
/*       On exit nsub is not changed. */

/*     ind is an long int array of dimension nsub. */
/*       On entry ind specifies the indices of subspace variables. */
/*       On exit ind is unchanged. */

/*     nenter is an long int variable. */
/*       On entry nenter is the number of variables entering the */
/*         free set. */
/*       On exit nenter is unchanged. */

/*     ileave is an long int variable. */
/*       On entry indx2(ileave),...,indx2(n) are the variables leaving */
/*         the free set. */
/*       On exit ileave is unchanged. */

/*     indx2 is an long int array of dimension n. */
/*       On entry indx2(1),...,indx2(nenter) are the variables entering */
/*         the free set, while indx2(ileave),...,indx2(n) are the */
/*         variables leaving the free set. */
/*       On exit indx2 is unchanged. */

/*     iupdat is an long int variable. */
/*       On entry iupdat is the total number of BFGS updates made so far. */
/*       On exit iupdat is unchanged. */

/*     updatd is a long int variable. */
/*       On entry 'updatd' is true if the L-BFGS matrix is updatd. */
/*       On exit 'updatd' is unchanged. */

/*     wn is a double precision array of dimension 2m x 2m. */
/*       On entry wn is unspecified. */
/*       On exit the upper triangle of wn stores the LEL^T factorization */
/*         of the 2*col x 2*col indefinite matrix */
/*                     [-D -Y'ZZ'Y/theta     L_a'-R_z'  ] */
/*                     [L_a -R_z           theta*S'AA'S ] */

/*     wn1 is a double precision array of dimension 2m x 2m. */
/*       On entry wn1 stores the lower triangular part of */
/*                     [Y' ZZ'Y   L_a'+R_z'] */
/*                     [L_a+R_z   S'AA'S   ] */
/*         in the previous iteration. */
/*       On exit wn1 stores the corresponding updated matrices. */
/*       The purpose of wn1 is just to store these inner products */
/*       so they can be easily updated and inserted into wn. */

/*     m is an long int variable. */
/*       On entry m is the maximum number of variable metric corrections */
/*         used to define the limited memory matrix. */
/*       On exit m is unchanged. */

/*     ws, wy, sy, and wtyy are double precision arrays; */
/*     theta is a double precision variable; */
/*     col is an long int variable; */
/*     head is an long int variable. */
/*       On entry they store the information defining the */
/*                                          limited memory BFGS matrix: */
/*         ws(n,m) stores S, a set of s-vectors; */
/*         wy(n,m) stores Y, a set of y-vectors; */
/*         sy(m,m) stores S'Y; */
/*         wtyy(m,m) stores the Cholesky factorization */
/*                                   of (theta*S'S+LD^(-1)L') */
/*         theta is the scaling factor specifying B_0 = theta I; */
/*         col is the number of variable metric corrections stored; */
/*         head is the location of the 1st s- (or y-) vector in S (or Y). */
/*       On exit they are unchanged. */

/*     info is an long int variable. */
/*       On entry info is unspecified. */
/*       On exit info =  0 for normal return; */
/*                    = -1 when the 1st Cholesky factorization failed; */
/*                    = -2 when the 2st Cholesky factorization failed. */

/*     Subprograms called: */

/*       Linpack ... dcopy, dpofa, dtrsl. */


/*     References: */
/*       [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited */
/*       memory algorithm for bound constrained optimization'', */
/*       SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208. */

/*       [2] C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, ``L-BFGS-B: a */
/*       limited memory FORTRAN code for solving bound constrained */
/*       optimization problems'', Tech. Report, NAM-11, EECS Department, */
/*       Northwestern University, 1994. */

/*       (Postscript files of these papers are available via anonymous */
/*        ftp to eecs.nwu.edu in the directory pub/lbfgs/lbfgs_bcm.) */

/*                           *  *  * */

/*     NEOS, November 1994. (Latest revision June 1996.) */
/*     Optimization Technology Center. */
/*     Argonne National Laboratory and Northwestern University. */
/*     Written by */
/*                        Ciyou Zhu */
/*     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal. */


/*     ************ */
/*     Form the lower triangular part of */
/*               WN1 = [Y' ZZ'Y   L_a'+R_z'] */
/*                     [L_a+R_z   S'AA'S   ] */
/*        where L_a is the strictly lower triangular part of S'AA'Y */
/*              R_z is the upper triangular part of S'ZZ'Y. */
    /* Parameter adjustments */
    --indx2;
    --ind;
    sy_dim1 = *m;
    sy_offset = 1 + sy_dim1 * 1;
    sy -= sy_offset;
    wy_dim1 = *n;
    wy_offset = 1 + wy_dim1 * 1;
    wy -= wy_offset;
    ws_dim1 = *n;
    ws_offset = 1 + ws_dim1 * 1;
    ws -= ws_offset;
    wn1_dim1 = 2 * *m;
    wn1_offset = 1 + wn1_dim1 * 1;
    wn1 -= wn1_offset;
    wn_dim1 = 2 * *m;
    wn_offset = 1 + wn_dim1 * 1;
    wn -= wn_offset;

    /* Function Body */
    if (*updatd) {
	if (*iupdat > *m) {
/*                                 shift old part of WN1. */
	    i__1 = *m - 1;
	    for (jy = 1; jy <= i__1; ++jy) {
		js = *m + jy;
		i__2 = *m - jy;
		dcopy_(&i__2, &wn1[jy + 1 + (jy + 1) * wn1_dim1], &c__1, &wn1[
			jy + jy * wn1_dim1], &c__1);
		i__2 = *m - jy;
		dcopy_(&i__2, &wn1[js + 1 + (js + 1) * wn1_dim1], &c__1, &wn1[
			js + js * wn1_dim1], &c__1);
		i__2 = *m - 1;
		dcopy_(&i__2, &wn1[*m + 2 + (jy + 1) * wn1_dim1], &c__1, &wn1[
			*m + 1 + jy * wn1_dim1], &c__1);
/* L10: */
	    }
	}
/*          put new rows in blocks (1,1), (2,1) and (2,2). */
	pbegin = 1;
	pend = *nsub;
	dbegin = *nsub + 1;
	dend = *n;
	iy = *col;
	is = *m + *col;
	ipntr = *head + *col - 1;
	if (ipntr > *m) {
	    ipntr -= *m;
	}
	jpntr = *head;
	i__1 = *col;
	for (jy = 1; jy <= i__1; ++jy) {
	    js = *m + jy;
	    temp1 = 0.;
	    temp2 = 0.;
	    temp3 = 0.;
/*             compute element jy of row 'col' of Y'ZZ'Y */
	    i__2 = pend;
	    for (k = pbegin; k <= i__2; ++k) {
		k1 = ind[k];
		temp1 += wy[k1 + ipntr * wy_dim1] * wy[k1 + jpntr * wy_dim1];
/* L15: */
	    }
/*             compute elements jy of row 'col' of L_a and S'AA'S */
	    i__2 = dend;
	    for (k = dbegin; k <= i__2; ++k) {
		k1 = ind[k];
		temp2 += ws[k1 + ipntr * ws_dim1] * ws[k1 + jpntr * ws_dim1];
		temp3 += ws[k1 + ipntr * ws_dim1] * wy[k1 + jpntr * wy_dim1];
/* L16: */
	    }
	    wn1[iy + jy * wn1_dim1] = temp1;
	    wn1[is + js * wn1_dim1] = temp2;
	    wn1[is + jy * wn1_dim1] = temp3;
	    jpntr = jpntr % *m + 1;
/* L20: */
	}
/*          put new column in block (2,1). */
	jy = *col;
	jpntr = *head + *col - 1;
	if (jpntr > *m) {
	    jpntr -= *m;
	}
	ipntr = *head;
	i__1 = *col;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    is = *m + i__;
	    temp3 = 0.;
/*             compute element i of column 'col' of R_z */
	    i__2 = pend;
	    for (k = pbegin; k <= i__2; ++k) {
		k1 = ind[k];
		temp3 += ws[k1 + ipntr * ws_dim1] * wy[k1 + jpntr * wy_dim1];
/* L25: */
	    }
	    ipntr = ipntr % *m + 1;
	    wn1[is + jy * wn1_dim1] = temp3;
/* L30: */
	}
	upcl = *col - 1;
    } else {
	upcl = *col;
    }
/*       modify the old parts in blocks (1,1) and (2,2) due to changes */
/*       in the set of free variables. */
    ipntr = *head;
    i__1 = upcl;
    for (iy = 1; iy <= i__1; ++iy) {
	is = *m + iy;
	jpntr = *head;
	i__2 = iy;
	for (jy = 1; jy <= i__2; ++jy) {
	    js = *m + jy;
	    temp1 = 0.;
	    temp2 = 0.;
	    temp3 = 0.;
	    temp4 = 0.;
	    i__3 = *nenter;
	    for (k = 1; k <= i__3; ++k) {
		k1 = indx2[k];
		temp1 += wy[k1 + ipntr * wy_dim1] * wy[k1 + jpntr * wy_dim1];
		temp2 += ws[k1 + ipntr * ws_dim1] * ws[k1 + jpntr * ws_dim1];
/* L35: */
	    }
	    i__3 = *n;
	    for (k = *ileave; k <= i__3; ++k) {
		k1 = indx2[k];
		temp3 += wy[k1 + ipntr * wy_dim1] * wy[k1 + jpntr * wy_dim1];
		temp4 += ws[k1 + ipntr * ws_dim1] * ws[k1 + jpntr * ws_dim1];
/* L36: */
	    }
	    wn1[iy + jy * wn1_dim1] = wn1[iy + jy * wn1_dim1] + temp1 - temp3;
	    wn1[is + js * wn1_dim1] = wn1[is + js * wn1_dim1] - temp2 + temp4;
	    jpntr = jpntr % *m + 1;
/* L40: */
	}
	ipntr = ipntr % *m + 1;
/* L45: */
    }
/*       modify the old parts in block (2,1). */
    ipntr = *head;
    i__1 = *m + upcl;
    for (is = *m + 1; is <= i__1; ++is) {
	jpntr = *head;
	i__2 = upcl;
	for (jy = 1; jy <= i__2; ++jy) {
	    temp1 = 0.;
	    temp3 = 0.;
	    i__3 = *nenter;
	    for (k = 1; k <= i__3; ++k) {
		k1 = indx2[k];
		temp1 += ws[k1 + ipntr * ws_dim1] * wy[k1 + jpntr * wy_dim1];
/* L50: */
	    }
	    i__3 = *n;
	    for (k = *ileave; k <= i__3; ++k) {
		k1 = indx2[k];
		temp3 += ws[k1 + ipntr * ws_dim1] * wy[k1 + jpntr * wy_dim1];
/* L51: */
	    }
	    if (is <= jy + *m) {
		wn1[is + jy * wn1_dim1] = wn1[is + jy * wn1_dim1] + temp1 - 
			temp3;
	    } else {
		wn1[is + jy * wn1_dim1] = wn1[is + jy * wn1_dim1] - temp1 + 
			temp3;
	    }
	    jpntr = jpntr % *m + 1;
/* L55: */
	}
	ipntr = ipntr % *m + 1;
/* L60: */
    }
/*     Form the upper triangle of WN = [D+Y' ZZ'Y/theta   -L_a'+R_z' ] */
/*                                     [-L_a +R_z        S'AA'S*theta] */
    m2 = *m << 1;
    i__1 = *col;
    for (iy = 1; iy <= i__1; ++iy) {
	is = *col + iy;
	is1 = *m + iy;
	i__2 = iy;
	for (jy = 1; jy <= i__2; ++jy) {
	    js = *col + jy;
	    js1 = *m + jy;
	    wn[jy + iy * wn_dim1] = wn1[iy + jy * wn1_dim1] / *theta;
	    wn[js + is * wn_dim1] = wn1[is1 + js1 * wn1_dim1] * *theta;
/* L65: */
	}
	i__2 = iy - 1;
	for (jy = 1; jy <= i__2; ++jy) {
	    wn[jy + is * wn_dim1] = -wn1[is1 + jy * wn1_dim1];
/* L66: */
	}
	i__2 = *col;
	for (jy = iy; jy <= i__2; ++jy) {
	    wn[jy + is * wn_dim1] = wn1[is1 + jy * wn1_dim1];
/* L67: */
	}
	wn[iy + iy * wn_dim1] += sy[iy + iy * sy_dim1];
/* L70: */
    }
/*     Form the upper triangle of WN= [  LL'            L^-1(-L_a'+R_z')] */
/*                                    [(-L_a +R_z)L'^-1   S'AA'S*theta  ] */
/*        first Cholesky factor (1,1) block of wn to get LL' */
/*                          with L' stored in the upper triangle of wn. */
    dpofa_(&wn[wn_offset], &m2, col, info);
    if (*info != 0) {
	*info = -1;
	return 0;
    }
/*        then form L^-1(-L_a'+R_z') in the (1,2) block. */
    col2 = *col << 1;
    i__1 = col2;
    for (js = *col + 1; js <= i__1; ++js) {
	dtrsl_(&wn[wn_offset], &m2, col, &wn[js * wn_dim1 + 1], &c__11, info);
/* L71: */
    }
/*     Form S'AA'S*theta + (L^-1(-L_a'+R_z'))'L^-1(-L_a'+R_z') in the */
/*        upper triangle of (2,2) block of wn. */
    i__1 = col2;
    for (is = *col + 1; is <= i__1; ++is) {
	i__2 = col2;
	for (js = is; js <= i__2; ++js) {
	    wn[is + js * wn_dim1] += ddot_(col, &wn[is * wn_dim1 + 1], &c__1, 
		    &wn[js * wn_dim1 + 1], &c__1);
/* L74: */
	}
/* L72: */
    }
/*     Cholesky factorization of (2,2) block of wn. */
    dpofa_(&wn[*col + 1 + (*col + 1) * wn_dim1], &m2, col, info);
    if (*info != 0) {
	*info = -2;
	return 0;
    }
    return 0;
} /* formk_ */

/* ======================= The end of formk ============================== */
/* Subroutine */ int formt_(long int *m, double *wt, double *sy, double *ss, long int *col, double *theta, long int *info)
/*long int *m;
double *wt, *sy, *ss;
long int *col;
double *theta;
long int *info;*/
{
    /* System generated locals */
    long int wt_dim1, wt_offset, sy_dim1, sy_offset, ss_dim1, ss_offset, i__1, 
	    i__2, i__3;

    /* Local variables */
    static double ddum;
    static long int i__, j, k;
    extern /* Subroutine */ int dpofa_(double *a, long int *lda, long int *n, long int *info);
    static long int k1;

/*     ************ */

/*     Subroutine formt */

/*       This subroutine forms the upper half of the pos. def. and symm. */
/*         T = theta*SS + L*D^(-1)*L', stores T in the upper triangle */
/*         of the array wt, and performs the Cholesky factorization of T */
/*         to produce J*J', with J' stored in the upper triangle of wt. */

/*     Subprograms called: */

/*       Linpack ... dpofa. */


/*                           *  *  * */

/*     NEOS, November 1994. (Latest revision June 1996.) */
/*     Optimization Technology Center. */
/*     Argonne National Laboratory and Northwestern University. */
/*     Written by */
/*                        Ciyou Zhu */
/*     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal. */


/*     ************ */
/*     Form the upper half of  T = theta*SS + L*D^(-1)*L', */
/*        store T in the upper triangle of the array wt. */
    /* Parameter adjustments */
    ss_dim1 = *m;
    ss_offset = 1 + ss_dim1 * 1;
    ss -= ss_offset;
    sy_dim1 = *m;
    sy_offset = 1 + sy_dim1 * 1;
    sy -= sy_offset;
    wt_dim1 = *m;
    wt_offset = 1 + wt_dim1 * 1;
    wt -= wt_offset;

    /* Function Body */
    i__1 = *col;
    for (j = 1; j <= i__1; ++j) {
	wt[j * wt_dim1 + 1] = *theta * ss[j * ss_dim1 + 1];
/* L52: */
    }
    i__1 = *col;
    for (i__ = 2; i__ <= i__1; ++i__) {
	i__2 = *col;
	for (j = i__; j <= i__2; ++j) {
	    k1 = min_(i__,j) - 1;
	    ddum = 0.;
	    i__3 = k1;
	    for (k = 1; k <= i__3; ++k) {
		ddum += sy[i__ + k * sy_dim1] * sy[j + k * sy_dim1] / sy[k + 
			k * sy_dim1];
/* L53: */
	    }
	    wt[i__ + j * wt_dim1] = ddum + *theta * ss[i__ + j * ss_dim1];
/* L54: */
	}
/* L55: */
    }
/*     Cholesky factorize T to J*J' with */
/*        J' stored in the upper triangle of wt. */
    dpofa_(&wt[wt_offset], m, col, info);
    if (*info != 0) {
	*info = -3;
    }
    return 0;
} /* formt_ */

/* ======================= The end of formt ============================== */
/* Subroutine */ int freev_(long int *n, long int *nfree, long int *index, long int *nenter, long int *ileave, long int *indx2, long int *iwhere, 
			long int *wrk, long int *updatd, long int *cnstnd, long int *iprint, long int *iter)
/*long int *n, *nfree, *index, *nenter, *ileave, *indx2, *iwhere;
long int *wrk, *updatd, *cnstnd;
long int *iprint, *iter;*/
{
    /* System generated locals */
    long int i__1;

    /* Builtin functions */
//    long int s_wsle(cilist *), do_lio(long int *, long int *, char *, long int), e_wsle();

    /* Local variables */
    static long int iact, i__, k;

    /* Fortran I/O blocks */
/*    static cilist io___168 = { 0, 6, 0, 0, 0 };
    static cilist io___169 = { 0, 6, 0, 0, 0 };
    static cilist io___170 = { 0, 6, 0, 0, 0 };
    static cilist io___172 = { 0, 6, 0, 0, 0 };*/


/*     ************ */

/*     Subroutine freev */

/*     This subroutine counts the entering and leaving variables when */
/*       iter > 0, and finds the index set of free and active variables */
/*       at the GCP. */

/*     cnstnd is a long int variable indicating whether bounds are present */

/*     index is an long int array of dimension n */
/*       for i=1,...,nfree, index(i) are the indices of free variables */
/*       for i=nfree+1,...,n, index(i) are the indices of bound variables */
/*       On entry after the first iteration, index gives */
/*         the free variables at the previous iteration. */
/*       On exit it gives the free variables based on the determination */
/*         in cauchy using the array iwhere. */

/*     indx2 is an long int array of dimension n */
/*       On entry indx2 is unspecified. */
/*       On exit with iter>0, indx2 indicates which variables */
/*          have changed status since the previous iteration. */
/*       For i= 1,...,nenter, indx2(i) have changed from bound to free. */
/*       For i= ileave+1,...,n, indx2(i) have changed from free to bound. */


/*                           *  *  * */

/*     NEOS, November 1994. (Latest revision June 1996.) */
/*     Optimization Technology Center. */
/*     Argonne National Laboratory and Northwestern University. */
/*     Written by */
/*                        Ciyou Zhu */
/*     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal. */


/*     ************ */
    /* Parameter adjustments */
    --iwhere;
    --indx2;
    --index;

    /* Function Body */
    *nenter = 0;
    *ileave = *n + 1;
    if (*iter > 0 && *cnstnd) {
/*                           count the entering and leaving variables. */
	i__1 = *nfree;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    k = index[i__];
	    if (iwhere[k] > 0) {
		--(*ileave);
		indx2[*ileave] = k;
		if (*iprint >= 100) {
/*		    s_wsle(&io___168);
		    do_lio(&c__9, &c__1, "Variable ", (long int)9);
		    do_lio(&c__3, &c__1, (char *)&k, (long int)sizeof(long int));
		    do_lio(&c__9, &c__1, " leaves the set of free variables", 
			    (long int)33);
		    e_wsle();*/
		}
	    }
/* L20: */
	}
	i__1 = *n;
	for (i__ = *nfree + 1; i__ <= i__1; ++i__) {
	    k = index[i__];
	    if (iwhere[k] <= 0) {
		++(*nenter);
		indx2[*nenter] = k;
		if (*iprint >= 100) {
/*		    s_wsle(&io___169);
		    do_lio(&c__9, &c__1, "Variable ", (long int)9);
		    do_lio(&c__3, &c__1, (char *)&k, (long int)sizeof(long int));
		    do_lio(&c__9, &c__1, " enters the set of free variables", 
			    (long int)33);
		    e_wsle();*/
		}
	    }
/* L22: */
	}
	if (*iprint >= 99) {
/*	    s_wsle(&io___170);
	    i__1 = *n + 1 - *ileave;
	    do_lio(&c__3, &c__1, (char *)&i__1, (long int)sizeof(long int));
	    do_lio(&c__9, &c__1, " variables leave; ", (long int)18);
	    do_lio(&c__3, &c__1, (char *)&(*nenter), (long int)sizeof(long int));
	    do_lio(&c__9, &c__1, " variables enter", (long int)16);
	    e_wsle();*/
	}
    }
    *wrk = *ileave < *n + 1 || *nenter > 0 || *updatd;
/*     Find the index set of free and active variables at the GCP. */
    *nfree = 0;
    iact = *n + 1;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (iwhere[i__] <= 0) {
	    ++(*nfree);
	    index[*nfree] = i__;
	} else {
	    --iact;
	    index[iact] = i__;
	}
/* L24: */
    }
    if (*iprint >= 99) {
/*	s_wsle(&io___172);
	do_lio(&c__3, &c__1, (char *)&(*nfree), (long int)sizeof(long int));
	do_lio(&c__9, &c__1, " variables are free at GCP ", (long int)27);
	i__1 = *iter + 1;
	do_lio(&c__3, &c__1, (char *)&i__1, (long int)sizeof(long int));
	e_wsle();*/
    }
    return 0;
} /* freev_ */

/* ======================= The end of freev ============================== */
/* Subroutine */ int hpsolb_(long int *n, double *t, long int *iorder, long int *iheap)
/*long int *n;
double *t;
long int *iorder, *iheap;*/
{
    /* System generated locals */
    long int i__1;

    /* Local variables */
    static double ddum;
    static long int i__, j, k, indxin, indxou;
    static double out;

/*     ************ */

/*     Subroutine hpsolb */

/*     This subroutine sorts out the least element of t, and puts the */
/*       remaining elements of t in a heap. */

/*     n is an long int variable. */
/*       On entry n is the dimension of the arrays t and iorder. */
/*       On exit n is unchanged. */

/*     t is a double precision array of dimension n. */
/*       On entry t stores the elements to be sorted, */
/*       On exit t(n) stores the least elements of t, and t(1) to t(n-1) */
/*         stores the remaining elements in the form of a heap. */

/*     iorder is an long int array of dimension n. */
/*       On entry iorder(i) is the index of t(i). */
/*       On exit iorder(i) is still the index of t(i), but iorder may be */
/*         permuted in accordance with t. */

/*     iheap is an long int variable specifying the task. */
/*       On entry iheap should be set as follows: */
/*         iheap .eq. 0 if t(1) to t(n) is not in the form of a heap, */
/*         iheap .ne. 0 if otherwise. */
/*       On exit iheap is unchanged. */


/*     References: */
/*       Algorithm 232 of CACM (J. W. J. Williams): HEAPSORT. */

/*                           *  *  * */

/*     NEOS, November 1994. (Latest revision June 1996.) */
/*     Optimization Technology Center. */
/*     Argonne National Laboratory and Northwestern University. */
/*     Written by */
/*                        Ciyou Zhu */
/*     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal. */

/*     ************ */
    /* Parameter adjustments */
    --iorder;
    --t;

    /* Function Body */
    if (*iheap == 0) {
/*        Rearrange the elements t(1) to t(n) to form a heap. */
	i__1 = *n;
	for (k = 2; k <= i__1; ++k) {
	    ddum = t[k];
	    indxin = iorder[k];
/*           Add ddum to the heap. */
	    i__ = k;
L10:
	    if (i__ > 1) {
		j = i__ / 2;
		if (ddum < t[j]) {
		    t[i__] = t[j];
		    iorder[i__] = iorder[j];
		    i__ = j;
		    goto L10;
		}
	    }
	    t[i__] = ddum;
	    iorder[i__] = indxin;
/* L20: */
	}
    }
/*     Assign to 'out' the value of t(1), the least member of the heap, */
/*        and rearrange the remaining members to form a heap as */
/*        elements 1 to n-1 of t. */
    if (*n > 1) {
	i__ = 1;
	out = t[1];
	indxou = iorder[1];
	ddum = t[*n];
	indxin = iorder[*n];
/*        Restore the heap */
L30:
	j = i__ + i__;
	if (j <= *n - 1) {
	    if (t[j + 1] < t[j]) {
		++j;
	    }
	    if (t[j] < ddum) {
		t[i__] = t[j];
		iorder[i__] = iorder[j];
		i__ = j;
		goto L30;
	    }
	}
	t[i__] = ddum;
	iorder[i__] = indxin;
/*     Put the least member in t(n). */
	t[*n] = out;
	iorder[*n] = indxou;
    }
    return 0;
} /* hpsolb_ */

/* ====================== The end of hpsolb ============================== */
/* Subroutine */ int lnsrlb_(long int *n, double *l, double *u, long int *nbd, double *x, double *f, double *fold,
			 double *gd, double *gdold, double *g, double *d__, double *r__, double *t, 
			 double *z__, double *stp, double *dnorm, double *dtd, double *xstep, double *stpmx,
			 long int *iter, long int *ifun, long int *iback, long int *nfgv, long int *info, char *task, long int *boxed, 
			 long int *cnstnd, char *csave, long int *isave, double *dsave, long int task_len, long int csave_len)
/*long int *n;
double *l, *u;
long int *nbd;
double *x, *f, *fold, *gd, *gdold, *g, *d__, *r__, *t, *z__, *stp, *dnorm,
	 *dtd, *xstep, *stpmx;
long int *iter, *ifun, *iback, *nfgv, *info;
char *task;
long int *boxed, *cnstnd;
char *csave;
long int *isave;
double *dsave;
long int task_len;
long int csave_len;*/
{
    /* System generated locals */
    long int i__1;
    double d__1;

    /* Builtin functions */
    long int s_cmp(char *, const char *const, long int, long int);
    //double sqrt();
    /* Subroutine */ int s_copy(char *, const char *const, long int, long int);

    /* Local variables */
    extern double ddot_(long int *n, double *dx, long int *incx, double *dy, long int *incy);
    static long int i__;
    extern /* Subroutine */ int dcopy_(long int *n, double *dx, long int *incx, double *dy, long int *incy);
    static double a1, a2;
    extern /* Subroutine */ int dcsrch_(double *f, double *g, double *stp, double *ftol, double *gtol, double *xtol,
			 double *stpmin, double *stpmax, char *task, long int *isave, double *dsave, long int task_len);

/*     ********** */

/*     Subroutine lnsrlb */

/*     This subroutine calls subroutine dcsrch from the Minpack2 library */
/*       to perform the line search.  Subroutine dscrch is safeguarded so */
/*       that all trial points lie within the feasible region. */

/*     Subprograms called: */

/*       Minpack2 Library ... dcsrch. */

/*       Linpack ... dtrsl, ddot. */


/*                           *  *  * */

/*     NEOS, November 1994. (Latest revision June 1996.) */
/*     Optimization Technology Center. */
/*     Argonne National Laboratory and Northwestern University. */
/*     Written by */
/*                        Ciyou Zhu */
/*     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal. */


/*     ********** */
    /* Parameter adjustments */
    --z__;
    --t;
    --r__;
    --d__;
    --g;
    --x;
    --nbd;
    --u;
    --l;
    --isave;
    --dsave;

    /* Function Body */
    if (s_cmp(task, "FG_LN", (long int)5, (long int)5) == 0) {
	goto L556;
    }
    *dtd = ddot_(n, &d__[1], &c__1, &d__[1], &c__1);
    *dnorm = sqrt(*dtd);
/*     Determine the maximum step length. */
    *stpmx = 1e10;
    if (*cnstnd) {
	if (*iter == 0) {
	    *stpmx = 1.;
	} else {
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		a1 = d__[i__];
		if (nbd[i__] != 0) {
		    if (a1 < 0. && nbd[i__] <= 2) {
			a2 = l[i__] - x[i__];
			if (a2 >= 0.) {
			    *stpmx = 0.;
			} else if (a1 * *stpmx < a2) {
			    *stpmx = a2 / a1;
			}
		    } else if (a1 > 0. && nbd[i__] >= 2) {
			a2 = u[i__] - x[i__];
			if (a2 <= 0.) {
			    *stpmx = 0.;
			} else if (a1 * *stpmx > a2) {
			    *stpmx = a2 / a1;
			}
		    }
		}
/* L43: */
	    }
	}
    }
    if (*iter == 0 && ! (*boxed)) {
/* Computing MIN */
	d__1 = 1. / *dnorm;
	*stp = min_(d__1,*stpmx);
    } else {
	*stp = 1.;
    }
    dcopy_(n, &x[1], &c__1, &t[1], &c__1);
    dcopy_(n, &g[1], &c__1, &r__[1], &c__1);
    *fold = *f;
    *ifun = 0;
    *iback = 0;
    s_copy(csave, "START", (long int)60, (long int)5);
L556:
    *gd = ddot_(n, &g[1], &c__1, &d__[1], &c__1);
    if (*ifun == 0) {
	*gdold = *gd;
	if (*gd >= 0.) {
/*                               the directional derivative >=0. */
/*                               Line search is impossible. */
	    *info = -4;
	    return 0;
	}
    }
    dcsrch_(f, gd, stp, &c_b275, &c_b276, &c_b277, &c_b9, stpmx, csave, &
	    isave[1], &dsave[1], (long int)60);
    *xstep = *stp * *dnorm;
    if (s_cmp(csave, "CONV", (long int)4, (long int)4) != 0 && s_cmp(csave, "WARN"
	    , (long int)4, (long int)4) != 0) {
	s_copy(task, "FG_LNSRCH", (long int)60, (long int)9);
	++(*ifun);
	++(*nfgv);
	*iback = *ifun - 1;
	if (*stp == 1.) {
	    dcopy_(n, &z__[1], &c__1, &x[1], &c__1);
	} else {
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		x[i__] = *stp * d__[i__] + t[i__];
/* L41: */
	    }
	}
    } else {
	s_copy(task, "NEW_X", (long int)60, (long int)5);
    }
    return 0;
} /* lnsrlb_ */

/* ======================= The end of lnsrlb ============================= */
/* Subroutine */ int matupd_(long int *n, long int *m, double *ws, double *wy, double *sy, double *ss,
			 double *d__, double *r__, long int *itail, long int *iupdat, long int *col, long int *head,
			 double *theta, double *rr, double *dr, double *stp, double *dtd)
/*long int *n, *m;
double *ws, *wy, *sy, *ss, *d__, *r__;
long int *itail, *iupdat, *col, *head;
double *theta, *rr, *dr, *stp, *dtd;*/
{
    /* System generated locals */
    long int ws_dim1, ws_offset, wy_dim1, wy_offset, sy_dim1, sy_offset, 
	    ss_dim1, ss_offset, i__1, i__2;

    /* Local variables */
    extern double ddot_(long int *n, double *dx, long int *incx, double *dy, long int *incy);
    static long int j;
    extern /* Subroutine */ int dcopy_(long int *n, double *dx, long int *incx, double *dy, long int *incy);
    static long int pointr;

/*     ************ */

/*     Subroutine matupd */

/*       This subroutine updates matrices WS and WY, and forms the */
/*         middle matrix in B. */

/*     Subprograms called: */

/*       Linpack ... dcopy, ddot. */


/*                           *  *  * */

/*     NEOS, November 1994. (Latest revision June 1996.) */
/*     Optimization Technology Center. */
/*     Argonne National Laboratory and Northwestern University. */
/*     Written by */
/*                        Ciyou Zhu */
/*     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal. */


/*     ************ */
/*     Set pointers for matrices WS and WY. */
    /* Parameter adjustments */
    --r__;
    --d__;
    ss_dim1 = *m;
    ss_offset = 1 + ss_dim1 * 1;
    ss -= ss_offset;
    sy_dim1 = *m;
    sy_offset = 1 + sy_dim1 * 1;
    sy -= sy_offset;
    wy_dim1 = *n;
    wy_offset = 1 + wy_dim1 * 1;
    wy -= wy_offset;
    ws_dim1 = *n;
    ws_offset = 1 + ws_dim1 * 1;
    ws -= ws_offset;

    /* Function Body */
    if (*iupdat <= *m) {
	*col = *iupdat;
	*itail = (*head + *iupdat - 2) % *m + 1;
    } else {
	*itail = *itail % *m + 1;
	*head = *head % *m + 1;
    }
/*     Update matrices WS and WY. */
    dcopy_(n, &d__[1], &c__1, &ws[*itail * ws_dim1 + 1], &c__1);
    dcopy_(n, &r__[1], &c__1, &wy[*itail * wy_dim1 + 1], &c__1);
/*     Set theta=yy/ys. */
    *theta = *rr / *dr;
/*     Form the middle matrix in B. */
/*        update the upper triangle of SS, */
/*                                         and the lower triangle of SY: */
    if (*iupdat > *m) {
/*                              move old information */
	i__1 = *col - 1;
	for (j = 1; j <= i__1; ++j) {
	    dcopy_(&j, &ss[(j + 1) * ss_dim1 + 2], &c__1, &ss[j * ss_dim1 + 1]
		    , &c__1);
	    i__2 = *col - j;
	    dcopy_(&i__2, &sy[j + 1 + (j + 1) * sy_dim1], &c__1, &sy[j + j * 
		    sy_dim1], &c__1);
/* L50: */
	}
    }
/*        add new information: the last row of SY */
/*                                             and the last column of SS: */
    pointr = *head;
    i__1 = *col - 1;
    for (j = 1; j <= i__1; ++j) {
	sy[*col + j * sy_dim1] = ddot_(n, &d__[1], &c__1, &wy[pointr * 
		wy_dim1 + 1], &c__1);
	ss[j + *col * ss_dim1] = ddot_(n, &ws[pointr * ws_dim1 + 1], &c__1, &
		d__[1], &c__1);
	pointr = pointr % *m + 1;
/* L51: */
    }
    if (*stp == 1.) {
	ss[*col + *col * ss_dim1] = *dtd;
    } else {
	ss[*col + *col * ss_dim1] = *stp * *stp * *dtd;
    }
    sy[*col + *col * sy_dim1] = *dr;
    return 0;
} /* matupd_ */

/* ======================= The end of matupd ============================= */
/* Subroutine */ int prn1lb_(long int *n, long int *m, double *l, double *u, double *x, long int *iprint, long int *itfile, double *epsmch)
/*long int *n, *m;
double *l, *u, *x;
long int *iprint, *itfile;
double *epsmch;*/
{
    /* Format strings *//*
    static char fmt_7001[] = "(\002RUNNING THE L-BFGS-B CODE\002,/,/,\002   \
        * * *\002,/,/,\002Machine precision =\002,1p,d10.3)";
    static char fmt_2001[] = "(\002RUNNING THE L-BFGS-B CODE\002,/,/,\002it \
   = iteration number\002,/,\002nf    = number of function evaluations\002,/,\
\002nint  = number of segments explored during the Cauchy search\002,/,\002n\
act  = number of active bounds at the generalized Cauchy point\002,/,\002sub\
   = manner in which the subspace minimization terminated:\002,/,\002       \
 con = converged, bnd = a bound was reached\002,/,\002itls  = number of iter\
ations performed in the line search\002,/,\002stepl = step length used\002,/,\
\002tstep = norm of the displacement (total step)\002,/,\002projg = norm of \
the projected gradient\002,/,\002f     = function value\002,/,/,\002        \
   * * *\002,/,/,\002Machine precision =\002,1p,d10.3)";
    static char fmt_9001[] = "(/,3x,\002it\002,3x,\002nf\002,2x,\002nint\002\
,2x,\002nact\002,2x,\002sub\002,2x,\002itls\002,2x,\002stepl\002,4x,\002tstep\
\002,5x,\002projg\002,8x,\002f\002)";
    static char fmt_1004[] = "(/,a4,1p,6(1x,d11.4),/,(4x,1p,6(1x,d11.4)))";*/

    /* System generated locals */
//    long int i__1;

    /* Builtin functions */
//    long int s_wsfe(cilist *), do_fio(long int*, char*, long int), e_wsfe(), s_wsle(cilist *), do_lio(long int *, long int *, char *, long int), e_wsle();

    /* Local variables */
//    static long int i__;

    /* Fortran I/O blocks */
/*    static cilist io___185 = { 0, 6, 0, fmt_7001, 0 };
    static cilist io___186 = { 0, 6, 0, 0, 0 };
    static cilist io___187 = { 0, 0, 0, fmt_2001, 0 };
    static cilist io___188 = { 0, 0, 0, 0, 0 };
    static cilist io___189 = { 0, 0, 0, fmt_9001, 0 };
    static cilist io___190 = { 0, 6, 0, fmt_1004, 0 };
    static cilist io___192 = { 0, 6, 0, fmt_1004, 0 };
    static cilist io___193 = { 0, 6, 0, fmt_1004, 0 }; */


/*     ************ */

/*     Subroutine prn1lb */

/*     This subroutine prints the input data, initial point, upper and */
/*       lower bounds of each variable, machine precision, as well as */
/*       the headings of the output. */


/*                           *  *  * */

/*     NEOS, November 1994. (Latest revision June 1996.) */
/*     Optimization Technology Center. */
/*     Argonne National Laboratory and Northwestern University. */
/*     Written by */
/*                        Ciyou Zhu */
/*     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal. */


/*     ************ */
    /* Parameter adjustments */
    --x;
    --u;
    --l;

    /* Function Body */
    if (*iprint >= 0) {
/*	s_wsfe(&io___185);
	do_fio(&c__1, (char *)&(*epsmch), (long int)sizeof(double));
	e_wsfe();
	s_wsle(&io___186);
	do_lio(&c__9, &c__1, "N = ", (long int)4);
	do_lio(&c__3, &c__1, (char *)&(*n), (long int)sizeof(long int));
	do_lio(&c__9, &c__1, "    M = ", (long int)8);
	do_lio(&c__3, &c__1, (char *)&(*m), (long int)sizeof(long int));
	e_wsle();
	if (*iprint >= 1) {
	    io___187.ciunit = *itfile;
	    s_wsfe(&io___187);
	    do_fio(&c__1, (char *)&(*epsmch), (long int)sizeof(double));
	    e_wsfe();
	    io___188.ciunit = *itfile;
	    s_wsle(&io___188);
	    do_lio(&c__9, &c__1, "N = ", (long int)4);
	    do_lio(&c__3, &c__1, (char *)&(*n), (long int)sizeof(long int));
	    do_lio(&c__9, &c__1, "    M = ", (long int)8);
	    do_lio(&c__3, &c__1, (char *)&(*m), (long int)sizeof(long int));
	    e_wsle();
	    io___189.ciunit = *itfile;
	    s_wsfe(&io___189);
	    e_wsfe();
	    if (*iprint > 100) {
		s_wsfe(&io___190);
		do_fio(&c__1, "L =", (long int)3);
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    do_fio(&c__1, (char *)&l[i__], (long int)sizeof(double))
			    ;
		}
		e_wsfe();
		s_wsfe(&io___192);
		do_fio(&c__1, "X0 =", (long int)4);
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    do_fio(&c__1, (char *)&x[i__], (long int)sizeof(double))
			    ;
		}
		e_wsfe();
		s_wsfe(&io___193);
		do_fio(&c__1, "U =", (long int)3);
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    do_fio(&c__1, (char *)&u[i__], (long int)sizeof(double))
			    ;
		}
		e_wsfe();
	    }
	}*/
    }
    return 0;
} /* prn1lb_ */

/* ======================= The end of prn1lb ============================= */
/* Subroutine */ int prn2lb_(long int *n, double *x, double *f, double *g, long int *iprint, long int *itfile,
			 long int *iter, long int *nfgv, long int *nact, double *sbgnrm, long int *nint, char *word,
			 long int *iword, long int *iback, double *stp, double *xstep, long int word_len)
/*long int *n;
double *x, *f, *g;
long int *iprint, *itfile, *iter, *nfgv, *nact;
double *sbgnrm;
long int *nint;
char *word;
long int *iword, *iback;
double *stp, *xstep;
long int word_len;*/
{
    /* Format strings *//*
    static char fmt_2001[] = "(/,\002At iterate\002,i5,4x,\002f= \002,1p,d12\
.5,4x,\002|proj g|= \002,1p,d12.5)";
    static char fmt_1004[] = "(/,a4,1p,6(1x,d11.4),/,(4x,1p,6(1x,d11.4)))";
    static char fmt_3001[] = "(2(1x,i4),2(1x,i5),2x,a3,1x,i4,1p,2(2x,d7.1),1\
p,2(1x,d10.3))";*/

    /* System generated locals */
//    long int i__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, const char *const, long int, long int);
 //   long int s_wsle(cilist *), do_lio(long int *, long int *, char *, long int), e_wsle(), s_wsfe(cilist *), do_fio(long int*, char*, long int), e_wsfe();

    /* Local variables */
//    static long int imod, i__;

    /* Fortran I/O blocks */
/*    static cilist io___194 = { 0, 6, 0, 0, 0 };
    static cilist io___195 = { 0, 6, 0, fmt_2001, 0 };
    static cilist io___196 = { 0, 6, 0, fmt_1004, 0 };
    static cilist io___198 = { 0, 6, 0, fmt_1004, 0 };
    static cilist io___200 = { 0, 6, 0, fmt_2001, 0 };
    static cilist io___201 = { 0, 0, 0, fmt_3001, 0 };*/


/*     ************ */

/*     Subroutine prn2lb */

/*     This subroutine prints out new information after a successful */
/*       line search. */


/*                           *  *  * */

/*     NEOS, November 1994. (Latest revision June 1996.) */
/*     Optimization Technology Center. */
/*     Argonne National Laboratory and Northwestern University. */
/*     Written by */
/*                        Ciyou Zhu */
/*     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal. */


/*     ************ */
/*           'word' records the status of subspace solutions. */
    /* Parameter adjustments */
    --g;
    --x;

    /* Function Body */
    if (*iword == 0) {
/*                            the subspace minimization converged. */
	s_copy(word, "con", (long int)3, (long int)3);
    } else if (*iword == 1) {
/*                          the subspace minimization stopped at a bound. */
	s_copy(word, "bnd", (long int)3, (long int)3);
    } else if (*iword == 5) {
/*                             the truncated Newton step has been used. */
	s_copy(word, "TNT", (long int)3, (long int)3);
    } else {
	s_copy(word, "---", (long int)3, (long int)3);
    }
    if (*iprint >= 99) {
/*	s_wsle(&io___194);
	do_lio(&c__9, &c__1, "LINE SEARCH", (long int)11);
	do_lio(&c__3, &c__1, (char *)&(*iback), (long int)sizeof(long int));
	do_lio(&c__9, &c__1, " times; norm of step = ", (long int)23);
	do_lio(&c__5, &c__1, (char *)&(*xstep), (long int)sizeof(double));
	e_wsle();
	s_wsfe(&io___195);
	do_fio(&c__1, (char *)&(*iter), (long int)sizeof(long int));
	do_fio(&c__1, (char *)&(*f), (long int)sizeof(double));
	do_fio(&c__1, (char *)&(*sbgnrm), (long int)sizeof(double));
	e_wsfe();
	if (*iprint > 100) {
	    s_wsfe(&io___196);
	    do_fio(&c__1, "X =", (long int)3);
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		do_fio(&c__1, (char *)&x[i__], (long int)sizeof(double));
	    }
	    e_wsfe();
	    s_wsfe(&io___198);
	    do_fio(&c__1, "G =", (long int)3);
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		do_fio(&c__1, (char *)&g[i__], (long int)sizeof(double));
	    }
	    e_wsfe();
	}*/
    } else if (*iprint > 0) {
/*	imod = *iter % *iprint;
	if (imod == 0) {
	    s_wsfe(&io___200);
	    do_fio(&c__1, (char *)&(*iter), (long int)sizeof(long int));
	    do_fio(&c__1, (char *)&(*f), (long int)sizeof(double));
	    do_fio(&c__1, (char *)&(*sbgnrm), (long int)sizeof(double));
	    e_wsfe();
	}*/
    }
    if (*iprint >= 1) {
/*	io___201.ciunit = *itfile;
	s_wsfe(&io___201);
	do_fio(&c__1, (char *)&(*iter), (long int)sizeof(long int));
	do_fio(&c__1, (char *)&(*nfgv), (long int)sizeof(long int));
	do_fio(&c__1, (char *)&(*nint), (long int)sizeof(long int));
	do_fio(&c__1, (char *)&(*nact), (long int)sizeof(long int));
	do_fio(&c__1, word, (long int)3);
	do_fio(&c__1, (char *)&(*iback), (long int)sizeof(long int));
	do_fio(&c__1, (char *)&(*stp), (long int)sizeof(double));
	do_fio(&c__1, (char *)&(*xstep), (long int)sizeof(double));
	do_fio(&c__1, (char *)&(*sbgnrm), (long int)sizeof(double));
	do_fio(&c__1, (char *)&(*f), (long int)sizeof(double));
	e_wsfe();*/
    }
    return 0;
} /* prn2lb_ */

/* ======================= The end of prn2lb ============================= */
/* Subroutine */ int prn3lb_(long int *n, double *x, double *f, char *task, long int *iprint, long int *info, 
			 long int *itfile, long int *iter, long int *nfgv, long int *nintol, long int *nskip, long int *nact,
			 double *sbgnrm, double *time, long int *nint, char *word, long int *iback, double *stp,
			 double *xstep, long int *k, double *cachyt, double *sbtime, double *lnscht,
			 long int task_len, long int word_len)
/*long int *n;
double *x, *f;
char *task;
long int *iprint, *info, *itfile, *iter, *nfgv, *nintol, *nskip, *nact;
double *sbgnrm, *time;
long int *nint;
char *word;
long int *iback;
double *stp, *xstep;
long int *k;
double *cachyt, *sbtime, *lnscht;
long int task_len;
long int word_len;*/
{
    /* Format strings *//*
    static char fmt_3003[] = "(/,\002           * * *\002,/,/,\002Tit   = to\
tal number of iterations\002,/,\002Tnf   = total number of function evaluati\
ons\002,/,\002Tnint = total number of segments explored during\002,\002 Cauc\
hy searches\002,/,\002Skip  = number of BFGS updates skipped\002,/,\002Nact \
 = number of active bounds at final generalized\002,\002 Cauchy point\002,/\
,\002Projg = norm of the final projected gradient\002,/,\002F     = final fu\
nction value\002,/,/,\002           * * *\002)";
    static char fmt_3004[] = "(/,3x,\002N\002,3x,\002Tit\002,2x,\002Tnf\002,\
2x,\002Tnint\002,2x,\002Skip\002,2x,\002Nact\002,5x,\002Projg\002,8x,\002\
F\002)";
    static char fmt_3005[] = "(i5,2(1x,i4),(1x,i6),(2x,i4),(1x,i5),1p,2(2x,d\
10.3))";
    static char fmt_1004[] = "(/,a4,1p,6(1x,d11.4),/,(4x,1p,6(1x,d11.4)))";
    static char fmt_3009[] = "(/,a60)";
    static char fmt_9011[] = "(/,\002 Matrix in 1st Cholesky factorization i\
n formk is not Pos. Def.\002)";
    static char fmt_9012[] = "(/,\002 Matrix in 2st Cholesky factorization i\
n formk is not Pos. Def.\002)";
    static char fmt_9013[] = "(/,\002 Matrix in the Cholesky factorization i\
n formt is not Pos. Def.\002)";
    static char fmt_9014[] = "(/,\002 Derivative >= 0, backtracking line sea\
rch impossible.\002,/,\002   Previous x, f and g restored.\002,/,\002 Possib\
le causes: 1 error in function or gradient evaluation;\002,/,\002           \
       2 rounding errors dominate computation.\002)";
    static char fmt_9015[] = "(/,\002 Warning:  more than 10 function and gr\
adient\002,/,\002   evaluations in the last line search.  Termination\002,/\
,\002   may possibly be caused by a bad search direction.\002)";
    static char fmt_9018[] = "(/,\002 The triangular system is singular.\002)"
	    ;
    static char fmt_9019[] = "(/,\002 Line search cannot locate an adequate \
point after 20 function\002,/,\002  and gradient evaluations.  Previous x, f\
 and g restored.\002,/,\002 Possible causes: 1 error in function or gradient\
 evaluation;\002,/,\002                  2 rounding error dominate computati\
on.\002)";
    static char fmt_3007[] = "(/,\002 Cauchy                time\002,1p,e10.\
3,\002 seconds.\002,/\002 Subspace minimization time\002,1p,e10.3,\002 secon\
ds.\002,/\002 Line search           time\002,1p,e10.3,\002 seconds.\002)";
    static char fmt_3008[] = "(/,\002 Total User time\002,1p,e10.3,\002 seco\
nds.\002,/)";
    static char fmt_3002[] = "(2(1x,i4),2(1x,i5),2x,a3,1x,i4,1p,2(2x,d7.1),6\
x,\002-\002,10x,\002-\002)";*/

    /* System generated locals */
//    long int i__1;

    /* Builtin functions */
    long int s_cmp(char *, const char *const, long int, long int);
//    long int s_wsfe(cilist *), e_wsfe(), do_fio(long int*, char*, long int), s_wsle(cilist *), do_lio(long int *, long int *, char *, long int), e_wsle();

    /* Local variables */
//    static long int i__;

    /* Fortran I/O blocks */
/*    static cilist io___202 = { 0, 6, 0, fmt_3003, 0 };
    static cilist io___203 = { 0, 6, 0, fmt_3004, 0 };
    static cilist io___204 = { 0, 6, 0, fmt_3005, 0 };
    static cilist io___205 = { 0, 6, 0, fmt_1004, 0 };
    static cilist io___207 = { 0, 6, 0, 0, 0 };
    static cilist io___208 = { 0, 6, 0, fmt_3009, 0 };
    static cilist io___209 = { 0, 6, 0, fmt_9011, 0 };
    static cilist io___210 = { 0, 6, 0, fmt_9012, 0 };
    static cilist io___211 = { 0, 6, 0, fmt_9013, 0 };
    static cilist io___212 = { 0, 6, 0, fmt_9014, 0 };
    static cilist io___213 = { 0, 6, 0, fmt_9015, 0 };
    static cilist io___214 = { 0, 6, 0, 0, 0 };
    static cilist io___215 = { 0, 6, 0, 0, 0 };
    static cilist io___216 = { 0, 6, 0, fmt_9018, 0 };
    static cilist io___217 = { 0, 6, 0, fmt_9019, 0 };
    static cilist io___218 = { 0, 6, 0, fmt_3007, 0 };
    static cilist io___219 = { 0, 6, 0, fmt_3008, 0 };
    static cilist io___220 = { 0, 0, 0, fmt_3002, 0 };
    static cilist io___221 = { 0, 0, 0, fmt_3009, 0 };
    static cilist io___222 = { 0, 0, 0, fmt_9011, 0 };
    static cilist io___223 = { 0, 0, 0, fmt_9012, 0 };
    static cilist io___224 = { 0, 0, 0, fmt_9013, 0 };
    static cilist io___225 = { 0, 0, 0, fmt_9014, 0 };
    static cilist io___226 = { 0, 0, 0, fmt_9015, 0 };
    static cilist io___227 = { 0, 0, 0, fmt_9018, 0 };
    static cilist io___228 = { 0, 0, 0, fmt_9019, 0 };
    static cilist io___229 = { 0, 0, 0, fmt_3008, 0 };*/


/*     ************ */

/*     Subroutine prn3lb */

/*     This subroutine prints out information when either a built-in */
/*       convergence test is satisfied or when an error message is */
/*       generated. */


/*                           *  *  * */

/*     NEOS, November 1994. (Latest revision June 1996.) */
/*     Optimization Technology Center. */
/*     Argonne National Laboratory and Northwestern University. */
/*     Written by */
/*                        Ciyou Zhu */
/*     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal. */


/*     ************ */
    /* Parameter adjustments */
    --x;

    /* Function Body */
    if (s_cmp(task, "ERROR", (long int)5, (long int)5) == 0) {
	goto L999;
    }
    if (*iprint >= 0) {
/*	s_wsfe(&io___202);
	e_wsfe();
	s_wsfe(&io___203);
	e_wsfe();
	s_wsfe(&io___204);
	do_fio(&c__1, (char *)&(*n), (long int)sizeof(long int));
	do_fio(&c__1, (char *)&(*iter), (long int)sizeof(long int));
	do_fio(&c__1, (char *)&(*nfgv), (long int)sizeof(long int));
	do_fio(&c__1, (char *)&(*nintol), (long int)sizeof(long int));
	do_fio(&c__1, (char *)&(*nskip), (long int)sizeof(long int));
	do_fio(&c__1, (char *)&(*nact), (long int)sizeof(long int));
	do_fio(&c__1, (char *)&(*sbgnrm), (long int)sizeof(double));
	do_fio(&c__1, (char *)&(*f), (long int)sizeof(double));
	e_wsfe();
	if (*iprint >= 100) {
	    s_wsfe(&io___205);
	    do_fio(&c__1, "X =", (long int)3);
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		do_fio(&c__1, (char *)&x[i__], (long int)sizeof(double));
	    }
	    e_wsfe();
	}
	if (*iprint >= 1) {
	    s_wsle(&io___207);
	    do_lio(&c__9, &c__1, " F =", (long int)4);
	    do_lio(&c__5, &c__1, (char *)&(*f), (long int)sizeof(double));
	    e_wsle();
	}*/
    }
L999:
    if (*iprint >= 0) {
/*	s_wsfe(&io___208);
	do_fio(&c__1, task, (long int)60);
	e_wsfe();
	if (*info != 0) {
	    if (*info == -1) {
		s_wsfe(&io___209);
		e_wsfe();
	    }
	    if (*info == -2) {
		s_wsfe(&io___210);
		e_wsfe();
	    }
	    if (*info == -3) {
		s_wsfe(&io___211);
		e_wsfe();
	    }
	    if (*info == -4) {
		s_wsfe(&io___212);
		e_wsfe();
	    }
	    if (*info == -5) {
		s_wsfe(&io___213);
		e_wsfe();
	    }
	    if (*info == -6) {
		s_wsle(&io___214);
		do_lio(&c__9, &c__1, " Input nbd(", (long int)11);
		do_lio(&c__3, &c__1, (char *)&(*k), (long int)sizeof(long int));
		do_lio(&c__9, &c__1, ") is invalid.", (long int)13);
		e_wsle();
	    }
	    if (*info == -7) {
		s_wsle(&io___215);
		do_lio(&c__9, &c__1, " l(", (long int)3);
		do_lio(&c__3, &c__1, (char *)&(*k), (long int)sizeof(long int));
		do_lio(&c__9, &c__1, ") > u(", (long int)6);
		do_lio(&c__3, &c__1, (char *)&(*k), (long int)sizeof(long int));
		do_lio(&c__9, &c__1, ").  No feasible solution.", (long int)25);
		e_wsle();
	    }
	    if (*info == -8) {
		s_wsfe(&io___216);
		e_wsfe();
	    }
	    if (*info == -9) {
		s_wsfe(&io___217);
		e_wsfe();
	    }
	}
	if (*iprint >= 1) {
	    s_wsfe(&io___218);
	    do_fio(&c__1, (char *)&(*cachyt), (long int)sizeof(double));
	    do_fio(&c__1, (char *)&(*sbtime), (long int)sizeof(double));
	    do_fio(&c__1, (char *)&(*lnscht), (long int)sizeof(double));
	    e_wsfe();
	}
	s_wsfe(&io___219);
	do_fio(&c__1, (char *)&(*time), (long int)sizeof(double));
	e_wsfe();
	if (*iprint >= 1) {
	    if (*info == -4 || *info == -9) {
		io___220.ciunit = *itfile;
		s_wsfe(&io___220);
		do_fio(&c__1, (char *)&(*iter), (long int)sizeof(long int));
		do_fio(&c__1, (char *)&(*nfgv), (long int)sizeof(long int));
		do_fio(&c__1, (char *)&(*nint), (long int)sizeof(long int));
		do_fio(&c__1, (char *)&(*nact), (long int)sizeof(long int));
		do_fio(&c__1, word, (long int)3);
		do_fio(&c__1, (char *)&(*iback), (long int)sizeof(long int));
		do_fio(&c__1, (char *)&(*stp), (long int)sizeof(double));
		do_fio(&c__1, (char *)&(*xstep), (long int)sizeof(double));
		e_wsfe();
	    }
	    io___221.ciunit = *itfile;
	    s_wsfe(&io___221);
	    do_fio(&c__1, task, (long int)60);
	    e_wsfe();
	    if (*info != 0) {
		if (*info == -1) {
		    io___222.ciunit = *itfile;
		    s_wsfe(&io___222);
		    e_wsfe();
		}
		if (*info == -2) {
		    io___223.ciunit = *itfile;
		    s_wsfe(&io___223);
		    e_wsfe();
		}
		if (*info == -3) {
		    io___224.ciunit = *itfile;
		    s_wsfe(&io___224);
		    e_wsfe();
		}
		if (*info == -4) {
		    io___225.ciunit = *itfile;
		    s_wsfe(&io___225);
		    e_wsfe();
		}
		if (*info == -5) {
		    io___226.ciunit = *itfile;
		    s_wsfe(&io___226);
		    e_wsfe();
		}
		if (*info == -8) {
		    io___227.ciunit = *itfile;
		    s_wsfe(&io___227);
		    e_wsfe();
		}
		if (*info == -9) {
		    io___228.ciunit = *itfile;
		    s_wsfe(&io___228);
		    e_wsfe();
		}
	    }
	    io___229.ciunit = *itfile;
	    s_wsfe(&io___229);
	    do_fio(&c__1, (char *)&(*time), (long int)sizeof(double));
	    e_wsfe();
	}*/
    }
/* L3006: */
    return 0;
} /* prn3lb_ */

/* ======================= The end of prn3lb ============================= */
/* Subroutine */ int projgr_(long int *n, double *l, double *u, long int *nbd, double *x, double *g, double *sbgnrm)
/*long int *n;
double *l, *u;
long int *nbd;
double *x, *g, *sbgnrm;*/
{
    /* System generated locals */
    long int i__1;
    double d__1, d__2;

    /* Local variables */
    static long int i__;
    static double gi;

/*     ************ */

/*     Subroutine projgr */

/*     This subroutine computes the infinity norm of the projected */
/*       gradient. */


/*                           *  *  * */

/*     NEOS, November 1994. (Latest revision June 1996.) */
/*     Optimization Technology Center. */
/*     Argonne National Laboratory and Northwestern University. */
/*     Written by */
/*                        Ciyou Zhu */
/*     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal. */


/*     ************ */
    /* Parameter adjustments */
    --g;
    --x;
    --nbd;
    --u;
    --l;

    /* Function Body */
    *sbgnrm = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	gi = g[i__];
	if (nbd[i__] != 0) {
	    if (gi < 0.) {
		if (nbd[i__] >= 2) {
/* Computing MAX */
		    d__1 = x[i__] - u[i__];
		    gi = max_(d__1,gi);
		}
	    } else {
		if (nbd[i__] <= 2) {
/* Computing MIN */
		    d__1 = x[i__] - l[i__];
		    gi = min_(d__1,gi);
		}
	    }
	}
/* Computing MAX */
	d__1 = *sbgnrm, d__2 = abs_(gi);
	*sbgnrm = max_(d__1,d__2);
/* L15: */
    }
    return 0;
} /* projgr_ */

/* ======================= The end of projgr ============================= */
/* Subroutine */ int subsm_(long int *n, long int *m, long int *nsub, long int *ind, double *l, double *u, long int *nbd,
			double *x, double *d__, double *ws, double *wy, double *theta, long int *col,
			long int *head, long int *iword, double *wv, double *wn, long int *iprint, long int *info)
/*long int *n, *m, *nsub, *ind;
double *l, *u;
long int *nbd;
double *x, *d__, *ws, *wy, *theta;
long int *col, *head, *iword;
double *wv, *wn;
long int *iprint, *info;*/
{
    /* Format strings *//*
    static char fmt_1001[] = "(/,\002----------------SUBSM entered----------\
-------\002,/)";
    static char fmt_1002[] = "(\002ALPHA = \002,f7.5,\002 backtrack to the B\
OX\002)";
    static char fmt_1003[] = "(\002Subspace solution X =  \002,/,(4x,1p,6(1x\
,d11.4)))";
    static char fmt_1004[] = "(/,\002----------------exit SUBSM ------------\
--------\002,/)";*/

    /* System generated locals */
    long int ws_dim1, ws_offset, wy_dim1, wy_offset, wn_dim1, wn_offset, i__1, 
	    i__2;

    /* Builtin functions */
//    long int s_wsfe(cilist *), e_wsfe(), do_fio(long int*, char*, long int), s_wsle(cilist *), do_lio(long int *, long int *, char *, long int), e_wsle();

    /* Local variables */
    static double temp1, temp2;
    static long int i__, j, k;
    static double alpha;
    //extern /* Subroutine */ int dtrsl_();
    extern int dtrsl_(double *t, long int *ldt, long int *n, double *b, long int *job, long int *info);
    static long int m2;
    static double dk;
    static long int js, jy, pointr, ibd, col2;

    /* Fortran I/O blocks */
/*    static cilist io___232 = { 0, 6, 0, fmt_1001, 0 };
    static cilist io___246 = { 0, 6, 0, fmt_1002, 0 };
    static cilist io___247 = { 0, 6, 0, 0, 0 };
    static cilist io___248 = { 0, 6, 0, fmt_1003, 0 };
    static cilist io___249 = { 0, 6, 0, fmt_1004, 0 }; */


/*     ************ */

/*     Subroutine subsm */

/*     Given xcp, l, u, r, an index set that specifies */
/* 	the active set at xcp, and an l-BFGS matrix B */
/* 	(in terms of WY, WS, SY, WT, head, col, and theta), */
/* 	this subroutine computes an approximate solution */
/* 	of the subspace problem */

/*     	(P)   min Q(x) = r'(x-xcp) + 1/2 (x-xcp)' B (x-xcp) */

/*             subject to l<=x<=u */
/* 	  	        x_i=xcp_i for all i in A(xcp) */

/* 	along the subspace unconstrained Newton direction */

/* 	   d = -(Z'BZ)^(-1) r. */

/*       The formula for the Newton direction, given the L-BFGS matrix */
/*       and the Sherman-Morrison formula, is */

/* 	   d = (1/theta)r + (1/theta*2) Z'WK^(-1)W'Z r. */

/*       where */
/*                 K = [-D -Y'ZZ'Y/theta     L_a'-R_z'  ] */
/*                     [L_a -R_z           theta*S'AA'S ] */

/*     Note that this procedure for computing d differs */
/*     from that described in [1]. One can show that the matrix K is */
/*     equal to the matrix M^[-1]N in that paper. */

/*     n is an long int variable. */
/*       On entry n is the dimension of the problem. */
/*       On exit n is unchanged. */

/*     m is an long int variable. */
/*       On entry m is the maximum number of variable metric corrections */
/*         used to define the limited memory matrix. */
/*       On exit m is unchanged. */

/*     nsub is an long int variable. */
/*       On entry nsub is the number of free variables. */
/*       On exit nsub is unchanged. */

/*     ind is an long int array of dimension nsub. */
/*       On entry ind specifies the coordinate indices of free variables. */
/*       On exit ind is unchanged. */

/*     l is a double precision array of dimension n. */
/*       On entry l is the lower bound of x. */
/*       On exit l is unchanged. */

/*     u is a double precision array of dimension n. */
/*       On entry u is the upper bound of x. */
/*       On exit u is unchanged. */

/*     nbd is a long int array of dimension n. */
/*       On entry nbd represents the type of bounds imposed on the */
/*         variables, and must be specified as follows: */
/*         nbd(i)=0 if x(i) is unbounded, */
/*                1 if x(i) has only a lower bound, */
/*                2 if x(i) has both lower and upper bounds, and */
/*                3 if x(i) has only an upper bound. */
/*       On exit nbd is unchanged. */

/*     x is a double precision array of dimension n. */
/*       On entry x specifies the Cauchy point xcp. */
/*       On exit x(i) is the minimizer of Q over the subspace of */
/*                                                        free variables. */

/*     d is a double precision array of dimension n. */
/*       On entry d is the reduced gradient of Q at xcp. */
/*       On exit d is the Newton direction of Q. */

/*     ws and wy are double precision arrays; */
/*     theta is a double precision variable; */
/*     col is an long int variable; */
/*     head is an long int variable. */
/*       On entry they store the information defining the */
/*                                          limited memory BFGS matrix: */
/*         ws(n,m) stores S, a set of s-vectors; */
/*         wy(n,m) stores Y, a set of y-vectors; */
/*         theta is the scaling factor specifying B_0 = theta I; */
/*         col is the number of variable metric corrections stored; */
/*         head is the location of the 1st s- (or y-) vector in S (or Y). */
/*       On exit they are unchanged. */

/*     iword is an long int variable. */
/*       On entry iword is unspecified. */
/*       On exit iword specifies the status of the subspace solution. */
/*         iword = 0 if the solution is in the box, */
/*                 1 if some bound is encountered. */

/*     wv is a double precision working array of dimension 2m. */

/*     wn is a double precision array of dimension 2m x 2m. */
/*       On entry the upper triangle of wn stores the LEL^T factorization */
/*         of the indefinite matrix */

/*              K = [-D -Y'ZZ'Y/theta     L_a'-R_z'  ] */
/*                  [L_a -R_z           theta*S'AA'S ] */
/*                                                    where E = [-I  0] */
/*                                                              [ 0  I] */
/*       On exit wn is unchanged. */

/*     iprint is an long int variable that must be set by the user. */
/*       It controls the frequency and type of output generated: */
/*        iprint<0    no output is generated; */
/*        iprint=0    print only one line at the last iteration; */
/*        0<iprint<99 print also f and |proj g| every iprint iterations; */
/*        iprint=99   print details of every iteration except n-vectors; */
/*        iprint=100  print also the changes of active set and final x; */
/*        iprint>100  print details of every iteration including x and g; */
/*       When iprint > 0, the file iterate.dat will be created to */
/*                        summarize the iteration. */

/*     info is an long int variable. */
/*       On entry info is unspecified. */
/*       On exit info = 0       for normal return, */
/*                    = nonzero for abnormal return */
/*                                  when the matrix K is ill-conditioned. */

/*     Subprograms called: */

/*       Linpack dtrsl. */


/*     References: */

/*       [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited */
/*       memory algorithm for bound constrained optimization'', */
/*       SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208. */



/*                           *  *  * */

/*     NEOS, November 1994. (Latest revision June 1996.) */
/*     Optimization Technology Center. */
/*     Argonne National Laboratory and Northwestern University. */
/*     Written by */
/*                        Ciyou Zhu */
/*     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal. */


/*     ************ */
    /* Parameter adjustments */
    --d__;
    --x;
    --nbd;
    --u;
    --l;
    wn_dim1 = 2 * *m;
    wn_offset = 1 + wn_dim1 * 1;
    wn -= wn_offset;
    --wv;
    wy_dim1 = *n;
    wy_offset = 1 + wy_dim1 * 1;
    wy -= wy_offset;
    ws_dim1 = *n;
    ws_offset = 1 + ws_dim1 * 1;
    ws -= ws_offset;
    --ind;

    /* Function Body */
    if (*nsub <= 0) {
	return 0;
    }
    if (*iprint >= 99) {
/*	s_wsfe(&io___232);
	e_wsfe();*/
    }
/*     Compute wv = W'Zd. */
    pointr = *head;
    i__1 = *col;
    for (i__ = 1; i__ <= i__1; ++i__) {
	temp1 = 0.;
	temp2 = 0.;
	i__2 = *nsub;
	for (j = 1; j <= i__2; ++j) {
	    k = ind[j];
	    temp1 += wy[k + pointr * wy_dim1] * d__[j];
	    temp2 += ws[k + pointr * ws_dim1] * d__[j];
/* L10: */
	}
	wv[i__] = temp1;
	wv[*col + i__] = *theta * temp2;
	pointr = pointr % *m + 1;
/* L20: */
    }
/*     Compute wv:=K^(-1)wv. */
    m2 = *m << 1;
    col2 = *col << 1;
    dtrsl_(&wn[wn_offset], &m2, &col2, &wv[1], &c__11, info);
    if (*info != 0) {
	return 0;
    }
    i__1 = *col;
    for (i__ = 1; i__ <= i__1; ++i__) {
	wv[i__] = -wv[i__];
/* L25: */
    }
    dtrsl_(&wn[wn_offset], &m2, &col2, &wv[1], &c__1, info);
    if (*info != 0) {
	return 0;
    }
/*     Compute d = (1/theta)d + (1/theta**2)Z'W wv. */
    pointr = *head;
    i__1 = *col;
    for (jy = 1; jy <= i__1; ++jy) {
	js = *col + jy;
	i__2 = *nsub;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    k = ind[i__];
	    d__[i__] = d__[i__] + wy[k + pointr * wy_dim1] * wv[jy] / *theta 
		    + ws[k + pointr * ws_dim1] * wv[js];
/* L30: */
	}
	pointr = pointr % *m + 1;
/* L40: */
    }
    i__1 = *nsub;
    for (i__ = 1; i__ <= i__1; ++i__) {
	d__[i__] /= *theta;
/* L50: */
    }
/*     Backtrack to the feasible region. */
    alpha = 1.;
    temp1 = alpha;
    i__1 = *nsub;
    for (i__ = 1; i__ <= i__1; ++i__) {
	k = ind[i__];
	dk = d__[i__];
	if (nbd[k] != 0) {
	    if (dk < 0. && nbd[k] <= 2) {
		temp2 = l[k] - x[k];
		if (temp2 >= 0.) {
		    temp1 = 0.;
		} else if (dk * alpha < temp2) {
		    temp1 = temp2 / dk;
		}
	    } else if (dk > 0. && nbd[k] >= 2) {
		temp2 = u[k] - x[k];
		if (temp2 <= 0.) {
		    temp1 = 0.;
		} else if (dk * alpha > temp2) {
		    temp1 = temp2 / dk;
		}
	    }
	    if (temp1 < alpha) {
		alpha = temp1;
		ibd = i__;
	    }
	}
/* L60: */
    }
    if (alpha < 1.) {
	dk = d__[ibd];
	k = ind[ibd];
	if (dk > 0.) {
	    x[k] = u[k];
	    d__[ibd] = 0.;
	} else if (dk < 0.) {
	    x[k] = l[k];
	    d__[ibd] = 0.;
	}
    }
    i__1 = *nsub;
    for (i__ = 1; i__ <= i__1; ++i__) {
	k = ind[i__];
	x[k] += alpha * d__[i__];
/* L70: */
    }
    if (*iprint >= 99) {
/*	if (alpha < 1.) {
	    s_wsfe(&io___246);
	    do_fio(&c__1, (char *)&alpha, (long int)sizeof(double));
	    e_wsfe();
	} else { 
	    s_wsle(&io___247);
	    do_lio(&c__9, &c__1, "SM solution inside the box", (long int)26);
	    e_wsle();
	}
	if (*iprint > 100) {
	    s_wsfe(&io___248);
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		do_fio(&c__1, (char *)&x[i__], (long int)sizeof(double));
	    }
	    e_wsfe();
	}*/
    }
    if (alpha < 1.) {
	*iword = 1;
    } else {
	*iword = 0;
    }
    if (*iprint >= 99) {
/*	s_wsfe(&io___249);
	e_wsfe();*/
    }
    return 0;
} /* subsm_ */

/* ====================== The end of subsm =============================== */
/* Subroutine */ int dcsrch_(double *f, double *g, double *stp, double *ftol, double *gtol, double *xtol,
			 double *stpmin, double *stpmax, char *task, long int *isave, double *dsave, long int task_len)
/*double *f, *g, *stp, *ftol, *gtol, *xtol, *stpmin, *stpmax;
char *task;
long int *isave;
double *dsave;
long int task_len;*/
{
    /* System generated locals */
    double d__1;

    /* Builtin functions */
    long int s_cmp(char *, const char *const, long int, long int);
    /* Subroutine */ int s_copy(char *, const char *const, long int, long int);

    /* Local variables */
    static long int stage;
    static double finit, ginit, width, ftest, gtest, stmin, stmax, width1,
	     fm, gm, fx, fy, gx, gy;
    static long int brackt;
    //extern /* Subroutine */ int dcstep_();
    extern int dcstep_(double *stx, double *fx, double *dx, double *sty, double *fy, double *dy,
		     double *stp, double *fp, double *dp, long int *brackt, double *stpmin, double *stpmax);

    static double fxm, fym, gxm, gym, stx, sty;

/*     ********** */

/*     Subroutine dcsrch */

/*     This subroutine finds a step that satisfies a sufficient */
/*     decrease condition and a curvature condition. */

/*     Each call of the subroutine updates an interval with */
/*     endpoints stx and sty. The interval is initially chosen */
/*     so that it contains a minimizer of the modified function */

/*           psi(stp) = f(stp) - f(0) - ftol*stp*f'(0). */

/*     If psi(stp) <= 0 and f'(stp) >= 0 for some step, then the */
/*     interval is chosen so that it contains a minimizer of f. */

/*     The algorithm is designed to find a step that satisfies */
/*     the sufficient decrease condition */

/*           f(stp) <= f(0) + ftol*stp*f'(0), */

/*     and the curvature condition */

/*           abs_(f'(stp)) <= gtol*abs_(f'(0)). */

/*     If ftol is less than gtol and if, for example, the function */
/*     is bounded below, then there is always a step which satisfies */
/*     both conditions. */

/*     If no step can be found that satisfies both conditions, then */
/*     the algorithm stops with a warning. In this case stp only */
/*     satisfies the sufficient decrease condition. */

/*     A typical invocation of dcsrch has the following outline: */

/*     task = 'START' */
/*  10 continue */
/*        call dcsrch( ... ) */
/*        if (task .eq. 'FG') then */
/*           Evaluate the function and the gradient at stp */
/*           goto 10 */
/*           end if */

/*     NOTE: The user must no alter work arrays between calls. */

/*     The subroutine statement is */

/*        subroutine dcsrch(f,g,stp,ftol,gtol,xtol,stpmin,stpmax, */
/*                          task,isave,dsave) */
/*     where */

/*       f is a double precision variable. */
/*         On initial entry f is the value of the function at 0. */
/*            On subsequent entries f is the value of the */
/*            function at stp. */
/*         On exit f is the value of the function at stp. */

/* 	g is a double precision variable. */
/*         On initial entry g is the derivative of the function at 0. */
/*            On subsequent entries g is the derivative of the */
/*            function at stp. */
/*         On exit g is the derivative of the function at stp. */

/* 	stp is a double precision variable. */
/*         On entry stp is the current estimate of a satisfactory */
/*            step. On initial entry, a positive initial estimate */
/*            must be provided. */
/*         On exit stp is the current estimate of a satisfactory step */
/*            if task = 'FG'. If task = 'CONV' then stp satisfies */
/*            the sufficient decrease and curvature condition. */

/*       ftol is a double precision variable. */
/*         On entry ftol specifies a nonnegative tolerance for the */
/*            sufficient decrease condition. */
/*         On exit ftol is unchanged. */

/*       gtol is a double precision variable. */
/*         On entry gtol specifies a nonnegative tolerance for the */
/*            curvature condition. */
/*         On exit gtol is unchanged. */

/* 	xtol is a double precision variable. */
/*         On entry xtol specifies a nonnegative relative tolerance */
/*            for an acceptable step. The subroutine exits with a */
/*            warning if the relative difference between sty and stx */
/*            is less than xtol. */
/*         On exit xtol is unchanged. */

/* 	stpmin is a double precision variable. */
/*         On entry stpmin is a nonnegative lower bound for the step. */
/*         On exit stpmin is unchanged. */

/* 	stpmax is a double precision variable. */
/*         On entry stpmax is a nonnegative upper bound for the step. */
/*         On exit stpmax is unchanged. */

/*       task is a character variable of length at least 60. */
/*         On initial entry task must be set to 'START'. */
/*         On exit task indicates the required action: */

/*            If task(1:2) = 'FG' then evaluate the function and */
/*            derivative at stp and call dcsrch again. */

/*            If task(1:4) = 'CONV' then the search is successful. */

/*            If task(1:4) = 'WARN' then the subroutine is not able */
/*            to satisfy the convergence conditions. The exit value of */
/*            stp contains the best point found during the search. */

/*            If task(1:5) = 'ERROR' then there is an error in the */
/*            input arguments. */

/*         On exit with convergence, a warning or an error, the */
/*            variable task contains additional information. */

/*       isave is an long int work array of dimension 2. */

/*       dsave is a double precision work array of dimension 13. */

/*     Subprograms called */

/* 	MINPACK-2 ... dcstep */

/*     MINPACK-1 Project. June 1983. */
/*     Argonne National Laboratory. */
/*     Jorge J. More' and David J. Thuente. */

/*     MINPACK-2 Project. October 1993. */
/*     Argonne National Laboratory and University of Minnesota. */
/*     Brett M. Averick, Richard G. Carter, and Jorge J. More'. */

/*     ********** */
/*     Initialization block. */
    /* Parameter adjustments */
    --dsave;
    --isave;

    /* Function Body */
    if (s_cmp(task, "START", (long int)5, (long int)5) == 0) {
/*        Check the input arguments for errors. */
	if (*stp < *stpmin) {
	    s_copy(task, "ERROR: STP .LT. STPMIN", task_len, (long int)22);
	}
	if (*stp > *stpmax) {
	    s_copy(task, "ERROR: STP .GT. STPMAX", task_len, (long int)22);
	}
	if (*g >= 0.) {
	    s_copy(task, "ERROR: INITIAL G .GE. ZERO", task_len, (long int)26);
	}
	if (*ftol < 0.) {
	    s_copy(task, "ERROR: FTOL .LT. ZERO", task_len, (long int)21);
	}
	if (*gtol < 0.) {
	    s_copy(task, "ERROR: GTOL .LT. ZERO", task_len, (long int)21);
	}
	if (*xtol < 0.) {
	    s_copy(task, "ERROR: XTOL .LT. ZERO", task_len, (long int)21);
	}
	if (*stpmin < 0.) {
	    s_copy(task, "ERROR: STPMIN .LT. ZERO", task_len, (long int)23);
	}
	if (*stpmax < *stpmin) {
	    s_copy(task, "ERROR: STPMAX .LT. STPMIN", task_len, (long int)25);
	}
/*        Exit if there are errors on input. */
	if (s_cmp(task, "ERROR", (long int)5, (long int)5) == 0) {
	    return 0;
	}
/*        Initialize local variables. */
	brackt = FALSE_;
	stage = 1;
	finit = *f;
	ginit = *g;
	gtest = *ftol * ginit;
	width = *stpmax - *stpmin;
	width1 = width / .5;
/*        The variables stx, fx, gx contain the values of the step, */
/*        function, and derivative at the best step. */
/*        The variables sty, fy, gy contain the value of the step, */
/*        function, and derivative at sty. */
/*        The variables stp, f, g contain the values of the step, */
/*        function, and derivative at stp. */
	stx = 0.;
	fx = finit;
	gx = ginit;
	sty = 0.;
	fy = finit;
	gy = ginit;
	stmin = 0.;
	stmax = *stp + *stp * 4.;
	s_copy(task, "FG", task_len, (long int)2);
	goto L1000;
    } else {
/*        Restore local variables. */
	if (isave[1] == 1) {
	    brackt = TRUE_;
	} else {
	    brackt = FALSE_;
	}
	stage = isave[2];
	ginit = dsave[1];
	gtest = dsave[2];
	gx = dsave[3];
	gy = dsave[4];
	finit = dsave[5];
	fx = dsave[6];
	fy = dsave[7];
	stx = dsave[8];
	sty = dsave[9];
	stmin = dsave[10];
	stmax = dsave[11];
	width = dsave[12];
	width1 = dsave[13];
    }
/*     If psi(stp) <= 0 and f'(stp) >= 0 for some step, then the */
/*     algorithm enters the second stage. */
    ftest = finit + *stp * gtest;
    if (stage == 1 && *f <= ftest && *g >= 0.) {
	stage = 2;
    }
/*     Test for warnings. */
    if (brackt && (*stp <= stmin || *stp >= stmax)) {
	s_copy(task, "WARNING: ROUNDING ERRORS PREVENT PROGRESS", task_len, (
		long int)41);
    }
    if (brackt && stmax - stmin <= *xtol * stmax) {
	s_copy(task, "WARNING: XTOL TEST SATISFIED", task_len, (long int)28);
    }
    if (*stp == *stpmax && *f <= ftest && *g <= gtest) {
	s_copy(task, "WARNING: STP = STPMAX", task_len, (long int)21);
    }
    if (*stp == *stpmin && (*f > ftest || *g >= gtest)) {
	s_copy(task, "WARNING: STP = STPMIN", task_len, (long int)21);
    }
/*     Test for convergence. */
    if (*f <= ftest && abs_(*g) <= *gtol * (-ginit)) {
	s_copy(task, "CONVERGENCE", task_len, (long int)11);
    }
/*     Test for termination. */
    if (s_cmp(task, "WARN", (long int)4, (long int)4) == 0 || s_cmp(task, "CONV", 
	    (long int)4, (long int)4) == 0) {
	goto L1000;
    }
/*     A modified function is used to predict the step during the */
/*     first stage if a lower function value has been obtained but */
/*     the decrease is not sufficient. */
    if (stage == 1 && *f <= fx && *f > ftest) {
/*        Define the modified function and derivative values. */
	fm = *f - *stp * gtest;
	fxm = fx - stx * gtest;
	fym = fy - sty * gtest;
	gm = *g - gtest;
	gxm = gx - gtest;
	gym = gy - gtest;
/*        Call dcstep to update stx, sty, and to compute the new step. */
	dcstep_(&stx, &fxm, &gxm, &sty, &fym, &gym, stp, &fm, &gm, &brackt, &
		stmin, &stmax);
/*        Reset the function and derivative values for f. */
	fx = fxm + stx * gtest;
	fy = fym + sty * gtest;
	gx = gxm + gtest;
	gy = gym + gtest;
    } else {
/*       Call dcstep to update stx, sty, and to compute the new step. */
	dcstep_(&stx, &fx, &gx, &sty, &fy, &gy, stp, f, g, &brackt, &stmin, &
		stmax);
    }
/*     Decide if a bisection step is needed. */
    if (brackt) {
	if ((d__1 = sty - stx, abs_(d__1)) >= width1 * .66) {
	    *stp = stx + (sty - stx) * .5;
	}
	width1 = width;
	width = (d__1 = sty - stx, abs_(d__1));
    }
/*     Set the minimum and maximum steps allowed for stp. */
    if (brackt) {
	stmin = min_(stx,sty);
	stmax = max_(stx,sty);
    } else {
	stmin = *stp + (*stp - stx) * 1.1;
	stmax = *stp + (*stp - stx) * 4.;
    }
/*     Force the step to be within the bounds stpmax and stpmin. */
    *stp = max_(*stp,*stpmin);
    *stp = min_(*stp,*stpmax);
/*     If further progress is not possible, let stp be the best */
/*     point obtained during the search. */
    if (brackt && (*stp <= stmin || *stp >= stmax) || brackt && stmax - stmin 
	    <= *xtol * stmax) {
	*stp = stx;
    }
/*     Obtain another function and derivative. */
    s_copy(task, "FG", task_len, (long int)2);
L1000:
/*     Save local variables. */
    if (brackt) {
	isave[1] = 1;
    } else {
	isave[1] = 0;
    }
    isave[2] = stage;
    dsave[1] = ginit;
    dsave[2] = gtest;
    dsave[3] = gx;
    dsave[4] = gy;
    dsave[5] = finit;
    dsave[6] = fx;
    dsave[7] = fy;
    dsave[8] = stx;
    dsave[9] = sty;
    dsave[10] = stmin;
    dsave[11] = stmax;
    dsave[12] = width;
    dsave[13] = width1;
    return 0;
} /* dcsrch_ */

/* ====================== The end of dcsrch ============================== */
/* Subroutine */ int dcstep_(double *stx, double *fx, double *dx, double *sty, double *fy, double *dy,
			 double *stp, double *fp, double *dp, long int *brackt, double *stpmin, double *stpmax)
/*double *stx, *fx, *dx, *sty, *fy, *dy, *stp, *fp, *dp;
long int *brackt;
double *stpmin, *stpmax;*/
{
    /* System generated locals */
    double d__1, d__2, d__3;

    /* Builtin functions */
    //double sqrt();

    /* Local variables */
    static double sgnd, stpc, stpf, stpq, p, q, gamma, r__, s, theta;

/*     ********** */

/*     Subroutine dcstep */

/*     This subroutine computes a safeguarded step for a search */
/*     procedure and updates an interval that contains a step that */
/*     satisfies a sufficient decrease and a curvature condition. */

/*     The parameter stx contains the step with the least function */
/*     value. If brackt is set to .true. then a minimizer has */
/*     been bracketed in an interval with endpoints stx and sty. */
/*     The parameter stp contains the current step. */
/*     The subroutine assumes that if brackt is set to .true. then */

/*           min_(stx,sty) < stp < max_(stx,sty), */

/*     and that the derivative at stx is negative in the direction */
/*     of the step. */

/*     The subroutine statement is */

/*       subroutine dcstep(stx,fx,dx,sty,fy,dy,stp,fp,dp,brackt, */
/*                         stpmin,stpmax) */

/*     where */

/*       stx is a double precision variable. */
/*         On entry stx is the best step obtained so far and is an */
/*            endpoint of the interval that contains the minimizer. */
/*         On exit stx is the updated best step. */

/*       fx is a double precision variable. */
/*         On entry fx is the function at stx. */
/*         On exit fx is the function at stx. */

/*       dx is a double precision variable. */
/*         On entry dx is the derivative of the function at */
/*            stx. The derivative must be negative in the direction of */
/*            the step, that is, dx and stp - stx must have opposite */
/*            signs. */
/*         On exit dx is the derivative of the function at stx. */

/*       sty is a double precision variable. */
/*         On entry sty is the second endpoint of the interval that */
/*            contains the minimizer. */
/*         On exit sty is the updated endpoint of the interval that */
/*            contains the minimizer. */

/*       fy is a double precision variable. */
/*         On entry fy is the function at sty. */
/*         On exit fy is the function at sty. */

/*       dy is a double precision variable. */
/*         On entry dy is the derivative of the function at sty. */
/*         On exit dy is the derivative of the function at the exit sty. */

/*       stp is a double precision variable. */
/*         On entry stp is the current step. If brackt is set to .true. */
/*            then on input stp must be between stx and sty. */
/*         On exit stp is a new trial step. */

/*       fp is a double precision variable. */
/*         On entry fp is the function at stp */
/*         On exit fp is unchanged. */

/*       dp is a double precision variable. */
/*         On entry dp is the the derivative of the function at stp. */
/*         On exit dp is unchanged. */

/*       brackt is an long int variable. */
/*         On entry brackt specifies if a minimizer has been bracketed. */
/*            Initially brackt must be set to .false. */
/*         On exit brackt specifies if a minimizer has been bracketed. */
/*            When a minimizer is bracketed brackt is set to .true. */

/*       stpmin is a double precision variable. */
/*         On entry stpmin is a lower bound for the step. */
/*         On exit stpmin is unchanged. */

/*       stpmax is a double precision variable. */
/*         On entry stpmax is an upper bound for the step. */
/*         On exit stpmax is unchanged. */

/*     MINPACK-1 Project. June 1983 */
/*     Argonne National Laboratory. */
/*     Jorge J. More' and David J. Thuente. */

/*     MINPACK-2 Project. October 1993. */
/*     Argonne National Laboratory and University of Minnesota. */
/*     Brett M. Averick and Jorge J. More'. */

/*     ********** */
    sgnd = *dp * (*dx / abs_(*dx));
/*     First case: A higher function value. The minimum is bracketed. */
/*     If the cubic step is closer to stx than the quadratic step, the */
/*     cubic step is taken, otherwise the average of the cubic and */
/*     quadratic steps is taken. */
    if (*fp > *fx) {
	theta = (*fx - *fp) * 3. / (*stp - *stx) + *dx + *dp;
/* Computing MAX */
	d__1 = abs_(theta), d__2 = abs_(*dx), d__1 = max_(d__1,d__2), d__2 = abs_(
		*dp);
	s = max_(d__1,d__2);
/* Computing 2nd power */
	d__1 = theta / s;
	gamma = s * sqrt(d__1 * d__1 - *dx / s * (*dp / s));
	if (*stp < *stx) {
	    gamma = -gamma;
	}
	p = gamma - *dx + theta;
	q = gamma - *dx + gamma + *dp;
	r__ = p / q;
	stpc = *stx + r__ * (*stp - *stx);
	stpq = *stx + *dx / ((*fx - *fp) / (*stp - *stx) + *dx) / 2. * (*stp 
		- *stx);
	if ((d__1 = stpc - *stx, abs_(d__1)) < (d__2 = stpq - *stx, abs_(d__2)))
		 {
	    stpf = stpc;
	} else {
	    stpf = stpc + (stpq - stpc) / 2.;
	}
	*brackt = TRUE_;
/*     Second case: A lower function value and derivatives of opposite */
/*     sign. The minimum is bracketed. If the cubic step is farther from */
/*     stp than the secant step, the cubic step is taken, otherwise the */
/*     secant step is taken. */
    } else if (sgnd < 0.) {
	theta = (*fx - *fp) * 3. / (*stp - *stx) + *dx + *dp;
/* Computing MAX */
	d__1 = abs_(theta), d__2 = abs_(*dx), d__1 = max_(d__1,d__2), d__2 = abs_(
		*dp);
	s = max_(d__1,d__2);
/* Computing 2nd power */
	d__1 = theta / s;
	gamma = s * sqrt(d__1 * d__1 - *dx / s * (*dp / s));
	if (*stp > *stx) {
	    gamma = -gamma;
	}
	p = gamma - *dp + theta;
	q = gamma - *dp + gamma + *dx;
	r__ = p / q;
	stpc = *stp + r__ * (*stx - *stp);
	stpq = *stp + *dp / (*dp - *dx) * (*stx - *stp);
	if ((d__1 = stpc - *stp, abs_(d__1)) > (d__2 = stpq - *stp, abs_(d__2)))
		 {
	    stpf = stpc;
	} else {
	    stpf = stpq;
	}
	*brackt = TRUE_;
/*     Third case: A lower function value, derivatives of the same sign, */
/*     and the magnitude of the derivative decreases. */
    } else if (abs_(*dp) < abs_(*dx)) {
/*        The cubic step is computed only if the cubic tends to infinity */
/*        in the direction of the step or if the minimum of the cubic */
/*        is beyond stp. Otherwise the cubic step is defined to be the */
/*        secant step. */
	theta = (*fx - *fp) * 3. / (*stp - *stx) + *dx + *dp;
/* Computing MAX */
	d__1 = abs_(theta), d__2 = abs_(*dx), d__1 = max_(d__1,d__2), d__2 = abs_(
		*dp);
	s = max_(d__1,d__2);
/*        The case gamma = 0 only arises if the cubic does not tend */
/*        to infinity in the direction of the step. */
/* Computing MAX */
/* Computing 2nd power */
	d__3 = theta / s;
	d__1 = 0., d__2 = d__3 * d__3 - *dx / s * (*dp / s);
	gamma = s * sqrt((max_(d__1,d__2)));
	if (*stp > *stx) {
	    gamma = -gamma;
	}
	p = gamma - *dp + theta;
	q = gamma + (*dx - *dp) + gamma;
	r__ = p / q;
	if (r__ < 0. && gamma != 0.) {
	    stpc = *stp + r__ * (*stx - *stp);
	} else if (*stp > *stx) {
	    stpc = *stpmax;
	} else {
	    stpc = *stpmin;
	}
	stpq = *stp + *dp / (*dp - *dx) * (*stx - *stp);
	if (*brackt) {
/*           A minimizer has been bracketed. If the cubic step is */
/*           closer to stp than the secant step, the cubic step is */
/*           taken, otherwise the secant step is taken. */
	    if ((d__1 = stpc - *stp, abs_(d__1)) < (d__2 = stpq - *stp, abs_(
		    d__2))) {
		stpf = stpc;
	    } else {
		stpf = stpq;
	    }
	    if (*stp > *stx) {
/* Computing MIN */
		d__1 = *stp + (*sty - *stp) * .66;
		stpf = min_(d__1,stpf);
	    } else {
/* Computing MAX */
		d__1 = *stp + (*sty - *stp) * .66;
		stpf = max_(d__1,stpf);
	    }
	} else {
/*           A minimizer has not been bracketed. If the cubic step is */
/*           farther from stp than the secant step, the cubic step is */
/*           taken, otherwise the secant step is taken. */
	    if ((d__1 = stpc - *stp, abs_(d__1)) > (d__2 = stpq - *stp, abs_(
		    d__2))) {
		stpf = stpc;
	    } else {
		stpf = stpq;
	    }
	    stpf = min_(*stpmax,stpf);
	    stpf = max_(*stpmin,stpf);
	}
/*     Fourth case: A lower function value, derivatives of the same sign, */
/*     and the magnitude of the derivative does not decrease. If the */
/*     minimum is not bracketed, the step is either stpmin or stpmax, */
/*     otherwise the cubic step is taken. */
    } else {
	if (*brackt) {
	    theta = (*fp - *fy) * 3. / (*sty - *stp) + *dy + *dp;
/* Computing MAX */
	    d__1 = abs_(theta), d__2 = abs_(*dy), d__1 = max_(d__1,d__2), d__2 = 
		    abs_(*dp);
	    s = max_(d__1,d__2);
/* Computing 2nd power */
	    d__1 = theta / s;
	    gamma = s * sqrt(d__1 * d__1 - *dy / s * (*dp / s));
	    if (*stp > *sty) {
		gamma = -gamma;
	    }
	    p = gamma - *dp + theta;
	    q = gamma - *dp + gamma + *dy;
	    r__ = p / q;
	    stpc = *stp + r__ * (*sty - *stp);
	    stpf = stpc;
	} else if (*stp > *stx) {
	    stpf = *stpmax;
	} else {
	    stpf = *stpmin;
	}
    }
/*     Update the interval which contains a minimizer. */
    if (*fp > *fx) {
	*sty = *stp;
	*fy = *fp;
	*dy = *dp;
    } else {
	if (sgnd < 0.) {
	    *sty = *stx;
	    *fy = *fx;
	    *dy = *dx;
	}
	*stx = *stp;
	*fx = *fp;
	*dx = *dp;
    }
/*     Compute the new step. */
    *stp = stpf;
    return 0;
} /* dcstep_ */

/* ====================== The end of dcstep ============================== */
/* Subroutine */ int timer_(double *ttime)
//double *ttime;
{
 //   static float temp;
 //    extern double etime_(float *);
    static float tarray[2];

/*     ********* */

/*     Subroutine timer */

/*     This subroutine is used to determine user time. In a typical */
/*     application, the user time for a code segment requires calls */
/*     to subroutine timer to determine the initial and final time. */

/*     The subroutine statement is */

/*       subroutine timer(ttime) */

/*     where */

/*       ttime is an output variable which specifies the user time. */

/*     Argonne National Laboratory and University of Minnesota. */
/*     MINPACK-2 Project. */

/*     Modified October 1990 by Brett M. Averick. */

/*     ********** */
/*     The first element of the array tarray specifies user time */
    //temp = etime_(tarray);
    tarray[0]=0.0;
    tarray[1]=1.0;
    *ttime = (double) tarray[0];
    return 0;
} /* timer_ */

/* ====================== The end of timer =============================== */
double dnrm2_(long int *n, double *x, long int *incx)
/*long int *n;
double *x;
long int *incx;*/
{
    /* System generated locals */
    long int i__1, i__2;
    double ret_val, d__1, d__2, d__3;

    /* Builtin functions */
    //double sqrt();

    /* Local variables */
    static long int i__;
    static double scale;

/*     ********** */

/*     Function dnrm2 */

/*     Given a vector x of length n, this function calculates the */
/*     Euclidean norm of x with stride incx. */

/*     The function statement is */

/*       double precision function dnrm2(n,x,incx) */

/*     where */

/*       n is a positive long int input variable. */

/*       x is an input array of length n. */

/*       incx is a positive long int variable that specifies the */
/*         stride of the vector. */

/*     Subprograms called */

/*       FORTRAN-supplied ... abs, max, sqrt */

/*     MINPACK-2 Project. February 1991. */
/*     Argonne National Laboratory. */
/*     Brett M. Averick. */

/*     ********** */
    /* Parameter adjustments */
    --x;

    /* Function Body */
    ret_val = 0.;
    scale = 0.;
    i__1 = *n;
    i__2 = *incx;
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MAX */
	d__2 = scale, d__3 = (d__1 = x[i__], abs_(d__1));
	scale = max_(d__2,d__3);
/* L10: */
    }
    if (scale == 0.) {
	return ret_val;
    }
    i__2 = *n;
    i__1 = *incx;
    for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1) {
/* Computing 2nd power */
	d__1 = x[i__] / scale;
	ret_val += d__1 * d__1;
/* L20: */
    }
    ret_val = scale * sqrt(ret_val);
    return ret_val;
} /* dnrm2_ */

/* ====================== The end of dnrm2 =============================== */
double dpmeps_()
{
    /* Initialized data */

    static double zero = 0.;
    static double one = 1.;
    static double two = 2.;

    /* System generated locals */
    long int i__1;
    double ret_val;

    /* Local variables */
    static double beta;
    static long int irnd;
    static volatile double temp, temp1, a, b;
    static long int i__;
    static double betah;
    static long int ibeta, negep;
    static double tempa;
    static long int itemp, it;
    static double betain;

/*     ********** */

/*     Subroutine dpeps */

/*     This subroutine computes the machine precision parameter */
/*     dpmeps as the smallest floating point number such that */
/*     1 + dpmeps differs from 1. */

/*     This subroutine is based on the subroutine machar described in */

/*     W. J. Cody, */
/*     MACHAR: A subroutine to dynamically determine machine parameters, */
/*     ACM Transactions on Mathematical Software, 14, 1988, pages 303-311. */

/*     The subroutine statement is: */

/*       subroutine dpeps(dpmeps) */

/*     where */

/*       dpmeps is a double precision variable. */
/*         On entry dpmeps need not be specified. */
/*         On exit dpmeps is the machine precision. */

/*     MINPACK-2 Project. February 1991. */
/*     Argonne National Laboratory and University of Minnesota. */
/*     Brett M. Averick. */

/*     ******* */
/*     determine ibeta, beta ala malcolm. */
    a = one;
    b = one;
L10:
    a += a;
    temp = a + one;
    temp1 = temp - a;
    if (temp1 - one == zero) {
	goto L10;
    }
L20:
    b += b;
    temp = a + b;
    itemp = (long int) (temp - a);
    if (itemp == 0) {
	goto L20;
    }
    ibeta = itemp;
    beta = (double) ibeta;
/*     determine it, irnd. */
    it = 0;
    b = one;
L30:
    ++it;
    b *= beta;
    temp = b + one;
    temp1 = temp - b;
    if (temp1 - one == zero) {
	goto L30;
    }
    irnd = 0;
    betah = beta / two;
    temp = a + betah;
    if (temp - a != zero) {
	irnd = 1;
    }
    tempa = a + beta;
    temp = tempa + betah;
    if (irnd == 0 && temp - tempa != zero) {
	irnd = 2;
    }
/*     determine dpmeps. */
    negep = it + 3;
    betain = one / beta;
    a = one;
    i__1 = negep;
    for (i__ = 1; i__ <= i__1; ++i__) {
	a *= betain;
/* L40: */
    }
L50:
    temp = one + a;
    if (temp - one != zero) {
	goto L60;
    }
    a *= beta;
    goto L50;
L60:
    ret_val = a;
    if (ibeta == 2 || irnd == 0) {
	goto L70;
    }
    a = a * (one + a) / two;
    temp = one + a;
    if (temp - one != zero) {
	ret_val = a;
    }
L70:
    return ret_val;
} /* dpmeps_ */

/* ====================== The end of dpmeps ============================== */
/* Subroutine */ int daxpy_(long int *n, double *da, double *dx, long int *incx, double *dy, long int *incy)
/*long int *n;
double *da, *dx;
long int *incx;
double *dy;
long int *incy;*/
{
    /* System generated locals */
    long int i__1;

    /* Local variables */
    static long int i__, m, ix, iy, mp1;


/*     constant times a vector plus a vector. */
/*     uses unrolled loops for increments equal to one. */
/*     jack dongarra, linpack, 3/11/78. */


    /* Parameter adjustments */
    --dy;
    --dx;

    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
    if (*da == 0.) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*        code for unequal increments or equal increments */
/*          not equal to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dy[iy] += *da * dx[ix];
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;

/*        code for both increments equal to 1 */


/*        clean-up loop */

L20:
    m = *n % 4;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dy[i__] += *da * dx[i__];
/* L30: */
    }
    if (*n < 4) {
	return 0;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 4) {
	dy[i__] += *da * dx[i__];
	dy[i__ + 1] += *da * dx[i__ + 1];
	dy[i__ + 2] += *da * dx[i__ + 2];
	dy[i__ + 3] += *da * dx[i__ + 3];
/* L50: */
    }
    return 0;
} /* daxpy_ */

/* ====================== The end of daxpy =============================== */
/* Subroutine */ int dcopy_(long int *n, double *dx, long int *incx, double *dy, long int *incy)
/*long int *n;
double *dx;
long int *incx;
double *dy;
long int *incy;*/
{
    /* System generated locals */
    long int i__1;

    /* Local variables */
    static long int i__, m, ix, iy, mp1;


/*     copies a vector, x, to a vector, y. */
/*     uses unrolled loops for increments equal to one. */
/*     jack dongarra, linpack, 3/11/78. */


    /* Parameter adjustments */
    --dy;
    --dx;

    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*        code for unequal increments or equal increments */
/*          not equal to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dy[iy] = dx[ix];
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;

/*        code for both increments equal to 1 */


/*        clean-up loop */

L20:
    m = *n % 7;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dy[i__] = dx[i__];
/* L30: */
    }
    if (*n < 7) {
	return 0;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 7) {
	dy[i__] = dx[i__];
	dy[i__ + 1] = dx[i__ + 1];
	dy[i__ + 2] = dx[i__ + 2];
	dy[i__ + 3] = dx[i__ + 3];
	dy[i__ + 4] = dx[i__ + 4];
	dy[i__ + 5] = dx[i__ + 5];
	dy[i__ + 6] = dx[i__ + 6];
/* L50: */
    }
    return 0;
} /* dcopy_ */

/* ====================== The end of dcopy =============================== */
double ddot_(long int *n, double *dx, long int *incx, double *dy, long int *incy)
/*long int *n;
double *dx;
long int *incx;
double *dy;
long int *incy;*/
{
    /* System generated locals */
    long int i__1;
    double ret_val;

    /* Local variables */
    static long int i__, m;
    static double dtemp;
    static long int ix, iy, mp1;


/*     forms the dot product of two vectors. */
/*     uses unrolled loops for increments equal to one. */
/*     jack dongarra, linpack, 3/11/78. */


    /* Parameter adjustments */
    --dy;
    --dx;

    /* Function Body */
    ret_val = 0.;
    dtemp = 0.;
    if (*n <= 0) {
	return ret_val;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*        code for unequal increments or equal increments */
/*          not equal to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dtemp += dx[ix] * dy[iy];
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    ret_val = dtemp;
    return ret_val;

/*        code for both increments equal to 1 */


/*        clean-up loop */

L20:
    m = *n % 5;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dtemp += dx[i__] * dy[i__];
/* L30: */
    }
    if (*n < 5) {
	goto L60;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 5) {
	dtemp = dtemp + dx[i__] * dy[i__] + dx[i__ + 1] * dy[i__ + 1] + dx[
		i__ + 2] * dy[i__ + 2] + dx[i__ + 3] * dy[i__ + 3] + dx[i__ + 
		4] * dy[i__ + 4];
/* L50: */
    }
L60:
    ret_val = dtemp;
    return ret_val;
} /* ddot_ */

/* ====================== The end of ddot ================================ */
/* Subroutine */ int dpofa_(double *a, long int *lda, long int *n, long int *info)
/*double *a;
long int *lda, *n, *info;*/
{
    /* System generated locals */
    long int a_dim1, a_offset, i__1, i__2, i__3;

    /* Builtin functions */
    //double sqrt();

    /* Local variables */
    //extern double ddot_();
    static long int j, k;
    static double s, t;
    static long int jm1;


/*     dpofa factors a double precision symmetric positive definite */
/*     matrix. */

/*     dpofa is usually called by dpoco, but it can be called */
/*     directly with a saving in time if  rcond  is not needed. */
/*     (time for dpoco) = (1 + 18/n)*(time for dpofa) . */

/*     on entry */

/*        a       double precision(lda, n) */
/*                the symmetric matrix to be factored.  only the */
/*                diagonal and upper triangle are used. */

/*        lda     long int */
/*                the leading dimension of the array  a . */

/*        n       long int */
/*                the order of the matrix  a . */

/*     on return */

/*        a       an upper triangular matrix  r  so that  a = trans(r)*r */
/*                where  trans(r)  is the transpose. */
/*                the strict lower triangle is unaltered. */
/*                if  info .ne. 0 , the factorization is not complete. */

/*        info    long int */
/*                = 0  for normal return. */
/*                = k  signals an error condition.  the leading minor */
/*                     of order  k  is not positive definite. */

/*     linpack.  this version dated 08/14/78 . */
/*     cleve moler, university of new mexico, argonne national lab. */

/*     subroutines and functions */

/*     blas ddot */
/*     fortran sqrt */

/*     internal variables */

/*     begin block with ...exits to 40 */


    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;

    /* Function Body */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	*info = j;
	s = 0.;
	jm1 = j - 1;
	if (jm1 < 1) {
	    goto L20;
	}
	i__2 = jm1;
	for (k = 1; k <= i__2; ++k) {
	    i__3 = k - 1;
	    t = a[k + j * a_dim1] - ddot_(&i__3, &a[k * a_dim1 + 1], &c__1, &
		    a[j * a_dim1 + 1], &c__1);
	    t /= a[k + k * a_dim1];
	    a[k + j * a_dim1] = t;
	    s += t * t;
/* L10: */
	}
L20:
	s = a[j + j * a_dim1] - s;
/*     ......exit */
	if (s <= 0.) {
	    goto L40;
	}
	a[j + j * a_dim1] = sqrt(s);
/* L30: */
    }
    *info = 0;
L40:
    return 0;
} /* dpofa_ */

/* ====================== The end of dpofa =============================== */
/* Subroutine */ int dscal_(long int *n, double *da, double *dx, long int *incx)
/*long int *n;
double *da, *dx;
long int *incx;*/
{
    /* System generated locals */
    long int i__1, i__2;

    /* Local variables */
    static long int i__, m, nincx, mp1;


/*     scales a vector by a constant. */
/*     uses unrolled loops for increment equal to one. */
/*     jack dongarra, linpack, 3/11/78. */
/*     modified 3/93 to return if incx .le. 0. */


    /* Parameter adjustments */
    --dx;

    /* Function Body */
    if (*n <= 0 || *incx <= 0) {
	return 0;
    }
    if (*incx == 1) {
	goto L20;
    }

/*        code for increment not equal to 1 */

    nincx = *n * *incx;
    i__1 = nincx;
    i__2 = *incx;
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
	dx[i__] = *da * dx[i__];
/* L10: */
    }
    return 0;

/*        code for increment equal to 1 */


/*        clean-up loop */

L20:
    m = *n % 5;
    if (m == 0) {
	goto L40;
    }
    i__2 = m;
    for (i__ = 1; i__ <= i__2; ++i__) {
	dx[i__] = *da * dx[i__];
/* L30: */
    }
    if (*n < 5) {
	return 0;
    }
L40:
    mp1 = m + 1;
    i__2 = *n;
    for (i__ = mp1; i__ <= i__2; i__ += 5) {
	dx[i__] = *da * dx[i__];
	dx[i__ + 1] = *da * dx[i__ + 1];
	dx[i__ + 2] = *da * dx[i__ + 2];
	dx[i__ + 3] = *da * dx[i__ + 3];
	dx[i__ + 4] = *da * dx[i__ + 4];
/* L50: */
    }
    return 0;
} /* dscal_ */

/* ====================== The end of dscal =============================== */
/* Subroutine */ int dtrsl_(double *t, long int *ldt, long int *n, double *b, long int *job, long int *info)
/*double *t;
long int *ldt, *n;
double *b;
long int *job, *info;*/
{
    /* System generated locals */
    long int t_dim1, t_offset, i__1, i__2;

    /* Local variables */
    static long int case__;
    //extern double ddot_(void);
    static double temp;
    static long int j;
    //extern /* Subroutine */ int daxpy_(double *a, long int *n, double *c, long int *d, long int *e);
    extern int daxpy_(long int *n, double *da, double *dx, long int *incx, double *dy, long int *incy);
    static long int jj;



/*     dtrsl solves systems of the form */

/*                   t * x = b */
/*     or */
/*                   trans(t) * x = b */

/*     where t is a triangular matrix of order n. here trans(t) */
/*     denotes the transpose of the matrix t. */

/*     on entry */

/*         t         double precision(ldt,n) */
/*                   t contains the matrix of the system. the zero */
/*                   elements of the matrix are not referenced, and */
/*                   the corresponding elements of the array can be */
/*                   used to store other information. */

/*         ldt       long int */
/*                   ldt is the leading dimension of the array t. */

/*         n         long int */
/*                   n is the order of the system. */

/*         b         double precision(n). */
/*                   b contains the right hand side of the system. */

/*         job       long int */
/*                   job specifies what kind of system is to be solved. */
/*                   if job is */

/*                        00   solve t*x=b, t lower triangular, */
/*                        01   solve t*x=b, t upper triangular, */
/*                        10   solve trans(t)*x=b, t lower triangular, */
/*                        11   solve trans(t)*x=b, t upper triangular. */

/*     on return */

/*         b         b contains the solution, if info .eq. 0. */
/*                   otherwise b is unaltered. */

/*         info      long int */
/*                   info contains zero if the system is nonsingular. */
/*                   otherwise info contains the index of */
/*                   the first zero diagonal element of t. */

/*     linpack. this version dated 08/14/78 . */
/*     g. w. stewart, university of maryland, argonne national lab. */

/*     subroutines and functions */

/*     blas daxpy,ddot */
/*     fortran mod */

/*     internal variables */


/*     begin block permitting ...exits to 150 */

/*        check for zero diagonal elements. */

    /* Parameter adjustments */
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1 * 1;
    t -= t_offset;
    --b;

    /* Function Body */
    i__1 = *n;
    for (*info = 1; *info <= i__1; ++(*info)) {
/*     ......exit */
	if (t[*info + *info * t_dim1] == 0.) {
	    goto L150;
	}
/* L10: */
    }
    *info = 0;

/*        determine the task and go to it. */

    case__ = 1;
    if (*job % 10 != 0) {
	case__ = 2;
    }
    if (*job % 100 / 10 != 0) {
	case__ += 2;
    }
    switch ((int)case__) {
	case 1:  goto L20;
	case 2:  goto L50;
	case 3:  goto L80;
	case 4:  goto L110;
    }

/*        solve t*x=b for t lower triangular */

L20:
    b[1] /= t[t_dim1 + 1];
    if (*n < 2) {
	goto L40;
    }
    i__1 = *n;
    for (j = 2; j <= i__1; ++j) {
	temp = -b[j - 1];
	i__2 = *n - j + 1;
	daxpy_(&i__2, &temp, &t[j + (j - 1) * t_dim1], &c__1, &b[j], &c__1);
	b[j] /= t[j + j * t_dim1];
/* L30: */
    }
L40:
    goto L140;

/*        solve t*x=b for t upper triangular. */

L50:
    b[*n] /= t[*n + *n * t_dim1];
    if (*n < 2) {
	goto L70;
    }
    i__1 = *n;
    for (jj = 2; jj <= i__1; ++jj) {
	j = *n - jj + 1;
	temp = -b[j + 1];
	daxpy_(&j, &temp, &t[(j + 1) * t_dim1 + 1], &c__1, &b[1], &c__1);
	b[j] /= t[j + j * t_dim1];
/* L60: */
    }
L70:
    goto L140;

/*        solve trans(t)*x=b for t lower triangular. */

L80:
    b[*n] /= t[*n + *n * t_dim1];
    if (*n < 2) {
	goto L100;
    }
    i__1 = *n;
    for (jj = 2; jj <= i__1; ++jj) {
	j = *n - jj + 1;
	i__2 = jj - 1;
	b[j] -= ddot_(&i__2, &t[j + 1 + j * t_dim1], &c__1, &b[j + 1], &c__1);
	b[j] /= t[j + j * t_dim1];
/* L90: */
    }
L100:
    goto L140;

/*        solve trans(t)*x=b for t upper triangular. */

L110:
    b[1] /= t[t_dim1 + 1];
    if (*n < 2) {
	goto L130;
    }
    i__1 = *n;
    for (j = 2; j <= i__1; ++j) {
	i__2 = j - 1;
	b[j] -= ddot_(&i__2, &t[j * t_dim1 + 1], &c__1, &b[1], &c__1);
	b[j] /= t[j + j * t_dim1];
/* L120: */
    }
L130:
L140:
L150:
    return 0;
} /* dtrsl_ */


