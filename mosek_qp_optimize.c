/* Wrapper for Mosek QP solver */
/* 6 November 2007 */

#include <stdio.h>
#include <assert.h>
#include "mosek.h"

static void MSKAPI printstr(void *handle, char str[]) {
  //printf("%s", str);
} /* printstr */

int mosek_qp_optimize(double** psiDiffs, double* delta, double* w, double* cur_slack, long k, double C, double *primal_obj, long w_size, long w1_size) {
  long i,j,t,l;
  double *c;
  MSKlidxt *aptrb;
  MSKlidxt *aptre;
  MSKidxt *asub=NULL;
  double *aval=NULL;
  MSKboundkeye bkc[k];
  double blc[k];
  double buc[k];
  MSKboundkeye *bkx;
  double *blx;
  double *bux;
  MSKidxt *qsubi,*qsubj;
  double *qval;

  MSKenv_t env;
  MSKtask_t task;
  MSKrescodee r;

  c = (double*) malloc(sizeof(double)*(w_size+1));
  assert(c!=NULL);
  aptrb = (MSKlidxt*) malloc(sizeof(MSKlidxt)*(w_size+1));
  assert(aptrb!=NULL);
  aptre = (MSKlidxt*) malloc(sizeof(MSKlidxt)*(w_size+1));
  assert(aptre!=NULL);

  bkx = (MSKboundkeye*) malloc(sizeof(MSKboundkeye)*(w_size+1));
  assert(bkx!=NULL);
  blx = (double*) malloc(sizeof(double)*(w_size+1));
  assert(blx!=NULL);
  bux = (double*) malloc(sizeof(double)*(w_size+1));
  assert(bux!=NULL);
  qsubi = (MSKidxt*) malloc(sizeof(MSKidxt)*w_size);
  assert(qsubi!=NULL);  
  qsubj = (MSKidxt*) malloc(sizeof(MSKidxt)*w_size);
  assert(qsubj!=NULL);  
  qval = (double*) malloc(sizeof(double)*w_size);
  assert(qval!=NULL);  

  double *solutionVec = (double *) malloc(sizeof(double)*(w_size+1));

  l=0;
  for(i=0; i<1; i++){
    c[i] = C;
    aptrb[i] = l;
    for(j=0;j<k;j++){
      asub = (MSKidxt*) realloc(asub, sizeof(MSKidxt)*(l+1));
      assert(asub!=NULL);
      aval = (double*) realloc(aval, sizeof(double)*(l+1));
      assert(aval!=NULL);
      asub[l] = j;
      aval[l] = 1;
      l++;
    }
    aptre[i] = l;
  }  
  
  for (i=1;i<(w_size+1);i++) {
		c[i] = 0;
    aptrb[i] = l;
    for(j=0;j<k;j++){
      if(psiDiffs[j][i-1] != 0){
          asub = (MSKidxt*) realloc(asub, sizeof(MSKidxt)*(l+1));
          assert(asub!=NULL);
          aval = (double*) realloc(aval, sizeof(double)*(l+1));
          assert(aval!=NULL);
          asub[l] = j;
          aval[l] = psiDiffs[j][i-1];
          l++;
      }
    }
		aptre[i] = l;
  }

  for(i=0; i<1; i++){
      bkx[i] = MSK_BK_LO;
      blx[i] = 0;
      bux[i] = MSK_INFINITY;
  }
 
  for (i=1;i<(w1_size+1);i++) {
      bkx[i] = MSK_BK_FR;
      blx[i] = -MSK_INFINITY;
      bux[i] = MSK_INFINITY;
  }
  for (i=(w1_size+1);i<(w_size+1);i++) {
      bkx[i] = MSK_BK_UP;
      blx[i] = -MSK_INFINITY;
      bux[i] = 0;
  }
  for (i=0; i<k; i++) {
      bkc[i] = MSK_BK_LO;
      blc[i] = delta[i];
      buc[i] = MSK_INFINITY;
  }
  
  /* create mosek environment */
  r = MSK_makeenv(&env, NULL, NULL, NULL, NULL);

  /* check return code */
  if (r==MSK_RES_OK) {
    /* directs output to printstr function */
    MSK_linkfunctoenvstream(env, MSK_STREAM_LOG, NULL, printstr);
  }

  /* initialize the environment */
  r = MSK_initenv(env);

  if (r==MSK_RES_OK) {
    /* create the optimization task */
    r = MSK_maketask(env,k,(w_size+1),&task);
	
    if (r==MSK_RES_OK) {
      r = MSK_linkfunctotaskstream(task, MSK_STREAM_LOG,NULL,printstr);
	  
      if (r==MSK_RES_OK) {
      	r = MSK_inputdata(task,
      			  k,(w_size+1),
      			  k,(w_size+1),
      			  c,0.0,
      			  aptrb,aptre,
      			  asub,aval,
      			  bkc,blc,buc,
      			  bkx,blx,bux);						  
      }
	  
      if (r==MSK_RES_OK) {
      	/* coefficients for the Gram matrix */
      	t = 0;
      	for (i=1;i<(w_size+1);i++) {
      	    qsubi[t] = i;
      	    qsubj[t] = i;
      			qval[t] = 1;
      	    t++;
	      }
	    
	       r = MSK_putqobj(task, w_size, qsubi,qsubj,qval);
      }
      

      /* DEBUG */
      /*
      printf("t: %ld\n", t);
      for (i=0;i<t;i++) {
	       printf("qsubi: %d, qsubj: %d, qval: %.4f\n", qsubi[i], qsubj[i], qval[i]);
      }
      fflush(stdout);
      */
      /* DEBUG */

      /* set relative tolerance gap (DEFAULT = 1E-8)*/
      //MSK_putdouparam(task, MSK_DPAR_INTPNT_TOL_REL_GAP, 1E-10);
      MSK_putdouparam(task, MSK_DPAR_INTPNT_TOL_REL_GAP, 1E-14);

      if (r==MSK_RES_OK) {
	       r = MSK_optimize(task);
      }
      
      if (r==MSK_RES_OK) {
      	MSK_getsolutionslice(task,
      			     MSK_SOL_ITR,
      			     MSK_SOL_ITEM_XX,
      			     0,
      			     (w_size+1),
      			     solutionVec);
        /* print out alphas */
      	/*
      	for (i=0;i<k;i++) {
      	  printf("alpha[%ld]: %.8f\n", i, alpha[i]); fflush(stdout);
      	}
      	*/
      	/* output the objective value */
      	MSK_getprimalobj(task, MSK_SOL_ITR, primal_obj);
      	//printf("ITER DUAL_OBJ %.8g\n", -(*dual_obj)); fflush(stdout);
      }
      MSK_deletetask(&task);
    }
    MSK_deleteenv(&env);
  }

  for(i=0; i<1; i++){
      cur_slack[i] = solutionVec[i];
  }

  for (i = 1; i < (w_size+1); i++)
  {
    w[i-k+1] = solutionVec[i];
  }
  
  
  /* free the memory */
  free(c);
  free(aptrb);
  free(aptre);
  free(asub);
  free(aval);
  free(bkx);
  free(blx);
  free(bux);
  free(qsubi);  
  free(qsubj);  
  free(qval);  

  free(solutionVec);
  
	if(r == MSK_RES_OK)
  	return(0);  
	else
		return(r);
}

