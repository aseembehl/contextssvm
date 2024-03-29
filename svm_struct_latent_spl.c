/************************************************************************/
/*                                                                      */
/*   svm_struct_latent_spl.c                                            */
/*                                                                      */
/*   Main Optimization Code for Latent SVM^struct using Self-Paced      */
/*   Learning. NOTE: This implementation modifies the CCCP code by      */
/*   Chun-Nam Yu, specifically the file svm_struct_latent_cccp.c,       */
/*   which is a part of the Latent SVM^struct package available on      */
/*   Chun-Nam Yu's webpage.                                             */
/*                                                                      */
/*   Authors: M. Pawan Kumar and Ben Packer                             */
/*                                                                      */
/*   This software is available for non-commercial use only. It must    */
/*   not be modified and distributed without prior permission of the    */
/*   author. The author is not responsible for implications from the    */
/*   use of this software.                                              */
/*                                                                      */
/************************************************************************/

#include <stdio.h>
#include <assert.h>
#include <float.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "svm_struct_latent_api.h"

	
	#include "./svm_light/svm_learn.h"



#define ALPHA_THRESHOLD 1E-14
#define IDLE_ITER 20
#define CLEANUP_CHECK 50
#define STOP_PREC 1E-2
#define UPDATE_BOUND 3
#define MAX_CURRICULUM_ITER 10

#define EQUALITY_EPSILON 1e-6
#define SLACK_EPSILON 1e-8

#define MAX(x,y) ((x) < (y) ? (y) : (x))
#define MIN(x,y) ((x) > (y) ? (y) : (x))

#define DEBUG_LEVEL 0

 //int mosek_qp_optimize(double**, double*, double*, long, double, double*, DOC **, int, int);
 //int mosek_qp_optimize(double**,double**, double*, double*, long, long, double, double*, double, double);
int mosek_qp_optimize(double**, double*, double*, double*, long, double, double*, long, long);
int mosek_qp_optimize_dual(double**, double**, double*, double*, long, long, double, double*, double, double);

void my_read_input_parameters(int argc, char* argv[], char *trainfile, char *modelfile, 
			      LEARN_PARM *learn_parm, KERNEL_PARM *kernel_parm, STRUCT_LEARN_PARM *struct_parm, 
						double *init_spl_weight, double *spl_factor);

void my_wait_any_key();

int resize_cleanup(int size_active, int **ptr_idle, double **ptr_cur_slack, double **ptr_delta, DOC ***ptr_dXc,
		double ***ptr_psiDiffs, int *mv_iter);

void approximate_to_psd(double **G, int size_active, double eps);

void Jacobi_Cyclic_Method(double eigenvalues[], double *eigenvectors, double *A, int n);

double sprod_nn(double *a, double *b, long n) {
  double ans=0.0;
  long i;
  for (i=1;i<n+1;i++) {
    ans+=a[i]*b[i];
  }
  return(ans);
}

void add_vector_nn(double *w, double *dense_x, long n, double factor) {
  long i;
  for (i=1;i<n+1;i++) {
    w[i]+=factor*dense_x[i];
  }
}

double* add_list_nn(SVECTOR *a, long totwords) 
     /* computes the linear combination of the SVECTOR list weighted
	by the factor of each SVECTOR. assumes that the number of
	features is small compared to the number of elements in the
	list */
{
    SVECTOR *f;
    long i;
    double *sum;

    sum=create_nvector(totwords);

    for(i=0;i<=totwords;i++) 
      sum[i]=0;

    for(f=a;f;f=f->next)  
      add_vector_ns(sum,f,f->factor);

    return(sum);
}


double current_obj_val(EXAMPLE *ex, SVECTOR **fycache, long m, STRUCTMODEL *sm, STRUCT_LEARN_PARM *sparm, double C, int *valid_examples) {

  long i, j;
  SVECTOR *f, *fy, *fybar, *lhs;
  LABEL       ybar;
  double lossval, margin;
  double *new_constraint;
	double obj = 0.0;

  /* find cutting plane */
  lhs = NULL;
  margin = 0;
  for (i=0;i<m;i++) {
		if(!valid_examples[i])
			continue;
    find_most_violated_constraint_marginrescaling(ex[i].x, ex[i].y, &ybar, sm, sparm);
    /* get difference vector */
    fy = copy_svector(fycache[i]);
    fybar = psi(ex[i].x,ybar,sm,sparm);
    lossval = loss(ex[i].y,ybar,sparm);

    /* scale difference vector */
    for (f=fy;f;f=f->next) {
      //f->factor*=1.0/m;
      f->factor*=ex[i].x.example_cost/m;
    }
    for (f=fybar;f;f=f->next) {
      //f->factor*=-1.0/m;
      f->factor*=-ex[i].x.example_cost/m;
    }
    /* add ybar to constraint */
    append_svector_list(fy,lhs);
    append_svector_list(fybar,fy);
    lhs = fybar;
    //margin+=lossval/m;
		margin += lossval*ex[i].x.example_cost/m;
  }

  /* compact the linear representation */
  new_constraint = add_list_nn(lhs, sm->sizePsi);
  free_svector(lhs);

	obj = margin;
	for(i = 1; i < sm->sizePsi+1; i++)
		obj -= new_constraint[i]*sm->w[i];
	if(obj < 0.0)
		obj = 0.0;
	obj *= C;
	for(i = 1; i < sm->sizePsi+1; i++)
		obj += 0.5*sm->w[i]*sm->w[i];
  free(new_constraint);

	return obj;
}

int compar(const void *a, const void *b)
{
  sortStruct *c = (sortStruct *) a;
  sortStruct *d = (sortStruct *) b;
  if(c->val < d->val)
    return -1;
  if(c->val > d->val)
    return 1;
  return 0;
}


SVECTOR* find_cutting_plane(EXAMPLE *ex, SVECTOR **fycache, double *margin, long m, STRUCTMODEL *sm, STRUCT_LEARN_PARM *sparm,
														int *valid_examples) {

  long i, j;
  SVECTOR *f, *fy, *fybar, *lhs;
  LABEL       ybar;
  double lossval;
  double *new_constraint;
	long valid_count = 0;

  long l,k;
  SVECTOR *fvec;
  WORD *words;  

  /* find cutting plane */
  lhs = NULL;
  *margin = 0;

	for (i=0;i<m;i++) {
		if (valid_examples[i]) {
			valid_count++;
		}
	}

  for (i=0;i<m;i++) {

		if (!valid_examples[i]) {
			continue;
		}

    find_most_violated_constraint_marginrescaling(ex[i].x, ex[i].y, &ybar, sm, sparm);
    /* get difference vector */
    fy = copy_svector(fycache[i]);
    fybar = psi(ex[i].x,ybar,sm,sparm);
    lossval = loss(ex[i].y,ybar,sparm);
    free_label(ybar);
		
    /* scale difference vector */
    for (f=fy;f;f=f->next) {
      //f->factor*=1.0/m;
      //f->factor*=ex[i].x.example_cost/m;
      f->factor*=ex[i].x.example_cost/valid_count;
    }
    for (f=fybar;f;f=f->next) {
      //f->factor*=-1.0/m;
      //f->factor*=-ex[i].x.example_cost/m;
      f->factor*=-ex[i].x.example_cost/valid_count;
    }
    /* add ybar to constraint */
    append_svector_list(fy,lhs);
    append_svector_list(fybar,fy);
    lhs = fybar;
    //*margin+=lossval/m;
    //*margin+=lossval*ex[i].x.example_cost/m;
    *margin+=lossval*ex[i].x.example_cost/valid_count;
  }

  /* compact the linear representation */
  new_constraint = add_list_nn(lhs, sm->sizePsi);
  free_svector(lhs);

  l=0;
  for (i=1;i<sm->sizePsi+1;i++) {
    if (fabs(new_constraint[i])>1E-10) l++; // non-zero
  }
  words = (WORD*)my_malloc(sizeof(WORD)*(l+1)); 
  assert(words!=NULL);
  k=0;
  for (i=1;i<sm->sizePsi+1;i++) {
    if (fabs(new_constraint[i])>1E-10) {
      words[k].wnum = i;
      words[k].weight = new_constraint[i]; 
      k++;
    }
  }
  words[k].wnum = 0;
  words[k].weight = 0.0;
  fvec = create_svector(words,"",1);

  free(words);
  free(new_constraint);

  return(fvec); 
}

double cutting_plane_algorithm(double *w, long m, int MAX_ITER, double C, double epsilon, SVECTOR **fycache, EXAMPLE *ex, 
															STRUCTMODEL *sm, STRUCT_LEARN_PARM *sparm, int *valid_examples) {
  long i,j,t;
  double *alpha;
  DOC **dXc; /* constraint matrix */
  double *delta; /* rhs of constraints */
  SVECTOR *new_constraint;
  int iter, size_active; 
  double value;
	double threshold = 0.0;
  double margin;
  double primal_obj, cur_obj;
	double *cur_slack = NULL;
	int mv_iter;
	int *idle = NULL;
	double **psiDiffs = NULL;
	SVECTOR *f;
	int r;
	long fnum, last_wnum;

  /* set parameters for hideo solver */
  LEARN_PARM lparm;
  KERNEL_PARM kparm;
  MODEL *svm_model=NULL;
  lparm.biased_hyperplane = 0;
  lparm.epsilon_crit = MIN(epsilon,0.001);
  lparm.svm_c = C;
  lparm.sharedslack = 1;
  kparm.kernel_type = LINEAR;

  lparm.remove_inconsistent=0;
  lparm.skip_final_opt_check=0;
  lparm.svm_maxqpsize=10;
  lparm.svm_newvarsinqp=0;
  lparm.svm_iter_to_shrink=-9999;
  lparm.maxiter=100000;
  lparm.kernel_cache_size=40;
  lparm.eps = epsilon; 
  lparm.transduction_posratio=-1.0;
  lparm.svm_costratio=1.0;
  lparm.svm_costratio_unlab=1.0;
  lparm.svm_unlabbound=1E-5;
  lparm.epsilon_a=1E-10;  /* changed from 1e-15 */
  lparm.compute_loo=0;
  lparm.rho=1.0;
  lparm.xa_depth=0;
  strcpy(lparm.alphafile,"");
  kparm.poly_degree=3;
  kparm.rbf_gamma=1.0;
  kparm.coef_lin=1;
  kparm.coef_const=1;
  strcpy(kparm.custom,"empty");
 
  iter = 0;
  size_active = 0;
  alpha = NULL;
  dXc = NULL;
  delta = NULL;

  printf("Running structural SVM solver: "); fflush(stdout); 

	new_constraint = find_cutting_plane(ex, fycache, &margin, m, sm, sparm, valid_examples);
 	value = margin - sprod_ns(w, new_constraint);
	while((value>threshold+epsilon)&&(iter<MAX_ITER)) {
		iter+=1;
		size_active+=1;

		printf("."); fflush(stdout); 


	    /* add  constraint */
	  	dXc = (DOC**)realloc(dXc, sizeof(DOC*)*size_active);
	   	assert(dXc!=NULL);
	   	dXc[size_active-1] = (DOC*)malloc(sizeof(DOC));
	   	dXc[size_active-1]->fvec = new_constraint; 
	   	dXc[size_active-1]->slackid = 1; // only one common slackid (one-slack)
	   	dXc[size_active-1]->costfactor = 1.0;

	   	delta = (double*)realloc(delta, sizeof(double)*size_active);
	   	assert(delta!=NULL);
	   	delta[size_active-1] = margin;

	   	/*alpha = (double*)realloc(alpha, sizeof(double)*size_active);
	   	assert(alpha!=NULL);
	   	alpha[size_active-1] = 0.0;*/

		/*idle = (int *) realloc(idle, sizeof(int)*size_active);
		assert(idle!=NULL);
		idle[size_active-1] = 0;*/

		/* update Gram matrix */
		psiDiffs = (double **) realloc(psiDiffs, sizeof(double *)*size_active);
		assert(psiDiffs!=NULL);
		psiDiffs[size_active-1] = NULL;
		psiDiffs[size_active-1] = (double *) realloc(psiDiffs[size_active-1], sizeof(double)*((sparm->phi1_size+sparm->phi2_size)*3));
		assert(psiDiffs[size_active-1]!=NULL);
		
		fnum = 0;
		last_wnum = 0;
		while(dXc[size_active-1]->fvec->words[fnum].wnum) {
			for (t = last_wnum+1; t < dXc[size_active-1]->fvec->words[fnum].wnum; t++)	{
				psiDiffs[size_active-1][t-1] = 0;
			}
			psiDiffs[size_active-1][dXc[size_active-1]->fvec->words[fnum].wnum-1] = dXc[size_active-1]->fvec->words[fnum].weight;
			/*if((psiDiffs[size_active-1][dXc[size_active-1]->fvec->words[fnum].wnum-1]<EQUALITY_EPSILON) && (psiDiffs[size_active-1][dXc[size_active-1]->fvec->words[fnum].wnum-1]>(-1*EQUALITY_EPSILON))){
				psiDiffs[size_active-1][dXc[size_active-1]->fvec->words[fnum].wnum-1] = 0;
			}*/
			last_wnum = dXc[size_active-1]->fvec->words[fnum].wnum;
			fnum++;
		}
		for (t = (last_wnum+1); t <= (sparm->phi1_size+sparm->phi2_size)*3; t++)	{
			psiDiffs[size_active-1][t-1] = 0;
		}			

   		/* solve QP to update w */
   		clear_nvector(w,sm->sizePsi);
   		//cur_slack = (double *) realloc(cur_slack,sizeof(double)*size_active);
   		cur_slack = (double *) realloc(cur_slack,sizeof(double));

		r = mosek_qp_optimize(psiDiffs, delta, w, cur_slack, (long) size_active, C, &cur_obj, (sparm->phi1_size+sparm->phi2_size)*3, (sparm->phi1_size+sparm->phi2_size)*2);

		if(r >= 1293 && r <= 1296)
		{
			printf("r:%d. G might not be psd due to numerical errors.\n",r);
			exit(1);
		}
		else if(r)
		{
			printf("Error %d in mosek_qp_optimize: Check ${MOSEKHOME}/${VERSION}/tools/platform/${PLATFORM}/h/mosek.h\n",r);
			exit(1);
		}

		for(j = 1; j <= (sparm->phi1_size+sparm->phi2_size)*3; j++) {
			if((w[j]<EQUALITY_EPSILON) && (w[j]>(-1*EQUALITY_EPSILON))){
	   			w[j] = 0;
   			}
		}

		/*for (j=0;j<size_active;j++) {
	     	if (cur_slack[j]>ALPHA_THRESHOLD) {
					idle[j] = 0;
	     	}
				else
					idle[j]++;
   		}*/

		/*mv_iter = 0;
		if(size_active > 1) {
			for(j = 0; j < size_active; j++) {
				if(cur_slack[j] >= cur_slack[mv_iter])
					mv_iter = j;
			}
		}*/

		if(size_active > 1)
			//threshold = cur_slack[mv_iter];
			threshold = cur_slack[0];
		else
			threshold = 0.0;

 		new_constraint = find_cutting_plane(ex, fycache, &margin, m, sm, sparm, valid_examples);
   		value = margin - sprod_ns(w, new_constraint);

		/*if((iter % CLEANUP_CHECK) == 0)
		{
			printf("+"); fflush(stdout);
			size_active = resize_cleanup(size_active, &idle, &cur_slack, &delta, &dXc, &psiDiffs, &mv_iter);
		}*/

 	} // end cutting plane while loop 

	primal_obj = current_obj_val(ex, fycache, m, sm, sparm, C, valid_examples);

  printf(" Inner loop optimization finished.\n"); fflush(stdout); 
      
  /* free memory */
  for (j=0;j<size_active;j++) {
		free(psiDiffs[j]);
    free_example(dXc[j],1);	
  }
	free(psiDiffs);
  free(dXc);
  //free(alpha);
  free(delta);
  free_svector(new_constraint);
	free(cur_slack);
	//free(idle);
  if (svm_model!=NULL) free_model(svm_model,0);

  return(primal_obj);
}

void cutting_plane_algorithm_dual(double *w, long m, int MAX_ITER, double C, double epsilon, SVECTOR **fycache, EXAMPLE *ex, 
															STRUCTMODEL *sm, STRUCT_LEARN_PARM *sparm, int *valid_examples) {
  long i,j;
  double *alpha;
  DOC **dXc; // constraint matrix 
  double *delta; // rhs of constraints 
  SVECTOR *new_constraint;
  int iter, size_active; 
  double value;
	double threshold = 0.0;
  double margin;
  double primal_obj, cur_obj;
	double *cur_slack = NULL;
	int mv_iter;
	int *idle = NULL;
	double **G = NULL;
	double **G2 = NULL;
	double **qmatrix = NULL;
	SVECTOR *f;
	int r;

  // set parameters for hideo solver 
  LEARN_PARM lparm;
  KERNEL_PARM kparm;
  MODEL *svm_model=NULL;
  lparm.biased_hyperplane = 0;
  lparm.epsilon_crit = MIN(epsilon,0.001);
  lparm.svm_c = C;
  lparm.sharedslack = 1;
  kparm.kernel_type = LINEAR;

  lparm.remove_inconsistent=0;
  lparm.skip_final_opt_check=0;
  lparm.svm_maxqpsize=10;
  lparm.svm_newvarsinqp=0;
  lparm.svm_iter_to_shrink=-9999;
  lparm.maxiter=100000;
  lparm.kernel_cache_size=40;
  lparm.eps = epsilon; 
  lparm.transduction_posratio=-1.0;
  lparm.svm_costratio=1.0;
  lparm.svm_costratio_unlab=1.0;
  lparm.svm_unlabbound=1E-5;
  lparm.epsilon_a=1E-10;  // changed from 1e-15 
  lparm.compute_loo=0;
  lparm.rho=1.0;
  lparm.xa_depth=0;
  strcpy(lparm.alphafile,"");
  kparm.poly_degree=3;
  kparm.rbf_gamma=1.0;
  kparm.coef_lin=1;
  kparm.coef_const=1;
  strcpy(kparm.custom,"empty");
 
  iter = 0;
  size_active = 0;
  alpha = NULL;
  dXc = NULL;
  delta = NULL;

  //qmatrix = (double **) malloc(sizeof(double *)*10);
  //assert(qmatrix!=NULL);

  printf("Running structural SVM solver: "); fflush(stdout); 
	new_constraint = find_cutting_plane(ex, fycache, &margin, m, sm, sparm, valid_examples);
 	value = margin - sprod_ns(w, new_constraint);
	while((value>threshold+epsilon)&&(iter<MAX_ITER)) {
		iter+=1;
		size_active+=1;

		printf("."); fflush(stdout); 


	    // add  constraint 
	  	dXc = (DOC**)realloc(dXc, sizeof(DOC*)*size_active);
	   	assert(dXc!=NULL);
	   	dXc[size_active-1] = (DOC*)malloc(sizeof(DOC));
	   	dXc[size_active-1]->fvec = new_constraint; 
	   	dXc[size_active-1]->slackid = 1; // only one common slackid (one-slack)
	   	dXc[size_active-1]->costfactor = 1.0;


	   	delta = (double*)realloc(delta, sizeof(double)*size_active);
	   	assert(delta!=NULL);
	   	delta[size_active-1] = margin;

	   	//alpha = (double*)malloc(sizeof(double)*(size_active+(sparm->phi1_size+sparm->phi2_size)));
	   	//assert(alpha!=NULL);
   		//for(j=0; j<(sparm->phi1_size+sparm->phi2_size)+size_active; j++){
   		//	alpha[j] = 0.0;
   		//}
   		alpha = (double*)realloc(alpha, sizeof(double)*(size_active+(sparm->phi1_size+sparm->phi2_size)));
	   	assert(alpha!=NULL);
	   	alpha[size_active-1] = 0.0;

		idle = (int *) realloc(idle, sizeof(int)*size_active);
		assert(idle!=NULL);
		idle[size_active-1] = 0;

		
		qmatrix = (double **) realloc(qmatrix, sizeof(double *)*size_active);
  		assert(qmatrix!=NULL);

		qmatrix[size_active-1] = malloc(sizeof(double)*(sparm->phi1_size+sparm->phi2_size));
		for(j = 0; j < (sparm->phi1_size+sparm->phi2_size); j++){
			qmatrix[size_active-1][j] = (-1)*returnWeightAtIndex(dXc[size_active-1]->fvec->words, ((sparm->phi1_size+sparm->phi2_size)*2+j+1));
		}

		// update Gram matrix 
		G = (double **) realloc(G, sizeof(double *)*size_active);
		assert(G!=NULL);
		G[size_active-1] = NULL;
		for(j = 0; j < size_active; j++) {
			G[j] = (double *) realloc(G[j], sizeof(double)*size_active);
			assert(G[j]!=NULL);
		}

		for(j = 0; j < size_active-1; j++) {
			G[size_active-1][j] = sprod_ss(dXc[size_active-1]->fvec, dXc[j]->fvec);
			G[size_active-1][j] = G[size_active-1][j]/2;
			G[j][size_active-1]  = G[size_active-1][j];
		}
		G[size_active-1][size_active-1] = sprod_ss(dXc[size_active-1]->fvec,dXc[size_active-1]->fvec);

		// hack: add a constant to the diagonal to make sure G is PSD 
		G[size_active-1][size_active-1] += 1e-6;

	   	// solve QP to update alpha 
		//r = mosek_qp_optimize(G, delta, alpha, (long) size_active, C, &cur_obj, dXc, (sparm->phi1_size+sparm->phi2_size)*2, (sparm->phi1_size+sparm->phi2_size));
		r = mosek_qp_optimize_dual(G, qmatrix, delta, alpha, (long) size_active, (long) (sparm->phi1_size+sparm->phi2_size), C, &cur_obj, 0, 0);
	    
		if(r >= 1293 && r <= 1296)
		{
			printf("r:%d. G might not be psd due to numerical errors.\n",r);
			fflush(stdout);
			//exit(1);
			while(r==1295) {
				printf("r:%d. G might not be psd due to numerical errors. Gram Reg=%0.7f\n",r, sparm->gram_regularization);
				fflush(stdout);
				for(i=0;i<size_active;i++) {
					G[i][i] += 10*sparm->gram_regularization-sparm->gram_regularization;
				}
				sparm->gram_regularization *= 10;
				r = mosek_qp_optimize_dual(G, qmatrix, delta, alpha, (long) size_active, (long) (sparm->phi1_size+sparm->phi2_size), C, &cur_obj, sparm->gram_regularization, sparm->gram_regularization*0.1);
			}
		}
		else if(r)
		{
			printf("Error %d in mosek_qp_optimize: Check ${MOSEKHOME}/${VERSION}/tools/platform/${PLATFORM}/h/mosek.h\n",r);
			exit(1);
		}

	   	clear_nvector(w,sm->sizePsi);
	   	for (j=0;j<size_active;j++) {
	     	if (alpha[j]>C*ALPHA_THRESHOLD) {
					add_vector_ns(w,dXc[j]->fvec,alpha[j]);
					idle[j] = 0;
	     	}
			else
				idle[j]++;
	   	}
	   	for(j=0; j<(sparm->phi1_size+sparm->phi2_size);j++){
	   		if (alpha[size_active+j] > EQUALITY_EPSILON){
	   			w[j+1+(sparm->phi1_size+sparm->phi2_size)*2] = w[j+1+(sparm->phi1_size+sparm->phi2_size)*2] - alpha[size_active+j];
	   		}	   		
	   	}

	   	for(j=1; j<=(sparm->phi1_size+sparm->phi2_size)*3;j++){
	   		if((w[j]<EQUALITY_EPSILON) && (w[j]>(-1*EQUALITY_EPSILON))){
	   			w[j] = 0;
	   		}
	   	}	   

	   	for(j=(sparm->phi1_size+sparm->phi2_size)*2+1; j<=(sparm->phi1_size+sparm->phi2_size)*3;j++){
	   		//assert(w[j] <= 0);
	   		if(w[j]>0){
	   			printf("j = %ld, w[j] = %0.6f\n", j, w[j]);
	   			fflush(stdout);
	   		}
	   		
	   	}	

		cur_slack = (double *) realloc(cur_slack,sizeof(double)*size_active);

		for(i = 0; i < size_active; i++) {
			cur_slack[i] = 0.0;
			for(f = dXc[i]->fvec; f; f = f->next) {
				j = 0;
				while(f->words[j].wnum) {
					cur_slack[i] += w[f->words[j].wnum]*f->words[j].weight;
					j++;
				}
			}
			if(cur_slack[i] >= delta[i])
				cur_slack[i] = 0.0;
			else
				cur_slack[i] = delta[i]-cur_slack[i];
		}

		mv_iter = 0;
		if(size_active > 1) {
			for(j = 0; j < size_active; j++) {
				if(cur_slack[j] >= cur_slack[mv_iter])
					mv_iter = j;
			}
		}

		if(size_active > 1)
			threshold = cur_slack[mv_iter];
		else
			threshold = 0.0;

 		new_constraint = find_cutting_plane(ex, fycache, &margin, m, sm, sparm, valid_examples);
   		value = margin - sprod_ns(w, new_constraint);

		if((iter % CLEANUP_CHECK) == 0)
		{
			printf("+"); fflush(stdout);
			size_active = resize_cleanup(size_active, &idle, &alpha, &delta, &dXc, &G, &mv_iter);
		}

		free(alpha);
		alpha=NULL;

 	} // end cutting plane while loop 

	//primal_obj = current_obj_val(ex, fycache, m, sm, sparm, C, valid_examples);

  printf(" Inner loop optimization finished.\n"); fflush(stdout); 
      
  // free memory
  for (j=0;j<size_active;j++) {
		free(G[j]);
    free_example(dXc[j],1);	
  }
	free(G);
  free(dXc);
  free(alpha);
  free(delta);
  free_svector(new_constraint);
	free(cur_slack);
	free(idle);
  if (svm_model!=NULL) free_model(svm_model,0);

  //return(primal_obj);
  return;
}

int check_acs_convergence(int *prev_valid_examples, int *valid_examples, long m)
{
	long i;
	int converged = 1;

	for (i=0;i<m;i++) {
		if (prev_valid_examples[i] != valid_examples[i]) {
			converged = 0;
			break;
		}
	}

	return converged;
}

int update_valid_examples(double *w, long m, double C, SVECTOR **fycache, EXAMPLE *ex, 
													STRUCTMODEL *sm, STRUCT_LEARN_PARM *sparm, int *valid_examples, double spl_weight) {

	long i, j;

	/* if self-paced learning weight is non-positive, all examples are valid */
	if(spl_weight <= 0.0) {
		for (i=0;i<m;i++)
			valid_examples[i] = 1;
		return (m);
	}

	sortStruct *slack = (sortStruct *) malloc(m*sizeof(sortStruct));
	LABEL ybar;
	SVECTOR *f, *fy, *fybar;
	double lossval;
	double penalty = 1.0/spl_weight;
	if(penalty < 0.0)
		penalty = DBL_MAX;

	for (i=0;i<m;i++) {
		find_most_violated_constraint_marginrescaling(ex[i].x, ex[i].y, &ybar, sm, sparm);
		fy = copy_svector(fycache[i]);
		fybar = psi(ex[i].x,ybar,sm,sparm);
		slack[i].index = i;
		slack[i].val = loss(ex[i].y,ybar,sparm);
		for (f=fy;f;f=f->next) {
			j = 0;
			while (1) {
				if(!f->words[j].wnum)
					break;
				slack[i].val -= sm->w[f->words[j].wnum]*f->words[j].weight;
				j++;
			}
		}
		for (f=fybar;f;f=f->next) {
			j = 0;
			while (1) {
				if(!f->words[j].wnum)
					break;
				slack[i].val += sm->w[f->words[j].wnum]*f->words[j].weight;
				j++;
			}
		}
		free_svector(fy);
		free_svector(fybar);
	}
	qsort(slack,m,sizeof(sortStruct),&compar);

	int nValid = 0;
	for (i=0;i<m;i++)
		valid_examples[i] = 0;
	for (i=0;i<m;i++) {
		if(slack[i].val*C/m > penalty)
			break;
		valid_examples[slack[i].index] = 1;
		nValid++;
	}

	free(slack);

	return nValid;
}

double alternate_convex_search(double *w, long m, int MAX_ITER, double C, double epsilon, SVECTOR **fycache, EXAMPLE *ex, 
                               STRUCTMODEL *sm, STRUCT_LEARN_PARM *sparm, int *valid_examples, double spl_weight) {

	long i;
	int iter = 0, converged, nValid;
	double last_relaxed_primal_obj = DBL_MAX, relaxed_primal_obj, decrement;

	int *prev_valid_examples = (int *) malloc(m*sizeof(int));
	double *best_w = (double *) malloc((sm->sizePsi+1)*sizeof(double));

	for (i=0;i<sm->sizePsi+1;i++)
		best_w[i] = w[i];
	nValid = update_valid_examples(w, m, C, fycache, ex, sm, sparm, valid_examples, spl_weight);
	//last_relaxed_primal_obj = current_obj_val(ex, fycache, m, sm, sparm, C, valid_examples);
	if(nValid < m)
		last_relaxed_primal_obj += (double)(m-nValid)/((double)spl_weight);

	for (i=0;i<m;i++) {
		prev_valid_examples[i] = 0;
	}

	for (iter=0;;iter++) {
		nValid = update_valid_examples(w, m, C, fycache, ex, sm, sparm, valid_examples, spl_weight);
		printf("ACS Iteration %d: number of examples = %d\n",iter,nValid); fflush(stdout);
		converged = check_acs_convergence(prev_valid_examples,valid_examples,m);
		if(converged)
			break;
		for (i=0;i<sm->sizePsi+1;i++)
			w[i] = 0.0;
		
    if(sparm->solve_dual){
        cutting_plane_algorithm_dual(w, m, MAX_ITER, C, epsilon, fycache, ex, sm, sparm, valid_examples);
    }
    else{
        cutting_plane_algorithm(w, m, MAX_ITER, C, epsilon, fycache, ex, sm, sparm, valid_examples);
    }
		/*relaxed_primal_obj = cutting_plane_algorithm(w, m, MAX_ITER, C, epsilon, fycache, ex, sm, sparm, valid_examples);
		if(nValid < m)
			relaxed_primal_obj += (double)(m-nValid)/((double)spl_weight);
		decrement = last_relaxed_primal_obj-relaxed_primal_obj;
    	printf("relaxed primal objective: %.4f\n", relaxed_primal_obj);
		if (iter) {
    	printf("decrement: %.4f\n", decrement); fflush(stdout);
		}
		else {
			printf("decrement: N/A\n"); fflush(stdout);
		}
		if (decrement>=0.0) {
			for (i=0;i<sm->sizePsi+1;i++) {
				best_w[i] = w[i];
			}
		}
		if (decrement <= C*epsilon) {
			break;
		}
		last_relaxed_primal_obj = relaxed_primal_obj;*/
		for (i=0;i<sm->sizePsi+1;i++) {
				best_w[i] = w[i];
		}
		for (i=0;i<m;i++) {
			prev_valid_examples[i] = valid_examples[i];
		}
	}

	for (i=0;i<m;i++) {
		prev_valid_examples[i] = 1;
	}

	if (iter) {
		for (i=0;i<sm->sizePsi+1;i++) {
			w[i] = best_w[i];
		}
	}

	double primal_obj;
	primal_obj = current_obj_val(ex, fycache, m, sm, sparm, C, prev_valid_examples);
	
	free(prev_valid_examples);
	free(best_w);

	//return;
	//return(relaxed_primal_obj);
	return(primal_obj);
}


long *randperm(long m)
{
  long *perm = (long *) malloc(sizeof(long)*m);
  long *map = (long *) malloc(sizeof(long)*m);
  long i, j;
  for(i = 0; i < m; i++)
    map[i] = i;
  srand(time(NULL));
  for(i = 0; i < m; i++)
  {
    int r = (int) (((double) m-i)*((double) rand())/(RAND_MAX+1.0));
    perm[i] = map[r];
    for(j = r; j < m-1; j++)
      map[j] = map[j+1];
  }
  free(map);
  return perm;
}

SAMPLE  generate_train_set(SAMPLE alldata, long *perm, int ntrain)
{
  SAMPLE  train;
  train.n = ntrain;
  long i;

  train.examples = (EXAMPLE *) malloc(train.n*sizeof(EXAMPLE));

  for(i = 0; i < train.n; i++)
  {
    train.examples[i] = alldata.examples[perm[i]];
  }

  return train;
}

SAMPLE  generate_validation_set(SAMPLE alldata, long *perm, int ntrain)
{
  SAMPLE  val;
  val.n = alldata.n - ntrain;
  long i;

  val.examples = (EXAMPLE *) malloc(val.n*sizeof(EXAMPLE));

  for(i = 0; i < val.n; i++)
  {
    val.examples[i] = alldata.examples[perm[ntrain+i]];
  }

  return val;
}

double compute_current_loss(SAMPLE val, STRUCTMODEL *sm, STRUCT_LEARN_PARM *sparm)
{
	long i;
	LABEL y;
	double cur_loss = 0.0;
	double store;
	for(i = 0; i < val.n; i++)
	{
		classify_struct_example(val.examples[i].x,&y,sm,sparm);
		store = loss(val.examples[i].y,y,sparm);
		cur_loss += store;
	}

	cur_loss /= (double) val.n;
	return cur_loss;
}


int main(int argc, char* argv[]) {

  double *w; /* weight vector */
  long m, i;
  double C, epsilon;
  LEARN_PARM learn_parm;
  KERNEL_PARM kernel_parm;
  char trainfile[1024];
  char modelfile[1024];
  int MAX_ITER;
  /* new struct variables */
  SVECTOR **fycache, *diff, *fy;
  EXAMPLE *ex;
	SAMPLE alldata;
  SAMPLE sample;
	SAMPLE val;
  STRUCT_LEARN_PARM sparm;
  STRUCTMODEL sm;
  
  double primal_obj;
  double stop_crit; 
	char itermodelfile[2000];

	/* self-paced learning variables */
	double init_spl_weight;
	double spl_weight;
	double spl_factor;
	int *valid_examples;
 

  /* read input parameters */
	my_read_input_parameters(argc, argv, trainfile, modelfile, &learn_parm, &kernel_parm, &sparm, 
													&init_spl_weight, &spl_factor); 

  epsilon = learn_parm.eps;
  C = learn_parm.svm_c;
  MAX_ITER = learn_parm.maxiter;

  /* read in examples */
  alldata = read_struct_examples(trainfile,&sparm);
  int ntrain = (int) round(1.0*alldata.n); /* no validation set */
	if(ntrain < alldata.n)
	{
 	 long *perm = randperm(alldata.n);
 	 sample = generate_train_set(alldata, perm, ntrain);
 	 val = generate_validation_set(alldata, perm, ntrain);
 	 free(perm);
	}
	else
	{
		sample = alldata;
	}
  ex = sample.examples;
  m = sample.n;
  
  /* initialization */
  init_struct_model(alldata,&sm,&sparm,&learn_parm,&kernel_parm); 

  w = create_nvector(sm.sizePsi);
  clear_nvector(w, sm.sizePsi);
  sm.w = w; /* establish link to w, as long as w does not change pointer */

  /* some training information */
  printf("C: %.8g\n", C);
	printf("spl weight: %.8g\n",init_spl_weight);
  printf("epsilon: %.8g\n", epsilon);
  printf("sample.n: %d\n", sample.n); 
  printf("sm.sizePsi: %ld\n", sm.sizePsi); fflush(stdout);


  /* prepare feature vector cache for correct labels with imputed latent variables */
  fycache = (SVECTOR**)malloc(m*sizeof(SVECTOR*));
  for (i=0;i<m;i++) {
    fy = psi(ex[i].x, ex[i].y, &sm, &sparm);
    diff = add_list_ss(fy);
    free_svector(fy);
    fy = diff;
    fycache[i] = fy;
  }

 	/* learn initial weight vector using all training examples */
	valid_examples = (int *) malloc(m*sizeof(int));     

  /* errors for validation set */

  double cur_loss, best_loss = DBL_MAX;
  int loss_iter;


	/* initializations */
	spl_weight = init_spl_weight;

	/* solve biconvex self-paced learning problem */
	primal_obj = alternate_convex_search(w, m, MAX_ITER, C, epsilon, fycache, ex, &sm, &sparm, valid_examples, spl_weight);
	printf("primal objective: %.4f\n", primal_obj);
	fflush(stdout);
	//alternate_convex_search(w, m, MAX_ITER, C, epsilon, fycache, ex, &sm, &sparm, valid_examples, spl_weight);
	int nValid = 0;
	for (i=0;i<m;i++) {
		if(valid_examples[i]) {
			nValid++;
		}
	}

		

	if(ntrain < alldata.n) {
		cur_loss = compute_current_loss(val,&sm,&sparm);
		printf("CURRENT LOSS: %f\n",cur_loss);
	}
  

  /* write structural model */
  write_struct_model(modelfile, &sm, &sparm);
  // skip testing for the moment  

  /* free memory */
  free_struct_sample(alldata);
	if(ntrain < alldata.n)
	{
		free(sample.examples);
		free(val.examples);
	}
  free_struct_model(sm, &sparm);
  for(i=0;i<m;i++) {
    free_svector(fycache[i]);
  }
  free(fycache);

	free(valid_examples);
   
  return(0); 
  
}



void my_read_input_parameters(int argc, char *argv[], char *trainfile, char* modelfile, 
			      LEARN_PARM *learn_parm, KERNEL_PARM *kernel_parm, STRUCT_LEARN_PARM *struct_parm,
						double *init_spl_weight, double *spl_factor) {
  
  long i;

  /* set default */
  learn_parm->maxiter=20000;
  learn_parm->svm_maxqpsize=100;
  learn_parm->svm_c=100.0;
  learn_parm->eps=0.001;
  learn_parm->biased_hyperplane=12345; /* store random seed */
  learn_parm->remove_inconsistent=10; 
  kernel_parm->kernel_type=0;
  kernel_parm->rbf_gamma=0.05;
  kernel_parm->coef_lin=1;
  kernel_parm->coef_const=1;
  kernel_parm->poly_degree=3;
	/* default: no self-paced learning */
	*init_spl_weight = 0.0;
	*spl_factor = 1.3;

	struct_parm->gram_regularization = 1E-7;
  struct_parm->solve_dual = 1;

  struct_parm->custom_argc=0;

  for(i=1;(i<argc) && ((argv[i])[0] == '-');i++) {
    switch ((argv[i])[1]) {
    case 'c': i++; learn_parm->svm_c=atof(argv[i]); break;
    case 'e': i++; learn_parm->eps=atof(argv[i]); break;
    case 's': i++; learn_parm->svm_maxqpsize=atol(argv[i]); break; 
    case 'g': i++; kernel_parm->rbf_gamma=atof(argv[i]); break;
    case 'd': i++; kernel_parm->poly_degree=atol(argv[i]); break;
    case 'r': i++; learn_parm->biased_hyperplane=atol(argv[i]); break; 
    case 't': i++; kernel_parm->kernel_type=atol(argv[i]); break;
    case 'n': i++; learn_parm->maxiter=atol(argv[i]); break;
    case 'p': i++; learn_parm->remove_inconsistent=atol(argv[i]); break; 
		case 'k': i++; *init_spl_weight = atof(argv[i]); break;
		case 'm': i++; *spl_factor = atof(argv[i]); break;
    case 'q': i++; struct_parm->solve_dual = atoi(argv[i]); break;
    case '-': strcpy(struct_parm->custom_argv[struct_parm->custom_argc++],argv[i]);i++; strcpy(struct_parm->custom_argv[struct_parm->custom_argc++],argv[i]);break; 
    default: printf("\nUnrecognized option %s!\n\n",argv[i]);
      exit(0);
    }

  }
	*init_spl_weight = (*init_spl_weight)/learn_parm->svm_c;

  if(i>=argc) {
    printf("\nNot enough input parameters!\n\n");
    my_wait_any_key();
    exit(0);
  }
  strcpy (trainfile, argv[i]);

  if((i+1)<argc) {
    strcpy (modelfile, argv[i+1]);
  }
	else {
		strcpy (modelfile, "lssvm.model");
	}

	/* self-paced learning weight should be non-negative */
	if(*init_spl_weight < 0.0)
		*init_spl_weight = 0.0;
	/* self-paced learning factor should be greater than 1.0 */
	if(*spl_factor < 1.0)
		*spl_factor = 1.1;

  
  parse_struct_parameters(struct_parm);
}


void my_wait_any_key()
{
  printf("\n(more)\n");
  (void)getc(stdin);
}

int resize_cleanup(int size_active, int **ptr_idle, double **ptr_cur_slack, double **ptr_delta, DOC ***ptr_dXc, 
		double ***ptr_psiDiffs, int *mv_iter) {
  int i,j, new_size_active;
  long k;

  int *idle=*ptr_idle;
  double *cur_slack=*ptr_cur_slack;
  double *delta=*ptr_delta;
	DOC	**dXc = *ptr_dXc;
	double **psiDiffs = *ptr_psiDiffs;
	int new_mv_iter;

  i=0;
  while ((i<size_active)&&(idle[i]<IDLE_ITER)) i++;
  j=i;
  while((j<size_active)&&(idle[j]>=IDLE_ITER)) j++;

  while (j<size_active) {
    // copying 
    cur_slack[i] = cur_slack[j];
    delta[i] = delta[j];
	free(psiDiffs[i]);
	psiDiffs[i] = psiDiffs[j];
	psiDiffs[j] = NULL;
    free_example(dXc[i],1);
    dXc[i] = dXc[j];
    dXc[j] = NULL;
	if(j == *mv_iter)
		new_mv_iter = i;

    i++;
    j++;
    while((j<size_active)&&(idle[j]>=IDLE_ITER)) j++;
  }
  for (k=i;k<size_active;k++) {
	if (psiDiffs[k]!=NULL) free(psiDiffs[k]);
    if (dXc[k]!=NULL) free_example(dXc[k],1);
  }
  *mv_iter = new_mv_iter;
  new_size_active = i;
  cur_slack = (double*)realloc(cur_slack, sizeof(double)*new_size_active);
  delta = (double*)realloc(delta, sizeof(double)*new_size_active);
  psiDiffs = (double **) realloc(psiDiffs, sizeof(double *)*new_size_active);
  dXc = (DOC**)realloc(dXc, sizeof(DOC*)*new_size_active);
  assert(dXc!=NULL);

  // resize idle 
  i=0;
  while ((i<size_active)&&(idle[i]<IDLE_ITER)) i++;
  j=i;
  while((j<size_active)&&(idle[j]>=IDLE_ITER)) j++;

  while (j<size_active) {
    idle[i] = idle[j];
    i++;
    j++;
    while((j<size_active)&&(idle[j]>=IDLE_ITER)) j++;
  }  
  idle = (int*)realloc(idle, sizeof(int)*new_size_active);

  *ptr_idle = idle;
  *ptr_cur_slack = cur_slack;
  *ptr_delta = delta;
  *ptr_psiDiffs = psiDiffs;
  *ptr_dXc = dXc;

  return(new_size_active);
}

/*int resize_cleanup(int size_active, int **ptr_idle, double **ptr_alpha, double **ptr_delta, DOC ***ptr_dXc, 
		double ***ptr_G, int *mv_iter) {
  int i,j, new_size_active;
  long k;

  int *idle=*ptr_idle;
  double *alpha=*ptr_alpha;
  double *delta=*ptr_delta;
	DOC	**dXc = *ptr_dXc;
	double **G = *ptr_G;
	int new_mv_iter;

  i=0;
  while ((i<size_active)&&(idle[i]<IDLE_ITER)) i++;
  j=i;
  while((j<size_active)&&(idle[j]>=IDLE_ITER)) j++;

  while (j<size_active) {
    // copying 
    alpha[i] = alpha[j];
    delta[i] = delta[j];
		free(G[i]);
		G[i] = G[j];
		G[j] = NULL;
    free_example(dXc[i],1);
    dXc[i] = dXc[j];
    dXc[j] = NULL;
		if(j == *mv_iter)
			new_mv_iter = i;

    i++;
    j++;
    while((j<size_active)&&(idle[j]>=IDLE_ITER)) j++;
  }
  for (k=i;k<size_active;k++) {
		if (G[k]!=NULL) free(G[k]);
    if (dXc[k]!=NULL) free_example(dXc[k],1);
  }
	*mv_iter = new_mv_iter;
  new_size_active = i;
  alpha = (double*)realloc(alpha, sizeof(double)*new_size_active);
  delta = (double*)realloc(delta, sizeof(double)*new_size_active);
	G = (double **) realloc(G, sizeof(double *)*new_size_active);
  dXc = (DOC**)realloc(dXc, sizeof(DOC*)*new_size_active);
  assert(dXc!=NULL);

  // resize idle 
  i=0;
  while ((i<size_active)&&(idle[i]<IDLE_ITER)) i++;
  j=i;
  while((j<size_active)&&(idle[j]>=IDLE_ITER)) j++;

  while (j<size_active) {
    idle[i] = idle[j];
		for (k=0;k<new_size_active;k++) {
			G[k][i] = G[k][j];
		}
    i++;
    j++;
    while((j<size_active)&&(idle[j]>=IDLE_ITER)) j++;
  }  
  idle = (int*)realloc(idle, sizeof(int)*new_size_active);
	for (k=0;k<new_size_active;k++) {
		G[k] = (double*)realloc(G[k], sizeof(double)*new_size_active);
	}

  *ptr_idle = idle;
  *ptr_alpha = alpha;
  *ptr_delta = delta;
	*ptr_G = G;
  *ptr_dXc = dXc;

  return(new_size_active);
}*/

void approximate_to_psd(double **G, int size_active, double eps)
{
	int i,j,k;
	double *copy_G = (double *)malloc(size_active*size_active*sizeof(double));
	double *eig_vec = (double *)malloc(size_active*size_active*sizeof(double));
	double *eig_val = (double *)malloc(size_active*sizeof(double));

	for(i = 0; i < size_active; i++)
		for(j = 0; j < size_active; j++)
			copy_G[i*size_active+j] = G[i][j];

	Jacobi_Cyclic_Method(eig_val,eig_vec,copy_G,size_active);

	for(i = 0; i < size_active; i++)
		for(j = 0; j < size_active; j++)
		{
			copy_G[i*size_active+j] = MAX(eig_val[i],eps)*eig_vec[j*size_active+i];
		}

	for(i = 0; i < size_active; i++)
		for(j = 0; j < size_active; j++)
		{
			G[i][j] = 0.0;
			for(k = 0; k < size_active; k++)
			{
				G[i][j] += eig_vec[i*size_active+k]*copy_G[k*size_active+j];
			}
		}

	free(copy_G);
	free(eig_vec);
	free(eig_val);
}

void Jacobi_Cyclic_Method(double eigenvalues[], double *eigenvectors, double *A, int n)
{
   int row, i, j, k, m;
   double *pAk, *pAm, *p_r, *p_e;
   double threshold_norm;
   double threshold;
   double tan_phi, sin_phi, cos_phi, tan2_phi, sin2_phi, cos2_phi;
   double sin_2phi, cos_2phi, cot_2phi;
   double dum1;
   double dum2;
   double dum3;
   double r;
   double max;

                  // Take care of trivial cases

   if ( n < 1) return;
   if ( n == 1) {
      eigenvalues[0] = *A;
      *eigenvectors = 1.0;
      return;
   }

          // Initialize the eigenvalues to the identity matrix.

   for (p_e = eigenvectors, i = 0; i < n; i++)
      for (j = 0; j < n; p_e++, j++)
         if (i == j) *p_e = 1.0; else *p_e = 0.0;
  
            // Calculate the threshold and threshold_norm.
 
   for (threshold = 0.0, pAk = A, i = 0; i < ( n - 1 ); pAk += n, i++) 
      for (j = i + 1; j < n; j++) threshold += *(pAk + j) * *(pAk + j);
   threshold = sqrt(threshold + threshold);
   threshold_norm = threshold * DBL_EPSILON;
   max = threshold + 1.0;
   while (threshold > threshold_norm) {
      threshold /= 10.0;
      if (max < threshold) continue;
      max = 0.0;
      for (pAk = A, k = 0; k < (n-1); pAk += n, k++) {
         for (pAm = pAk + n, m = k + 1; m < n; pAm += n, m++) {
            if ( fabs(*(pAk + m)) < threshold ) continue;

                 // Calculate the sin and cos of the rotation angle which
                 // annihilates A[k][m].

            cot_2phi = 0.5 * ( *(pAk + k) - *(pAm + m) ) / *(pAk + m);
            dum1 = sqrt( cot_2phi * cot_2phi + 1.0);
            if (cot_2phi < 0.0) dum1 = -dum1;
            tan_phi = -cot_2phi + dum1;
            tan2_phi = tan_phi * tan_phi;
            sin2_phi = tan2_phi / (1.0 + tan2_phi);
            cos2_phi = 1.0 - sin2_phi;
            sin_phi = sqrt(sin2_phi);
            if (tan_phi < 0.0) sin_phi = - sin_phi;
            cos_phi = sqrt(cos2_phi); 
            sin_2phi = 2.0 * sin_phi * cos_phi;
            cos_2phi = cos2_phi - sin2_phi;

                     // Rotate columns k and m for both the matrix A 
                     //     and the matrix of eigenvectors.

            p_r = A;
            dum1 = *(pAk + k);
            dum2 = *(pAm + m);
            dum3 = *(pAk + m);
            *(pAk + k) = dum1 * cos2_phi + dum2 * sin2_phi + dum3 * sin_2phi;
            *(pAm + m) = dum1 * sin2_phi + dum2 * cos2_phi - dum3 * sin_2phi;
            *(pAk + m) = 0.0;
            *(pAm + k) = 0.0;
            for (i = 0; i < n; p_r += n, i++) {
               if ( (i == k) || (i == m) ) continue;
               if ( i < k ) dum1 = *(p_r + k); else dum1 = *(pAk + i);
               if ( i < m ) dum2 = *(p_r + m); else dum2 = *(pAm + i);
               dum3 = dum1 * cos_phi + dum2 * sin_phi;
               if ( i < k ) *(p_r + k) = dum3; else *(pAk + i) = dum3;
               dum3 = - dum1 * sin_phi + dum2 * cos_phi;
               if ( i < m ) *(p_r + m) = dum3; else *(pAm + i) = dum3;
            }
            for (p_e = eigenvectors, i = 0; i < n; p_e += n, i++) {
               dum1 = *(p_e + k);
               dum2 = *(p_e + m);
               *(p_e + k) = dum1 * cos_phi + dum2 * sin_phi;
               *(p_e + m) = - dum1 * sin_phi + dum2 * cos_phi;
            }
         }
         for (i = 0; i < n; i++)
            if ( i == k ) continue;
            else if ( max < fabs(*(pAk + i))) max = fabs(*(pAk + i));
      }
   }
   for (pAk = A, k = 0; k < n; pAk += n, k++) eigenvalues[k] = *(pAk + k); 
}
