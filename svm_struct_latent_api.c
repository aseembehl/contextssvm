/************************************************************************/
/*                                                                      */
/*   svm_struct_latent_api.c                                            */
/*                                                                      */
/*   API function definitions for Latent SVM^struct                     */
/*                                                                      */
/*   Author: Chun-Nam Yu                                                */
/*   Date: 17.Dec.08                                                    */
/*                                                                      */
/*   This software is available for non-commercial use only. It must    */
/*   not be modified and distributed without prior permission of the    */
/*   author. The author is not responsible for implications from the    */
/*   use of this software.                                              */
/*                                                                      */
/************************************************************************/

#include <stdio.h>
#include <assert.h>
#include <string.h>

  #include "svm_struct_latent_api_types.h"

#include <errno.h>

#include "maxflowwrap.hpp"


#define MAX_INPUT_LINE_LENGTH 10000

void die(const char *message)
{
  if(errno) {
      perror(message); 
  } else {
      printf("ERROR: %s\n", message);
  }
  exit(1);
}

SVECTOR *read_sparse_vector(char *file_name, int object_id, STRUCT_LEARN_PARM *sparm){
    
    int scanned;
    WORD *words = NULL;
    char feature_file[1000];
    sprintf(feature_file, "%s_%d.feature", file_name, object_id);
    FILE *fp = fopen(feature_file, "r");
    
    int length = 0;
    while(!feof(fp)){
        length++;
        words = (WORD *) realloc(words, length*sizeof(WORD));
        if(!words) die("Memory error."); 
        scanned = fscanf(fp, " %d:%f", &words[length-1].wnum, &words[length-1].weight);
        if(scanned < 2) {
            words[length-1].wnum = 0;
            words[length-1].weight = 0.0;
        }
    }
    fclose(fp);

  SVECTOR *fvec = create_svector(words,"",1);
  free(words);

  return fvec;
}

SAMPLE read_struct_examples(char *file, STRUCT_LEARN_PARM *sparm) {
/*
  Read input examples {(x_1,y_1),...,(x_n,y_n)} from file.
  The type of pattern x and label y has to follow the definition in 
  svm_struct_latent_api_types.h.  
*/
  SAMPLE sample;

  int i;

  // open the file containing candidate bounding box dimensions/labels/featurePath and image label
  FILE *fp = fopen(file, "r");
  if(fp==NULL){
      printf("Error: Cannot open input file %s\n",file);
      exit(1);
  }

  sample.n = 1;  
  sample.examples = (EXAMPLE *) malloc(sample.n*sizeof(EXAMPLE));
  if(!sample.examples) die("Memory error.");
  sample.examples[0].x.n_pos = 0;
  sample.examples[0].x.n_neg = 0;

  fscanf(fp,"%d", &sample.examples[0].n_imgs);
    
  // Initialise pattern 
  sample.examples[0].x.example_cost = 1;

  sample.examples[0].x.x_is = (SUB_PATTERN *) malloc(sample.examples[0].n_imgs*sizeof(SUB_PATTERN));
  if(!sample.examples[0].x.x_is) die("Memory error.");
  sample.examples[0].y.labels = (int *) malloc(sample.examples[0].n_imgs*sizeof(int));
  if(!sample.examples[0].y.labels) die("Memory error.");

  for(i = 0; i < sample.examples[0].n_imgs; i++){  
      fscanf(fp,"%s",sample.examples[0].x.x_is[i].phi1_file_name);
      fscanf(fp,"%s",sample.examples[0].x.x_is[i].phi2_file_name);
      fscanf(fp, "%d", &sample.examples[0].x.x_is[i].id);
      fscanf(fp, "%d", &sample.examples[0].y.labels[i]);

      sample.examples[0].x.x_is[i].phi1_pos = read_sparse_vector(sample.examples[0].x.x_is[i].phi1_file_name, sample.examples[0].x.x_is[i].id, sparm);
      sample.examples[0].x.x_is[i].phi1_neg = create_svector_with_index(sample.examples[0].x.x_is[i].phi1_pos->words, "", 1, sparm->phi1_size);
      sample.examples[0].x.x_is[i].phi2 = read_sparse_vector(sample.examples[0].x.x_is[i].phi2_file_name, sample.examples[0].x.x_is[i].id, sparm);

      if(sample.examples[0].y.labels[i] == 1) {
          sample.examples[0].x.n_pos++;
      } 
      else{
          sample.examples[0].x.n_neg++;
      }
  }
  sample.examples[0].y.n_pos = sample.examples[0].x.n_pos;
  sample.examples[0].y.n_neg = sample.examples[0].x.n_neg;

  return(sample);
}

void init_struct_model(SAMPLE sample, STRUCTMODEL *sm, STRUCT_LEARN_PARM *sparm, LEARN_PARM *lparm, KERNEL_PARM *kparm) {
/*
  Initialize parameters in STRUCTMODEL sm. Set the diminension 
  of the feature space sm->sizePsi. Can also initialize your own
  variables in sm here. 
*/

	sm->n = sample.n;
  // \psi is concatanation of \psi1 and \psi2. Dimension is the sum of dimensions of \psi1 and \psi2
  sm->sizePsi = sparm->phi1_size*2 + (sparm->phi1_size+sparm->phi2_size);

}

SVECTOR *psi(PATTERN x, LABEL y, STRUCTMODEL *sm, STRUCT_LEARN_PARM *sparm) {
/*
  Creates the feature vector \Psi(x,y) and return a pointer to 
  sparse vector SVECTOR in SVM^light format. The dimension of the 
  feature vector returned has to agree with the dimension in sm->sizePsi. 
*/
  SVECTOR *fvec, *psi1, *psi2_1, *psi2_2=NULL;
  SVECTOR *temp_psi, *temp_sub, *temp_fvec=NULL;


  WORD *words = NULL;
  words = (WORD *) malloc(sizeof(WORD));
  if(!words) die("Memory error."); 
  words[0].wnum = 0;
  words[0].weight = 0;
  fvec = create_svector(words,"",1);
  psi1 = create_svector(words,"",1);
  psi2_1 = create_svector(words,"",1);
  psi2_2 = create_svector(words,"",1);
  free(words);

  int i,j = 0;
  
  for (i = 0; i < (x.n_pos+x.n_neg); i++){
      if(y.labels[i] == 1){
          temp_psi = add_ss(psi1, x.x_is[i].phi1_pos);
      }  
      else{
          temp_psi = add_ss(psi1, x.x_is[i].phi1_neg);
      }
      free_svector(psi1);
      psi1 = temp_psi;

      for (j= 0; j < (x.n_pos+x.n_neg); j++){
          if(y.labels[i] != y.labels[j]){
              temp_sub = sub_ss_abs(x.x_is[i].phi1_pos, x.x_is[j].phi1_pos);
              temp_psi = add_ss(psi2_1, temp_sub);
              free_svector(temp_sub);
              free_svector(psi2_1);
              psi2_1 = temp_psi;
              
              temp_sub = sub_ss_abs(x.x_is[i].phi2, x.x_is[j].phi2);
              temp_psi = add_ss(psi2_2, temp_sub);
              free_svector(temp_sub);              
              free_svector(psi2_2);              
              psi2_2 = temp_psi;
          }
      }
  }
  
  // scale w1 by 1/n
  temp_psi = smult_s(psi1, (float)1/(float)(x.n_pos+x.n_neg));
  free_svector(psi1);
  psi1 = temp_psi;
  
  // scale w2_2 by 1/n^2
  temp_psi = smult_s(psi2_1, (float)1/(float)((x.n_pos+x.n_neg)*(x.n_pos+x.n_neg)));
  free_svector(psi2_1);
  psi2_1 = temp_psi;
  
  // scale w2_2 by 1/n^2
  temp_psi = smult_s(psi2_2, (float)1/(float)((x.n_pos+x.n_neg)*(x.n_pos+x.n_neg)));
  free_svector(psi2_2);
  psi2_2 = temp_psi;
  
  // concatenate psi1, psi2_1 and psi2_2
  temp_psi = create_svector_with_index(psi2_1->words, "", 1, sparm->phi1_size*2);
  free_svector(psi2_1);
  temp_fvec = add_ss(psi1, temp_psi);
  free_svector(temp_psi);
  free_svector(psi1);
  temp_psi = create_svector_with_index(psi2_2->words, "", 1, (sparm->phi1_size*2 + sparm->phi1_size));
  free_svector(psi2_2);
  fvec = add_ss(temp_fvec, temp_psi);
  free_svector(temp_psi);
  free_svector(temp_fvec);  
  
  return(fvec);
}

void classify_struct_example(PATTERN x, LABEL *y, STRUCTMODEL *sm, STRUCT_LEARN_PARM *sparm) {
/*
  Makes prediction with input pattern x with weight vector in sm->w,
  i.e., computing argmax_{(y)} <w,psi(x,y)>. 
  Output pair (y) is stored at location pointed to by 
  pointers *y. 
*/
  int i,j;

  SVECTOR *temp_sub, *temp_sub_shifted=NULL;

  double *unary_pos = (double*)malloc((x.n_pos+x.n_neg)*sizeof(double));
  double *unary_neg =  (double*)malloc((x.n_pos+x.n_neg)*sizeof(double));
  double **binary =  (double**)malloc((x.n_pos+x.n_neg)*sizeof(double *));


  for (i = 0; i < (x.n_pos+x.n_neg); i++){
      binary[i] =  (double*)malloc((x.n_pos+x.n_neg)*sizeof(double));
      // compute unary potential for ybar.labels[i] == 1 
      unary_pos[i] = sprod_ns(sm->w, x.x_is[i].phi1_pos);
      if(unary_pos[i] != 0){
        unary_pos[i] = (float)(-1*unary_pos[i])/(float)(x.n_pos+x.n_neg);
      }
      // compute unary potential for ybar.labels[i] == -1
      unary_neg[i] = sprod_ns(sm->w, x.x_is[i].phi1_neg);
      if(unary_neg[i] != 0){
        unary_neg[i] = (float)(-1*unary_neg[i])/(float)(x.n_pos+x.n_neg);
      }
      for (j = (i+1); j < (x.n_pos+x.n_neg); j++){
          temp_sub = sub_ss_abs(x.x_is[i].phi1_pos, x.x_is[j].phi1_pos);
          temp_sub_shifted = create_svector_with_index(temp_sub->words, "", 1, sparm->phi1_size*2);
          binary[i][j] = sprod_ns(sm->w, temp_sub_shifted);
          free_svector(temp_sub);
          free_svector(temp_sub_shifted);

          temp_sub = sub_ss_abs(x.x_is[i].phi2, x.x_is[j].phi2);
          temp_sub_shifted = create_svector_with_index(temp_sub->words, "", 1, sparm->phi1_size*3);
          binary[i][j] += sprod_ns(sm->w, temp_sub_shifted);
          free_svector(temp_sub);
          free_svector(temp_sub_shifted);

          binary[i][j] = (double)(-1*binary[i][j])/(double)((x.n_pos+x.n_neg)*(x.n_pos+x.n_neg));
      }
  }

  y->labels = maxflowwrapper(unary_pos, unary_neg, binary, x.n_pos, x.n_neg);

  free(unary_pos);
  free(unary_neg);
  for (i = 0; i < (x.n_pos+x.n_neg); i++){
    free(binary[i]);
  }
  free(binary);
  
	return;

}

void find_most_violated_constraint_marginrescaling(PATTERN x, LABEL y, LABEL *ybar, STRUCTMODEL *sm, STRUCT_LEARN_PARM *sparm) {
/*
  Finds the most violated constraint (loss-augmented inference), i.e.,
  computing argmax_{(ybar,hbar)} [<w,psi(x,ybar,hbar)> + loss(y,ybar,hbar)].
  The output (ybar,hbar) are stored at location pointed by 
  pointers *ybar and *hbar. 
*/
  int i, j;

  SVECTOR *temp_sub, *temp_sub_shifted=NULL;

  double *unary_pos = (double*)malloc((x.n_pos+x.n_neg)*sizeof(double));
  double *unary_neg =  (double*)malloc((x.n_pos+x.n_neg)*sizeof(double));
  double **binary =  (double**)malloc((x.n_pos+x.n_neg)*sizeof(double *));


  for (i = 0; i < (x.n_pos+x.n_neg); i++){
      binary[i] =  (double*)malloc((x.n_pos+x.n_neg)*sizeof(double));
      // compute unary potential for ybar.labels[i] == 1 
      unary_pos[i] = sprod_ns(sm->w, x.x_is[i].phi1_pos);
      if(unary_pos[i] != 0){
        unary_pos[i] = (float)(-1*unary_pos[i])/(float)(x.n_pos+x.n_neg);
      }
      // compute unary potential for ybar.labels[i] == -1
      unary_neg[i] = sprod_ns(sm->w, x.x_is[i].phi1_neg);
      if(unary_neg[i] != 0){
        unary_neg[i] = (float)(-1*unary_neg[i])/(float)(x.n_pos+x.n_neg);
      }

      if(y.labels[i] == 1){
          // add 1/n to 'ybar == -1' unary term
          unary_neg[i] -= (float)1/(float)(x.n_pos+x.n_neg);
      }
      else{
          // add 1/n to 'ybar == 1' unary term
          unary_pos[i] -= (float)1/(float)(x.n_pos+x.n_neg);
      }

      for (j = (i+1); j < (x.n_pos+x.n_neg); j++){
          temp_sub = sub_ss_abs(x.x_is[i].phi1_pos, x.x_is[j].phi1_pos);
          temp_sub_shifted = create_svector_with_index(temp_sub->words, "", 1, sparm->phi1_size*2);
          binary[i][j] = sprod_ns(sm->w, temp_sub_shifted);
          free_svector(temp_sub);
          free_svector(temp_sub_shifted);

          temp_sub = sub_ss_abs(x.x_is[i].phi2, x.x_is[j].phi2);
          temp_sub_shifted = create_svector_with_index(temp_sub->words, "", 1, sparm->phi1_size*3);
          binary[i][j] += sprod_ns(sm->w, temp_sub_shifted);
          free_svector(temp_sub);
          free_svector(temp_sub_shifted);

          if(binary[i][j] != 0){
              binary[i][j] = (double)(-1*binary[i][j])/(double)((x.n_pos+x.n_neg)*(x.n_pos+x.n_neg));
          }
      }
  }

  ybar->labels = maxflowwrapper(unary_pos, unary_neg, binary, x.n_pos, x.n_neg);

  free(unary_pos);
  free(unary_neg);
  for (i = 0; i < (x.n_pos+x.n_neg); i++){
    free(binary[i]);
  }
  free(binary);

	return;

}

double loss(LABEL y, LABEL ybar, STRUCT_LEARN_PARM *sparm) {
/*
  Computes the loss of prediction (ybar,hbar) against the
  correct label y. 
*/ 
	double l  = 0;

  int i;
  for (i = 0; i < (y.n_pos+y.n_neg); i++){
    if (y.labels[i] != ybar.labels[i]){
      l++;
    }
  }
  l = l/(double)(y.n_pos+y.n_neg);

	return(l);

}

void write_struct_model(char *file, STRUCTMODEL *sm, STRUCT_LEARN_PARM *sparm) {
/*
  Writes the learned weight vector sm->w to file after training. 
*/
  FILE *modelfl;
  int i;
  
  modelfl = fopen(file,"w");
  if (modelfl==NULL) {
    printf("Cannot open model file %s for output!", file);
		exit(1);
  }
  
  for (i=1;i<sm->sizePsi+1;i++) {
    fprintf(modelfl, "%d:%.16g\n", i, sm->w[i]);
  }
  fclose(modelfl);
 
}

STRUCTMODEL read_struct_model(char *file, STRUCT_LEARN_PARM *sparm) {
/*
  Reads in the learned model parameters from file into STRUCTMODEL sm.
  The input file format has to agree with the format in write_struct_model().
*/
  STRUCTMODEL sm;

  FILE *modelfl;
  int sizePsi,i, fnum;
  double fweight;
  
  modelfl = fopen(file,"r");
  if (modelfl==NULL) {
    printf("Cannot open model file %s for input!", file);
	exit(1);
  }

	sizePsi = 1;
	sm.w = (double*)malloc((sizePsi+1)*sizeof(double));
	for (i=0;i<sizePsi+1;i++) {
		sm.w[i] = 0.0;
	}
	while (!feof(modelfl)) {
		fscanf(modelfl, "%d:%lf", &fnum, &fweight);
		if(fnum > sizePsi) {
			sizePsi = fnum;
			sm.w = (double *)realloc(sm.w,(sizePsi+1)*sizeof(double));
		}
		sm.w[fnum] = fweight;
	}

	fclose(modelfl);

	sm.sizePsi = sizePsi;

  return(sm);

}

void free_struct_model(STRUCTMODEL sm, STRUCT_LEARN_PARM *sparm) {
/*
  Free any memory malloc'ed in STRUCTMODEL sm after training. 
*/

  free(sm.w);

}

void free_pattern(PATTERN x) {
/*
  Free any memory malloc'ed when creating pattern x. 
*/

}

void free_label(LABEL y) {
/*
  Free any memory malloc'ed when creating label y. 
*/
  free(y.labels);

} 

void free_struct_sample(SAMPLE s) {
/*
  Free the whole training sample. 
*/
  int i;
  for (i=0;i<s.n;i++) {
    free_pattern(s.examples[i].x);
    free_label(s.examples[i].y);
  }
  free(s.examples);

}

void parse_struct_parameters(STRUCT_LEARN_PARM *sparm) {
/*
  Parse parameters for structured output learning passed 
  via the command line. 
*/
  int i;
  
  /* set default */
  sparm->phi1_size=24004;
  sparm->phi2_size=512;
  
  for (i=0;(i<sparm->custom_argc)&&((sparm->custom_argv[i])[0]=='-');i++) {
    switch ((sparm->custom_argv[i])[2]) {
      /* your code here */
      default: printf("\nUnrecognized option %s!\n\n", sparm->custom_argv[i]); exit(0);
    }
  }

}

void copy_label(LABEL l1, LABEL *l2)
{
}

void print_label(LABEL l, FILE	*flabel)
{
}
