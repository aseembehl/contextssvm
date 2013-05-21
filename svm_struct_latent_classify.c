/************************************************************************/
/*                                                                      */
/*   svm_struct_latent_classify.c                                       */
/*                                                                      */
/*   Classification Code for Latent SVM^struct                          */
/*                                                                      */
/*   Author: Chun-Nam Yu                                                */
/*   Date: 9.Nov.08                                                     */
/*                                                                      */
/*   This software is available for non-commercial use only. It must    */
/*   not be modified and distributed without prior permission of the    */
/*   author. The author is not responsible for implications from the    */
/*   use of this software.                                              */
/*                                                                      */
/************************************************************************/

#include <stdio.h>
#include "svm_struct_latent_api.h"


void read_input_parameters(int argc, char **argv, char *testfile, char *modelfile, char *outfile, STRUCT_LEARN_PARM *sparm);


int main(int argc, char* argv[]) {
  double avgloss,l;
  long i, correct;

  char testfile[1024];
  char modelfile[1024];
	char outfile[1024];
	FILE	*foutfile;

  STRUCTMODEL model;
  STRUCT_LEARN_PARM sparm;
  LEARN_PARM lparm;
  KERNEL_PARM kparm;

  SAMPLE testsample;
  LABEL y;

  /* read input parameters */
  read_input_parameters(argc,argv,testfile,modelfile,outfile,&sparm);
	foutfile = fopen(outfile,"w");

  /* read model file */
  printf("Reading model..."); fflush(stdout);
  model = read_struct_model(modelfile, &sparm);
  printf("done.\n"); 

  /* read test examples */
	printf("Reading test examples..."); fflush(stdout);
  testsample = read_struct_examples(testfile,&sparm);
	printf("done.\n");

  init_struct_model(testsample,&model,&sparm,&lparm,&kparm);
  
  avgloss = 0.0;
  correct = 0;
  for (i=0;i<testsample.n;i++) {
    classify_struct_example(testsample.examples[i].x,&y,&model,&sparm);
    l = loss(testsample.examples[i].y,y,&sparm);
    avgloss += l;
    if (l==0) correct++;

		//print_label(y,flabel);
		//fprintf(foutfile,"\n"); 

    free_label(y);
  }
  fprintf(foutfile, "%0.4f\n", avgloss/testsample.n);
  fflush(foutfile);
	fclose(foutfile);

  printf("Average loss on test set: %.4f\n", avgloss/testsample.n);
  
  //printf("Zero/one error on test set: %.4f\n", 1.0 - ((float) correct)/testsample.n);

  free_struct_sample(testsample);
  free_struct_model(model,&sparm);

  return(0);

}


void read_input_parameters(int argc, char **argv, char *testfile, char *modelfile, char *outfile, STRUCT_LEARN_PARM *sparm) {

  long i;
  
  /* set default */
  strcpy(modelfile, "lssvm_model");
  strcpy(outfile, "lssvm_outfile");
  sparm->custom_argc = 0;

  for (i=1;(i<argc)&&((argv[i])[0]=='-');i++) {
    switch ((argv[i])[1]) {
      case '-': strcpy(sparm->custom_argv[sparm->custom_argc++],argv[i]);i++; strcpy(sparm->custom_argv[sparm->custom_argc++],argv[i]);break;  
      default: printf("\nUnrecognized option %s!\n\n",argv[i]); exit(0);    
    }
  }

  if (i>=argc) {
    printf("\nNot enough input parameters!\n\n");
    exit(0);
  }

  strcpy(testfile, argv[i]);
	if(i+1<argc)
  	strcpy(modelfile, argv[i+1]);
	if(i+2<argc)
		strcpy(outfile,argv[i+2]);

  parse_struct_parameters(sparm);

}
