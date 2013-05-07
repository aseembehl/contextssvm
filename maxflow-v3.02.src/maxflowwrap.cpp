#include <stdlib.h>

#include "graph.h"

int scaleFactor = 1e6;
extern "C" int *maxflowwrapper(double *unary_pos, double *unary_neg, double **binary, int n_pos, int n_neg){
	int i,j;
	int *labels;
	typedef Graph<int,int,int> GraphType;

	GraphType *g = new GraphType( (n_pos+n_neg), (n_pos+n_neg)*(n_pos+n_neg-1)/2);

	for (i = 0; i < (n_pos+n_neg); i++){
		g -> add_node();
	}
	for (i = 0; i < (n_pos+n_neg); i++){
		g -> add_tweights( i, unary_pos[i]*scaleFactor, unary_neg[i]*scaleFactor);
			for (j = (i+1); j < (n_pos+n_neg); j++){
		          g -> add_edge( i, j,  binary[i][j]*scaleFactor, binary[i][j]*scaleFactor );
	       }
	}
	int flow = g -> maxflow();
	labels = (int *) malloc((n_pos+n_neg)*sizeof(int));
	for (i = 0; i < (n_pos+n_neg); i++){
	  if (g->what_segment(i) == GraphType::SOURCE){
	      labels[i] = -1;
	  }
	  else{
	      labels[i] = 1;
	  }
	}
	return labels;
}