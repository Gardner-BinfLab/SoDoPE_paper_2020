#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>



double weight(char c){
	switch(c){

case 'A' : return 0.8356956599678218; break;
case 'C' : return 0.5219207324456876; break;
case 'E' : return 0.9868660417547442; break;
case 'D' : return 0.9075983546378998; break;
case 'G' : return 0.8003827946673535; break;
case 'F' : return 0.5821934635876957; break;
case 'I' : return 0.6790449304566072; break;
case 'H' : return 0.8963977585570367; break;
case 'K' : return 0.9259165090012061; break;
case 'M' : return 0.6299964100098959; break;
case 'L' : return 0.6546922237065839; break;
case 'N' : return 0.8604957042204235; break;
case 'Q' : return 0.7895650031998229; break;
case 'P' : return 0.822104415564934; break;
case 'S' : return 0.7442464390120463; break;
case 'R' : return 0.771055152304471; break;
case 'T' : return 0.8098670971949234; break;
case 'W' : return 0.6386931894494416; break;
case 'V' : return 0.7344952876686051; break;
case 'Y' : return 0.6125581495225544; break;
default :  return 0; break; //unknown residues
}

}

double prob( float x){
	double a = 81.1496;
	double b = -62.8379;
	double p = 1/(1+ exp(-(a*x + b)));
	return p;

}

int main (int argc, char *argv[]) {
	FILE * fp;
	char * line = NULL;
	size_t len = 0;
	ssize_t read;

	if ((argc > 2) || (argc == 1)){
		printf("Usage: ./swi <sequence.fa>\n");
		exit(EXIT_FAILURE);
	}

	fp = fopen(argv[1], "r");
	if (fp == NULL)
		exit(EXIT_FAILURE);


	printf("\nAccession\tSWI\t\tProb.(Solubility)\n");
	printf("=========\t==========\t=================\n");

	int state = 0;
	while ((read = getline(&line, &len, fp)) != -1) {
		if (line[read - 1] == '\n')
			line[read - 1] = 0;
		if (line[0] == '>') {
			if (state == 1)
				printf("\n");
			printf("%s\t", line+1);
			state = 1;
		} else {
			char unknown_aa_msg[35] = {"\0"}; //len of message
			int length = (int)strlen(line);
			int total_known_aa = 0;
			double sum = 0.0;
			for (int i = 0; i < length; i++){
				float curr_wt = weight(line[i]);
				sum += curr_wt ;
				total_known_aa += 1;
				if (curr_wt == 0){
					strcpy(unknown_aa_msg, " Ignoring non-standard residue(s).");
					total_known_aa -= 1; //Don't count unknown aa
				}
			}
			double swi = sum/total_known_aa;
			printf("%.8f\t%.8f%s", swi, prob(swi), unknown_aa_msg);
		}
	}
	printf("\n\n");

	fclose(fp);
	if (line)
		free(line);
	exit(EXIT_SUCCESS);
}

