#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


double weight(char c) {
  switch (c) {

  case 'A':
    return 0.8356471476582918;
    break;
  case 'C':
    return 0.5208088354857734;
    break;
  case 'E':
    return 0.9876987431418378;
    break;
  case 'D':
    return 0.9079044671339564;
    break;
  case 'G':
    return 0.7997168496420723;
    break;
  case 'F':
    return 0.5849790194237692;
    break;
  case 'I':
    return 0.6784124413866582;
    break;
  case 'H':
    return 0.8947913996466419;
    break;
  case 'K':
    return 0.9267104557513497;
    break;
  case 'M':
    return 0.6296623675420369;
    break;
  case 'L':
    return 0.6554221515081433;
    break;
  case 'N':
    return 0.8597433107431216;
    break;
  case 'Q':
    return 0.789434648348208;
    break;
  case 'P':
    return 0.8235328714705341;
    break;
  case 'S':
    return 0.7440908318492778;
    break;
  case 'R':
    return 0.7712466317693457;
    break;
  case 'T':
    return 0.8096922697856334;
    break;
  case 'W':
    return 0.6374678690957594;
    break;
  case 'V':
    return 0.7357837119163659;
    break;
  case 'Y':
    return 0.6112801822947587;
    break;
  default:
    return 0;
    break; //unknown residues
  }

}

double prob(double x) {
  float a = 81.0581;
  float b = -62.7775;
  double p = 1 / (1 + exp(-(a * x + b)));
  return p;

}

int main(int argc, char * argv[]) {
  FILE * fp;
  char * line = NULL;
  size_t len = 0;
  ssize_t read;

  if ((argc > 2) || (argc == 1)) {
    printf("Usage: ./swi <protein_sequence.fa>\n");
    exit(EXIT_FAILURE);
  }

  fp = fopen(argv[1], "r");
  if (fp == NULL)
    exit(EXIT_FAILURE);

  printf("\nAccession\tSWI\t\tProb.(Solubility)\n");
  printf("=========\t==========\t=================\n");

  int state = 0;
  while ((read = getline( & line, & len, fp)) != -1) {
    if (line[read - 1] == '\n')
      line[read - 1] = 0;
    if (line[0] == '>') {
      if (state == 1)
        printf("\n");
      printf("%s\t", line + 1);
      state = 1;
    } else {
      char unknown_aa_msg[35] = {
        "\0"
      }; //len of message
      int length = (int) strlen(line);
      int total_known_aa = 0;
      double sum = 0.0;
      for (int i = 0; i < length; i++) {
        float curr_wt = weight(line[i]);
        sum += curr_wt;
        total_known_aa += 1;
        if (curr_wt == 0) {
          strcpy(unknown_aa_msg, " Ignoring non-standard residue(s).");
          total_known_aa -= 1; //Don't count unknown aa
        }
      }
      double swi = sum / total_known_aa;
      printf("%.8f\t%.8f%s", swi, prob(swi), unknown_aa_msg);
    }
  }
  printf("\n\n");

  fclose(fp);
  if (line)
    free(line);
  exit(EXIT_SUCCESS);
}
