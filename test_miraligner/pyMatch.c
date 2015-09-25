#include <stdio.h>
#include <string.h>
#include <malloc.h>


int Match(char *strl, char *strs, int err, int err2) {

	int lenl,lens;
	int error;
	int i,j;
	int pos = -1;
	
	lenl = strlen (strl);
	lens = strlen (strs);
	
	error = -1;
	
	i = 0;
	while (error < err && (i + lens) <= lenl) {
		j = 0;
		
		while (error < err && j < lens-err2+err) {
			if (strl[i+j] != strs[j]) error++;
			j++;
		}
		
		if (error < err) {
			pos = i;
			break;
		}
		
		error = -1;
		
		i++;
	}
	
	return pos;
}



int MMatch (int len, char **ls, char *ss, int err, int err2) {
	int i;
	int j;
	
	for (i = 0; i<len; i++) {
		j = Match (ls[i], ss, err, err2);
		
		if (j != -1) break;
	}

	if (i == len) i = -1;
	
	return i;
}


int Miraligner (char *lfn, char *sfn, char *rfn, int err, int err2) {
	FILE *plf;
	FILE *psf;
	FILE *prf;

	long size;

	char *pls;
	char *pss;

	char **ppln;
	char **ppls;
	char **ppsn;
	char **ppss;

	
	int i;
	int j;
	int k;
	int l;
	
	int nbLSeqs;
	int nbSSeqs;


	////////////////////////////////////////////////////////////////////////
	if ((plf = fopen(lfn, "r")) == NULL) {
		printf("\nError opening file %s", lfn);
		return -1;
	}

	fseek(plf, 0, SEEK_END); // seek to end of file
	size = ftell(plf) + 1;   // get current file pointer
	fseek(plf, 0, SEEK_SET); // seek back to beginning of file

	pls = (char *) malloc(size*sizeof(char));

	i = 0;
	j = 0;
	while (!feof(plf) && i < size) {
		pls[i] = fgetc(plf);
		
		if (pls[i] == '>') j++;
		
		i++;
	}
	
	fclose(plf);
	
	
	ppln = (char **) malloc(j*sizeof(char *));
	ppls = (char **) malloc(j*sizeof(char *));
	
	
	i = 0;
	j = 0;
	while (i < size) {
		if (pls[i] == '>') {
			ppln[j] = &pls[i+1];
			while (pls[i] != ' ' && pls[i] != '\n' && i < size) i++;
			if (pls[i] == '\n') {
				pls[i] = '\0';
			} else {
				pls[i] = '\0';
				while (pls[i] != '\n' && i < size) i++;
			}
			ppls[j] = &pls[i+1];
			while (pls[i] != '\n' && i < size) i++;
			pls[i] = '\0';
			
			j ++;
		}
		i++;
	}
	
	
	nbLSeqs = j;
	////////////////////////////////////////////////////////////////////////
	
	
	
	////////////////////////////////////////////////////////////////////////
	if ((psf = fopen(sfn, "r")) == NULL) {
		printf("\nError opening file %s", sfn);
		
		free(pls);	
		free(ppln);	
		free(ppls);	
		
		return -1;
	}

	fseek(psf, 0, SEEK_END); // seek to end of file
	size = ftell(psf) + 1;   // get current file pointer
	fseek(psf, 0, SEEK_SET); // seek back to beginning of file

	pss = (char *) malloc(size*sizeof(char));

	i = 0;
	j = 0;
	while (!feof(psf) && i < size) {
		pss[i] = fgetc(psf);
		
		if (pss[i] == '>') j++;
		
		i++;
	}
	
	fclose(psf);
	
	
	ppsn = (char **) malloc(j*sizeof(char *));
	ppss = (char **) malloc(j*sizeof(char *));
	
	
	i = 0;
	j = 0;
	while (i < size) {
		if (pss[i] == '>') {
			ppsn[j] = &pss[i+1];
			while (pss[i] != ' ' && pss[i] != '\n' && i < size) i++;
			if (pss[i] == '\n') {
				pss[i] = '\0';
			} else {
				pss[i] = '\0';
				while (pss[i] != '\n' && i < size) i++;
			}
			ppss[j] = &pss[i+1];
			while (pss[i] != '\n' && i < size) i++;
			pss[i] = '\0';
			
			j ++;
		}
		i++;
	}
	
	
	nbSSeqs = j;
	////////////////////////////////////////////////////////////////////////
	
	
	
	//printf("\nfsize %ld nbLSeqs %d", size, nbLSeqs);
	//printf("\nseq 100 name %s", ppln[100]);
	//printf("\nseq 100 seqs %s", ppls[100]);
	
	//printf("\nfsize %ld nbSSeqs %d", size, nbSSeqs);
	//printf("\nseq 100 name %s", ppsn[100]);
	//printf("\nseq 100 seqs %s", ppss[100]);
	
	//printf("\n\n");
	
	
	
	////////////////////////////////////////////////////////////////////////
	if ((prf = fopen(rfn, "w")) == NULL) {
		printf("\nError opening file %s", rfn);
		
		free(pls);	
		free(ppln);	
		free(ppls);	
		free(pss);	
		free(ppsn);	
		free(ppss);	
		
		return -1;
	}

	for (i = 0; i<nbSSeqs; i++) {
		l = strlen(ppss[i]);
		
		for (j = 0; j<nbLSeqs; j++) {
			if ((k = Match(ppls[j], ppss[i], err, err2)) != -1) {
				fprintf(prf ,"%s %s %d %d\n", ppsn[i], ppln[j], k, k+l-1);
			}
			//printf("\n%d", k);
		}		
	}
	
	
	fclose(prf);
	
	
	
	
	
	
	
	
	
	
	
	
	
	////////////////////////////////////////////////////////////////////////
	free(pls);	
	free(ppln);	
	free(ppls);	
	free(pss);	
	free(ppsn);	
	free(ppss);	
	
	
	

	return 1;
}



int main(int argc, char *argv) {
	//printf("\nResult: %d", Match("Hola Como Estas", "Como", 0));
	//printf("\nResult: %d", Match("Hola Como Estas", "Camo", 1));
	//printf("\nResult: %d", Match("Hola Como Estas", "Estas", 0));
	//printf("\n");
	
	//char *ls[3];
	
	//char *l0 = "TGGGATGAGGTAGTAGGTTGTATAGTTTTAGGGTCACACCCACCACTGGGAGATAACTATACAATCTACTGTCTTTCCTA";
	//char *l1 = "AGGTTGAGGTAGTAGGTTGTATAGTTTAGAATTACATCAAGGGAGATAACTGTACAGCCTCCTAGCTTTCCT";
	//char *l2 = "GGGTGAGGTAGTAGGTTGTATAGTTTGGGGCTCTGCCCTGCTATGGGATAACTATACAATCTACTGTCTTTCCT";

	//ls[0] = l0;
	//ls[1] = l1;
	//ls[2] = l2;
	
	//char *s0 = "GGGATGAGGTAGTAGGTTGTATA";
	//char *s1 = "GGGATGAGGTCGTAGGTTGTATA";
	
	//char *s3 = "GGGTGAGGTAGTAGGTTGTATAG";

	//printf("\nResult: %d", MMatch(3, ls,s0,0));
	//printf("\nResult: %d", MMatch(3, ls,s1,0));
	//printf("\nResult: %d", MMatch(3, ls,s3,0));
	//printf("\n");
	
	
	
	
	
	
	Miraligner("hairpin.hsa.fa", "sim.fasta", "results.res", 1, 3);
	
	
	
	
	
	return 1;
}