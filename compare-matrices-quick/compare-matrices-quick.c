/*  ================================================================
compmat.c: 25/09/2012 Sebastien Jaeger
 
compmat.c is a simplified C translation of compare_matrices developed by Jacques Van Helden.
	
Compilation: gcc compmat_quick_2modes.c -o compmat -lm
Exucution example: ./compmat -o result.tab -file1 DemoCompMat.txt -file2 JasparCore_Transfac.txt -lth_ncor 0.7 -lth_w 5

Synopsis: compmat -o <ouput_file> -file1 <ref_matrices_file> -file2 <query_matrices_file> \
			[-lth_w <min_aligned_columns>]	min treshold on matrices overlap (default 5)
			[-lth_ncor <min_ncor>]			min threshold on normalized correlation (default 0.4)
			[-lth_ncor1 <min_ncor1>]		min threshold on correlation normalized on ref. matrices (default 0.)
			[-lth_ncor2 <min_ncor2>]		min threshold on correlation normalized on query matrices (default 0.)

=================================================================== */
#include <stdio.h> 
#include <stdlib.h>  
#include <math.h>  
#include <string.h>
#include <time.h>
#include <sys/types.h>
#include <regex.h>

#define min(x,y) x<y ? x:y
#define max(x,y) x>y ? x:y
#define flipflap(x) (x+1)%2


//==================================================================
//=============== GLOBAL VARIABLES AND STRUCTURES ==================
//==================================================================

// --------------- User-specified parameters: ---------------------
int lth_w = 5;				// lower trheshold on the number of aligned columns 
double lth_cor = 0.7;		// lower threshold on correlation
double lth_ncor = 0.4;		// lower threshold on normalized correlation
double lth_ncor1 = 0.;		// lower threshold normalized on reference matrices
double lth_ncor2 = 0.;		// lower threshold normalized on query matrices
char *Qfile;				// Query matrices file 
char *Rfile;				// Reference conserved matrices file 
char *outfile="compmat_out.tab";	// Output file
int verbose = 0;
char detect_palindromes = 0;
char *mode="scan";

//-----------------------------------------------------------------

//static FILE *pf,*pf2;


typedef struct				// Matrix structure
{	
	char ID[128];
	char name[128];
	int width;
	int nrow;
	float **mat;
} pssm;

typedef struct				// Match coordinates structure
{	
	char ID[128];
	int offset;
} mat_coods;


typedef struct				// Correlations structure
{	
	double cor;
	double Ncor;
	double Ncor1;
	double Ncor2;
    int w;
    int W;
    double wr;
    double wr1;
    double wr2;
} correls;

typedef struct				// alignment structure (corresponds to output line format)
{	
	char *id1;
	char *id2;
	char *name1;
	char *name2;
	double cor;
	double Ncor;
	double Ncor1;
	double Ncor2;
	int w1;
	int w2;
    int w;
    int W;
    double wr;
    double wr1;
    double wr2;
	int offset;
	char strand;
	int fake_matches_count;
} match;




//==================================================================
//============================ PROTOTYPES ==========================
//==================================================================
static void read_arg(int argc, char *argv[]);		// reads main arguments
pssm *readmat(FILE *fp,int *matnum);				// loads ref an query matrices files
pssm reverse_matrix(pssm matrix);					// computes reverse-complement matrix
correls calc_corr(int offset, pssm M1, pssm M2);	// computes all correlations
double crude_freq(pssm m,int r,int c);				// computes crude frequences of one cell in matrix m
void print_help(void);								// print the help


//==================================================================
//============================== MAIN ==============================
//==================================================================
int main(int argc, char *argv[]){
	int i,j,k,Rnum=0,Qnum=0; //l
	FILE *fp;
	//char currline[64];
	//char test[20];
	pssm *Rmatab=NULL;
	pssm *Qmatab=NULL,*Qrevtab=NULL;
	long *last_match=NULL;
	correls **cor_tab[2];
	short int best_correl;
	//char ID[30];
	int min_offset,max_offset;
	match *res_tab=NULL;
	long res_count=0;
	float exec_time;
	clock_t t1,t2;
	
	t1 = clock();
	
	read_arg(argc, argv);							// Set user-defined parameters
		
	fp = fopen(Rfile,"r");							// loading of the input files
	Rmatab=readmat(fp,&Rnum);
	if (verbose > 0) {  printf("-> %d reference matrices retrieved from '%s'\n",Rnum,Rfile); }
	fclose(fp);
	fp = fopen(Qfile,"r");
	Qmatab=readmat(fp,&Qnum);
	if (verbose > 0) { printf("-> %d query matrices retrieved from '%s'\n",Qnum,Qfile); }
	fclose(fp);
	last_match=(long *)malloc(Qnum*sizeof(long));
	for (i=0; i<Qnum; i++) {
		last_match[i]=-1;
	}
	
	if (verbose > 0) {
		printf("Reference matrices file:\n");
		for (i=0; i<Rnum; i++) {
			printf("\tID: %s\tname: %s\t%dbp\n",Rmatab[i].ID,Rmatab[i].name,Rmatab[i].width);
		}
		printf("Query matrices file:\n");
		for (i=0; i<Qnum; i++) {
			printf("\tID: %s\tname: %s\t%dbp\n",Qmatab[i].ID,Qmatab[i].name,Qmatab[i].width);
		}
	}
	
	for (i=0; i<2; i++) {							// allocation of correlations tab
		cor_tab[i]=(correls **)malloc(Rnum*sizeof(correls *));
		for (j=0; j<Rnum; j++) {
			cor_tab[i][j]=(correls *)malloc(Qnum*sizeof(correls));
		}
	}
	
	Qrevtab=(pssm *)malloc(Qnum*sizeof(pssm));		// Computes reverse-complement query matrices
	for (i=0; i<Qnum; i++) {
		Qrevtab[i]=reverse_matrix(Qmatab[i]);
	}
	
		
	for (i=0; i<Rnum; i++) {						// Computes all correlations
		for (j=0; j<Qnum; j++) {
			min_offset = lth_w - Qmatab[j].width;
			max_offset = Rmatab[i].width - lth_w;
			for (k=min_offset; k<=max_offset; k++) {
				cor_tab[0][i][j] = calc_corr(k,Rmatab[i],Qmatab[j]);
				cor_tab[1][i][j] = calc_corr(k,Rmatab[i],Qrevtab[j]);
				(cor_tab[0][i][j].cor >= cor_tab[1][i][j].cor) ? (best_correl = 0) : (best_correl = 1);
				if ((cor_tab[best_correl][i][j].cor >= lth_cor) && (cor_tab[best_correl][i][j].Ncor >= lth_ncor) && (cor_tab[best_correl][i][j].Ncor1 >= lth_ncor1) && (cor_tab[best_correl][i][j].Ncor2 >= lth_ncor2)) {	// best strand cases
                    if (strcmp(mode,"scan") == 0) {
                        if ((last_match[j] >= 0) && (Rmatab[i].ID == res_tab[last_match[j]].id1) && ((k-res_tab[last_match[j]].offset)<res_tab[last_match[j]].w2)) {
                            if (cor_tab[best_correl][i][j].Ncor > res_tab[last_match[j]].Ncor) {
                                res_tab[last_match[j]].cor=cor_tab[best_correl][i][j].cor;
                                res_tab[last_match[j]].Ncor=cor_tab[best_correl][i][j].Ncor;
                                res_tab[last_match[j]].Ncor1=cor_tab[best_correl][i][j].Ncor1;
                                res_tab[last_match[j]].Ncor2=cor_tab[best_correl][i][j].Ncor2;
                                //res_tab[last_match[j]].w=cor_tab[best_correl][i][j].w;
                                res_tab[last_match[j]].offset=k;
                                (best_correl == 0) ? (res_tab[last_match[j]].strand='D') : (res_tab[last_match[j]].strand='R');
                            }
                            if ((cor_tab[flipflap(best_correl)][i][j].cor >= lth_cor) && (cor_tab[flipflap(best_correl)][i][j].Ncor >= lth_ncor) && (cor_tab[flipflap(best_correl)][i][j].Ncor1 >= lth_ncor1) && (cor_tab[flipflap(best_correl)][i][j].Ncor2 >= lth_ncor2)) {	// if other strand also matches...
                                res_tab[last_match[j]].fake_matches_count++;
                            }
                            res_tab[last_match[j]].fake_matches_count++;
                        }
                        else {
                            res_tab=(match *)realloc(res_tab,(res_count+1)*sizeof(match));
                            res_tab[res_count].fake_matches_count=0;
                            res_tab[res_count].id1=Rmatab[i].ID;
                            res_tab[res_count].id2=Qmatab[j].ID;
                            res_tab[res_count].name1=Rmatab[i].name;
                            res_tab[res_count].name2=Qmatab[j].name;
                            res_tab[res_count].cor=cor_tab[best_correl][i][j].cor;
                            res_tab[res_count].Ncor=cor_tab[best_correl][i][j].Ncor;
                            res_tab[res_count].Ncor1=cor_tab[best_correl][i][j].Ncor1;
                            res_tab[res_count].Ncor2=cor_tab[best_correl][i][j].Ncor2;
                            res_tab[res_count].w1=Rmatab[i].width;
                            res_tab[res_count].w2=Qmatab[j].width;
                            //res_tab[res_count].w=cor_tab[best_correl][i][j].w;
                            res_tab[res_count].offset=k;
                            (best_correl == 0) ? (res_tab[res_count].strand='D') : (res_tab[res_count].strand='R');
                            last_match[j]=res_count;
                            if ((cor_tab[flipflap(best_correl)][i][j].cor >= lth_cor) && (cor_tab[flipflap(best_correl)][i][j].Ncor >= lth_ncor) && (cor_tab[flipflap(best_correl)][i][j].Ncor1 >= lth_ncor1) && (cor_tab[flipflap(best_correl)][i][j].Ncor2 >=  lth_ncor2)) {	// if other strand also matches...
							res_tab[res_count].fake_matches_count++;
                            }
                            res_count++;
                        }
                    }
                    else if (strcmp(mode,"matches") == 0) { // mode modified by Morgane in August 2015
                        if ((last_match[j] >= 0) && (Rmatab[i].ID == res_tab[last_match[j]].id1)) {
                            if ((cor_tab[best_correl][i][j].Ncor > res_tab[last_match[j]].Ncor) || (cor_tab[best_correl][i][j].Ncor == res_tab[last_match[j]].Ncor && cor_tab[best_correl][i][j].w > res_tab[last_match[j]].w)) {
                                res_tab[last_match[j]].cor=cor_tab[best_correl][i][j].cor;
                                res_tab[last_match[j]].Ncor=cor_tab[best_correl][i][j].Ncor;
                                res_tab[last_match[j]].Ncor1=cor_tab[best_correl][i][j].Ncor1;
                                res_tab[last_match[j]].Ncor2=cor_tab[best_correl][i][j].Ncor2;
                                res_tab[last_match[j]].w=cor_tab[best_correl][i][j].w;
                                res_tab[last_match[j]].W=cor_tab[best_correl][i][j].W;
                            	res_tab[last_match[j]].wr=cor_tab[best_correl][i][j].wr;
                            	res_tab[last_match[j]].wr1=cor_tab[best_correl][i][j].wr1;
                            	res_tab[last_match[j]].wr2=cor_tab[best_correl][i][j].wr2;
                                res_tab[last_match[j]].offset=k;
                                (best_correl == 0) ? (res_tab[last_match[j]].strand='D') : (res_tab[last_match[j]].strand='R');
                            }
                            if ((cor_tab[flipflap(best_correl)][i][j].cor >= lth_cor) && (cor_tab[flipflap(best_correl)][i][j].Ncor >= lth_ncor) && (cor_tab[flipflap(best_correl)][i][j].Ncor1 >= lth_ncor1) && (cor_tab[flipflap(best_correl)][i][j].Ncor2 >= lth_ncor2)) {	// if other strand also matches...
                                res_tab[last_match[j]].fake_matches_count++;
                            }
                            res_tab[last_match[j]].fake_matches_count++;
                        }
                        else {
                            res_tab=(match *)realloc(res_tab,(res_count+1)*sizeof(match));
                            res_tab[res_count].fake_matches_count=0;
                            res_tab[res_count].id1=Rmatab[i].ID;
                            res_tab[res_count].id2=Qmatab[j].ID;
                            res_tab[res_count].name1=Rmatab[i].name;
                            res_tab[res_count].name2=Qmatab[j].name;
                            res_tab[res_count].cor=cor_tab[best_correl][i][j].cor;
                            res_tab[res_count].Ncor=cor_tab[best_correl][i][j].Ncor;
                            res_tab[res_count].Ncor1=cor_tab[best_correl][i][j].Ncor1;
                            res_tab[res_count].Ncor2=cor_tab[best_correl][i][j].Ncor2;
                            res_tab[res_count].w1=Rmatab[i].width;
                            res_tab[res_count].w2=Qmatab[j].width;
                            res_tab[res_count].w=cor_tab[best_correl][i][j].w;
                            res_tab[res_count].W=cor_tab[best_correl][i][j].W;
                            res_tab[res_count].wr=cor_tab[best_correl][i][j].wr;
                            res_tab[res_count].wr1=cor_tab[best_correl][i][j].wr1;
                            res_tab[res_count].wr2=cor_tab[best_correl][i][j].wr2;
                            res_tab[res_count].offset=k;
                            (best_correl == 0) ? (res_tab[res_count].strand='D') : (res_tab[res_count].strand='R');
                            last_match[j]=res_count;
                            if ((cor_tab[flipflap(best_correl)][i][j].cor >= lth_cor) && (cor_tab[flipflap(best_correl)][i][j].Ncor >= lth_ncor) && (cor_tab[flipflap(best_correl)][i][j].Ncor1 >= lth_ncor1) && (cor_tab[flipflap(best_correl)][i][j].Ncor2 >= lth_ncor2)) {	// if other strand also matches...
                                res_tab[res_count].fake_matches_count++;
                            }
                            res_count++;

                        }
                    }
                    else {printf("No valide mode set, analysis failed\n"); exit(0);}
                }
			}
		}
	}
	
	fp = fopen(outfile,"w");						// output printing... 
	fprintf(fp,";mode: %s\tthresholds:\tcor=%f\tncor=%f\tw=%d\tncor1=%f\tncor2=%f\n",mode,lth_cor,lth_ncor,lth_w,lth_ncor1,lth_ncor2);
	fprintf(fp,"#id1\tid2\tname1\tname2\tcor\tNcor\tNcor1\tNcor2\tw1\tw2\tw\tW\tWr\twr1\twr2\tstrand\toffset\tuncounted\n");
	for (i=0; i<res_count; i++) {
		fprintf(fp,"%s\t%s\t%s\t%s\t%lf\t%lf\t%lf\t%lf\t%d\t%d\t%d\t%d\t%lf\t%lf\t%lf\t%c\t%d\t%d\n",
		res_tab[i].id1,res_tab[i].id2,res_tab[i].name1,res_tab[i].name2,
		res_tab[i].cor,res_tab[i].Ncor,res_tab[i].Ncor1,res_tab[i].Ncor2,
		res_tab[i].w1,res_tab[i].w2,res_tab[i].w,
		res_tab[i].W,res_tab[i].wr,res_tab[i].wr1,res_tab[i].wr2,
		res_tab[i].strand,res_tab[i].offset,res_tab[i].fake_matches_count);
	}
	
	free(Rmatab);									// memory desallocations...
	free(Qmatab);
	free(Qrevtab);
	free(cor_tab[0]);
	free(cor_tab[1]);
	free(res_tab);
	
	if (verbose > 0) { printf("%ld matches found -> '%s'\n",res_count,outfile); }
	
	t2 = clock();
	exec_time = (float)(t2-t1)/CLOCKS_PER_SEC;
	fprintf(fp,";Analysis performed in %fs\n",exec_time);
	if (verbose > 0) { printf("Analysis performed in %fs\n",exec_time); }
	
	fclose(fp);
	
	return 0;
}

//==================================================================
//================ Read user specified arguments ===================
//==================================================================
void read_arg(int argc, char *argv[]){
  int  i;
  if (argc>=3) {
	  for (i = 1; i < argc; i++){
		  if (argv[i][0] == '-'){
			  if(strcmp(argv[i],"-o") == 0)  outfile  = argv[++i];
			  if(strcmp(argv[i],"-v") == 0)  verbose=0;
			  if(strcmp(argv[i],"-h") == 0)  {print_help(); exit(0);}
			  if(strcmp(argv[i],"-file1") == 0)  Rfile  = argv[++i];
			  if(strcmp(argv[i],"-file2") == 0)  Qfile  = argv[++i];
			  if(strcmp(argv[i],"-lth_w") == 0)  lth_w  = atoi(argv[++i]);
			  if(strcmp(argv[i],"-lth_cor") == 0)  lth_cor  = atof(argv[++i]);
			  if(strcmp(argv[i],"-lth_ncor") == 0)  lth_ncor  = atof(argv[++i]);
			  if(strcmp(argv[i],"-lth_ncor1") == 0)  lth_ncor1  = atof(argv[++i]);
			  if(strcmp(argv[i],"-lth_ncor2") == 0)  lth_ncor2  = atof(argv[++i]);
			  if(strcmp(argv[i],"-detect_palindromes") == 0)  detect_palindromes  = 1;
              if(strcmp(argv[i],"-mode") == 0)  mode  = argv[++i];
		  }
	  }
  }
  else {
	print_help();
	exit(0);
  }
}

//==================================================================
//====================== Read matrices files =======================
//==================================================================
pssm *readmat(FILE *fp,int *matnum){
	char currline[128];
	int width; //i, count=0;
	float a,c,g,t;
	pssm *matab;
	regex_t pregAC,pregID;
	//regmatch_t *pmatch = NULL;
	
	regcomp(&pregAC,"^AC[[:space:]]+[^[:space:]]+",REG_EXTENDED| REG_NOSUB);
	regcomp(&pregID,"^ID[[:space:]]+[^[:space:]]+",REG_EXTENDED| REG_NOSUB);
	
	*matnum=0;
	matab=NULL;
	while (fgets(currline,128,fp) != NULL) {
		width=0;
		if (regexec(&pregAC,currline,0,NULL,0)==0) {
			matab=(pssm *)realloc(matab,(++*matnum)*sizeof(pssm));
			matab[*matnum-1].mat=NULL;
			sscanf(currline,"AC %s",matab[*matnum-1].ID);
			fgets(currline,128,fp);
			while (strstr(currline,"//")==NULL) {
				if (regexec(&pregID,currline,0,NULL,0)==0) {
					sscanf(currline,"ID %s",matab[*matnum-1].name);
				}
				else if (sscanf(currline,"%*d %f %f %f %f",&a,&c,&g,&t)==4){
					width++;
					matab[*matnum-1].mat=(float **)realloc(matab[*matnum-1].mat, (width)*sizeof(float *));
					matab[*matnum-1].mat[width-1]=(float *)malloc(4*sizeof(float));
					matab[*matnum-1].mat[width-1][0]=a;
					matab[*matnum-1].mat[width-1][1]=c;
					matab[*matnum-1].mat[width-1][2]=g;
					matab[*matnum-1].mat[width-1][3]=t;
				}
				fgets(currline,128,fp);
			}
			matab[*matnum-1].width=width;
		}
	}	
	return matab;
}

//==================================================================
//============= Compute reverse complement matrix ==================
//==================================================================
pssm reverse_matrix(pssm matrix){
	int i,j;
	pssm rev_matrix;
	
	strcpy(rev_matrix.ID,matrix.ID);
	strcpy(rev_matrix.name,matrix.name);
	//rev_matrix.ID=matrix.ID;
	//rev_matrix.name=matrix.name;
	rev_matrix.width=matrix.width;

	rev_matrix.mat=(float **)malloc(matrix.width*sizeof(float *));
	for (i=0; i<matrix.width; i++) {
		rev_matrix.mat[i]=(float *)malloc(4*sizeof(float));
	}
	
	
	for (i=0; i<matrix.width; i++) {
		for (j=0; j<4; j++) {
			rev_matrix.mat[i][j]=matrix.mat[matrix.width-i-1][3-j];
		}
	}
	return rev_matrix;
}

//==========================================================================
//= Compute the normalized correlation between two matrices with an offset =
//==========================================================================
correls calc_corr(int offset, pssm M1, pssm M2) {
	double v1 = 0;			// Variance of aligned columns in matrix 1 
	double v2 = 0;			// Variance of aligned columns in matrix 2 
	double cov = 0;			// Covariance of aligned columns 
	//double cor = 0;			// Coefficient of correlation for aligned columns 
	double sum_f1 = 0;		// Sum of residue frequencies for aligned columns of matrix 1 
	double sum_f2 = 0;		// Sum of residue frequencies for aligned columns of matrix 2 
	double sum_sq_f1 = 0;	// sum of squared residue frequencies for aligned columns of matrix 1 (used to compute v1)  
	double sum_sq_f2 = 0;	// sum of squared residue frequencies for aligned columns of matrix 2 (used to compute v2) 
	double sum_f1f2 = 0;	// Sum of residue frequencies between the matrices (used to compute cov) 
	int end2,start2,w;
	int pos;				// Position of the current column relative to the number of aligned columns 
	//int col1,col2;			// Current column in M1 (query) and M2 (ref.) respectively 
	int r; //i,j,k
	int n;					// Total number of cells in the aligned matrix
	double f1,f2;			// frenquencies
	double norm_factor = 1;
	pssm m1,m2;				// aligned submatrices of M1 and M2, sizes should be w1=w2=w  
	correls Cor;			// Output structure containing cor, Ncor, Ncor1, Ncor2
	
	
	start2 = max(0, -offset);					// Compute w = the number of aligned columns 
	end2 = min(M2.width, M1.width-offset);
	w = end2-start2;
	

	m1.mat = M1.mat + (max(0,offset));			// !! not addition but pointer shift
	m1.width = w;
	m2.mat = M2.mat + start2; 
	m2.width = w;
	
    if (strcmp(mode,"scan")==0) { // scan mode :  strcmp function: if Return value = 0 then it indicates str1 is equal to str2
        norm_factor = (double)w / (M1.width+M2.width-w);
    }
    else if (strcmp(mode,"matches")==0) { // matches mode 
         norm_factor = (double)w / (M1.width+M2.width-w);
    }
	
	for (pos=0; pos<w; pos++) {					// Compute correlations for each pos 
		for (r=0;r<4;r++) {							// Iterate over residues (matrix rows) 
			f1 = crude_freq(m1,r,pos);				// Extract the frequency of residue r in the first matrix
			f2 = crude_freq(m2,r,pos);				// Extract the frequency of residue r in the second matrix
												// Update sums for correlations
			sum_f1 += f1;							// Update sum of residues for matrix 1
			sum_f2 += f2;							// Update sum of residues for matrix 2
			sum_sq_f1 += f1*f1;
			sum_sq_f2 += f2*f2;
			sum_f1f2 += f1*f2;
			}
	}	
	
												// Compute covariance and coefficient of correlation
	n = 4*w;										// Number of matrix cells in the alignment
	v1 = sum_sq_f1/n - (sum_f1/n)*(sum_f1/n);		// Variance of the aligned columns of first matrix
	v2 = sum_sq_f2/n - (sum_f2/n)*(sum_f2/n);		// Variance of the aligned columns of second matrix
	cov = sum_f1f2/n - sum_f1*sum_f2/(n*n);			// Covariance of aligned columns between first and second matrices
	if ((cov == 0) || (v1*v2==0)) {
			Cor.cor = 0;
	} 
	else {
		Cor.cor = cov/sqrt(v1*v2);
	}
	Cor.Ncor = Cor.cor * norm_factor; 
	Cor.Ncor1 = Cor.cor * (double)w / M1.width;
	Cor.Ncor2 = Cor.cor * (double)w / M2.width;
    Cor.w = w;
    Cor.W = M1.width+M2.width-w;
    Cor.wr = norm_factor; 
    Cor.wr1 = (double)w / M1.width ; // in Perl:  my $wr1 = $w / $w1; $Ncor1 = $cor * $w / $w1 ;
	Cor.wr2 = (double)w / M2.width ;		
	
	return Cor;
}

//=======================================================================
//== Return the crude fréquency of residue r in column c of matrix mat ==
//=======================================================================
double crude_freq(pssm m,int r,int c) {
	int i;
	double sum_col=0;
	
	for (i=0; i<4; i++) {
		sum_col += m.mat[c][i];
	}
	if (sum_col == 0) return 0.;
	else return (m.mat[c][r]/sum_col);
}
//=======================================================================
//====================== Print the help/manual  =========================
//=======================================================================
void print_help(void){
	printf("compare-matrices-quick is a simplified C translation of compare_matrices developed by Jacques Van Helden. ");
	printf("It takes as input a reference matrices file and a query matrices file (both transfac formated) and computes correlations for each possible couple of matrices.\n");
	printf("The output consists in a tab-delimited file containing one line per match (i.e. alignments with correlations scores above user-defined thresholds) following the format:\n\n");

	printf("<Query matrix ID>  <Ref. Matrix ID>  <Query name>  <Ref name>  <cor> <Ncor>  <Ncor1>  <Ncor2>  <Query matrix length>  <Ref length>  <offset>  <Strand>\n\n");

	printf("where:\t<cor>\t\tdenotes the coefficient of correlation for aligned columns,\n");
	printf("\t<Ncor>\t\tdenotes the coefficient of correlation normalized on alignment length,\n");
	printf("\t<Ncor1>\t\tdenotes the coefficient of correlation normalized on reference matrix length,\n");
	printf("\t<Ncor2>\t\tdenotes the coefficient of correlation normalized on query matrix length,\n");
	printf("\t<offset>\tis the relative position off the query matrix relatively to the reference matrix\n\t\t\t(-> zero if the first column of the query is aligned with the first column of the ref.,\n\t\t\tand can be negative if the query is shifted leftward),\n");
	printf("\t<Strand>\tcan be Direct ('D') or Reverse ('R') if the comparison considered the reverse-complement matrix of the query.\n\n");
	printf("\t<uncounted>\tindicates the number of matches that were not counted because of a beter match in a close vicinity or on the other strand\n\n");

	printf("Synopsis: compmat -o <ouput_file> -file1 <ref_matrices_file> -file2 <query_matrices_file> \\\n");
	printf("\t\t[-lth_w <min_aligned_columns>]\t\tmin treshold on matrices overlap (default 5)\n");
	printf("\t\t[-lth_cor <min_aligned_columns>]\t\tmin threshold on correlation (default 0.7)\n\n");
	printf("\t\t[-lth_ncor <min_ncor>]\t\t\tmin threshold on normalized correlation (default 0.4)\n");
	printf("\t\t[-lth_ncor1 <min_ncor1>]\t\tmin threshold on correlation normalized on ref. matrices (default 0)\n");
	printf("\t\t[-lth_ncor2 <min_ncor2>]\t\tmin threshold on correlation normalized on query matrices (default 0)\n");
	printf("\t\t[-h]\t\t\t\t\tprint this help\n");
	printf("\t\t[-v]\t\t\t\t\tVerbose mode (debuging...)\n");
    printf("\t\t[-mode <mode>]\t\t\t\t\twhere <mode> can be either \"scan\" or \"matches\" (default = scan)\n\n");
    printf("\t\"scan\" mode: reports all matching positions between matrix 2 (reference) and matrix 1 (query) that pass the thresholds on the metrics""\n\n");
    printf("\t\"matches\" mode: For each pair of matrices (one from file1 and one from file2), the program tests all possible offsets, and reports only\n\t\t the best Ncor matching position (if passing the Ncor threshold)\n\n");

	
	printf("Exucution example: ./compare-matrices-quick -o result.tab -file1 DemoCompMat.txt -file2 JasparCore_Transfac.txt -lth_ncor2 0.7 -lth_w 5 -mode matches\n");
}
	

