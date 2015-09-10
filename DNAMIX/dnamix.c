#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<string.h>

/**************************************************************/
/*                                                            */
/* This program will perform calculations for the methods and */
/* formulas presented in:                                     */
/* Weir BS, Triggs CM, Starling L, Stowell LI, Walsh KAJ,     */
/* Buckleton J. Interpreting DNA Mixtures. J Forensic Sci     */
/* 1997; 42(2):213-222.                                       */
/*                                                            */
/* Questions and/or comments ahould be sent to:               */
/* storey@statgen.ncsu.edu                                    */
/*                                                            */
/**************************************************************/

#define MAX 12
#define MAX2 924

int max(int a, int b)
{
	if(a>b)
		return a;
	return b;
}

int min(int a, int b)
{
	if(a<b)
		return a;
	return b;
}

void banner()
/**********************************************************************
* Displays a banner giving a reference to the paper from which the 
* formulas come
*********************************************************************/
{
  printf("\n******************************************************");
  printf("\n* This program performs calculations for the methods *");
  printf("\n* and formulas presented in:                         *");
  printf("\n*                                                    *");  
  printf("\n* Weir BS, Triggs CM, Starling L, Stowell LI,        *");
  printf("\n* Walsh KAJ, Buckleton J. Interpreting DNA           *");
  printf("\n* Mixtures. J Forensic Sci 1997; 42(2):213-222.      *");
  printf("\n******************************************************");
}


void sort(double freq[][MAX], int *us, int unk_al, int nbas)
/*********************************************************************
* Places the alleles whose sources are unknown in the front of the 
* array of allele frequencies
*********************************************************************/
{
	int q, tmp;
	int i,j,k;
	double temp;

	for(i=0;i<unk_al;i++)
		for(j=0;j<unk_al-1;j++)
			if(us[j]>us[j+1])
			{
				tmp=us[j];
				us[j]=us[j+1];
				us[j+1]=tmp;
			}
	
	for(i=0;i<nbas;i++)
		for(j=0;j<unk_al;j++)
		{
			q=us[j]-1;
			temp=freq[q][i];
			for(k=us[j]-2;k>j-1;k--)
				freq[k+1][i]=freq[k][i];
			freq[j][i]=temp;
		}
}

void comb(int n, int k, int **cm, int *c)
/*********************************************************************
* Computes all possible combinations of unknown alleles to be
* included in the various partial sums
*********************************************************************/
{
	int arr[MAX];
	int m, j, l;

	*c=1;

	for(m=0;m<MAX;m++)
		arr[m]=m+1;

	while(1)
	{
		for(j=0;j<k;j++)
			cm[*c-1][j]=arr[j];
  
		(*c)++;

		for(j=k;j>=1;j--)
		{
			arr[j-1]++;
			if(arr[j-1]<=j+n-k)
				break;
		}

		if(arr[0]>n-k+1)
			break;

		  
		for(l=j;l<=k;l++)
			arr[l]=arr[l-1]+1;
		  
	}
}

void calc(double f[][MAX], int allele, int unk_al, int x, double *s, int nbas)
/*********************************************************************
* Performs the calculations for the "p^2" method
*********************************************************************/
{
	double sum, ext, par;
	int **cm, c, n;
	int i,j,k,l,m;

	cm=(int **)malloc(MAX2*sizeof(int *));
	for(i=0;i<MAX2;i++)
		cm[i]=(int *)malloc(MAX*sizeof(int));
			
	ext=0;
	sum=0;
	par=0;

	for(i=unk_al;i<allele;i++)
		ext+=f[i][nbas-1];
      
	for(j=0;j<allele;j++)
		par+=f[j][nbas-1];
      
	sum+=pow(par,2*x);
	n=1;

	for(k=unk_al-1;k>=1;k--)
	{
		comb(unk_al, k, cm, &c);
		for(l=0;l<c-1;l++)
		{
			par=0;
			for(m=0;m<k;m++)
				par+=f[cm[l][m]-1][nbas-1];
 
			sum+=pow(-1,n)*pow(par+ext,2*x);
 		}
		n++;
 	}

	if(unk_al!=0)
		sum+=pow(-1,n)*pow(ext,2*x);

	for(i=0;i<MAX2;i++)
		free(cm[i]);
	free(cm);

	*s=sum;
}

void calc2p(double f[][MAX], int allele, int unk_al, int x, double *sm, int nbas)
/*********************************************************************
* Performs the calculations for the "2p" method
*********************************************************************/
{
	double sum, par, g[MAX+1];
	int **cm, c, n;
	int i,k,l,m,p,q,r,s,t,u,v;

	cm=(int **)malloc(MAX2*sizeof(int *));
	for(i=0;i<MAX2;i++)
		cm[i]=(int *)malloc(MAX*sizeof(int));
	
	par=0;
        sum=0;
	n=0;
	
	for(k=unk_al;k>=1;k--)
	{
		comb(unk_al, k, cm, &c);
		for(l=0;l<c-1;l++)
		{
			par=0;
			for(m=0;m<k;m++)
			  g[m]=f[cm[l][m]-1][nbas-1];

			for(p=k;p<k+allele-unk_al;p++)
				g[p]=f[p+unk_al-k][nbas-1];

			for(q=0;q<k+allele-unk_al;q++)
				par+=g[q];

			for(r=0;r<k+allele-unk_al-1;r++)
				for(s=r+1;s<k+allele-unk_al;s++)
					par+=g[r]*g[s];
					
            sum+=pow(-1,n)*pow(par,x);
 		}
		n++;
 	}

	par=0;
	for(t=unk_al;t<allele;t++)
		par+=f[t][nbas-1];
      
	for(u=unk_al;u<allele-1;u++)
		for(v=u+1;v<allele;v++)
			par+=f[u][nbas-1]*f[v][nbas-1];

	sum+=pow(-1,n)*pow(par,x);
	sum*=pow(2,x);

	for(i=0;i<MAX2;i++)
		free(cm[i]);
	free(cm);

	*sm=sum;
}

void copy(char temp[30], char a[9])
{
  int i=0;
  for( ; i<8; i++)
    a[i] = temp[i];
  a[8] = '\0';

  i = 0;
  for( ; i<30 ; i++)
    temp[i] = '\0';
}

main()
{
        int allele, form, unk_al, cont_l, cont_u, ans;

        int us[MAX], smry, nbas, cont_l2, cont_u2, us1[MAX];
	int unk_al2;

	double freq[MAX][MAX], fq[MAX][MAX], lkh_nm[MAX][MAX+1], lkh_dn[MAX][MAX+1];
	double check;

	char name[MAX][9], fil[9], base[MAX][9], locus[9];
	char paws, ch, temp[30];
	FILE *fp;

	int i,j, k, l, a;

	/*********************************************************************/
	/* Gets the data, i.e. allele names and frequencies present in the   */
	/* evidence sample and the databases to be used to calculate         */
	/* likelihoods                                                       */
	/*********************************************************************/

	banner();

        ans = 0;

        printf("\n");
	allele=-1;
	

	do {
		printf("\nEnter the number of alleles in the mixture (max. 12): ");
		scanf("%d", &allele);

		if(allele>MAX||allele<1)
			printf("Not a valid entry!\n");
	}while(allele>MAX||allele<1);

	while((ch = getchar()) != '\n');

	printf("\n");
	for(i=0;i<allele;i++) {
	printf("Enter the name of allele %d: ", i+1);
	scanf("%s", temp);
	copy(temp, name[i]);
	} 

	nbas=-1;
        while(nbas>MAX||nbas<1)
	{
		printf("\nEnter the number of databases to be used (max. 12): ");
		scanf("%d",&nbas);
		if(nbas>MAX||nbas<1)
			printf("Not a valid entry!\n");
	}

	while((ch = getchar()) != '\n');

	for(j=0;j<nbas;j++)
	{
		printf("\nEnter the name for database %d: ",j+1);
		scanf("%s", temp);
		copy(temp, base[j]);
		for(i=0;i<allele;i++)
		{
			printf("Enter the frequency of %s: ",name[i]);
			scanf("%lf",&freq[i][j]);
		}

		for(l=0;l<allele;l++)
			fq[l][j]=freq[l][j];
		check = 0;
		for(k=0;k<allele;k++)
			check+=freq[k][j];
		if((int)(1000*check)>1005)
		   printf("\nWarning: The sum of the allele frequencies may be greater than 1.\n");
	        while((ch = getchar()) != '\n');
	}


        do { 
	printf("\nDo you want to use the p^2 formulae or the 2p modification? (p^2=0, 2p=1): ");
        	scanf("%d",&form);
		if(form!=0&&form!=1)
			printf("Not a valid entry!\n");
	}while(form!=0&&form!=1);

	do {
		printf("\nDo you want a summary file? (No=0, Yes=1): ");
		scanf("%d",&smry);
		if(smry!=0&&smry!=1)
			printf("Not a valid entry!\n");
	}while(smry!=0&&smry!=1);
	   if(smry==1) {
		printf("\nEnter the name of the file: ");
		scanf("%s", temp);
		copy(temp, fil);
		fp=fopen(fil,"w");
		printf("\nEnter a title for the summary (such as a locus name): ");
		scanf("%s", temp);
		copy(temp, locus);
	   }


	/*********************************************************************/
	/* Gets the alleles whose sources are unknown and the possible number*/ 
	/* of unknown contributors for the numerator of the likelihood ratio */
	/*********************************************************************/

        do {
	printf("\n\n************************************************************************");
	printf("\n    THE FOLLOWING IS INPUT FOR THE NUMERATOR OF THE LIKELIHOOD RATIO");
	printf("\n************************************************************************");

	printf("\n");
	do {
	printf("\nEnter the number of alleles whose sources are unknown: ");
		scanf("%d",&unk_al);
		if(unk_al>allele||unk_al<0)
			printf("Not a valid entry!\n");
	}while(unk_al>allele||unk_al<0);

	for(i=0;i<allele;i++)
		us[i]=i+1;

	if(unk_al!=0&&unk_al!=allele)
	{
		printf("\nEnter the number(s) corresponding to the unknown allele(s):\n");
		for(i=0;i<allele;i++)
			printf("%2d = %s\n",i+1,name[i]);
		for(i=0;i<unk_al;i++)
		{
			do
			{
			  printf("Enter unknown allele %d: ", i+1);
			  scanf("%d",&us[i]);
			  if(us[i]>allele||us[i]<1)
			  printf("Not a valid entry!\n");
			  for(j=0;j<=i-1;j++) {
				if(us[i]==us[j])
				{
				     printf("This allele has already been entered!\n");
				     us[i]=-1;
				}
			  }
			} while(us[i]>allele||us[i]<1);
		}
		
		for(i=0;i<nbas;i++)
			for(l=0;l<allele;l++)
				freq[l][i]=fq[l][i];
     
        sort(freq, us, unk_al, nbas);
	}

 	do {
		printf("\nEnter the lower bound on the number of unknown contributors");
		printf(" (maximum of 12): ");
		scanf("%d",&cont_l);

		if(cont_l<0||unk_al>cont_l*2||cont_l>MAX)
   				printf("Not a valid entry!\n");
	}while(cont_l<0||unk_al>cont_l*2||cont_l>MAX);

	do {
		printf("\nEnter the upper bound on the number of unknown contributors");
		printf(" (maximum of 12): ");
		scanf("%d",&cont_u);
		if(cont_u<cont_l||cont_u>MAX)
   				printf("Not a valid entry!\n");
	}while(cont_u<cont_l||cont_u>MAX);

	printf("\n");
	for(j=cont_l;j<=cont_u;j++)
	{
	  printf("\nUnknown Contributors: %d\n", j);
	  printf("Database       Numerator Probability\n");
	  printf("------------------------------------\n");
	  for(i=0;i<nbas;i++)
	  {
	    if(form==0)
	      calc(freq, allele, unk_al, j, &lkh_nm[i][j], i+1);
	    else
	      calc2p(freq, allele, unk_al, j, &lkh_nm[i][j], i+1);
      
	    if((int)(1000*lkh_nm[i][j])>1000)
	      lkh_nm[i][j]=1;
     
	    printf("%s", base[i]);
	    for(k=15; k>strlen(base[i]); k--)
		  printf(" ");
	    printf("%16.6E\n", lkh_nm[i][j]);
	  }
	  printf("\n");
	}


/*********************************************************************
* Gets the alleles whose sources are unknown and the possible number 
* of unknown contributors for the denominator of the likelihood ratio
*********************************************************************/


	printf("\n************************************************************************");
	printf("\n   THE FOLLOWING IS INPUT FOR THE DENOMINATOR OF THE LIKELIHOOD RATIO");
	printf("\n************************************************************************");

	printf("\n");
	do {
	printf("\nEnter the number of alleles whose sources are unknown: ");
		scanf("%d",&unk_al2);
		if(unk_al2>allele||unk_al2<0)
			printf("Not a valid entry!\n");
	}while(unk_al2>allele||unk_al2<0);

	for(i=0;i<allele;i++)
		us1[i]=i+1;

	if(unk_al2!=0&&unk_al2!=allele)
	{
		printf("\nEnter the number(s) corresponding to the unknown allele(s):\n");
		for(i=0;i<allele;i++)
			printf("%2d = %s\n",i+1,name[i]);
		for(i=0;i<unk_al2;i++)
		{
			do
			{
			  printf("Enter unknown allele %d: ", i+1);
			  scanf("%d",&us1[i]);
			  if(us1[i]>allele||us1[i]<1)
			  printf("Not a valid entry!\n");
			  for(j=0;j<=i-1;j++) {
				if(us1[i]==us1[j])
				{
				     printf("This allele has already been entered!\n");
				     us1[i]=-1;
				}
			  }
			} while(us1[i]>allele||us1[i]<1);
		}
		
		for(i=0;i<nbas;i++)
			for(l=0;l<allele;l++)
				freq[l][i]=fq[l][i];
     
        sort(freq, us1, unk_al2, nbas);
	}

 	do {
		printf("\nEnter the lower bound on the number of unknown contributors");
		printf(" (maximum of 12): ");
		scanf("%d",&cont_l2);

		if(cont_l2<0||unk_al2>cont_l2*2||cont_l2>MAX)
   				printf("Not a valid entry!\n");
	}while(cont_l2<0||unk_al2>cont_l2*2||cont_l2>MAX);

	do {
		printf("\nEnter the upper bound on the number of unknown contributors");
		printf(" (maximum of 12): ");
		scanf("%d",&cont_u2);
		if(cont_u2<cont_l2||cont_u2>MAX)
   				printf("Not a valid entry!\n");
	}while(cont_u2<cont_l2||cont_u2>MAX);

	printf("\n");
	for(j=cont_l2;j<=cont_u2;j++)
	{
		printf("\nUnknown Contributors: %d\n", j);
		printf("Database     Denominator Probability\n");
		printf("------------------------------------\n");

		for(i=0;i<nbas;i++)
		{
			if(form==0)
				calc(freq, allele, unk_al2, j, &lkh_dn[i][j], i+1);
			else
				calc2p(freq, allele, unk_al2, j, &lkh_dn[i][j], i+1);
      
		  if((int)(1000*lkh_dn[i][j])>1000)
			  lkh_dn[i][j]=1;
      
		  printf("%s", base[i]);
		  for(k=15; k>strlen(base[i]); k--)
		    printf(" ");
		  printf("%16.6E\n", lkh_dn[i][j]);
		}
		printf("\n");
	}


/*********************************************************************
* Displays the likelihood ratios according to databases and number of 
* unknown contributors in the numerator and denominator
*********************************************************************/

	printf("\n************************************************************************");
	printf("\n       THE FOLLOWING ARE THE LIKELIHOOD RATIOS FOR EACH DATABASE");
	printf("\n************************************************************************");

	printf("\n\n(Press return after each database)");
	while((ch = getchar()) != '\n');

	printf("\n");
	printf("\n\n           Numerator Unknown    Denominator Unknown");
	printf("\nDatabase      Contributors         Contributors         Likelihood Ratio");
	printf("\n------------------------------------------------------------------------\n");
     
	for(i=0;i<nbas;i++)
	{
	  for(j=cont_l;j<=cont_u;j++)
	  for(k=cont_l2;k<=cont_u2;k++) {
	    printf("%s", base[i]);
	    for(a=12; a>strlen(base[i]); a--)
	       printf(" ");
	    printf("%8d %20d %30.2lf\n", j, k, lkh_nm[i][j]/lkh_dn[i][k]);
	  }
	    while((ch = getchar()) != '\n');
	}
 
/*********************************************************************
* Writes a summary of all the input and output performed so far
*********************************************************************/

	if(smry==1) {
	if(ans==0)
	{
	  fprintf(fp,"\nTitle: %s\n", locus);
	  if(form==0)
	    fprintf(fp,"Form : p^2\n\n");
	  else
	    fprintf(fp,"Form : 2p\n\n");			       
	  for(l=0;l<nbas;l+=3)
	  {
	    for(j=l;j<=min(l+2,nbas-1);j++) {
		fprintf(fp, "Database: %s",base[j]);
		for(k=13;k>strlen(base[j]);k--)
		  fprintf(fp, " ");
	    }
	     fprintf(fp,"\n");				
	     for(i=l;i<=min(l+2,nbas-1);i++)
	        fprintf(fp,"Allele        Freq     ");
	     fprintf(fp,"\n");
	     for(k=l;k<=min(l+2,nbas-1);k++)
		fprintf(fp,"------------------     ");
	     fprintf(fp,"\n");
	     for(a=0; a<allele; a++)
	     {
	       for(i=l;i<=min(l+2,nbas-1);i++) {
		   fprintf(fp, "%s", name[a]);
		   for(k=9;k>strlen(name[a]);k--)
		     fprintf(fp, " ");
	           fprintf(fp, "%9.6lf     ",fq[a][i]);
	     }
	     fprintf(fp,"\n");
	   }
	   for(k=l;k<=min(l+2,nbas-1);k++)
	     fprintf(fp,"------------------     ");
	   fprintf(fp, "\n\n");
	}
	fprintf(fp,"****************************************************************************\n");
	}
      
	fprintf(fp,"\nUnknown Alleles in the Numerator:\n");
	for(a=1;a<=unk_al;a++)
	  fprintf(fp,"%s \n",name[us[a-1]-1]);
	if(unk_al==0)
	fprintf(fp,"None\n");
	for(j=cont_l;j<=cont_u;j++)
	{
	  fprintf(fp,"\nUnknown Contributors: %d\n", j);
	  fprintf(fp,"Database       Numerator Probability\n");
	  fprintf(fp,"------------------------------------\n");
	  for(i=1;i<=nbas;i++) {
	    fprintf(fp,"%s", base[i-1]);
	  for(a=15;a>strlen(base[i-1]);a--)
	  fprintf(fp, " ");
	  fprintf(fp,"%16.6E\n",lkh_nm[i-1][j]);
	  }
	}
	fprintf(fp,"\n");

	fprintf(fp,"\nUnknown Alleles in the Denominator:\n");
	for(a=1;a<=unk_al2;a++)
	  fprintf(fp,"%s\n",name[us1[a-1]-1]);
	if(unk_al2==0)
	  fprintf(fp,"None\n");

	for(j=cont_l2;j<=cont_u2;j++)
	{
	  fprintf(fp,"\nUnknown Contributors: %d\n", j);
	  fprintf(fp,"Database     Denominator Probability\n");
	  fprintf(fp,"------------------------------------\n");
	  for(i=1;i<=nbas;i++) {
	    fprintf(fp,"%s", base[i-1]);
	  for(a=15;a>strlen(base[i-1]);a--)
	    fprintf(fp, " ");
	  fprintf(fp,"%16.6E\n",lkh_dn[i-1][j]);
	  }
	}
	fprintf(fp,"\n");
	fprintf(fp,"\n           Numerator Unknown    Denominator Unknown");
	fprintf(fp,"\nDatabase      Contributors         Contributors         Likelihood Ratio");
	fprintf(fp,"\n------------------------------------------------------------------------\n");
	for(i=0;i<nbas;i++)
	{
         for(j=cont_l;j<=cont_u;j++)
	   for(k=cont_l2;k<=cont_u2;k++) {
	     fprintf(fp,"%s", base[i]);
	     for(a=12; a>strlen(base[i]); a--)
	       fprintf(fp," ");
	     fprintf(fp,"%8d %20d %30.2lf\n", j, k, lkh_nm[i][j]/lkh_dn[i][k]);
	   }
	 fprintf(fp,"\n");
        }

	fprintf(fp,"\n****************************************************************************\n");
	}
      
	do {
	printf("\nDo you want to compute another likelihood ratio? (No=0, Yes=1): ");
	scanf("%d", &ans);
	if (ans != 0 && ans != 1)
	  printf("Not a valid entry!\n");
	}while(ans != 0 && ans != 1);

	}while(ans != 0);
	printf("\n\n           Program written by John Storey on July 7, 1997.\n");
	printf("            Questions and/or comments should be sent to:  \n");
	printf("                     storey@statgen.ncsu.edu              \n\n");
}

