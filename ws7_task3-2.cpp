#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#define N 5
#define M 6
#define NNZ 8
#define N_BIG 7000
#define NNZ_BIG 20000

struct coordinate_form_type_struct {
    double *AA;
    int *IA;
    int *JA;
    int length;
};

typedef struct coordinate_form_type_struct coordinate_form_type;

struct csr_type_struct {
    double *AA;
    int *JA;
    int *IA;
    int length_i;
    int length_j;
};

typedef struct csr_type_struct csr_type;

struct csr_extraction_type_struct {
    double *AA;
    int *JA;
    int length_j;
    int length_diag;
};

typedef struct csr_extraction_type_struct csr_extraction_type;

int min(int a, int b) {
    return a < b ? a : b;
}


int random_int(int max) {
    int i;

    while ((i = (int) (((double) rand()) / ((double) RAND_MAX) * (max + 1))) 
	   >= max);
    
    return i;
}


double** make_matrix(int n, int m) {
    double** A;
    int i, k;

    A = (double**) malloc(n * sizeof(double*));
    for (i=0;i<n;i++) 
	A[i] = (double*) malloc(m * sizeof(double));

    for (i=0;i<n;i++) 
	for (k=0;k<m;k++) 
	    A[i][k] = 0.0;

    return A;
}

void free_matrix(double** A, int n, int m) {
    int i;

    for (i=0;i<n; i++) 
	free(A[i]);
    free(A);
    
}

double **fill_random_sparse_matrix(double **A, int n, int m, int nnz) {
    int i, k, l;

    for (i=0;i<nnz/3;i++) {
	k = random_int(min(n,m));
	A[k][k] = ((double) rand()) / ((double) RAND_MAX); 
    }

    for (i=0;i<nnz - (nnz / 3);i++) {
	k = random_int(n);
	l = random_int(m);
	A[k][l] = ((double) rand()) / ((double) RAND_MAX); 
    }

    return A;
}

double **fill_random_matrix(double **A, int n, int m) {
    int i,k;

    for (i=0;i<n;i++) 
	for (k=0;k<m;k++) 
	    A[i][k] = ((double) rand()) / ((double) RAND_MAX);

    return A;
}

void print_matrix(double** A, int n, int m) {
    int i,k;

    for (i=0; i<n; i++) {
	for (k=0; k<m; k++) {
	    printf("%e ",A[i][k]);
	}
	printf("\n");
    }
    printf("\n\n");
}

double *make_vector(int n) {
    return (double*) malloc(n * sizeof(double));
}

void free_vector(double *v) {
    free(v);
}

double *fill_random_vector(double *v, int n) {
    int i;
    for (i=0;i<n;i++) 
	v[i] = ((double) rand()) / ((double) RAND_MAX);
    return v;
}

void print_vector(double *v, int n) {
    int i;
    for (i=0;i<n;i++) 
	printf("%e\n",v[i]);
    printf("\n");
}

double *multiply_matrix_vector(double **A, double *v, double *w, int n, int m) {
    int i, k;

    for (i=0;i<n;i++) 
	w[i] = 0.0;

    for (i=0;i<n;i++) {
	for(k=0;k<m;k++) {
	    w[i] += A[i][k] * v[k];
	}
    }

    return w;
}

int get_number_of_nonzeros(double **A, int n, int m) {
    int i,k,j;

    j = 0;
    for (i=0;i<n;i++) 
	for (k=0;k<m;k++) 
	    if (A[i][k] != 0.0) j++;

    return j;
}

coordinate_form_type *make_coordinate_form(int nnz) {
    coordinate_form_type *A;

    A = (coordinate_form_type*) malloc(sizeof(coordinate_form_type));

    A->length = nnz;
    A->AA = (double *) malloc(nnz * sizeof(double));
    A->JA = (int *) malloc(nnz * sizeof(int));
    A->IA = (int *) malloc(nnz * sizeof(int));

    return A;
}

void free_coordinate_form(coordinate_form_type *A) {
    free(A->IA);
    free(A->JA);
    free(A->AA);
    free(A);
}

void print_coordinate_form(coordinate_form_type *A) {
    int i;
    printf("AA:\n");
    for (i=0;i<A->length;i++) {
	printf("%e ",A->AA[i]);
    }
    printf("\nIA:\n");
    for (i=0;i<A->length;i++) {
	printf("%d ",A->IA[i]);
    }
    printf("\nJA:\n");
    for (i=0;i<A->length;i++) {
	printf("%d ",A->JA[i]);
    }
    printf("\n\n");
}


coordinate_form_type *matrix_to_coordinate_form(double **A, int n, int m, coordinate_form_type *A_coor) {
    int i,k,l;

    l = 0;
    for (i=0;i<n;i++) 
	for (k=0;k<M;k++) 
	    if (A[i][k] != 0) {
		A_coor->AA[l] = A[i][k];
		A_coor->JA[l] = k;
		A_coor->IA[l++] = i;
	    }
    
    return A_coor;
}

double **coordinate_form_to_matrix(coordinate_form_type *A_coor, double **A) {
    int l;
    
    for (l=0;l<A_coor->length;l++) 
	A[A_coor->IA[l]][A_coor->JA[l]] = A_coor->AA[l];

    return A;
}

double *multiply_matrix_coordinate_form_with_vector(coordinate_form_type *A_coor, double *v, double *w, int n) {
    int i,k,l;

    for (i=0;i<n;i++) 
	w[i] = 0.0;
    
    for (i=0;i<A_coor->length;i++) {
	w[A_coor->IA[i]] += A_coor->AA[i] * v[A_coor->JA[i]];
    }
    
    return w;
}

csr_type *make_csr(int length_i, int length_j) {
    csr_type *A;
    
    A = (csr_type*) malloc(sizeof(csr_type));

    A->AA = (double*) malloc(length_i * sizeof(double));
    A->JA = (int*) malloc(length_i * sizeof(int));
    A->IA = (int*) malloc((length_j + 1) * sizeof(int));
    A->length_i = length_i;
    A->length_j = length_j;

    return A;
}

void free_csr(csr_type *A) {
    free(A->AA);
    free(A->IA);
    free(A->JA);
    free(A);
}

void print_csr(csr_type *A) {
    int i;
    printf("AA:\n");
    for(i=0;i<A->length_i;i++) {
	printf("%e ",A->AA[i]);
    }
    printf("\nIA:\n");
    for(i=0;i<A->length_j+1;i++) {
	printf("%d ",A->IA[i]);
    }
    printf("\nJA:\n");
    for(i=0;i<A->length_i;i++) {
	printf("%d ",A->JA[i]);
    }
    printf("\n\n");
}

csr_type *matrix_to_csr(double **A, int n, int m, csr_type *A_csr) {
    int i,k,l,r;

    l = 0;
    A_csr->IA[0] = l;
    for (i=0; i<n; i++) {
	for (k=0; k<m; k++) {
	    if (A[i][k] != 0.0) {
		A_csr->AA[l] = A[i][k];
		A_csr->JA[l++] = k;
	    }
	}
	A_csr->IA[i+1] = l;
    }


    return A_csr;
}

double **csr_to_matrix(csr_type *A_csr, double** A) {
    int i,k,l;

    for (i=0;i<A_csr->length_j;i++) {
	for (k=A_csr->IA[i];k<A_csr->IA[i+1];k++) {
	    A[i][A_csr->JA[k]] = A_csr->AA[k];
	}
    }

    return A;
}

double *multiply_matrix_csr_with_vector(csr_type *A_csr, double *v, double *w, int n) {
    int i,k,l;

    for (i=0;i<n;i++) 
	w[i] = 0.0;
    
    for (i=0;i<A_csr->length_j;i++) {
	for (k=A_csr->IA[i];k<A_csr->IA[i+1];k++) {
	    w[i] += A_csr->AA[k] * v[A_csr->JA[k]];
	}
    }
    
    return w;
}



csr_extraction_type *make_csr_extraction(int length_j, int length_diag) {
    csr_extraction_type *A;

    A = (csr_extraction_type *) malloc(sizeof(csr_extraction_type));
    A->AA = (double *) malloc((length_j + 1) * sizeof(double));
    A->JA = (int *) malloc((length_j + 1) * sizeof(double));
    A->length_j = length_j;
    A->length_diag = length_diag;
    
    return A;
}

void free_csr_extraction(csr_extraction_type *A) {
    free(A->AA);
    free(A->JA);
    free(A);
}

void print_csr_extraction(csr_extraction_type *A) {
    int i;

    printf("AA:\n");
    for (i=0;i<A->length_j+1;i++) 
	printf("%e ",A->AA[i]);
    printf("\nJA:\n");
    for (i=0;i<A->length_j+1;i++) 
	printf("%d ",A->JA[i]);
    printf("\n\n");
}

int get_number_of_nonzeros_out_of_diagonal(double **A, int n, int m) {
    int i,k,l;

    l = 0;
    for (i=0;i<n;i++) {
	for (k=0;k<m;k++) {
	    if ((i != k) && (A[i][k] != 0.0)) l++;
	}
    }

    return l;
}

csr_extraction_type *matrix_to_csr_extraction(double **A, int n, int m, csr_extraction_type *A_csr) {
    int i,k,l,d;

    d = min(n,m);
    for (i=0;i<d;i++) 
	A_csr->AA[i] = A[i][i];

    l = d+1;
    A_csr->JA[0] = l;
    for (i=0; i<n; i++) {
	for (k=0; k<m; k++) {
	    if ((i != k) && (A[i][k] != 0.0)) {
		A_csr->AA[l] = A[i][k];
		A_csr->JA[l++] = k;
	    }
	}
	A_csr->JA[i+1] = l;
    }

    return A_csr;
}

double **csr_extraction_to_matrix(csr_extraction_type *A_csr, double** A) {
    int i,k,l;

    for (i=0;i<A_csr->length_diag; i++) {
	A[i][i] = A_csr->AA[i];
    }

    for (i=0;i<A_csr->length_diag;i++) {
	for (k=A_csr->JA[i];k<A_csr->JA[i+1];k++) {
	    A[i][A_csr->JA[k]] = A_csr->AA[k];
	}
    }

    return A;
}


double *multiply_matrix_csr_extraction_with_vector(csr_extraction_type *A_csr, double *v, double *w) {
    int i,k,l;

    for (i=0;i<A_csr->length_diag;i++) {
	w[i] = A_csr->AA[i] * v[i];
	for (k=A_csr->JA[i];k<A_csr->JA[i+1];k++) {
	    w[i] += A_csr->AA[k] * v[A_csr->JA[k]];
	}
    }
    
    return w;
}


int main(int argc, char** argv) {
    double **A, **A_coor_back, **A_csr_back, **A_csr_extr_back;
    double *b, *c, *c_coor, *c_csr, *c_csr_extr;
    coordinate_form_type *A_coor;
    csr_type *A_csr;
    csr_extraction_type *A_csr_extr;
    int nnz, nnzood, diag;
    clock_t start, end;


    A = make_matrix(N,M);
    A_coor_back = make_matrix(N,M);
    A_csr_back = make_matrix(N,M);
    A_csr_extr_back = make_matrix(N,M);
    b = make_vector(M);
    c = make_vector(N);
    c_coor = make_vector(N);
    c_csr = make_vector(N);
    c_csr_extr = make_vector(N);

    fill_random_sparse_matrix(A,N,M,NNZ);
    fill_random_vector(b,M);
   
    multiply_matrix_vector(A,b,c,N,M);
    
    printf("Matrix:\n");
    print_matrix(A,N,M);
    printf("Vector to be multiplied:\n");
    print_vector(b,M);
    printf("Multiplication result:\n");
    print_vector(c,N);

    nnz = get_number_of_nonzeros(A,N,M);
    nnzood = get_number_of_nonzeros_out_of_diagonal(A,N,M);
    diag = min(N,M);

    A_coor = make_coordinate_form(nnz);
    A_csr = make_csr(nnz,N);
    A_csr_extr = make_csr_extraction(nnzood + diag, diag);
    
    matrix_to_coordinate_form(A,N,M,A_coor);
    matrix_to_csr(A,N,M,A_csr);
    matrix_to_csr_extraction(A,N,M,A_csr_extr);


    printf("Matrix in coordinate form:\n");
    print_coordinate_form(A_coor);
    printf("Matrix in CSR form:\n");
    print_csr(A_csr);
    printf("Matrix in CSR form with extraction of the diagonal:\n");
    print_csr_extraction(A_csr_extr);

    coordinate_form_to_matrix(A_coor,A_coor_back);
    csr_to_matrix(A_csr,A_csr_back);
    csr_extraction_to_matrix(A_csr_extr,A_csr_extr_back);

    printf("Recovered coordinate form matrix:\n");
    print_matrix(A_coor_back,N,M);
    printf("Recovered CSR form matrix:\n");
    print_matrix(A_csr_back,N,M);
    printf("Recovered CSR with extraction of the diagonal form matrix:\n");
    print_matrix(A_csr_extr_back,N,M);

    multiply_matrix_coordinate_form_with_vector(A_coor,b,c_coor,N);
    multiply_matrix_csr_with_vector(A_csr,b,c_csr,N);
    multiply_matrix_csr_extraction_with_vector(A_csr_extr,b,c_csr_extr);

    printf("Multiplication result using the coordinate form:\n");
    print_vector(c_coor,N);
    printf("Multiplication result using the CSR form:\n");
    print_vector(c_csr,N);
    printf("Multiplication result using the CSR form with extraction of the diagonal:\n");
    print_vector(c_csr_extr,N);

    free_coordinate_form(A_coor);
    free_csr(A_csr);
    free_csr_extraction(A_csr_extr);
    free_matrix(A,N,M);
    free_matrix(A_coor_back,N,M);
    free_matrix(A_csr_back,N,M);
    free_matrix(A_csr_extr_back,N,M);
    free_vector(b);
    free_vector(c);
    free_vector(c_coor);
    free_vector(c_csr);
    free_vector(c_csr_extr);

    printf("Multiplying now a %dx%d matrix (%d nonzeros) with a %d vector using the different methods:\n",
	   N_BIG,N_BIG,NNZ_BIG,N_BIG);
    A = make_matrix(N_BIG,N_BIG);
    b = make_vector(N_BIG);
    fill_random_sparse_matrix(A,N_BIG,N_BIG,NNZ_BIG);
    fill_random_vector(b,N_BIG);


    nnz = get_number_of_nonzeros(A,N_BIG,N_BIG);
    nnzood = get_number_of_nonzeros_out_of_diagonal(A,N_BIG,N_BIG);
    diag = N_BIG;

    A_coor = make_coordinate_form(nnz);
    A_csr = make_csr(nnz,N_BIG);
    A_csr_extr = make_csr_extraction(nnzood + diag, diag);

    matrix_to_coordinate_form(A,N_BIG,N_BIG,A_coor);
    matrix_to_csr(A,N_BIG,N_BIG,A_csr);
    matrix_to_csr_extraction(A,N_BIG,N_BIG,A_csr_extr);

    c = make_vector(N_BIG);
    start = clock();
    multiply_matrix_vector(A,b,c,N_BIG,N_BIG);
    end = clock();
    printf("Directly stored matrix and vector multiplication: %d clock ticks\n",end - start);

    c = make_vector(N_BIG);
    start = clock();
    multiply_matrix_coordinate_form_with_vector(A_coor,b,c,N_BIG);
    end = clock();
    printf("Coordinate form matrix and vector multiplication: %d clock ticks\n",end - start);
    free_vector(c);

    c = make_vector(N_BIG);
    start = clock();
    multiply_matrix_csr_with_vector(A_csr,b,c,N_BIG);
    end = clock();
    printf("CSR matrix and vector multiplication: %d clock ticks\n",end - start);
    free_vector(c);

    c = make_vector(N_BIG);
    start = clock();
    multiply_matrix_csr_extraction_with_vector(A_csr_extr,b,c);
    end = clock();
    printf("CSR with extraction of the diagonal - stored matrix and vector multiplication: %d clock ticks\n",
	   end - start);
    free_vector(c);
    
    free_vector(b); 
    free_matrix(A,N_BIG,N_BIG);
    free_coordinate_form(A_coor);
    free_csr(A_csr);
    free_csr_extraction(A_csr_extr);

    return 0;
}
