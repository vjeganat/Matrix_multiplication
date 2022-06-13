/* 
    Please include compiler name below (you may also include any other modules you would like to be loaded)

COMPILER= gnu

module load gcc/4.9.3 mkl

OPT = -O3
CFLAGS = -Wall -std=gnu99 $(OPT) -m64 -I${MKLROOT}/include -mavx -ftree-vectorize  -ffast-math -fopt-info-vec-optimized  -fpeel-loops 

Above are the only lines changed in the Makefile.

    Please include All compiler flags and libraries as you want them run. You can simply copy this over from the Makefile's first few lines

*/


/*************************************************************************************
 *
 *
 *
 *
 *
 *
 ***************************************************************************************/


const char* dgemm_desc = "Simple blocked dgemm.";

#if !defined(BLOCK_SIZE)
#define BLOCK_SIZE 128 
	// block size of 104 works best for first option. 128 seems to work best for our ijk+=2 code	
	
#endif

#define min(a,b) (((a)<(b))?(a):(b))



// This is our best code. 
// ijk+=2 loop unrolling.
// comment out below and uncomment the indicated next section to run the other version.
//
// block size of 128 works best here. 
//

static void do_block_two_x_two_x_two_tiles (const int lda, const int M, const int N, const int K,  double* A, double* restrict B,  double* restrict C)
{
if(M%2==0&&N%2==0&&K%2==0){//if all M, N and K are even
	int i =0;
	int j =0;
	int k = 0;
  	// For each column of B and column of C 
  	for ( j = 0; j < M-1; j+=2){
    	   // For each row of B and column of A  
    	   for ( i = 0; i < N-1; i+=2){
		double bij = B[i+j*lda];
                double biij = B[i+1+j*lda];
                double bijj = B[i+(j+1)*lda];//loop unrolling in j
                double biijj = B[i+1+(j+1)*lda];
		// For each row of A and row of C
     	 	for ( k = 0; k < K; k+=2){
		  double aki= A[k+i*lda];
		  double akki= A[k+1+i*lda];//register blocking
		  double akii= A[k+(i+1)*lda];//loop unrolling in i
		  double akkii= A[k+1+(i+1)*lda];//loop unrolling in k
    		   C[k+j*lda] += aki * bij + akii *biij;
		   C[k+1+j*lda] += akki * bij + akkii *biij;
      		   C[k+(j+1)*lda] += aki * bijj+ akii * biijj;
		   C[k+1+(j+1)*lda] +=  akki * bijj + akkii * biijj;
		}

    	   }
 	}
}
else if(M%2==0&&N%2==0){//if M and N are even but K is odd
	int i =0;
	int j =0;
	int k = 0;
  	for (j = 0; j < M-1; j+=2){
    	   for (i = 0; i < N-1; i+=2){
		double bij = B[i+j*lda];
                double biij = B[i+1+j*lda];
                double bijj = B[i+(j+1)*lda];
                double biijj = B[i+1+(j+1)*lda];
     	 	for (k = 0; k < K-1; k+=2){
		  double aki= A[k+i*lda];
		  double akki= A[k+1+i*lda];
		  double akii= A[k+(i+1)*lda];
		  double akkii= A[k+1+(i+1)*lda];
    		   C[k+j*lda] += aki * bij + akii *biij;
		   C[k+1+j*lda] += akki * bij + akkii *biij;
      		   C[k+(j+1)*lda] += aki * bijj+ akii * biijj;
		   C[k+1+(j+1)*lda] +=  akki * bijj + akkii * biijj;
		}
		   double aki= A[k+i*lda];
		   double akii=A[k+(i+1)*lda];
                   C[k+j*lda] += aki * bij + akii *biij;
                   C[k+(j+1)*lda] += akii * biijj + aki * bijj;

    	   }
 	}
}
else if(M%2==0&&K%2==0){//if M and K are even but N is odd 
	int i =0;
	int j =0;
	int k = 0;
  	for (j = 0; j < M-1; j+=2){
    	   for (i = 0; i < N-1; i+=2){
		double bij = B[i+j*lda];
                double biij = B[i+1+j*lda];
                double bijj = B[i+(j+1)*lda];
                double biijj = B[i+1+(j+1)*lda];
     	 	for (k = 0; k < K; k+=2){
		  double aki= A[k+i*lda];
		  double akki= A[k+1+i*lda];
		  double akii= A[k+(i+1)*lda];
		  double akkii= A[k+1+(i+1)*lda];
    		   C[k+j*lda] += aki * bij + akii *biij;
		   C[k+1+j*lda] += akki * bij + akkii *biij;
      		   C[k+(j+1)*lda] += aki * bijj+ akii * biijj;
		   C[k+1+(j+1)*lda] +=  akki * bijj + akkii * biijj;
		}

    	   }
		double bij = B[i+j*lda];
                double bijj = B[i+(j+1)*lda];
     	 	for (k = 0; k < K; k+=2){
		   double aki= A[k+i*lda];
		   double akki= A[k+1+i*lda];
     		   C[k+j*lda] += aki * bij;
		   C[k+1+j*lda] += akki * bij;
      		   C[k+(j+1)*lda] += aki * bijj;
		   C[k+1+(j+1)*lda] += akki * bijj;
		}

 	}
}
else if(M%2==0){//if M is even and N, K odd
	int i =0;
	int j =0;
	int k = 0;
  	for (j = 0; j < M-1; j+=2){
    	   for (i = 0; i < N-1; i+=2){
		double bij = B[i+j*lda];
                double biij = B[i+1+j*lda];
                double bijj = B[i+(j+1)*lda];
                double biijj = B[i+1+(j+1)*lda];
     	 	for (k = 0; k < K-1; k+=2){
		  double aki= A[k+i*lda];
		  double akki= A[k+1+i*lda];
		  double akii= A[k+(i+1)*lda];
		  double akkii= A[k+1+(i+1)*lda];
    		   C[k+j*lda] += aki * bij + akii *biij;
		   C[k+1+j*lda] += akki * bij + akkii *biij;
      		   C[k+(j+1)*lda] += aki * bijj+ akii * biijj;
		   C[k+1+(j+1)*lda] +=  akki * bijj + akkii * biijj;
		}
		   double aki= A[k+i*lda];
		   double akii= A[k+(i+1)*lda];
                   C[k+j*lda] += aki * bij + akii *biij;
                   C[k+(j+1)*lda] += aki * bijj + akii * biijj ;

    	   }
		double bij = B[i+j*lda];
                double bijj = B[i+(j+1)*lda];
                for (k = 0; k < K-1; k+=2){
                   double aki= A[k+i*lda];
                   double akki= A[k+1+i*lda];
                   C[k+j*lda] += aki * bij;
                   C[k+1+j*lda] += akki * bij;
                   C[k+(j+1)*lda] += aki * bijj;
                   C[k+1+(j+1)*lda] += akki * bijj;
                }
                   double aki= A[k+i*lda];
                   C[k+j*lda] += aki * bij;
                   C[k+(j+1)*lda] += aki * bijj;

 	}
}
else if(N%2==0&&K%2==0){//if N,K even and M is odd
	int i =0;
	int j =0;
	int k = 0;
  	for (j = 0; j < M-1; j+=2){
    	   for (i = 0; i < N-1; i+=2){
		double bij = B[i+j*lda];
                double biij = B[i+1+j*lda];
                double bijj = B[i+(j+1)*lda];
                double biijj = B[i+1+(j+1)*lda];
     	 	for (k = 0; k < K; k+=2){
		  double aki= A[k+i*lda];
		  double akki= A[k+1+i*lda];
		  double akii= A[k+(i+1)*lda];
		  double akkii= A[k+1+(i+1)*lda];
    		   C[k+j*lda] += aki * bij + akii *biij;
		   C[k+1+j*lda] += akki * bij + akkii *biij;
      		   C[k+(j+1)*lda] += aki * bijj+ akii * biijj;
		   C[k+1+(j+1)*lda] +=  akki * bijj + akkii * biijj;
		}

    	   }
 	}
  	   for (i = 0; i < N-1; i+=2){
                double bij = B[i+j*lda];
                double biij = B[i+1+j*lda];
                for (k = 0; k < K; k+=2){
		   double aki= A[k+i*lda];
		   double akki=A[k+1+i*lda];
		   double akii=A[k+(i+1)*lda];
		   double akkii=A[k+1+(i+1)*lda];
                   C[k+j*lda] += aki * bij + akii*biij;
                   C[k+1+j*lda] += akkii*biij+akki * bij ;
                }

           }

}
else if(N%2==0){//if N is even and M, K are odd
	int i =0;
	int j =0;
	int k = 0;
  	for (j = 0; j < M-1; j+=2){
    	   for (i = 0; i < N-1; i+=2){
		double bij = B[i+j*lda];
                double biij = B[i+1+j*lda];
                double bijj = B[i+(j+1)*lda];
                double biijj = B[i+1+(j+1)*lda];
     	 	for (k = 0; k < K-1; k+=2){
		  double aki= A[k+i*lda];
		  double akki= A[k+1+i*lda];
		  double akii= A[k+(i+1)*lda];
		  double akkii= A[k+1+(i+1)*lda];
    		   C[k+j*lda] += aki * bij + akii *biij;
		   C[k+1+j*lda] += akki * bij + akkii *biij;
      		   C[k+(j+1)*lda] += aki * bijj+ akii * biijj;
		   C[k+1+(j+1)*lda] +=  akki * bijj + akkii * biijj;
		}
		   double aki= A[k+i*lda];
		   double akii=A[k+(i+1)*lda];
                   C[k+j*lda] += aki * bij + akii *biij;
                   C[k+(j+1)*lda] +=aki * bijj + akii * biijj ;

    	   }
 	}
	  for (i = 0; i < N-1; i+=2){
                double bij = B[i+j*lda];
                double biij = B[i+1+j*lda];
                for (k = 0; k < K-1; k+=2){
		   double aki=A[k+i*lda];
		   double akki=A[k+1+i*lda];
		   double akii=A[k+(i+1)*lda];
		   double akkii=A[k+1+(i+1)*lda];
                   C[k+j*lda] += aki * bij + akii *biij;
                   C[k+1+j*lda] += akki * bij + akkii *biij;
                }
                   C[k+j*lda] += A[k+i*lda] * bij + A[k+(i+1)*lda] *biij;

	 }	
}else if(K%2==0){//if M and N are odd but K is even
	int i =0;
	int j =0;
	int k = 0;
  	for (j = 0; j < M-1; j+=2){
    	   for (i = 0; i < N-1; i+=2){
		double bij = B[i+j*lda];
                double biij = B[i+1+j*lda];
                double bijj = B[i+(j+1)*lda];
                double biijj = B[i+1+(j+1)*lda];
     	 	for (k = 0; k < K; k+=2){
		  double aki= A[k+i*lda];
		  double akki= A[k+1+i*lda];
		  double akii= A[k+(i+1)*lda];
		  double akkii= A[k+1+(i+1)*lda];
    		   C[k+j*lda] += aki * bij + akii *biij;
		   C[k+1+j*lda] += akki * bij + akkii *biij;
      		   C[k+(j+1)*lda] += aki * bijj+ akii * biijj;
		   C[k+1+(j+1)*lda] +=  akki * bijj + akkii * biijj;
		}

    	   }
		double bij = B[i+j*lda];
                double bijj = B[i+(j+1)*lda];
     	 	for (k = 0; k < K; k+=2){
		   double aki= A[k+i*lda];
		   double akki= A[k+1+i*lda];
     		   C[k+j*lda] += aki * bij;
		   C[k+1+j*lda] += akki * bij;
      		   C[k+(j+1)*lda] += aki * bijj;
		   C[k+1+(j+1)*lda] += akki * bijj;
		}

 	}
           for (i = 0; i < N-1; i+=2){
                double bij = B[i+j*lda];
                double biij = B[i+1+j*lda];
                for (k = 0; k < K; k+=2){
		   double aki=A[k+i*lda];
		   double akki=A[k+1+i*lda];
		   double akii=A[k+(i+1)*lda];
		   double akkii=A[k+1+(i+1)*lda];
                   C[k+j*lda] += aki * bij + akii *biij;
                   C[k+1+j*lda] += akki* bij + akkii *biij;
                }

           }
                double bij = B[i+j*lda];
                for (k = 0; k < K; k+=2){
                   C[k+j*lda] += A[k+i*lda] * bij;
                   C[k+1+j*lda] += A[k+1+i*lda] * bij;
                }

}else{//if M, N and K are all odd
	int i =0;
	int j =0;
	int k = 0;
	
  	for (j = 0; j < M-1; j+=2){
    	   for (i = 0; i < N-1; i+=2){
		double bij = B[i+j*lda];
                double biij = B[i+1+j*lda];
                double bijj = B[i+(j+1)*lda];
                double biijj = B[i+1+(j+1)*lda];
     	 	for (k = 0; k < K-1; k+=2){
		  double aki= A[k+i*lda];
		  double akki= A[k+1+i*lda];
		  double akii= A[k+(i+1)*lda];
		  double akkii= A[k+1+(i+1)*lda];
    		   C[k+j*lda] += aki * bij + akii *biij;
		   C[k+1+j*lda] += akki * bij + akkii *biij;
      		   C[k+(j+1)*lda] += aki * bijj+ akii * biijj;
		   C[k+1+(j+1)*lda] +=  akki * bijj + akkii * biijj;
		}
		   double aki= A[k+i*lda];
		   double akii=A[k+(i+1)*lda];
                   C[k+j*lda] += aki * bij + akii *biij;
                   C[k+(j+1)*lda] += aki * bijj + akii * biijj;

    	   }
		double bij = B[i+j*lda];
                double bijj = B[i+(j+1)*lda];
                for (k = 0; k < K-1; k+=2){
                   double aki= A[k+i*lda];
                   double akki= A[k+1+i*lda];
                   C[k+j*lda] += aki * bij;
                   C[k+1+j*lda] += akki * bij;
                   C[k+(j+1)*lda] += aki * bijj;
                   C[k+1+(j+1)*lda] += akki * bijj;
                }
                   double aki= A[k+i*lda];
                   C[k+j*lda] += aki * bij;
                   C[k+(j+1)*lda] += aki * bijj;

 	}
           for (i = 0; i < N-1; i+=2){
                double bij = B[i+j*lda];
                double biij = B[i+1+j*lda];
                for (k = 0; k < K-1; k+=2){
		   double aki=A[k+i*lda];
		   double akki=A[k+1+i*lda];
		   double akii=A[k+(i+1)*lda];
		   double akkii=A[k+1+(i+1)*lda];
                   C[k+j*lda] += aki * bij + akii *biij;
                   C[k+1+j*lda] += akki * bij + akkii *biij;
                }
                   C[k+j*lda] += A[k+i*lda] * bij + A[k+(i+1)*lda] *biij;

           }
                double bij = B[i+j*lda];
                for (k = 0; k < K-1; k+=2){
                   C[k+j*lda] += A[k+i*lda] * bij;
                   C[k+1+j*lda] += A[k+1+i*lda] * bij;
                }
                   C[k+j*lda] += A[k+i*lda] * bij;

}
}

void square_dgemm (int lda, double* A, double* B, double* C)
{

/* For each column of B and each row of A */ 
  for (int j = 0; j < lda; j += BLOCK_SIZE)
    /* For each column of B and each column of C */
    for (int i = 0; i < lda; i += BLOCK_SIZE)
      /* For each row of A and each row of C */
      for (int k = 0; k < lda; k += BLOCK_SIZE)
      {
	/* Correct block dimensions if block "goes off edge of" the matrix */
	int M = min (BLOCK_SIZE, lda-j);//Number of columns is B and number of rows in A
	int N = min (BLOCK_SIZE, lda-i);//Number of column in B and number of columns in C
	int K = min (BLOCK_SIZE, lda-k);//Number of rows in B and number of row in C
	
	// Perform individual block dgemm
	//Holding B and iterating through A and C blocks
	do_block_two_x_two_x_two_tiles (lda, M, N, K, A + k + i*lda, B + i + j*lda, C + k + j*lda);

	}
}


//======================================================================
//======================================================================

/*

//
//This is the our second one, and not the best.
//jik inner loop ordering, kji block ordering.
// i,j +=2 loop unrolling
// k++
// efficiency around 28%
//
//Comment out this whole section if not using
//


static void do_block_two_x_two_tiles (const int lda, const int M, const int N, const int K, const double* A, const double* B,  double* restrict  C)
{

// IF THE BLOCK IS EVEN IN M AND N
if(M%2==0&&N%2==0){
//	int i =0;
//	int j =0;
  	for (int j = 0; j < N-1; j+=2)
    	   for (int i = 0; i < M-1; i+=2){
     	 	double cij = C[i+j*lda];
     	 	double ciij = C[i+1+j*lda];
     	 	double cijj = C[i+(j+1)*lda];
     	 	double ciijj = C[i+1+(j+1)*lda];
    	 	for (int k = 0; k < K; k++){			//auto-vectorized?
  		   cij += A[i+k*lda] * B[k+j*lda];
 		   cijj += A[i+k*lda] * B[k+(j+1)*lda];
		   ciijj += A[i+1+k*lda] * B[k+(j+1)*lda];
		   ciij += A[i+1 +k*lda] *B[k+j*lda];
	 	}
     	 	C[i+j*lda] = cij;
    	 	C[i+1+j*lda] = ciij;
     	 	C[i+(j+1)*lda]= cijj;
     	 	C[i+1+(j+1)*lda] = ciijj;	
    	   }

// IF THE BLOCK HAS ODD NUMBER OF COLUMNS 
}else if(M%2==0){  //N-odd only
//	int i =0;
	int j =0;
//	int k = 0;
  	for (j = 0; j < N-1; j+=2){
    	   for (int i = 0; i < M-1; i+=2){
		double cij = C[i+j*lda];
     	 	double ciij = C[i+1+j*lda];
     	 	double cijj = C[i+(j+1)*lda];
     	 	double ciijj = C[i+1+(j+1)*lda];
     	 	for (int k = 0; k < K; k++){		//auto-vectorized?
     		   cij += A[i+k*lda] * B[k+j*lda];
      		   cijj += A[i+k*lda] * B[k+(j+1)*lda];
      		   ciijj += A[i+1+k*lda] * B[k+(j+1)*lda];
      		   ciij += A[i+1 +k*lda] *B[k+j*lda];
	 	}
     	 	C[i+j*lda] = cij;
     	 	C[i+1+j*lda] = ciij;
     	 	C[i+(j+1)*lda]= cijj;
     	 	C[i+1+(j+1)*lda] = ciijj;	
    	   }

   	}
	
	for(int i =0;i< M-1;i+=2){
		double cij = C[i+j*lda];
		double ciij = C[i+1+j*lda]; 	
		for(int k=0;k<K;k++){			//auto-vectorized?
		cij += A[i+k*lda] *B[k+j*lda];
		ciij += A[i+1+k*lda]*B[k+j*lda];
		}
		C[i+j*lda] = cij;
		C[i+1+j*lda] = ciij;
	}

// IF THE BLOCK HAS AN ODD NUMBER OF ROWS
}else if(N%2==0){ //M-odd only
	int i =0;
	int j =0;
	int k = 0;
 	for (j = 0; j < N-1; j+=2){
  	   for (i = 0; i < M-1; i+=2){
     	 	double cij = C[i+j*lda];
     	 	double ciij = C[i+1+j*lda];
     	 	double cijj = C[i+(j+1)*lda];
     	 	double ciijj = C[i+1+(j+1)*lda];
     	 	for (k = 0; k < K; ++k){		//auto-vectorized?
     		   cij += A[i+k*lda] * B[k+j*lda];
      		   cijj += A[i+k*lda] * B[k+(j+1)*lda];
      		   ciijj += A[i+1+k*lda] * B[k+(j+1)*lda];
      		   ciij += A[i+1+k*lda] *B[k+j*lda];
	 	}
     	 	C[i+j*lda] = cij;
     	 	C[i+1+j*lda] = ciij;
     	 	C[i+(j+1)*lda]= cijj;
     	 	C[i+1+(j+1)*lda] = ciijj;	
    	   }
	double cij = C[i+j*lda];
	double cijj = C[i+(j+1)*lda];
	 for(k=0;k<K;k++){				//auto-vectorized?
		cij += A[i+k*lda] * B[k+j*lda];
		cijj += A[i+k*lda] * B[k+(j+1)*lda];
	   }
		C[i+j*lda] = cij;
		C[i+(j+1)*lda] = cijj;	
	}

// Below is pretty much a combination of the two loops above, 
// and ends up being called only once, if at all. 
// (just the bottom-right block)


// IF BOTH M AND N ARE ODD
}else{ //both M and N are odd
	int i =0;
	int j =0;
	int k = 0;
  	// For each row i of A 
  	for (j = 0; j < N-1; j+=2){
    	   // For each column j of B  
    	   for (i = 0; i < M-1; i+=2){
		// Compute C(i,j) 
     	 	double cij = C[i+j*lda];
     	 	double ciij = C[i+1+j*lda];
     	 	double cijj = C[i+(j+1)*lda];
     	 	double ciijj = C[i+1+(j+1)*lda];
     	 	for (k = 0; k < K; k++){	//auto-vectorized?
     		   cij += A[i+k*lda] * B[k+j*lda];
      		   cijj += A[i+k*lda] * B[k+(j+1)*lda];
      		   ciijj += A[i+1+k*lda] * B[k+(j+1)*lda];
      		   ciij += A[i+1 +k*lda] *B[k+j*lda];
	 	}
     	 	C[i+j*lda] = cij;
     	 	C[i+1+j*lda] = ciij;
     	 	C[i+(j+1)*lda]= cijj;
     	 	C[i+1+(j+1)*lda] = ciijj;	
   	   }

     	   double cij = C[i+j*lda];
     	   double cijj = C[i+(j+1)*lda];
     	 	for (k = 0; k < K; k++){	//auto-vectorized?
     		   cij += A[i+k*lda] * B[k+j*lda];
      		   cijj += A[i+k*lda] * B[k+(j+1)*lda];
	 	}
     	 	C[i+j*lda] = cij;
     	 	C[i+(j+1)*lda]= cijj;
  	}

	for(i=0;i<M-1;i+=2){
	 double cij = C[i+j*lda];
	   double ciij = C[i+1+(j)*lda]; 
	for(k=0;k<K;k++){			//auto-vectorized?
		cij += A[i+k*lda] * B[k+(j)*lda];
		ciij += A[i+1+k*lda] * B[k+(j)*lda];
	   }
		C[i+j*lda] = cij;
		C[i+1+(j)*lda] = ciij;
    	}

	double cij = C[i+j*lda];
	for(k=0;k<K;k++){			//auto-vectorized?
	   cij += A[i+k*lda] * B[k+j*lda];
	}
	C[i+j*lda] = cij;
}

} //end of method

// 

 



//Uncomment this square_dgemm() if you want to run the above code.
//Comment out the below square_dgemm(). 




void square_dgemm (const int lda, const double* A, const double* B, double* C)
{

// For each block-row of A 
  for (int k = 0; k < lda; k += BLOCK_SIZE)
    // For each block-column of B 
    for (int j = 0; j < lda; j += BLOCK_SIZE)
      // Accumulate block dgemms into block of C 
      for (int i = 0; i < lda;  i+= BLOCK_SIZE)
      {
	// Correct block dimensions if block "goes off edge of" the matrix 
	int M = min (BLOCK_SIZE, lda-i);
	int N = min (BLOCK_SIZE, lda-j);
	int K = min (BLOCK_SIZE, lda-k);
	
	// Perform individual block dgemm 
	
	do_block_two_x_two_tiles (lda, M, N, K, A + i + k*lda, B + k + j*lda, C + i + j*lda);
	}
}
// */








static void do_block_two_x_two_x_two_tiles (const int lda, const int M, const int N, const int K,  double* A, double* restrict B,  double* restrict C)
{
if(M%2==0&&N%2==0&&K%2==0){//if all M, N and K are even
	int i =0;
	int j =0;
	int k = 0;
  	// For each column of B and column of C 
  	for ( j = 0; j < M-1; j+=2){
    	   // For each row of B and column of A  
    	   for ( i = 0; i < N-1; i+=2){
		double bij = B[i+j*lda];
                double biij = B[i+1+j*lda];
                double bijj = B[i+(j+1)*lda];//loop unrolling in j
                double biijj = B[i+1+(j+1)*lda];
		// For each row of A and row of C
     	 	for ( k = 0; k < K; k+=2){
		  double aki= A[k+i*lda];
		  double akki= A[k+1+i*lda];//register blocking
		  double akii= A[k+(i+1)*lda];//loop unrolling in i
		  double akkii= A[k+1+(i+1)*lda];//loop unrolling in k
    		   C[k+j*lda] += aki * bij + akii *biij;
		   C[k+1+j*lda] += akki * bij + akkii *biij;
      		   C[k+(j+1)*lda] += aki * bijj+ akii * biijj;
		   C[k+1+(j+1)*lda] +=  akki * bijj + akkii * biijj;
		}

    	   }
 	}
}
else if(M%2==0&&N%2==0){//if M and N are even but K is odd
	int i =0;
	int j =0;
	int k = 0;
  	for (j = 0; j < M-1; j+=2){
    	   for (i = 0; i < N-1; i+=2){
		double bij = B[i+j*lda];
                double biij = B[i+1+j*lda];
                double bijj = B[i+(j+1)*lda];
                double biijj = B[i+1+(j+1)*lda];
     	 	for (k = 0; k < K-1; k+=2){
		  double aki= A[k+i*lda];
		  double akki= A[k+1+i*lda];
		  double akii= A[k+(i+1)*lda];
		  double akkii= A[k+1+(i+1)*lda];
    		   C[k+j*lda] += aki * bij + akii *biij;
		   C[k+1+j*lda] += akki * bij + akkii *biij;
      		   C[k+(j+1)*lda] += aki * bijj+ akii * biijj;
		   C[k+1+(j+1)*lda] +=  akki * bijj + akkii * biijj;
		}
		   double aki= A[k+i*lda];
		   double akii=A[k+(i+1)*lda];
                   C[k+j*lda] += aki * bij + akii *biij;
                   C[k+(j+1)*lda] += akii * biijj + aki * bijj;

    	   }
 	}
}
else if(M%2==0&&K%2==0){//if M and K are even but N is odd 
	int i =0;
	int j =0;
	int k = 0;
  	for (j = 0; j < M-1; j+=2){
    	   for (i = 0; i < N-1; i+=2){
		double bij = B[i+j*lda];
                double biij = B[i+1+j*lda];
                double bijj = B[i+(j+1)*lda];
                double biijj = B[i+1+(j+1)*lda];
     	 	for (k = 0; k < K; k+=2){
		  double aki= A[k+i*lda];
		  double akki= A[k+1+i*lda];
		  double akii= A[k+(i+1)*lda];
		  double akkii= A[k+1+(i+1)*lda];
    		   C[k+j*lda] += aki * bij + akii *biij;
		   C[k+1+j*lda] += akki * bij + akkii *biij;
      		   C[k+(j+1)*lda] += aki * bijj+ akii * biijj;
		   C[k+1+(j+1)*lda] +=  akki * bijj + akkii * biijj;
		}

    	   }
		double bij = B[i+j*lda];
                double bijj = B[i+(j+1)*lda];
     	 	for (k = 0; k < K; k+=2){
		   double aki= A[k+i*lda];
		   double akki= A[k+1+i*lda];
     		   C[k+j*lda] += aki * bij;
		   C[k+1+j*lda] += akki * bij;
      		   C[k+(j+1)*lda] += aki * bijj;
		   C[k+1+(j+1)*lda] += akki * bijj;
		}

 	}
}
else if(M%2==0){//if M is even and N, K odd
	int i =0;
	int j =0;
	int k = 0;
  	for (j = 0; j < M-1; j+=2){
    	   for (i = 0; i < N-1; i+=2){
		double bij = B[i+j*lda];
                double biij = B[i+1+j*lda];
                double bijj = B[i+(j+1)*lda];
                double biijj = B[i+1+(j+1)*lda];
     	 	for (k = 0; k < K-1; k+=2){
		  double aki= A[k+i*lda];
		  double akki= A[k+1+i*lda];
		  double akii= A[k+(i+1)*lda];
		  double akkii= A[k+1+(i+1)*lda];
    		   C[k+j*lda] += aki * bij + akii *biij;
		   C[k+1+j*lda] += akki * bij + akkii *biij;
      		   C[k+(j+1)*lda] += aki * bijj+ akii * biijj;
		   C[k+1+(j+1)*lda] +=  akki * bijj + akkii * biijj;
		}
		   double aki= A[k+i*lda];
		   double akii= A[k+(i+1)*lda];
                   C[k+j*lda] += aki * bij + akii *biij;
                   C[k+(j+1)*lda] += aki * bijj + akii * biijj ;

    	   }
		double bij = B[i+j*lda];
                double bijj = B[i+(j+1)*lda];
                for (k = 0; k < K-1; k+=2){
                   double aki= A[k+i*lda];
                   double akki= A[k+1+i*lda];
                   C[k+j*lda] += aki * bij;
                   C[k+1+j*lda] += akki * bij;
                   C[k+(j+1)*lda] += aki * bijj;
                   C[k+1+(j+1)*lda] += akki * bijj;
                }
                   double aki= A[k+i*lda];
                   C[k+j*lda] += aki * bij;
                   C[k+(j+1)*lda] += aki * bijj;

 	}
}
else if(N%2==0&&K%2==0){//if N,K even and M is odd
	int i =0;
	int j =0;
	int k = 0;
  	for (j = 0; j < M-1; j+=2){
    	   for (i = 0; i < N-1; i+=2){
		double bij = B[i+j*lda];
                double biij = B[i+1+j*lda];
                double bijj = B[i+(j+1)*lda];
                double biijj = B[i+1+(j+1)*lda];
     	 	for (k = 0; k < K; k+=2){
		  double aki= A[k+i*lda];
		  double akki= A[k+1+i*lda];
		  double akii= A[k+(i+1)*lda];
		  double akkii= A[k+1+(i+1)*lda];
    		   C[k+j*lda] += aki * bij + akii *biij;
		   C[k+1+j*lda] += akki * bij + akkii *biij;
      		   C[k+(j+1)*lda] += aki * bijj+ akii * biijj;
		   C[k+1+(j+1)*lda] +=  akki * bijj + akkii * biijj;
		}

    	   }
 	}
  	   for (i = 0; i < N-1; i+=2){
                double bij = B[i+j*lda];
                double biij = B[i+1+j*lda];
                for (k = 0; k < K; k+=2){
		   double aki= A[k+i*lda];
		   double akki=A[k+1+i*lda];
		   double akii=A[k+(i+1)*lda];
		   double akkii=A[k+1+(i+1)*lda];
                   C[k+j*lda] += aki * bij + akii*biij;
                   C[k+1+j*lda] += akkii*biij+akki * bij ;
                }

           }

}
else if(N%2==0){//if N is even and M, K are odd
	int i =0;
	int j =0;
	int k = 0;
  	for (j = 0; j < M-1; j+=2){
    	   for (i = 0; i < N-1; i+=2){
		double bij = B[i+j*lda];
                double biij = B[i+1+j*lda];
                double bijj = B[i+(j+1)*lda];
                double biijj = B[i+1+(j+1)*lda];
     	 	for (k = 0; k < K-1; k+=2){
		  double aki= A[k+i*lda];
		  double akki= A[k+1+i*lda];
		  double akii= A[k+(i+1)*lda];
		  double akkii= A[k+1+(i+1)*lda];
    		   C[k+j*lda] += aki * bij + akii *biij;
		   C[k+1+j*lda] += akki * bij + akkii *biij;
      		   C[k+(j+1)*lda] += aki * bijj+ akii * biijj;
		   C[k+1+(j+1)*lda] +=  akki * bijj + akkii * biijj;
		}
		   double aki= A[k+i*lda];
		   double akii=A[k+(i+1)*lda];
                   C[k+j*lda] += aki * bij + akii *biij;
                   C[k+(j+1)*lda] +=aki * bijj + akii * biijj ;

    	   }
 	}
	  for (i = 0; i < N-1; i+=2){
                double bij = B[i+j*lda];
                double biij = B[i+1+j*lda];
                for (k = 0; k < K-1; k+=2){
		   double aki=A[k+i*lda];
		   double akki=A[k+1+i*lda];
		   double akii=A[k+(i+1)*lda];
		   double akkii=A[k+1+(i+1)*lda];
                   C[k+j*lda] += aki * bij + akii *biij;
                   C[k+1+j*lda] += akki * bij + akkii *biij;
                }
                   C[k+j*lda] += A[k+i*lda] * bij + A[k+(i+1)*lda] *biij;

	 }	
}else if(K%2==0){//if M and N are odd but K is even
	int i =0;
	int j =0;
	int k = 0;
  	for (j = 0; j < M-1; j+=2){
    	   for (i = 0; i < N-1; i+=2){
		double bij = B[i+j*lda];
                double biij = B[i+1+j*lda];
                double bijj = B[i+(j+1)*lda];
                double biijj = B[i+1+(j+1)*lda];
     	 	for (k = 0; k < K; k+=2){
		  double aki= A[k+i*lda];
		  double akki= A[k+1+i*lda];
		  double akii= A[k+(i+1)*lda];
		  double akkii= A[k+1+(i+1)*lda];
    		   C[k+j*lda] += aki * bij + akii *biij;
		   C[k+1+j*lda] += akki * bij + akkii *biij;
      		   C[k+(j+1)*lda] += aki * bijj+ akii * biijj;
		   C[k+1+(j+1)*lda] +=  akki * bijj + akkii * biijj;
		}

    	   }
		double bij = B[i+j*lda];
                double bijj = B[i+(j+1)*lda];
     	 	for (k = 0; k < K; k+=2){
		   double aki= A[k+i*lda];
		   double akki= A[k+1+i*lda];
     		   C[k+j*lda] += aki * bij;
		   C[k+1+j*lda] += akki * bij;
      		   C[k+(j+1)*lda] += aki * bijj;
		   C[k+1+(j+1)*lda] += akki * bijj;
		}

 	}
           for (i = 0; i < N-1; i+=2){
                double bij = B[i+j*lda];
                double biij = B[i+1+j*lda];
                for (k = 0; k < K; k+=2){
		   double aki=A[k+i*lda];
		   double akki=A[k+1+i*lda];
		   double akii=A[k+(i+1)*lda];
		   double akkii=A[k+1+(i+1)*lda];
                   C[k+j*lda] += aki * bij + akii *biij;
                   C[k+1+j*lda] += akki* bij + akkii *biij;
                }

           }
                double bij = B[i+j*lda];
                for (k = 0; k < K; k+=2){
                   C[k+j*lda] += A[k+i*lda] * bij;
                   C[k+1+j*lda] += A[k+1+i*lda] * bij;
                }

}else{//if M, N and K are all odd
	int i =0;
	int j =0;
	int k = 0;
	
  	for (j = 0; j < M-1; j+=2){
    	   for (i = 0; i < N-1; i+=2){
		double bij = B[i+j*lda];
                double biij = B[i+1+j*lda];
                double bijj = B[i+(j+1)*lda];
                double biijj = B[i+1+(j+1)*lda];
     	 	for (k = 0; k < K-1; k+=2){
		  double aki= A[k+i*lda];
		  double akki= A[k+1+i*lda];
		  double akii= A[k+(i+1)*lda];
		  double akkii= A[k+1+(i+1)*lda];
    		   C[k+j*lda] += aki * bij + akii *biij;
		   C[k+1+j*lda] += akki * bij + akkii *biij;
      		   C[k+(j+1)*lda] += aki * bijj+ akii * biijj;
		   C[k+1+(j+1)*lda] +=  akki * bijj + akkii * biijj;
		}
		   double aki= A[k+i*lda];
		   double akii=A[k+(i+1)*lda];
                   C[k+j*lda] += aki * bij + akii *biij;
                   C[k+(j+1)*lda] += aki * bijj + akii * biijj;

    	   }
		double bij = B[i+j*lda];
                double bijj = B[i+(j+1)*lda];
                for (k = 0; k < K-1; k+=2){
                   double aki= A[k+i*lda];
                   double akki= A[k+1+i*lda];
                   C[k+j*lda] += aki * bij;
                   C[k+1+j*lda] += akki * bij;
                   C[k+(j+1)*lda] += aki * bijj;
                   C[k+1+(j+1)*lda] += akki * bijj;
                }
                   double aki= A[k+i*lda];
                   C[k+j*lda] += aki * bij;
                   C[k+(j+1)*lda] += aki * bijj;

 	}
           for (i = 0; i < N-1; i+=2){
                double bij = B[i+j*lda];
                double biij = B[i+1+j*lda];
                for (k = 0; k < K-1; k+=2){
		   double aki=A[k+i*lda];
		   double akki=A[k+1+i*lda];
		   double akii=A[k+(i+1)*lda];
		   double akkii=A[k+1+(i+1)*lda];
                   C[k+j*lda] += aki * bij + akii *biij;
                   C[k+1+j*lda] += akki * bij + akkii *biij;
                }
                   C[k+j*lda] += A[k+i*lda] * bij + A[k+(i+1)*lda] *biij;

           }
                double bij = B[i+j*lda];
                for (k = 0; k < K-1; k+=2){
                   C[k+j*lda] += A[k+i*lda] * bij;
                   C[k+1+j*lda] += A[k+1+i*lda] * bij;
                }
                   C[k+j*lda] += A[k+i*lda] * bij;

}
}

void square_dgemm (int lda, double* A, double* B, double* C)
{

/* For each column of B and each row of A */ 
  for (int j = 0; j < lda; j += BLOCK_SIZE)
    /* For each column of B and each column of C */
    for (int i = 0; i < lda; i += BLOCK_SIZE)
      /* For each row of A and each row of C */
      for (int k = 0; k < lda; k += BLOCK_SIZE)
      {
	/* Correct block dimensions if block "goes off edge of" the matrix */
	int M = min (BLOCK_SIZE, lda-j);//Number of columns is B and number of rows in A
	int N = min (BLOCK_SIZE, lda-i);//Number of column in B and number of columns in C
	int K = min (BLOCK_SIZE, lda-k);//Number of rows in B and number of row in C
	
	// Perform individual block dgemm
	//Holding B and iterating through A and C blocks
	do_block_two_x_two_x_two_tiles (lda, M, N, K, A + k + i*lda, B + i + j*lda, C + k + j*lda);

	}
}
