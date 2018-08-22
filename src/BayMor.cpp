#include <RcppArmadillo.h>
#include <omp.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]

			

double exp_rs(double a, double b)
{
  double  z, w, hr;
  hr = 1/a;
   z = R::rexp(hr);
   while(z > (b-a)) z = R::rexp(hr);
   w = R::runif(0.0, 1.0);

   while( log(w) > (-0.5*z*z))
   {
      z = R::rexp(hr);
      while(z > (b-a)) z = R::rexp(hr);
      w = R::runif(0.0,1.0);
   }
   return(z+a);
}


// [[Rcpp::export]]
double nors(double a, double b)
{
   double  x;
   x = Rf_rnorm(0.0, 1.0);
   while( (x < a) || (x > b) ) x = norm_rand();
   return x;
}

// [[Rcpp::export]]
double urs(double a, double b)
{
   double ar, lor, x, lu;
   if(a <= 0.0) ar = 0.0;
   else ar = a;
   lor = R::dnorm(ar, 0.0, 1.0, 1.0);

   x = R::runif(a, b);
   lu = log(R::runif(0.0, 1.0));
   while( lu > (R::dnorm(x, 0.0, 1.0,1.0) - lor))
   {
      x = R::runif(a, b);
      lu = log(R::runif(0.0, 1.0));
   }
   return x;
}

// [[Rcpp::export]]
double hnors(double a, double b)
{
   double   x;
   x = fabs(norm_rand());
   while( (x<a) || (x>b) ) x = fabs(norm_rand());
   return x;
}


// [[Rcpp::export]]
double rtruncnorm(double mu, double sigma, double lower, double upper)
{
int mr;
 double a, b;
 double po = log(0.150), logt2 = log(2.18), t3 = 0.725;
 double z, tmp, loat;

 mr = 0;
 a = (lower - mu)/sigma;
 b = (upper - mu)/sigma;

 if( (a == R_NegInf) || (b == R_PosInf))
   {
     if(a == R_NegInf)
       {
     mr = 1;
     a = -b;
     b = R_PosInf;
       }
     if(a <= 0.45) z = nors(a, b);
     else z = exp_rs(a, b);
     if(mr) z = -z;
   }
 else if((a * b) <= 0.0)
   {
     if((R::dnorm(a, 0.0, 1.0,1.0) <= po) || (R::dnorm(b, 0.0, 1.0, 1.0) <= po))
       {
     z = nors(a, b);
       }
     else z = urs(a,b);
   }
 else
   {
     if(b < 0)
       {
     tmp = b; b = -a; a = -tmp; mr = 1;
       }
     loat = R::dnorm(a, 0.0, 1.0, 1.0) - R::dnorm(b, 0.0, 1.0, 1.0);
     if(loat <= logt2) z = urs(a,b);
     else if((loat > po) && (a < t3)) z = hnors(a,b);
     else z = exp_rs(a,b);
     if(mr) z = -z;
   }
   double output;
   output = sigma*z + mu;
 return (output);
}



#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat mor(arma::mat Y, arma::rowvec l, Rcpp::List X, arma::rowvec nx, int Nit) 
{
			
                       
             
			int N=Y.n_rows ; 
			
			
			
			int m= Y.n_cols ;
			

			
			double Ml=max(l);   

			arma::mat a(m,(Ml+1),arma::fill::zeros);
	
			for (int i=0; i<m; i++){
		
				a(i,0)=-10;
				a(i,1)=1.5;
				a(i,l(i))= 10;
				a(i,l(i)-1)=l(i)-0.5;
		
					if(l(i)>3) { for ( int k=2;k<=l(i)-2;k++) {
		                   a(i,k)=k+0.5; }

								}
									}
			
		
			
			l.insert_cols(0,1);
			
			arma::mat Idy(sum(l), N, arma::fill::zeros); 
			 
			arma::rowvec nid(sum(l), arma::fill::zeros);   
			
			 
			
	    for (int i=0 ; i<m; i++)  {   
		    
			for (int k=1; k<=l(i+1);k++ ){
	 

				arma::uvec  rr= find(Y.col(i)==k);   
				arma::rowvec   r= arma::conv_to<arma::rowvec>::from(rr);
				
				int lr= r.n_elem;   
				nid.col((sum(l.cols(0,i))+k-1))=lr; 
				
				if(lr>0) {
				Idy.submat(sum(l.cols(0,i))+k-1,0,sum(l.cols(0,i))+k-1,lr-1)=r; 
                          }

					}
			}	
			
	
			
			
			    nx.for_each( [](arma::rowvec::elem_type& val) { val += 1; } ); 
			    nx.insert_cols(0,1); 
			
			arma::mat D((m*N), sum(nx),arma::fill::zeros);    
			for (int j=1;j <= m; j++) {
				arma::mat dd=X[j-1];
				arma::mat d= arma::join_rows(arma::ones(N,1),dd);
				D(arma::span( (N*(j-1)),((N*j)-1)), arma::span((accu(nx.cols(0,j-1))), (accu(nx.cols(0,j))-1)))= d;    
								   	}
			
	double  ss, up=1, lo=1;

			
			arma::mat Beta(Nit,D.n_cols,arma::fill::zeros);   
			
			
			arma::mat Z(N,m); 
			Z=Y; 
			
		

            arma::mat bhat ;
            arma::mat bv;
            arma::mat b;
            arma::mat E;
            arma::mat Iprcs =arma::eye(m,m);
			arma::mat iN = arma::eye(N,N);	

			
	for (int it=1; it<Nit; it++){		
			
			arma::mat Ii=arma::kron (Iprcs, iN);
			
			bv = inv(D.t()*Ii* D)  ;  	 	
            bhat = bv* D.t() * Ii * arma:: vectorise(Z)   ;
            b = mvnrnd(bhat, bv);
            arma::rowvec  bt= b.t();
            Beta.row(it) =  bt ;  
			
			
            arma::mat db=  D*b;
			db.set_size(N,m); 
            E = Z - db;  
            Iprcs = arma::wishrnd(inv(E.t()*E),N);    

			for (int j=0; j<m; j++) { 
				for (int i=0; i<N; i++) { 

					arma::mat sIp = Iprcs;
					arma::mat sZ= Z;
					arma::mat sD= D;
					sIp.shed_col(j); 
					sZ.shed_col(j);				
					sD.shed_rows(N*j,(N*(j+1))-1); 
			arma::mat	uj= sD* b;  
			            uj.set_size(N, m-1); 
			arma::mat	zu =( (sZ.row(i)) - (uj.row(i))); 
					arma::mat kj= sIp.row(j);
			       
   	        arma::mat    mu = (D.row(i+(j*N))*b) -(( 1/Iprcs(j,j))* (kj *(zu.t())));
					      ss= 1/Iprcs(j,j) ;
				 
			double muu=mu(0,0);
		
						for (int k=1; k<=l(j+1); k++) { 
						     if(Y(i,j)==k) { up =a(j,k);  lo=a(j,(k-1)); break; }
							
													 }
							 
					Z(i,j)= rtruncnorm(muu, ss, lo, up);
						
  		                           }
		                        }

		for (int i=0 ; i<m; i++) {
	
				if (l(i+1) >3){
		
					for (int j=2 ; j<l(i+1)-1; j++) {
		            arma::mat vMax; double zMax;
						if ( (nid(sum(l.cols(0,i))+j-1)) && (nid(sum(l.cols(0,i))+j))!=0){
					
						    vMax= Idy.submat(sum(l.cols(0,i))+j-1,0,sum(l.cols(0,i))+j-1, nid(sum(l.cols(0,i))+j-1)-1);
							arma::uvec   VMa = arma::conv_to<arma::uvec>::from(vMax);
							arma::mat dma=Z.col(i);
							arma::vec ddma=dma.rows(VMa);
							zMax=max(ddma);
							 
  
							arma::mat vMin= Idy.submat(sum(l.cols(0,i))+j,0,sum(l.cols(0,i))+j, nid(sum(l.cols(0,i))+j)-1);
							arma::uvec   VMi = arma::conv_to<arma::uvec>::from(vMin);
							arma::mat dmi=Z.col(i);
							arma::vec ddmi=dmi.rows(VMi);
							double zmin=min(ddmi); 
   
					a(i,j)=R::runif(zMax,zmin);                                            
																						 }		
													}
							}
								} 	
			
			
			
		


								}

return(Beta);
}



