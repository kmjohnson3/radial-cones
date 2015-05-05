/* 

   Recon Code for Cartesian and Non-Cartesian Data
   (in process of reorganizing for template + 4D recons)


 */
			 

#include "omp.h"

#include <iostream>
#include <fstream>
#include <armadillo>

#include "gridFFT.h"
#include "ArrayTemplates.hpp"
#include "tictoc.hpp"

using namespace std;
using arma::mat;
using arma::vec;
using namespace NDarray;

double gettime(void);
void temp_recon(int argc, char **argv);
mat get_rotation_matrix(float ux,float uy, float uz,float omega);
float randn(void);

int main(int argc, char **argv){

	// ------------------------------------
	// Setup Recon
	// ------------------------------------
	cout<< "Start your engines" << endl << flush;
	FILE *fid = fopen("RadialCones.struct","r");

	int N = 0;
	fread(&N,sizeof(int),1,fid);
	cout << "K-space Uses " << N << " points" << endl;
	float *kx = new float[N];
	float *ky = new float[N];
	float *kz = new float[N];
	float *kw_line = new float[N];
	fread(kx,sizeof(float),N,fid);	
	fread(ky,sizeof(float),N,fid);	
	fread(kz,sizeof(float),N,fid);	
	fread(kw_line,sizeof(float),N,fid);	

	int shots = 0;
	fread(&shots,sizeof(int),1,fid);
	cout << "K-space Uses " << shots << " shots" << endl;
	float *xdir = new float[shots];
	float *ydir = new float[shots];
	float *zdir = new float[shots];
	float *omega = new float[shots];
	float *omega_test = new float[shots];
	fread(xdir,sizeof(float),shots,fid);	
	fread(ydir,sizeof(float),shots,fid);	
	fread(zdir,sizeof(float),shots,fid);	
	fclose(fid);

	// Final Image Solution
	Array< float, 3> kxspace(N,shots,1,ColumnMajorArray<3>());
	Array< float, 3> kyspace(N,shots,1,ColumnMajorArray<3>());
	Array< float, 3> kzspace(N,shots,1,ColumnMajorArray<3>());
	Array< float, 3> kw(N,shots,1,ColumnMajorArray<3>());
	Array< float, 3> kspace( N,shots,1,ColumnMajorArray<3>());

	float initial_energy=0.0;
	float best_error;

	// Setup Gridding + FFT Structure
	gridFFT gridding;
	gridding.kernel_type =TRIANGLE_KERNEL;
	gridding.dwinX = 1;
	gridding.dwinY = 1;
	gridding.dwinZ = 1;
	gridding.grid_x = 1;
	gridding.grid_y = 1;
	gridding.grid_z = 1;
	gridding.grid_in_x = 1;
	gridding.grid_in_y = 1;
	gridding.grid_in_z = 1;
	gridding.precalc_kernel(); 

	// Get K-space Trajectory
	for(int shot=0; shot<shots; shot++){
		double ux = xdir[shot];
		double uy = ydir[shot];
		double uz = zdir[shot];

		omega_test[shot] = 2.0*3.14156*( (float)rand() / (float)RAND_MAX);

		mat Rnet(3,3);
		Rnet = get_rotation_matrix(ux,uy,uz,omega_test[shot]);

		for(int pos = 0; pos<N; pos++){
			vec k(3);
			k(0) = kx[pos];
			k(1) = ky[pos];
			k(2) = kz[pos];
			vec kR = Rnet*k;					

			kxspace(pos,shot,0)=kR(0);
			kyspace(pos,shot,0)=kR(1);
			kzspace(pos,shot,0)=kR(2);
			kw(pos,shot,0) = kw_line[pos];
		}
	}	

	int sizeX = 8*(int)(0.5+ (max(kxspace) - min(kxspace) ) / 8 );
	int sizeY = 8*(int)(0.5+ (max(kyspace) - min(kyspace) ) / 8 );
	int sizeZ = 8*(int)(0.5+ (max(kzspace) - min(kzspace) ) / 8 );
	cout << "Size is " << sizeX << " x " << sizeY << " x " << sizeZ << endl;
	Array< float, 3> X(sizeX, sizeY, sizeZ,ColumnMajorArray<3>());
		

	// ---- Optimization 
	int test_shots = 1000;
	Array< float, 3> kxTest(N,test_shots,1,ColumnMajorArray<3>());
	Array< float, 3> kyTest(N,test_shots,1,ColumnMajorArray<3>());
	Array< float, 3> kzTest(N,test_shots,1,ColumnMajorArray<3>());
	Array< float, 3> kTest(N,test_shots,1,ColumnMajorArray<3>());
	Array< float, 3> kwTest(N,test_shots,1,ColumnMajorArray<3>());


	// ALlocated arrays for single shot
	Array< float,3>kx_shot(N,1,1,ColumnMajorArray<3>());
	Array< float,3>ky_shot(N,1,1,ColumnMajorArray<3>());
	Array< float,3>kz_shot(N,1,1,ColumnMajorArray<3>());
	Array< float,3>kw_shot(N,1,1,ColumnMajorArray<3>());

	// Set Densities
	kw_shot(Range::all(),0,0) = kw(Range::all(),0,0);
	for( int shot = 0; shot < test_shots; shot++){
		kwTest(Range::all(),shot,0) = kw(Range::all(),0,0);
	}

	// These are data values
	Array< float, 3> kOnes(1,shots,N,ColumnMajorArray<3>());
	Array< float, 3> kPhantom(1,shots,N,ColumnMajorArray<3>());
	kPhantom = 1.0;
	kOnes = 1.0;


	cout << "Start optimization" << endl;
	
	// Get Initial Density
	X =0;
	gridding.grid_forward( X, kw, kxspace, kyspace,kzspace);
	gridding.grid_backward( X, kspace, kxspace, kyspace,kzspace);
	
	bool check_psf = true;
	tictoc T;
	arma::uvec indices;
	for(int iter=0; iter< shots*20; iter++){
		
		//cout << "Iter = " << iter << endl;
		
		// Get cost
		if( iter%100==0){
			omp_set_num_threads(omp_get_max_threads());
		
			T.tic();
			X =0;
			gridding.grid_forward( X, kw, kxspace, kyspace,kzspace);
			
			if(iter==0){
				initial_energy = sum(sqr(X));
			}

			// Get Density estimate
			kspace = 0.0;
			gridding.grid_backward( X, kspace, kxspace, kyspace,kzspace);

			// Look for worst shots
			vec cost(shots);
			for(int shot=0; shot<shots; shot++){
				cost(shot) = 0.0;
				for(int pos =0; pos<N; pos++){
					double temp = (double)kspace(pos,shot,0);
					cost(shot) += temp*temp;
				}
			}
		 	indices = sort_index( cost,"descend");


			if( check_psf){
				Array< complex<float>,3>PSF(sizeX,sizeY,sizeZ,ColumnMajorArray<3>());
				for(int k=0; k< PSF.length(thirdDim); k++){
					for(int j=0; j< PSF.length(secondDim); j++){
						for(int i=0; i< PSF.length(firstDim); i++){
							PSF(i,j,k) = X(i,j,k);
						}
					}
				}
				//fftshift(PSF);
				fft(PSF);
				{
					Array< complex<float>,2> Pslice = PSF(Range::all(),Range::all(),(int)(PSF.length(thirdDim)/2)); 
					ArrayWriteMagAppend(Pslice,"PSF_Slice.dat");
				}
			}

			cout << "Iter " << iter << "::Image Energy = " << (sum(sqr(X))/initial_energy) << endl;
		}


		int shot = indices(iter%100);

		
		// Degrid the shot
		T.tic();
		{
			kx_shot(Range::all(),0,0) = kxspace(Range::all(),shot,0);
			ky_shot(Range::all(),0,0) = kyspace(Range::all(),shot,0);
			kz_shot(Range::all(),0,0) = kzspace(Range::all(),shot,0);
			kw_shot(Range::all(),0,0) = kw(Range::all(),shot,0);
			kw_shot *= -1;
			omp_set_num_threads(1);
			gridding.grid_forward( X, kw_shot, kx_shot, ky_shot,kz_shot);
		}
		//cout << "\tSubtraction took " << T << endl;


		
		T.tic();
		// Setup Test Angles		
		double ux = xdir[shot];
		double uy = ydir[shot];
		double uz = zdir[shot];

		for(int test=0; test<test_shots; test++){
			mat Rnet(3,3);
			vec k(3);
			vec kR(3);

			double test_omega = omega[shot] + (double)test/((double)test_shots-1.0)*2.0*PI;
			Rnet = get_rotation_matrix(ux,uy,uz,test_omega);

			for(int pos = 0; pos<N; pos++){
				k(0) = kx[pos];
				k(1) = ky[pos];
				k(2) = kz[pos];
				
				kR = Rnet*k;					

				kxTest(pos,test,0)=kR(0);
				kyTest(pos,test,0)=kR(1);
				kzTest(pos,test,0)=kR(2);
			}
		}
		//cout << "\tRotate = " << T << endl;


		// Grid backwards
		T.tic();
		omp_set_num_threads(omp_get_max_threads());
		gridding.grid_backward(X, kTest, kxTest, kyTest,kzTest);
		//cout << "\tBackwards = " << T << endl;

		//Find Min point
		T.tic();
		double best_cost = 0;
		int best_omega;
		for(int test=0; test<test_shots; test++){

			double current_cost = 0.0;
			for(int pos = 0; pos<N; pos++){
				current_cost += abs(kTest(pos,test,0));
			}

			if( (current_cost < best_cost) || (test==0)){
				best_cost  =current_cost;
				best_omega = test;
			}
		}
		
		// Update Shot
		omega[shot] = omega[shot] + (double)best_omega/((double)test_shots-1.0)*2.0*PI;
		for(int pos = 0; pos<N; pos++){
			kxspace(pos,shot,0)=kxTest(pos,best_omega,0);
			kyspace(pos,shot,0)=kyTest(pos,best_omega,0);
			kzspace(pos,shot,0)=kzTest(pos,best_omega,0);
		}
		//cout << "\tMax = " << T << endl;
		
		// cout << " Best omega = " <<  best_omega << endl;

		// Grid back onto 
		T.tic();
		omp_set_num_threads(1);
		{
			kx_shot(Range::all(),0,0) = kxspace(Range::all(),shot,0);
			ky_shot(Range::all(),0,0) = kyspace(Range::all(),shot,0);
			kz_shot(Range::all(),0,0) = kzspace(Range::all(),shot,0);
			kw_shot(Range::all(),0,0) = kw(Range::all(),shot,0);
			gridding.grid_forward( X, kw_shot, kx_shot, ky_shot,kz_shot);
		}
	
		//gridding.grid_forward( kOnes[0][shot], kxspace[0][shot], kyspace[0][shot],kzspace[0][shot],kw[0][shot],N);
		//gridding.grid_backward( kspace[0][shot], kxspace[0][shot], kyspace[0][shot],kzspace[0][shot],kw[0][shot],N);
		//cout << "\tRegrid = " << T << endl;				

		{
			Array< float,2> Xslice = X(Range::all(),Range::all(),(int)(X.length(thirdDim)/2)); 
			ArrayWriteAppend(Xslice,"X.dat");
		}


	}	

	
	
	fid = fopen("Optimal.rotation","w");
	Array< float,3>Roptimal(3,3,shots,ColumnMajorArray<3>());
	for(int shot=0; shot<shots; shot++){
		
		// Setup Test Angles		
		double ux = xdir[shot];
		double uy = ydir[shot];
		double uz = zdir[shot];
		double om = omega[shot];
		
		mat Rnet(3,3);
		Rnet = get_rotation_matrix(ux,uy,uz,omega[shot]);
		for(int i=0; i<3; i++){
			for(int j=0; j<3; j++){
				float temp = Rnet(i,j);
				fwrite(&temp,sizeof(float),1,fid); 
			}
		}
	
	}
	fclose(fid);


	return(0);
}

float randn(void){
	int counts = 20;

	float val=0.0;
	for(float pos=0; pos< counts; pos++){
		val += (float)rand() - ((float)RAND_MAX)/2.0;
	}	
	val /=  ( ((float)RAND_MAX)/sqrt(12.0) )*sqrt( (float)counts);
	return(val);
}

mat get_rotation_matrix(float ux,float uy, float uz,float omega){

	mat Rphi(3,3);
	mat Rtheta(3,3);
	mat Rutu(3,3);
	mat Rusk(3,3);
	mat Raxis(3,3);
	mat Rnet(3,3);
	mat I;
	I.eye(3,3);

	double phi  = acos( uz /sqrt( uz*uz + uy*uy + ux*ux)); 
	double theta= PI/2 + atan2( uy,ux);

	// Rotation about X
	Rphi(0,0)= 1.0; 	Rphi(0,1)= 0.0;			Rphi(0,2)= 0.0;
	Rphi(1,0)= 0.0; 	Rphi(1,1)= cos(phi);	Rphi(1,2)= -sin(phi);				
	Rphi(2,0)= 0.0; 	Rphi(2,1)= sin(phi);	Rphi(2,2)=  cos(phi);				

	// About Z 		
	Rtheta(0,0)= cos(theta); 	Rtheta(0,1)= -sin(theta);		Rtheta(0,2)= 0.0;
	Rtheta(1,0)= sin(theta); 	Rtheta(1,1)=  cos(theta);		Rtheta(1,2)= 0.0;				
	Rtheta(2,0)= 0.0; 			Rtheta(2,1)= 0.0;				Rtheta(2,2)= 1.0;				

	// About Axis
	Rutu(0,0)= ux*ux; 	Rutu(0,1)= uy*ux;		Rutu(0,2)= uz*ux;
	Rutu(1,0)= ux*uy; 	Rutu(1,1)= uy*uy;		Rutu(1,2)= uz*uy;				
	Rutu(2,0)= ux*uz; 	Rutu(2,1)= uy*uz;		Rutu(2,2)= uz*uz;				

	Rusk(0,0)= 0; 	Rusk(0,1)= -uz;		Rusk(0,2)= uy;
	Rusk(1,0)= uz; 	Rusk(1,1)= 0;		Rusk(1,2)= -ux;				
	Rusk(2,0)= -uy; Rusk(2,1)= ux;		Rusk(2,2)= 0;				

	// cout << "Omega = " << omega << endl;
	Raxis = Rutu + cos(omega)*( I - Rutu) + sin(omega)*( Rusk );		

	Rnet = Raxis*( Rtheta*Rphi );

	return(Rnet);
}


