/* 

   Recon Code for Cartesian and Non-Cartesian Data
   (in process of reorganizing for template + 4D recons)


 */
			 

#include "omp.h"

#include <iostream>
#include <fstream>
#include <recon_lib.h>
#include <armadillo>

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
	Array< float, 3> kspace(N,shots,1,ColumnMajorArray<3>());
	Array< float, 3> kw(N,shots,1,ColumnMajorArray<3>());
	kw = 1.0;
	
	Array< float, 3> kwN(N,shots,1,ColumnMajorArray<3>());
	kwN = 1.0;
	
	

	float initial_energy;
	float best_error;

	// Setup Gridding + FFT Structure
	gridFFT gridding;
	gridding.read_commandline(argc,argv);
	gridding.precalc_gridding(256,256,256,THREED,THREEDNONCARTESIAN);

	Array< float, 3> X(256,256,256,ColumnMajorArray<3>());
	
	for(int iter =0; iter< 1; iter++){
		
		// Get K-space Trajectory
		for(int shot=0; shot<shots; shot++){
			double ux = xdir[shot];
			double uy = ydir[shot];
			double uz = zdir[shot];

			double phi  = acos( uz /sqrt( uz*uz + uy*uy + ux*ux)); 

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
				kspace(pos,shot,0) = 1.0;
				kw(pos,shot,0) = kw_line[pos];
				kwN(pos,shot,0) = -kw_line[pos];

			}
		}	
		
		// Get density
		X =0;
		gridding.grid_forward( X, kspace, kxspace, kyspace,kzspace);
				
		//gridding.grid_backward(X, kspace, kw, kxspace, kyspace,kzspace);

		if(iter==0){ 
			initial_energy = sum( abs(X) );
		}
		
		float current_error = sum(abs(X));
		if( (current_error < best_error) || (iter==0)){
			best_error = current_error;
			for(int shot=0; shot<shots; shot++){
				omega[shot] = omega_test[shot];
			}
			cout << "Iter " << iter << "::Energy = " << (best_error/initial_energy) << endl;
		}
	}

#ifdef KJHKJH


	// ---- Optimization 
	int test_shots = 1000;

	array3D< float >kxTest ;
	kxTest.alloc(1,test_shots,N);
	array3D< float >kyTest;
	kyTest.alloc(1,test_shots,N);
	array3D< float >kzTest;
	kzTest.alloc(1,test_shots,N);

	array3D< float >kTest;
	kTest.alloc(1,test_shots,N);


	array3D< float >kOnes;
	kOnes.alloc(1,shots,N);
	kOnes = 1.0;

	array3D< float >kPhantom;
	kPhantom.alloc(1,shots,N);
	kPhantom = 1.0;

	cout << "Start optimization" << endl;
	gridding.k3d_grid.zero();
	gridding.grid_forward( kOnes[0][0], kxspace[0][0], kyspace[0][0],kzspace[0][0],kw[0][0],N*shots);
	X = gridding.k3d_grid;

	if(test_shots > shots){
		kw.freeArray();
		kwN.freeArray();
		kw.alloc(1,test_shots,N);
		kwN.alloc(1,test_shots,N);
		for(int shot=0; shot<test_shots; shot++){
			for(int pos = 0; pos<N; pos++){
				kw[0][shot][pos] = kw_line[pos];
				kwN[0][shot][pos] = -kw_line[pos];
			}}		
	}

	
	tictoc T;
	initial_energy = X.Menergy();
	for(int iter=0; iter< shots*20; iter++){
		
		//cout << "Iter = " << iter << endl;
		
		// Find worst shot
		// int shot = rand()%shots; // rand()%shots;
		int shot;
		double worst_cost = 0;
		for(int test=0; test<shots; test++){
			double current_cost = 0.0;
			for(int pos = 0; pos<N; pos++){
				current_cost += abs(kspace[0][test][pos]);
			}
			if( (current_cost > worst_cost) || (test==0)){
				worst_cost  =current_cost;
				shot = test;
			}
		}
		
		// Degrid the shot
		T.tic();
		omp_set_num_threads(16);
		gridding.grid_forward( kOnes[0][shot], kxspace[0][shot], kyspace[0][shot],kzspace[0][shot],kwN[0][shot],N);
		//cout << "\n\tForward = " << T << endl;

		T.tic();
		// Setup Test Angles		
		double ux = xdir[shot];
		double uy = ydir[shot];
		double uz = zdir[shot];

		mat Rnet(3,3);
		vec k(3);
		vec kR(3);
		for(int test=0; test<test_shots; test++){
			double test_omega = omega[shot] + (double)test/((double)test_shots-1.0)*2.0*PI;
			Rnet = get_rotation_matrix(ux,uy,uz,test_omega);

			for(int pos = 0; pos<N; pos++){
				k(0) = kx[pos];
				k(1) = ky[pos];
				k(2) = kz[pos];
				kR = Rnet*k;					

				kxTest[0][test][pos]=kR(0);
				kyTest[0][test][pos]=kR(1);
				kzTest[0][test][pos]=kR(2);
			}
		}
		//cout << "\tRotate = " << T << endl;


		// Grid backwards
		T.tic();
		omp_set_num_threads(1);
		gridding.grid_backward( kTest[0][0], kxTest[0][0], kyTest[0][0],kzTest[0][0],kw[0][0],test_shots*N);
		//cout << "\tBackwards = " << T << endl;

		//Find Min point
		T.tic();
		double best_cost = 0;
		int best_omega;
		for(int test=0; test<test_shots; test++){

			double current_cost = 0.0;
			for(int pos = 0; pos<N; pos++){
				current_cost += abs(kTest[0][test][pos]);
			}

			if( (current_cost < best_cost) || (test==0)){
				best_cost  =current_cost;
				best_omega = test;
			}
		}
		omega[shot] = omega[shot] + (double)best_omega/((double)test_shots-1.0)*2.0*PI;
		for(int pos = 0; pos<N; pos++){
			kxspace[0][shot][pos]=kxTest[0][best_omega][pos];
			kyspace[0][shot][pos]=kyTest[0][best_omega][pos];
			kzspace[0][shot][pos]=kzTest[0][best_omega][pos];
		}
		//cout << "\tMax = " << T << endl;
		
		// cout << " Best omega = " <<  best_omega << endl;

		// Grid back onto 
		T.tic();
		omp_set_num_threads(1);
		gridding.grid_forward( kOnes[0][shot], kxspace[0][shot], kyspace[0][shot],kzspace[0][shot],kw[0][shot],N);
		gridding.grid_backward( kspace[0][shot], kxspace[0][shot], kyspace[0][shot],kzspace[0][shot],kw[0][shot],N);
		//cout << "\tRegrid = " << T << endl;				
		
		// Get cost
		if( iter%100==0){
			omp_set_num_threads(16);

			T.tic();
			gridding.k3d_grid.zero();
			gridding.grid_forward( kOnes[0][0], kxspace[0][0], kyspace[0][0],kzspace[0][0],kw[0][0],N*shots);
			gridding.grid_backward( kspace[0][0], kxspace[0][0], kyspace[0][0],kzspace[0][0],kw[0][0],N*shots);
			X = gridding.k3d_grid;

			cout << "Iter " << iter << "::Image Energy = " << (X.Menergy()/initial_energy) << endl;
			cout << "Iter " << iter << "::Kspace Energy = " << (kspace.Menergy()/initial_energy) << endl;
			X.write_mag("Kspace_Slice.dat",X.Nz/2,"a+");

			kxspace.write("KMAPX_VD_0.dat");
			kyspace.write("KMAPY_VD_0.dat");
			kzspace.write("KMAPZ_VD_0.dat");
			kw.write("KWEIGHT.dat");
			
			{
				FILE *fid;
				fid = fopen("OmegaOpt.dat","w");
				fwrite( omega,shots, sizeof(float), fid);
     			fclose(fid); 
				
			}
			
		}
	}	



	kspace.write_mag("kw.dat");
#endif
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


