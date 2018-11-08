#include <iostream>
#include <random>
#include <chrono>
using namespace std;

#include <blitz/array.h>
typedef blitz::Array<double,4> MatD4;
typedef blitz::Array<double,3> MatD3;
typedef blitz::Array<double,2> MatD2;
typedef blitz::Array<double,1> VecD;


int main()
{
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator (seed);
	std::uniform_real_distribution<double> distribution(0.0,1.0);

	blitz::firstIndex  i;
	blitz::secondIndex j;
	blitz::thirdIndex  k;
	blitz::fourthIndex l;
	
	MatD2 B(3,3);
	MatD2 Iden(3,3);
	MatD4 BTensor(3,3,3,3);
	MatD4 BTensor2(3,3,3,3);
	
	Iden = 1,0,0,
	       0,1,0,
		   0,0,1;
	cout << Iden << endl;	   
	for ( int i = 0; i < 3; ++i ) {
		for ( int j = 0; j < 3; ++j ) {
			B(i,j) = distribution(generator);			
		}
	}
	
	BTensor = Iden(i,k)*B(j,l)+Iden(j,k)*B(i,l);
	
	// compute B_ijkl
	// B_ijkl = delta_ik B^e,tr_jl + delta_jk B^e,tr_il (see Eq. 14.102, pg. 598 in complas peric)
	for ( int i = 0; i < 3; ++i ) {
		for ( int j = 0; j < 3; ++j ) {
			for ( int k = 0; k < 3; ++k ) {
				for ( int l = 0; l < 3; ++l ) {
					BTensor2(i,j,k,l) = Iden(i,k)*B(j,l)+Iden(j,k)*B(i,l);
				}
			}
		}
	}
	
	cout << "BTensor  = " << endl << BTensor << endl;
	cout << "BTensor2 = " << endl << BTensor2 << endl;
	
	MatD2 test1(3,3);
	MatD2 test2(3,3);
	test1 = blitz::sum(B(i,k)*Iden(k,j),k);
	
	for ( int i = 0; i < 3; ++i ) {
		for ( int j = 0; j < 3; ++j ) {
			for ( int k = 0; k < 3; ++k ) {
				test2(i,j) += B(i,k)*Iden(k,j);
			}
		}
	}
	cout << "test1 = " << endl << test1 << endl;
	cout << "test2 = " << endl << test2 << endl;
	
	
}

