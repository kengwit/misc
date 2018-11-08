#include <iostream>
#include <vector>
using namespace std;
#include <Eigen/Eigenvalues>

struct eigen { 
    int index;
	double value;
	double vector[3];
    bool operator<(eigen const &other) const { 
        return value < other.value;
    }
};

int main()
{
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es;
	Eigen::Matrix3d X = Eigen::Matrix3d::Random(3,3);
	Eigen::Matrix3d A = X + X.transpose();
	es.compute(A);

	vector<eigen> EigArray(3);
	
	cout << "es.eigenvectors().col(0)[2] = " << endl << es.eigenvectors().col(0)[2] << endl;
	
	for ( int i = 0; i < 3; ++i )
	{
		EigArray[i].index = i;
		EigArray[i].value = es.eigenvalues()[i];
		for ( int j = 0; j < 3; ++j )
		{
			EigArray[i].vector[j] = es.eigenvectors().col(i)[j];
		}
	}
	
	//EigArray[0].value = -1;
	//EigArray[1].value = 3;
	//EigArray[2].value = 1;
	
	
	// sort based on ascending order of eigenvalues
	sort(EigArray.begin(),EigArray.end());
	
	// print using sorted indices
	cout << "sorted eigenvalues: " << endl;
	for ( int i = 0; i < 3; ++i )
	{
		cout << EigArray[i].index << ", " << EigArray[i].value << ", eigevec: (" 
		     << EigArray[i].vector[0] << "," <<  EigArray[i].vector[1] << "," << EigArray[i].vector[2] << ")\n";
	}
}

