#include <iostream>
using namespace std;

#include <blitz/array.h>
BZ_USING_NAMESPACE(blitz)

typedef Array<double,2> MatD2
typedef Array<double,1> VecD

class Element
{
public:
	VecD coord;
	
}
	
void onem(Array<double,2> A)
{
	A = 1,0,0,
	0,1,0,
	0,0,1;
	
}

void modify_subarray(Array<double,2> A)
{
	A(Range::all(),0) = 0.5;
	A(Range::all(),1) = -0.5;

}

int main()
{
	int nen = 4;
	int ndm = 2;
	Array<double,3> Bmat(3,nen,ndm);
	Array<double,2> F(3,3);
	Array<double,1> vec(3);
	
	vec = 1.333333333;
	vec(0)=-3;
	cout << "vec = " << endl;
	cout << vec << endl;
	
	Bmat = 1.0;
	F = 2.0;
	
	cout << "F before\n";
	cout << F << endl;
	
	onem(F);
	
	cout << "F after\n";
	cout << F << endl;
	
	cout << "Bmat before\n";
	cout << Bmat << endl;
	
	Array<double,2> BmatView = Bmat(2,Range::all(),Range::all());
	modify_subarray(BmatView);
	
	cout << "Bmat after\n";
	cout << Bmat << endl;
		
		
}

