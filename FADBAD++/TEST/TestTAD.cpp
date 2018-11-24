#define BaseType double
#include "tadiff.h"

#include "TestTAD.h"

#include <iostream>

using namespace std;

template <class F, int ORDER>
class DIntegral
{
	F m_func;                       // Function to integrate
public:
	DIntegral():m_func(){}
	DIntegral(const F& f):m_func(f){}
	template <class aType>
	aType operator ()(const aType a, const aType b)
	{
		T< aType > x,fx;
		aType h2((b-a)/2),res;
		x=(a+b)/2;              // Point of evaluation
		x[1]=1;                 // dx/dx=1
		fx=m_func(x);           // Record function
		fx.eval(ORDER);         // Compute up to ORDER coefficients
		res=0;
		// Evaluate the integral from a to b of the Taylor-polynomial.
		for(int i=ORDER;i>=0;i--)
			res+=fx[i]*(pow(h2,i+1)-pow(-h2,i+1))/(i+1);
		return res;
	}
};

template <class F, int ORDER=5, int N=10>
class Integral
{
	DIntegral<F,ORDER> m_dInt;
public:
	Integral():m_dInt(){}
	Integral(const F& f):m_dInt(f){}
	template <class aType>
	aType operator ()(const aType a, const aType b)
	{
		int i,j;
		aType res=0;
		for(i=0;i<N;i++)
		{
			j=i+1;
			res+=m_dInt((a*(N-i)+b*i)/N,(a*(N-j)+b*j)/N);
		}
		return res;
	}

};

class F1
{
public:
	template <class aType>
	aType operator ()(const aType& x)
	{
		aType exx(exp(x*x));
		return 2*x*exx*sin(exx);
	}
};

class F2
{
	double m_a;
public:
	F2(const double& a):m_a(a){}
	template <class aType>
	aType operator ()(const aType& x)
	{
		return 1/sqrt(m_a*m_a-x*x);
	}
};

class F3
{
	double m_a;
public:
	F3(const double& a):m_a(a){}
	template <class aType>
	aType operator ()(const aType& x)
	{
		return sqrt(m_a*m_a-x*x);
	}
};

class F4
{
	double m_a;
public:
	F4(const double& a):m_a(a){}
	template <class aType>
	aType operator ()(const aType& x)
	{
		return 1/(1-2*m_a*cos(x)+m_a*m_a);
	}
};

class F5
{
public:
	template <class aType>
	aType operator ()(const aType& x)
	{
		return log(1-x)/x;
	}
};

void TestTAD::run(IReportLog& rlog)
{
	double numerical, analytical;

	numerical=Integral<F1,6,160>()(0.0,2.0);
	analytical=cos(1.0)-cos(exp(4.0));
//	cout << "Numerical:  " << numerical << endl;
//	cout << "Analytical: " << analytical << endl;
	if (fabs(numerical-analytical)>0.01)
		rlog.failed() << "Integral of F1 failed" << endl;
	else
		rlog.succeeded() << "Integral of F1 succeeded" << endl;


	double a=5;
	F2 f2(a);
	numerical=Integral<F2,10,500>(f2)(0.0,a);
	analytical=3.141592653589793/2;
//	cout << "Numerical:  " << numerical << endl;
//	cout << "Analytical: " << analytical << endl;
	if (fabs(numerical-analytical)>0.01)
		rlog.failed() << "Integral of F2 failed" << endl;
	else
		rlog.succeeded() << "Integral of F2 succeeded" << endl;

	F3 f3(a);
	numerical=Integral<F3,6,160>(f3)(0.0,a);
	analytical=3.141592653589793*a*a/4;
//	cout << "Numerical:  " << numerical << endl;
//	cout << "Analytical: " << analytical << endl;
	if (fabs(numerical-analytical)>0.01)
		rlog.failed() << "Integral of F3 failed" << endl;
	else
		rlog.succeeded() << "Integral of F3 succeeded" << endl;

	a=0.4;
	F4 f4(a);
	numerical=Integral<F4,6,160>(f4)(0.0,2*3.141592653589793);
	analytical=2*3.141592653589793/(1-a*a);
//	cout << "Numerical:  " << numerical << endl;
//	cout << "Analytical: " << analytical << endl;
	if (fabs(numerical-analytical)>0.01)
		rlog.failed() << "Integral of F4 failed" << endl;
	else
		rlog.succeeded() << "Integral of F4 succeeded" << endl;

	numerical=Integral<F5,6,160>()(0.0,1.0);
	analytical=-pow(3.141592653589793,2)/6;
//	cout << "Numerical:  " << numerical << endl;
//	cout << "Analytical: " << analytical << endl;
	if (fabs(numerical-analytical)>0.01)
		rlog.failed() << "Integral of F5 failed" << endl;
	else
		rlog.succeeded() << "Integral of F5 succeeded" << endl;


}







