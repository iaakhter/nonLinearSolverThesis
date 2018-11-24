#define BaseType double
#include "badiff.h"
#include "fadiff.h"

#include <vector>
#include "TestFADBAD.h"


// Functions to differentiate:

class Add
{
public:
	static const char* Name(){ return "simple addition"; }
	template <class T>
	std::vector<T> operator()(const std::vector<T>& v)
	{
		std::vector<T> out(1);
		out[0]=v[0]+v[1];
		return out;
	}
	static std::vector<double> point()
	{
		std::vector< double > in(2);
		in[0]=1.5;
		in[1]=-1.5;
		return in;
	}
};
class Sub
{
public:
	static const char* Name(){ return "simple subtraction"; }
	template <class T>
	std::vector<T> operator()(const std::vector<T>& v)
	{
		std::vector<T> out(1);
		out[0]=v[0]-v[1];
		return out;
	}
	static std::vector<double> point()
	{
		std::vector< double > in(2);
		in[0]=1.5;
		in[1]=-1.5;
		return in;
	}
};
class Mul
{
public:
	static const char* Name(){ return "simple multiplication"; }
	template <class T>
	std::vector<T> operator()(const std::vector<T>& v)
	{
		std::vector<T> out(1);
		out[0]=v[0]*v[1];
		return out;
	}
	static std::vector<double> point()
	{
		std::vector< double > in(2);
		in[0]=1.5;
		in[1]=-1.5;
		return in;
	}
};
class Div
{
public:
	static const char* Name(){ return "simple division"; }
	template <class T>
	std::vector<T> operator()(const std::vector<T>& v)
	{
		std::vector<T> out(1);
		out[0]=v[0]/v[1];
		return out;
	}
	static std::vector<double> point()
	{
		std::vector< double > in(2);
		in[0]=1.5;
		in[1]=-1.5;
		return in;
	}
};
class Pow
{
public:
	static const char* Name(){ return "simple power"; }
	template <class T>
	std::vector<T> operator()(const std::vector<T>& v)
	{
		std::vector<T> out(1);
		out[0]=pow(v[0],3);
		return out;
	}
	static std::vector<double> point()
	{
		std::vector< double > in(1);
		in[0]=1.5;
		return in;
	}
};
class Pow2
{
public:
	static const char* Name(){ return "power"; }
	template <class T>
	std::vector<T> operator()(const std::vector<T>& v)
	{
		std::vector<T> out(1);
		out[0]=pow(v[0],v[1]);
		return out;
	}
	static std::vector<double> point()
	{
		std::vector< double > in(2);
		in[0]=1.5;
		in[0]=3.3;
		return in;
	}
};
class Sqr
{
public:
	static const char* Name(){ return "simple square"; }
	template <class T>
	std::vector<T> operator()(const std::vector<T>& v)
	{
		std::vector<T> out(1);
		out[0]=sqr(v[0]);
		return out;
	}
	static std::vector<double> point()
	{
		std::vector< double > in(1);
		in[0]=1.5;
		return in;
	}
};
class Exp
{
public:
	static const char* Name(){ return "simple exponential"; }
	template <class T>
	std::vector<T> operator()(const std::vector<T>& v)
	{
		std::vector<T> out(1);
		out[0]=exp(v[0]);
		return out;
	}
	static std::vector<double> point()
	{
		std::vector< double > in(1);
		in[0]=1.5;
		return in;
	}
};
class Log
{
public:
	static const char* Name(){ return "simple logarithmic"; }
	template <class T>
	std::vector<T> operator()(const std::vector<T>& v)
	{
		std::vector<T> out(1);
		out[0]=log(v[0]);
		return out;
	}
	static std::vector<double> point()
	{
		std::vector< double > in(1);
		in[0]=1.5;
		return in;
	}
};
class Sqrt
{
public:
	static const char* Name(){ return "simple square-root"; }
	template <class T>
	std::vector<T> operator()(const std::vector<T>& v)
	{
		std::vector<T> out(1);
		out[0]=sqrt(v[0]);
		return out;
	}
	static std::vector<double> point()
	{
		std::vector< double > in(1);
		in[0]=1.5;
		return in;
	}
};
class Sin
{
public:
	static const char* Name(){ return "simple sine"; }
	template <class T>
	std::vector<T> operator()(const std::vector<T>& v)
	{
		std::vector<T> out(1);
		out[0]=sin(v[0]);
		return out;
	}
	static std::vector<double> point()
	{
		std::vector< double > in(1);
		in[0]=1.5;
		return in;
	}
};
class Cos
{
public:
	static const char* Name(){ return "simple cosine"; }
	template <class T>
	std::vector<T> operator()(const std::vector<T>& v)
	{
		std::vector<T> out(1);
		out[0]=cos(v[0]);
		return out;
	}
	static std::vector<double> point()
	{
		std::vector< double > in(1);
		in[0]=1.5;
		return in;
	}
};
class Tan
{
public:
	static const char* Name(){ return "simple tangent"; }
	template <class T>
	std::vector<T> operator()(const std::vector<T>& v)
	{
		std::vector<T> out(1);
		out[0]=tan(v[0]);
		return out;
	}
	static std::vector<double> point()
	{
		std::vector< double > in(1);
		in[0]=1.5;
		return in;
	}
};
class ASin
{
public:
	static const char* Name(){ return "simple arc sine"; }
	template <class T>
	std::vector<T> operator()(const std::vector<T>& v)
	{
		std::vector<T> out(1);
		out[0]=asin(v[0]);
		return out;
	}
	static std::vector<double> point()
	{
		std::vector< double > in(1);
		in[0]=1.5;
		return in;
	}
};
class ACos
{
public:
	static const char* Name(){ return "simple arc cosine"; }
	template <class T>
	std::vector<T> operator()(const std::vector<T>& v)
	{
		std::vector<T> out(1);
		out[0]=acos(v[0]);
		return out;
	}
	static std::vector<double> point()
	{
		std::vector< double > in(1);
		in[0]=1.5;
		return in;
	}
};
class ATan
{
public:
	static const char* Name(){ return "simple arc tangent"; }
	template <class T>
	std::vector<T> operator()(const std::vector<T>& v)
	{
		std::vector<T> out(1);
		out[0]=atan(v[0]);
		return out;
	}
	static std::vector<double> point()
	{
		std::vector< double > in(1);
		in[0]=1.5;
		return in;
	}
};



class csmap // The cos-sin map.
{
public:
	static const char* Name(){ return "Cos-Sine map"; }
	template <class T>
	std::vector<T> operator()(const std::vector<T>& v)
	{
		std::vector<T> out(2);
		out[0]=cos(v[0]+4*v[1]);
		out[1]=sin(4*v[0]+v[1]);
		return out;
	}
	static std::vector<double> point()
	{
		std::vector< double > in(2);
		in[0]=1.5;
		in[1]=-1.5;
		return in;
	}
};
class func1
{
public:
	static const char* Name(){ return "atan[(a sin x)/(1-a cos x)]"; }
	template <class T>
	std::vector<T> operator()(const std::vector<T>& v)
	{
		std::vector<T> out(1);
		T x(v[0]);
		T a(v[1]);
		out[0]=atan((a*sin(x))/(1-a*cos(x)));
		return out;
	}
	static std::vector<double> point()
	{
		std::vector< double > in(2);
		in[0]=1.5;
		in[1]=-.456;
		return in;
	}
};



// The classes Fdiff<C> and Bdiff are capable of differentiating a function
// f = [f0,f1,..,f(m-1)] : T^n->T^m by using either the forward or the backward
// modes of automatic differentiation. Where T is real, interval, etc.
// The argument when using the evaluating operation on Fdiff is an
// n-dimensional std::vector x in which the function should be evaluated
// and differentiated. The return value is a std::vector that should
// be interpreted as follows:
//		Vout[0..m-1] is the function value.
//		Vout[m..m+n-1] is the value of [df1/dx1, ..., df1/dxn]
//		Vout[m+i*n..m+(i+1)*n-1] is the value of [dfi/dx1, ..., dfi/dxn].

template <class C>
class Fdiff
{
public:
	template <class T>
	std::vector<T> operator()(const std::vector<T>& v)
	{
		int i,j,inSize=v.size();;
		std::vector< F<T> > out;
		std::vector< F<T> > in(inSize);
		for(i=0;i<inSize;++i)
		{
			in[i]=v[i];
			in[i].diff(i,inSize);
		}
		out=C()(in);
		int outSize=out.size();
		std::vector<T> retval(outSize*(1+inSize));
		for(j=0;j<outSize;++j)
		{
			retval[j]=out[j].x();
			for(i=0;i<inSize;++i)
			{
				retval[outSize+j*inSize+i]=out[j].d(i);
			}
		}
		return retval;
	}
};

template <class C>
class Bdiff
{
public:
	template <class T>
	std::vector<T> operator()(const std::vector<T>& v)
	{
		int i,j,inSize=v.size();
		std::vector< B<T> > out;
		std::vector< B<T> > in(inSize);
		for(i=0;i<inSize;++i)
		{
			in[i]=v[i];
		}
		out=C()(in);
		int outSize=out.size();
		for(j=0;j<outSize;++j)
		{
			out[j].diff(j,outSize);
		}
		std::vector<T> retval(outSize*(1+inSize));
		for(j=0;j<outSize;++j)
		{
			retval[j]=out[j].x();
			for(i=0;i<inSize;++i)
			{
				retval[outSize+j*inSize+i]=in[i].d(j);
			}
		}
		return retval;
	}
};

template <class T>
std::ostream& operator << (std::ostream& os, const std::vector<T>& v)
{
	os << "[";
	for(int i=0;i<v.size();i++)
	{
		os << v[i];
		if (i<v.size()-1) os << ",";
	}
	os << "]";
	return os;
}

template <class T>
class check
{
	public:
	void run(IReportLog& log);
};

bool diff(const std::vector<double> v1, const std::vector<double> v2)
{
	if (v1.size()!=v2.size()) return true;
	for(int i=0;i<v1.size();++i)
		if (fabs(v1[i]-v2[i])>1e-6) return true;
	return false;
}

template <class T>
void check<T>::run(IReportLog& log)
{
	std::vector<double> v1;
	std::vector<double> v2;

	{
		Fdiff< T > Fd;
		v1=Fd(T::point());
	}

	{
		Bdiff< T > Bd;
		v2=Bd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 1. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 1. order derivative succeeded for " << T::Name() << std::endl;

	{
		Fdiff< Fdiff < T > > FFd;
		v1=FFd(T::point());
	}

	{
		Bdiff< Fdiff < T > > BFd;
		v2=BFd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 2. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 2. order derivative succeeded for " << T::Name() << std::endl;

	{
		Fdiff< Bdiff < T > > FBd;
		v2=FBd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 2. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 2. order derivative succeeded for " << T::Name() << std::endl;

	{
		Bdiff< Bdiff < T > > BBd;
		v2=BBd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 2. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 2. order derivative succeeded for " << T::Name() << std::endl;

	{
		Fdiff < Fdiff< Fdiff < T > > > FFFd;
		v1=FFFd(T::point());
	}

	{
		Bdiff < Fdiff< Fdiff < T > > > BFFd;
		v2=BFFd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 3. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 3. order derivative succeeded for " << T::Name() << std::endl;

	{
		Fdiff < Bdiff< Fdiff < T > > > FBFd;
		v2=FBFd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 3. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 3. order derivative succeeded for " << T::Name() << std::endl;

	{
		Bdiff < Bdiff< Fdiff < T > > > BBFd;
		v2=BBFd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 3. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 3. order derivative succeeded for " << T::Name() << std::endl;

	{
		Fdiff < Fdiff< Bdiff < T > > > FFBd;
		v2=FFBd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 3. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 3. order derivative succeeded for " << T::Name() << std::endl;

	{
		Bdiff < Fdiff< Bdiff < T > > > BFBd;
		v2=BFBd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 3. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 3. order derivative succeeded for " << T::Name() << std::endl;

	{
		Fdiff < Bdiff< Bdiff < T > > > FBBd;
		v2=FBBd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 3. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 3. order derivative succeeded for " << T::Name() << std::endl;

	{
		Bdiff < Bdiff< Bdiff < T > > > BBBd;
		v2=BBBd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 3. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 3. order derivative succeeded for " << T::Name() << std::endl;
/* Uncomment this if you want to wait even longer for compilation to finish
	{
		Fdiff < Fdiff < Fdiff< Fdiff < T > > > > FFFFd;
		v1=FFFFd(T::point());
	}

	{
		Bdiff < Fdiff < Fdiff< Fdiff < T > > > > BFFFd;
		v2=BFFFd(T::point());
	}
	if (diff(v1,v2)) log.failed() << "Testing 4. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 4. order derivative succeeded for " << T::Name() << std::endl;

	{
		Fdiff < Bdiff < Fdiff< Fdiff < T > > > > FBFFd;
		v2=FBFFd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 4. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 4. order derivative succeeded for " << T::Name() << std::endl;

	{
		Bdiff < Bdiff < Fdiff< Fdiff < T > > > > BBFFd;
		v2=BBFFd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 4. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 4. order derivative succeeded for " << T::Name() << std::endl;

	{
		Fdiff < Fdiff < Bdiff< Fdiff < T > > > > FFBFd;
		v2=FFBFd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 4. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 4. order derivative succeeded for " << T::Name() << std::endl;

	{
		Bdiff < Fdiff < Bdiff< Fdiff < T > > > > BFBFd;
		v2=BFBFd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 4. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 4. order derivative succeeded for " << T::Name() << std::endl;

	{
		Fdiff < Bdiff < Bdiff< Fdiff < T > > > > FBBFd;
		v2=FBBFd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 4. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 4. order derivative succeeded for " << T::Name() << std::endl;

	{
		Bdiff < Bdiff < Bdiff< Fdiff < T > > > > BBBFd;
		v2=BBBFd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 4. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 4. order derivative succeeded for " << T::Name() << std::endl;

	{
		Fdiff < Fdiff < Fdiff< Bdiff < T > > > > FFFBd;
		v2=FFFBd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 4. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 4. order derivative succeeded for " << T::Name() << std::endl;

	{
		Bdiff < Fdiff < Fdiff< Bdiff < T > > > > BFFBd;
		v2=BFFBd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 4. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 4. order derivative succeeded for " << T::Name() << std::endl;

	{
		Fdiff < Bdiff < Fdiff< Bdiff < T > > > > FBFBd;
		v2=FBFBd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 4. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 4. order derivative succeeded for " << T::Name() << std::endl;

	{
		Bdiff < Bdiff < Fdiff< Bdiff < T > > > > BBFBd;
		v2=BBFBd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 4. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 4. order derivative succeeded for " << T::Name() << std::endl;

	{
		Fdiff < Fdiff < Bdiff< Bdiff < T > > > > FFBBd;
		v2=FFBBd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 4. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 4. order derivative succeeded for " << T::Name() << std::endl;

	{
		Bdiff < Fdiff < Bdiff< Bdiff < T > > > > BFBBd;
		v2=BFBBd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 4. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 4. order derivative succeeded for " << T::Name() << std::endl;

	{
		Fdiff < Bdiff < Bdiff< Bdiff < T > > > > FBBBd;
		v2=FBBBd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 4. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 4. order derivative succeeded for " << T::Name() << std::endl;

	{
		Bdiff < Bdiff < Bdiff< Bdiff < T > > > > BBBBd;
		v2=BBBBd(T::point());
	}

	if (diff(v1,v2)) log.failed() << "Testing 4. order derivative failed for " << T::Name() << std::endl;
	else log.succeeded() << "Testing 4. order derivative succeeded for " << T::Name() << std::endl;
*/
}


void TestFADBAD::run(IReportLog& rlog)
{
	{ check<Add> chk; chk.run(rlog); }
	{ check<Sub> chk; chk.run(rlog); }
	{ check<Mul> chk; chk.run(rlog); }
	{ check<Div> chk; chk.run(rlog); }
	{ check<Pow> chk; chk.run(rlog); }
	{ check<Pow2> chk; chk.run(rlog); }
	{ check<Sqr> chk; chk.run(rlog); }
	{ check<Exp> chk; chk.run(rlog); }
	{ check<Log> chk; chk.run(rlog); }
	{ check<Sqrt> chk; chk.run(rlog); }
	{ check<Sin> chk; chk.run(rlog); }
	{ check<Cos> chk; chk.run(rlog); }
	{ check<Tan> chk; chk.run(rlog); }
	{ check<ASin> chk; chk.run(rlog); }
	{ check<ACos> chk; chk.run(rlog); }
	{ check<ATan> chk; chk.run(rlog); }
	{ check<csmap> chk; chk.run(rlog); }
	{ check<func1> chk; chk.run(rlog); }
}
