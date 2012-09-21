template<int dim>
class Base
{
public:
	Base();
	~Base();
	virtual int return_int(){return 321;};
};

template<int dim>
class Derived : public Base<dim>
{
public:
	Derived();
	~Derived();
	virtual int return_int(){ return 123;}
};

template<int dim>
class DerivedDerived : public Derived<dim>
{
public:
	DerivedDerived();
	~DerivedDerived();
	virtual int return_int(){ return 456;}
};