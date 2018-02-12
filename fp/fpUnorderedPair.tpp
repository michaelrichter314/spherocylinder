namespace fp
{

	template<typename T>
	auto UnorderedPair<T>::has(T p) const
	{
		return ( a == p || b == p );
	}


	template<typename T>
	UnorderedPair<T>::UnorderedPair(T a_, T b_) : a(a_), b(b_) {}

	template<typename T>
	auto UnorderedPair<T>::operator==(const UnorderedPair<T>& p)
	{
		return p.has(a) && p.has(b);
	}

	template<typename T>
	auto UnorderedPair<T>::operator!=(const UnorderedPair<T>& p)
	{
		return !(p.has(a)) || !(p.has(b));
	}
	
	template<typename T>
	T UnorderedPair<T>::other(T p) const
	{
		if (p == a) {return b;}
		if (p == b) {return a;}
		throw 100;
	}
	
	template<typename T>
	T UnorderedPair<T>::getA() const {return a;}
	
	template<typename T>
	T UnorderedPair<T>::getB() const {return b;}
	
	
}
