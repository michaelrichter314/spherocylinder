#ifndef fpunorderedpair_h
#define fpunorderedpair_h

namespace fp
{
	template<typename T>
	class UnorderedPair
	{
		public:
			UnorderedPair(T, T);
			
			auto has(T) const;
			T other(T) const;
			
			T getA() const;
			T getB() const;
			
			auto operator==(const UnorderedPair<T>& p);
			auto operator!=(const UnorderedPair<T>& p);
			
			
		
		protected:
			T a;
			T b;
	};
}

#include "fpUnorderedPair.tpp"

#endif
