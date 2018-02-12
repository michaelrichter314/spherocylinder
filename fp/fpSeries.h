#ifndef fpseries_h
#define fpseries_h

/*
	Enhanced array, optimized for speed. Fixed length.
*/

#include "fpMemory.h"

#ifdef debug
	#include <iostream>
#endif

namespace fp
{
	template <class type>
	class series
	{
		public:
			series(int l);
			~series();
			
			size_t length();
			type& operator[](const int &pos);
			//int& operator[] (const int nIndex);
			
		private:
			size_t N; //length
			type* data;
	};
	
	// TODO: this should be in cpp-file
	template <class type>
	series<type>::series(int l) : N(l)
	{
		#ifdef debug
			if (N < 1)
			{
				std::cerr << "fpSeries: Length of series has to be strictly positive.\n";
				return;
			}
		#endif
		
		data = (type*)std::malloc(N * sizeof(type));
		
		#ifdef debug
			if (data == NULL)
			{
				std::cerr << "fpSeries: Could not allocate memory (probably out of RAM).\n";
			}
		#endif
	}
	
	template <class type>
	series<type>::~series()
	{
		std::free(data);
	}
	
	template <class type>
	size_t series<type>::length()
	{
		return N;
	}


	template <class type>
	type& series<type>::operator[](const int &pos)
	{
			#ifdef debug
				if (pos < 0 || pos >= N)
				{
					std::cerr << "fpSeries: Could not find element #" << pos << ": Has to be between 0 and " << N-1 << ".\n";
				}
			#endif
		
		return data[pos];
	}
}
#endif
