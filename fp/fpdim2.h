#ifndef fpdim2h
#define fpdim2h

#include <iostream>

#include "fpmath.h"

namespace fp
{

    template <class type>
    class dim2
    {
    public:
        type x, y;
    };
    
    template <typename type>
	bool operator==(const dim2<type>& lhs, const dim2<type>& rhs);
    
	template <typename type>
	bool operator!=(const dim2<type>& lhs, const dim2<type>& rhs);
	
	template <typename type>
	bool operator<(const dim2<type>& lhs, const dim2<type>& rhs);
	
	template <typename type>
	bool operator>(const dim2<type>& lhs, const dim2<type>& rhs);
	
	template <typename type>
	bool operator<=(const dim2<type>& lhs, const dim2<type>& rhs);
	
	template <typename type>
	bool operator>=(const dim2<type>& lhs, const dim2<type>& rhs);
	
	template <typename type>
	type operator*(const dim2<type>& lhs, const dim2<type>& rhs);
	
	template <typename type>
	dim2<type> operator+(const dim2<type>& lhs, const dim2<type>& rhs);
	
	template <typename type>
	dim2<type> operator*(const dim2<type>& lhs, const type& rhs);
	
	template <typename type>
	dim2<type> operator*(const type& lhs, const dim2<type>& rhs);
	
	template <typename type>
	dim2<type> operator*=(dim2<type>& lhs, const type& rhs);
	
	template <typename type>
	dim2<double> operator/(const dim2<type>& lhs, const double& rhs);
	
	template <typename type>
	dim2<type> operator/=(dim2<type>& lhs, const double& rhs);
	
	template <typename type>
	dim2<type> operator-(const dim2<type>& lhs, const dim2<type>& rhs);
	
	template <typename type>
	dim2<type> operator-(const dim2<type>& rhs);
	
	template <typename type>
	dim2<type> operator+=(dim2<type>& lhs, const dim2<type>& rhs);
	
	template <typename type>
	dim2<type> operator-=(dim2<type>& lhs, const dim2<type>& rhs);
	
	template <typename type>
	double norm(dim2<type> vec);
	
	template <typename type>
	type norm2(dim2<type> vec);
	
	template <typename type>
	dim2<type> operator/(const dim2<type>& lhs, const dim2<type>& rhs);
	
	template <typename type>
	dim2<type> operator/=(dim2<type>& lhs, const dim2<type>& rhs);
	
	template <typename type>
	dim2<type> operator%(const dim2<type>& lhs, const dim2<type>& rhs);
	
	template <typename type>
	dim2<type> operator%=(dim2<type>& lhs, const dim2<type>& rhs);

	template <typename type>
    std::ostream &operator<<(std::ostream &os, const dim2<type> &vec);
	
	template <typename type>
	double dist(dim2<type> vec1, dim2<type> vec2);
	
	template <typename type>
	type dist2(dim2<type> vec1, dim2<type> vec2);
	

    typedef dim2<int> int2;
    typedef dim2<double> double2;
    
}

#include "fpdim2.cpp"

#endif