#include "fpdim2.h"

namespace fp
{
    template <typename type>
    bool operator==(const dim2<type>& lhs, const dim2<type>& rhs)
	{
		return (lhs.x == rhs.x) && (lhs.y == rhs.y);
	}
	
	template <typename type>
	bool operator!=(const dim2<type>& lhs, const dim2<type>& rhs)
	{
		return (lhs.x != rhs.x) || (lhs.y != rhs.y);
	}
	
	template <typename type>
	bool operator<(const dim2<type>& lhs, const dim2<type>& rhs)
	{
		return lhs.x < rhs.x && lhs.y < rhs.y;
	}
	
	template <typename type>
	bool operator>(const dim2<type>& lhs, const dim2<type>& rhs)
	{
		return lhs.x > rhs.x && lhs.y > rhs.y;
	}
	
	template <typename type>
	bool operator<=(const dim2<type>& lhs, const dim2<type>& rhs)
	{
		return lhs.x <= rhs.x && lhs.y <= rhs.y;
	}
	
	template <typename type>
	bool operator>=(const dim2<type>& lhs, const dim2<type>& rhs)
	{
		return lhs.x >= rhs.x && lhs.y >= rhs.y;
	}
	
	template <typename type>
	type operator*(const dim2<type>& lhs, const dim2<type>& rhs)
	{
		return lhs.x*rhs.x + lhs.y*rhs.y;
	}
	
	template <typename type>
	dim2<type> operator+(const dim2<type>& lhs, const dim2<type>& rhs)
	{
		dim2<type> toReturn;
		toReturn.x = lhs.x + rhs.x;
		toReturn.y = lhs.y + rhs.y;
		return toReturn;
	}
	
	template <typename type>
	dim2<type> operator*(const dim2<type>& lhs, const type& rhs)
	{
		dim2<type> toReturn;
		toReturn.x = lhs.x * rhs;
		toReturn.y = lhs.y * rhs;
		return toReturn;
	}
	
	template <typename type>
	dim2<type> operator*(const type& lhs, const dim2<type>& rhs)
	{
		dim2<type> toReturn;
		toReturn.x = rhs.x * lhs;
		toReturn.y = rhs.y * lhs;
		return toReturn;
	}
	
	template <typename type>
	dim2<type> operator*=(dim2<type>& lhs, const type& rhs)
	{
		lhs.x *= rhs;
		lhs.y *= rhs;
		return lhs;
	}
	
	template <typename type>
	dim2<double> operator/(const dim2<type>& lhs, const double& rhs)
	{
		dim2<type> toReturn;
		toReturn.x = lhs.x / rhs;
		toReturn.y = lhs.y / rhs;
		return toReturn;
	}
	
	template <typename type>
	dim2<type> operator/=(dim2<type>& lhs, const double& rhs)
	{
		lhs.x /= rhs;
		lhs.y /= rhs;
		return lhs;
	}
	
	template <typename type>
	dim2<type> operator-(const dim2<type>& lhs, const dim2<type>& rhs)
	{
		dim2<type> toReturn;
		toReturn.x = lhs.x - rhs.x;
		toReturn.y = lhs.y - rhs.y;
		return toReturn;
	}
	
	template <typename type>
	dim2<type> operator-(const dim2<type>& rhs)
	{
		dim2<type> toReturn;
		toReturn.x = -rhs.x;
		toReturn.y = -rhs.y;
		return toReturn;
	}
	
	template <typename type>
	dim2<type> operator+=(dim2<type>& lhs, const dim2<type>& rhs)
	{
		lhs.x += rhs.x;
		lhs.y += rhs.y;
		return lhs;
	}
	
	template <typename type>
	dim2<type> operator-=(dim2<type>& lhs, const dim2<type>& rhs)
	{
		lhs.x -= rhs.x;
		lhs.y -= rhs.y;
		return lhs;
	}
	
	template <typename type>
	double norm(dim2<type> vec)
	{
        return std::sqrt(vec.x*vec.x + vec.y*vec.y);
	}
	
	template <typename type>
	// square of norm
	type norm2(dim2<type> vec)
	{
		return vec.x*vec.x + vec.y*vec.y;
	}
	
	template <typename type>
	dim2<type> operator/(const dim2<type>& lhs, const dim2<type>& rhs)
	{
		dim2<type> toReturn;
		double normRhs = norm(rhs);
		toReturn = lhs/normRhs;
		return toReturn;
	}
	
	template <typename type>
	dim2<type> operator/=(dim2<type>& lhs, const dim2<type>& rhs)
	{
		double normRhs = norm(rhs);
		lhs/=normRhs;
		return lhs;
	}
	
	template <typename type>
	dim2<type> operator%(const dim2<type>& lhs, const dim2<type>& rhs)
	{
		dim2<type> toReturn = lhs;
		toReturn.x %= rhs.x;
		toReturn.y %= rhs.y;
		return toReturn;
	}
	
	template <typename type>
	dim2<type> operator%=(dim2<type>& lhs, const dim2<type>& rhs)
	{
		lhs.x %= rhs.x;
		lhs.y %= rhs.y;
		return lhs;
	}

	template <typename type>
    std::ostream &operator<<(std::ostream &os, const dim2<type> &vec)
    {
        return os << "(" << vec.x << ", " << vec.y << ")";
    }
	
	template <typename type>
	double dist(dim2<type> vec1, dim2<type> vec2)
	{
		return norm(vec2 - vec1);
	}
	
	template <typename type>
	type dist2(dim2<type> vec1, dim2<type> vec2)
	{
		return norm2(vec2-vec1);
	}
}