#ifndef fpdim3h
#define fpdim3h

#include <iostream>
#include <vector>

#include "fpmath.h"

namespace fp
{

    template <class type>
    class dim3
    {
    public:
        type x, y, z;
    };
    
    template <typename type>
	bool operator==(const dim3<type>& lhs, const dim3<type>& rhs);
    
	template <typename type>
	bool operator!=(const dim3<type>& lhs, const dim3<type>& rhs);
	
	template <typename type>
	bool operator<(const dim3<type>& lhs, const dim3<type>& rhs);
	
	template <typename type>
	bool operator>(const dim3<type>& lhs, const dim3<type>& rhs);
	
	template <typename type>
	bool operator<=(const dim3<type>& lhs, const dim3<type>& rhs);
	
	template <typename type>
	bool operator>=(const dim3<type>& lhs, const dim3<type>& rhs);
	
	template <typename type>
	type operator*(const dim3<type>& lhs, const dim3<type>& rhs);
	
	template <typename type>
	dim3<type> operator+(const dim3<type>& lhs, const dim3<type>& rhs);
	
	template <typename type>
	dim3<type> operator*(const dim3<type>& lhs, const type& rhs);
	
	template <typename type>
	dim3<type> operator*(const type& lhs, const dim3<type>& rhs);
	
	template <typename type>
	dim3<type> operator*=(dim3<type>& lhs, const type& rhs);
	
	template <typename type>
	dim3<double> operator/(const dim3<type>& lhs, const double& rhs);
	
	template <typename type>
	dim3<type> operator/=(dim3<type>& lhs, const double& rhs);
	
	template <typename type>
	dim3<type> operator-(const dim3<type>& lhs, const dim3<type>& rhs);
	
	template <typename type>
	dim3<type> operator-(const dim3<type>& rhs);
	
	template <typename type>
	dim3<type> operator+=(dim3<type>& lhs, const dim3<type>& rhs);
	
	template <typename type>
	dim3<type> operator-=(dim3<type>& lhs, const dim3<type>& rhs);
	
	template <typename type>
	double norm(dim3<type> vec);
	
	template <typename type>
	type norm2(dim3<type> vec);
	
	template <typename type>
	dim3<type> operator/(const dim3<type>& lhs, const dim3<type>& rhs);
	
	template <typename type>
	dim3<type> operator/=(dim3<type>& lhs, const dim3<type>& rhs);
	
	template <typename type>
	dim3<type> operator%(const dim3<type>& lhs, const dim3<type>& rhs);
	
	template <typename type>
	dim3<type> operator%=(dim3<type>& lhs, const dim3<type>& rhs);

	template <typename type>
    std::ostream &operator<<(std::ostream &os, const dim3<type> &vec);
	
	template <typename type>
	double dist(dim3<type> vec1, dim3<type> vec2);
	
	template <typename type>
	type dist2(dim3<type> vec1, dim3<type> vec2);
	
	template <typename type>
	dim3<type> cross(const dim3<type>& a, const dim3<type>& b);
	

    typedef dim3<int> int3;
    typedef dim3<double> double3;
	
    template <typename type>
    bool operator==(const dim3<type>& lhs, const dim3<type>& rhs)
	{
		return (lhs.x == rhs.x) && (lhs.y == rhs.y) && (lhs.z == rhs.z);
	}
	
	template <typename type>
	bool operator!=(const dim3<type>& lhs, const dim3<type>& rhs)
	{
		return (lhs.x != rhs.x) || (lhs.y != rhs.y) || (lhs.z != rhs.z);
	}
	
	template <typename type>
	bool operator<(const dim3<type>& lhs, const dim3<type>& rhs)
	{
		return lhs.x < rhs.x && lhs.y < rhs.y && lhs.z < rhs.z;
	}
	
	template <typename type>
	bool operator>(const dim3<type>& lhs, const dim3<type>& rhs)
	{
		return lhs.x > rhs.x && lhs.y > rhs.y && lhs.z > rhs.z;
	}
	
	template <typename type>
	bool operator<=(const dim3<type>& lhs, const dim3<type>& rhs)
	{
		return lhs.x <= rhs.x && lhs.y <= rhs.y && lhs.z <= rhs.z;
	}
	
	template <typename type>
	bool operator>=(const dim3<type>& lhs, const dim3<type>& rhs)
	{
		return lhs.x >= rhs.x && lhs.y >= rhs.y && lhs.z >= rhs.z;
	}
	
	template <typename type>
	type operator*(const dim3<type>& lhs, const dim3<type>& rhs)
	{
		return lhs.x*rhs.x + lhs.y*rhs.y + lhs.z*rhs.z;
	}
	
	template <typename type>
	dim3<type> operator+(const dim3<type>& lhs, const dim3<type>& rhs)
	{
		dim3<type> toReturn;
		toReturn.x = lhs.x + rhs.x;
		toReturn.y = lhs.y + rhs.y;
        toReturn.z = lhs.z + rhs.z;
		return toReturn;
	}
	
	template <typename type>
	dim3<type> operator*(const dim3<type>& lhs, const type& rhs)
	{
		dim3<type> toReturn;
		toReturn.x = lhs.x * rhs;
		toReturn.y = lhs.y * rhs;
        toReturn.z = lhs.z * rhs;
		return toReturn;
	}
	
	template <typename type>
	dim3<type> operator*(const type& lhs, const dim3<type>& rhs)
	{
		dim3<type> toReturn;
		toReturn.x = rhs.x * lhs;
		toReturn.y = rhs.y * lhs;
        toReturn.z = rhs.z * lhs;
		return toReturn;
	}
	
	template <typename type>
	dim3<type> operator*=(dim3<type>& lhs, const type& rhs)
	{
		lhs.x *= rhs;
		lhs.y *= rhs;
        lhs.z *= rhs;
		return lhs;
	}
	
	template <typename type>
	dim3<double> operator/(const dim3<type>& lhs, const double& rhs)
	{
		dim3<type> toReturn;
		toReturn.x = lhs.x / rhs;
		toReturn.y = lhs.y / rhs;
        toReturn.z = lhs.z / rhs;
		return toReturn;
	}
	
	template <typename type>
	dim3<type> operator/=(dim3<type>& lhs, const double& rhs)
	{
		lhs.x /= rhs;
		lhs.y /= rhs;
        lhs.z /= rhs;
		return lhs;
	}
	
	template <typename type>
	dim3<type> operator-(const dim3<type>& lhs, const dim3<type>& rhs)
	{
		dim3<type> toReturn;
		toReturn.x = lhs.x - rhs.x;
		toReturn.y = lhs.y - rhs.y;
        toReturn.z = lhs.z - rhs.z;
		return toReturn;
	}
	
	template <typename type>
	dim3<type> operator-(const dim3<type>& rhs)
	{
		dim3<type> toReturn;
		toReturn.x = -rhs.x;
		toReturn.y = -rhs.y;
        toReturn.z = -rhs.z;
		return toReturn;
	}
	
	template <typename type>
	dim3<type> operator+=(dim3<type>& lhs, const dim3<type>& rhs)
	{
		lhs.x += rhs.x;
		lhs.y += rhs.y;
        lhs.z += rhs.z;
		return lhs;
	}
	
	template <typename type>
	dim3<type> operator-=(dim3<type>& lhs, const dim3<type>& rhs)
	{
		lhs.x -= rhs.x;
		lhs.y -= rhs.y;
        lhs.z -= rhs.z;
		return lhs;
	}
	
	template <typename type>
	double norm(dim3<type> vec)
	{
        return std::sqrt(vec.x*vec.x + vec.y*vec.y + vec.z*vec.z);
	}
	
	template <typename type>
	// square of norm
	type norm2(dim3<type> vec)
	{
		return vec.x*vec.x + vec.y*vec.y + vec.z*vec.z;
	}
	
	template <typename type>
	dim3<type> operator/(const dim3<type>& lhs, const dim3<type>& rhs)
	{
		dim3<type> toReturn;
		double normRhs = norm(rhs);
		toReturn = lhs/normRhs;
		return toReturn;
	}
	
	template <typename type>
	dim3<type> operator/=(dim3<type>& lhs, const dim3<type>& rhs)
	{
		double normRhs = norm(rhs);
		lhs/=normRhs;
		return lhs;
	}
	
	template <typename type>
	dim3<type> operator%(const dim3<type>& lhs, const dim3<type>& rhs)
	{
		dim3<type> toReturn = lhs;
		toReturn.x %= rhs.x;
		toReturn.y %= rhs.y;
        toReturn.z %= rhs.z;
		return toReturn;
	}
	
	template <typename type>
	dim3<type> operator%=(dim3<type>& lhs, const dim3<type>& rhs)
	{
		lhs.x %= rhs.x;
		lhs.y %= rhs.y;
        lhs.z %= rhs.z;
		return lhs;
	}

	template <typename type>
    std::ostream &operator<<(std::ostream &os, const dim3<type> &vec)
    {
        return os << vec.x << "\t" << vec.y << "\t" << vec.z;
    }
	
	template <typename type>
	double dist(dim3<type> vec1, dim3<type> vec2)
	{
		return norm(vec2 - vec1);
	}
	
	template <typename type>
	type dist2(dim3<type> vec1, dim3<type> vec2)
	{
		return norm2(vec2-vec1);
	}
	
	template <typename type>
	dim3<type> cross(const dim3<type>& a, const dim3<type>& b)
	{
		dim3<type> toReturn;
		toReturn.x = a.y*b.z - a.z*b.y;
		toReturn.y = a.z*b.x - a.x*b.z;
		toReturn.z = a.x*b.y - a.y*b.x;
		
		return toReturn;
	}

	double3 operator+ (const double3& lhs, const int3& rhs);
	double3 operator+ (const int3& lhs, const double3& rhs);
	double3 operator += (double3& lhs, const int3& rhs);

	double3 operator*(const int3&, const double&);
	double3 operator*(const double&, const int3&);

	double3 normalize(const double3&);
	double angleBetween(const double3&, const double3&);
	double radBetween(const double3&, const double3&);
	double3 rotate(const double3&, const double3&, double);

	int3 round(const double3&);

	template <typename type>
		class vector3
		{
			public:
				vector3();
				vector3(int3);
				void resize(int3);
				type& at(const int3&);
				type ac(const int3&) const;
				int3 size() const;
				void clear();

			private:
				int3 isize;
				int xy;
				std::vector<type> v;
		};

	template<typename type>
		void vector3<type>::resize(int3 s)
		{
			isize = s;
			xy = isize.x*isize.y;

			if (s.x > 0 && s.y > 0 && s.z > 0)
			{
				v.resize(s.x*s.y*s.z);
			}
			else
			{
				v.resize(0);
			}
		}

	template<typename type>
		vector3<type>::vector3()
		{
			resize(int3{0,0,0});
		}

	template<typename type>
		int3 vector3<type>::size() const
		{
			return isize;
		}

	template<typename type>
		vector3<type>::vector3(int3 s)
		{
			resize(s);
		}

	template<typename type>
		type& vector3<type>::at(const int3& p)
		{
			int ip = p.x + p.y*isize.x + p.z*xy;

#ifdef debug
			if (ip >= v.size())
			{
				std::cerr << "vector3<>.at(): requesting " << p << " outside range " << isize << std::endl;
			}
#endif
			return v[ip];
		}

	template<typename type>
		type vector3<type>::ac(const int3& p) const
		{
			int ip = p.x + p.y*isize.x + p.z*xy;

#ifdef debug
			if (ip >= v.size())
			{
				std::cerr << "vector3<>.ac(): requesting " << p << " outside range " << isize << std::endl;
			}
#endif
			return v[ip];
		}

	template<typename type>
		void vector3<type>::clear()
		{
			v.clear();
		}
}


#endif
