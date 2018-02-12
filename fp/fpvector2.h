#ifndef fpvector2h
#define fpvector2h

#include <vector>

namespace fp
{

template<typename type>
class vector2
{
	public:
		vector2();
		vector2(int, int);

		void resize(int, int);

		type& at(int, int);
		type ac(int, int) const;
		int sizeX() const;
		int sizeY() const;

	private:
		int sizex;
		int sizey;
		std::vector<type> v;
};

	template<typename type>
void vector2<type>::resize(int x, int y)
{
	sizex = x;
	sizey = y;

	if (x > 0 && y > 0)
	{
		v.resize(x*y);
	}
	else
	{
		v.resize(0);
	}
}

	template<typename type>
vector2<type>::vector2()
{
	resize(0,0);
}

template<typename type>
int vector2<type>::sizeX() const
{
	return sizex;
}
template<typename type>
int vector2<type>::sizeY() const
{
	return sizey;
}

	template<typename type>
vector2<type>::vector2(int x, int y)
{
	resize(x,y);
}

	template<typename type>
type& vector2<type>::at(int x, int y)
{
	int ip = x + y*sizex;

#ifdef debug
	if (ip >= v.size())
	{
		std::cerr << "vector2<>.at(): requesting " << x << ", " << y << " outside range " << sizex << "," << sizey << std::endl;
	}
#endif
	return v[ip];
}

template<typename type>
type vector2<type>::ac(int x, int y) const
{
	int ip = x + y*sizex;

#ifdef debug
	if (ip >= v.size())
	{
		std::cerr << "vector2<>.at(): requesting " << x << "," << y << " outside range " << sizex << "," << sizey << std::endl;
	}
#endif
	return v[ip];
}
}

#endif
