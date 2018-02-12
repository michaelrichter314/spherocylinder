#pragma once

#include <iostream>
#include <vector>
#include <ostream>

using namespace std;

namespace fp
{

	template<typename type>
		type max(const vector<type>&);

	template<typename type>
		type min(const vector<type>&);

	template<class type>
		class Histogram
		{
			public:
				Histogram(int, const vector<type>&);
				Histogram(int, type, type);

				void insert(type);
				void insert(const vector<type>&);
				unsigned long getCount(int);



			protected:
				vector<unsigned long> count;
				vector<type> tic;
				int gaps;
				int tics;

				int underVal = 0;
				int overVal = 0;

				friend ostream& operator<<(ostream& os, const Histogram& h)
				{
					os << "#lTic\trTic\tmTic\tcount\n";
					os << "#out of range: " << h.underVal + h.overVal << " (" << h.underVal << " too small, " << h.overVal << " too large)\n";

					for (int i = 0; i < h.gaps; i++)
					{
						os << h.tic[i] << "\t" << h.tic[i+1] << "\t" << h.tic[i] + (h.tic[i+1]-h.tic[i])/2 << "\t" << h.count[i] << endl;
					}
				}
		};

	template<typename type>
		Histogram<type>::Histogram(int N, const vector<type>& v)
		{
			if (N <= 1)
			{
				throw 0;
			}

			gaps = N;
			tics = N + 1;

			count.resize(gaps);
			tic.resize(tics);

			for (auto& value : count) {value = 0;}

			// set tics
			// get lowest and highest center value
			type minSampleValue = fp::min(v);
			type maxSampleValue = fp::max(v);

			// calculate end points a,b
			type delta = (maxSampleValue - minSampleValue) / (2 * ( gaps-1 ) );
			type a = minSampleValue - delta;
			type b = maxSampleValue + delta;

			// set tic positions
			for (int i = 0; i < tics; i++)
			{
				tic[i] = a + ( (b-a)*i )/gaps;
			}

			// fill histogram with life
			insert(v);
		}

	template<typename type>
		Histogram<type>::Histogram(int N, type minvalue, type maxvalue) : Histogram(N, vector<type>{minvalue, maxvalue})
		{
			// fix minvalue and maxvalue
			count[0] = 0;
			count[gaps - 1] = 0;
		}

	template<typename type>
		void Histogram<type>::insert(type v)
		{
			int upperBoundId = tics - 1;
			int lowerBoundId = 0;

			if (v < tic[lowerBoundId]) {underVal++; return;}
			if (v > tic[upperBoundId]) {overVal++; return;}

			while (upperBoundId - lowerBoundId > 1)
			{
				int dist = upperBoundId - lowerBoundId;
				int midId = lowerBoundId + dist/2;

				if (v <= tic[midId])
				{
					upperBoundId = midId;
				}
				else
				{
					lowerBoundId = midId;
				}
			}

			//DEBUG
			if (upperBoundId - lowerBoundId != 1)
			{
				cerr << "lower bound ID: " << lowerBoundId << ", upper bound ID: " << upperBoundId << endl;
			}

			count[lowerBoundId]++;
		}

	template<typename type>
		void Histogram<type>::insert(const vector<type>& vec)
		{
			for (const auto& v : vec)
			{
				insert(v);
			}
		}

/*	template<typename type>
		ostream& operator<<(ostream& os, Histogram<type>& h)
		{
			for (int i = 0; i < h.gaps; i++)
			{
				os << h.tic[i] << " " << h.tic[i+1] << " " << h.tic[i] + (h.tic[i+1]-h.tic[i])/2 << " " << h.count[i] << endl;
			}
		}
*/

	template<typename type>
		unsigned long Histogram<type>::getCount(int i)
		{
			return count.at(i);
		}


	template<typename type>
		type max(const vector<type>& in)
		{
			if (in.size() == 0)
			{
				throw 0;
			}

			type toReturn = in[0];

			for (const auto& v : in)
			{
				if (v > toReturn) {toReturn = v;}
			}

			return toReturn;
		}

	template<typename type>
		type min(const vector<type>& in)
		{
			if (in.size() == 0)
			{
				throw 0;
			}

			type toReturn = in[0];

			for (const auto& v : in)
			{
				if (v < toReturn) {toReturn = v;}
			}

			return toReturn;
		}

} //end namespace fp
