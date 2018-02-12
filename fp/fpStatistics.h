#ifndef fpstatistics_h
#define fpstatistics_h

#include <fstream>
#include <string>
#include "fpmath.h"
#include "fpMemory.h"
#include "fpSeries.h"
#include "fpSample.h"

namespace fp
{
	template <typename type>
	type max(type* data, int N)
	{
		if (N > 0)
		{
			type max = data[0];
			for (int i = 1; i < N; i++)
				if (data[i] > max) {max = data[i];}
			
			return max;
		}
	}
	
	template <typename type>
	type max(sample<type>& data)
	{
		int N = data.getN();
		if (N > 0)
		{
			type max = data[0];
			for (int i = 1; i < N; i++)
				if (data[i] > max) {max = data[i];}
			
			return max;
		}
	}
	
	
	
	template <typename type>
	type min(type* data, int N)
	{
		if (N > 0)
		{
			type min = data[0];
			for (int i = 1; i < N; i++)
				if (data[i] < min) {min = data[i];}
			
			return min;
		}
	}
	
	template <typename type>
	type min(sample<type>& data)
	{
		int N = data.getN();
		if (N > 0)
		{
			type min = data[0];
			for (int i = 1; i < N; i++)
				if (data[i] < min) {min = data[i];}
			
			return min;
		}
	}
	
	template <typename type>
	type mean(type* data, int N)
	{
		if (N > 0)
		{
			type average = data[0];
			for (int i = 1; i < N; i++)
				average = average + data[i];
			
			return average / N;
			//return average;
		}
	}
	
	template <typename type>
	type mean(sample<type>& data)
	{
		int N = data.getN();
		if (N > 0)
		{
			type average = data[0];
			for (int i = 1; i < N; i++)
				average = average + data[i];
			
			return average / N;
			//return average;
		}
	}
	
	template <typename type>
	double meanSquare(type* data, int N)
	{
		if (N > 0)
		{
			double average;
			for (int i = 0; i < N; i++)
				average += data[i]*data[i];
			
			return average / N;
		}
	}
	
	template <typename type>
	double meanSquare(sample<type>& data)
	{
		int N = data.getN();
		if (N > 0)
		{
			double average;
			for (int i = 0; i < N; i++)
				average += data[i]*data[i];
			
			return average / N;
		}
	}
	
	template <typename type>
	double RMS(type* data, int N)
	{
		return std::sqrt(meanSquare(data, N));
	}
	
	template <typename type>
	double RMS(sample<type>& data)
	{
		return std::sqrt(meanSquare(data));
	}
	
	template <typename type>
	double variance(type* data, int N)
	{
		if (N <= 0) {return 0.0;}
		
		double sum = 0.0;
		type Mean = mean(data, N);
		type temp;
		
		for (int i = 0; i < N; i++)
		{
			temp = data[i] - Mean;
			sum += temp * temp;
		}
		
		return sum / N;
	}
	
	template <typename type>
	double variance(sample<type>& data)
	{
		double N = data.getN();
		if (N <= 0) {return 0.0;}
		
		double sum = 0.0;
		type Mean = mean(data);
		type temp;
		
		for (int i = 0; i < N; i++)
		{
			temp = data[i] - Mean;
			sum += temp * temp;
		}
		
		return sum / N;
	}
	
	template <typename type>
	double stdev(sample<type>& data)
	{
		return std::sqrt(variance(data));
	}
	
	template <typename type>
	double stdevMean(sample<type>& data)
	{
		return std::sqrt(variance(data)) / data.getN();
	}
	
	
	
	template <typename type>
	double correlation(type* data1, type* data2, int N)
	{
		type mean1 = mean(data1, N);
		type mean2 = mean(data2, N);

		double nominator = 0.0;
		double sqsum1 = 0.0;
		double sqsum2 = 0.0;

		for (int i = 0; i < N; i++)
		{
			nominator += (data1[i] - mean1) * (data2[i] - mean2);
			sqsum1 += (data1[i] - mean1) * (data1[i] - mean1);
			sqsum2 += (data2[i] - mean2) * (data2[i] - mean2);
		}

		double denominator = std::sqrt(sqsum1 * sqsum2);

		return nominator / denominator;
	}
	
	struct hp
	{
		int count;
		double pos;
	};
	
	sample<hp> histogram(sample<double> &input, int sections);
	sample<hp> histogram(sample<double> &input, int sections, double, double);
	
	struct hpint
	{
        int count;
        int pos;
    };
    
	sample<hpint> histogram(sample<int> &input);

	template<>
	void sample<hp>::print();

	template<>
	void sample<hp>::save(std::string);

	template<>
	void sample<hpint>::print();
}

#endif
