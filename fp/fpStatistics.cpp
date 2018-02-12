#include <fpStatistics.h>

namespace fp
{

	sample<hp> histogram(sample<double> &input, int sections)
	{
		double minValue = min(input);
		double maxValue = max(input);
		return histogram(input, sections, minValue, maxValue);
	}

	sample<hp> histogram(sample<double> &input, int sections, double minValue, double maxValue)
	{
		if (minValue == maxValue)
		{
			sample<hp> toReturn;
			hp p;
			p.count = input.getN();
			p.pos = minValue;
			toReturn.add(p);
			
			return toReturn;
		}
		
		sample<hp> toReturn;
		
		int* count = (int*)std::calloc(sections, sizeof(int));
		
		for (int i = 0; i < input.getN(); i++)
		{
			count[(int)((sections-1)*((input[i]-minValue)/(maxValue - minValue)))]++;
		}
		
		for (int i = 0; i < sections; i++)
		{
			hp toAdd;
			toAdd.count = count[i];
			toAdd.pos = minValue + i * (maxValue - minValue)/(sections-1);
			toReturn.add(toAdd);
		}
		
		free(count);
		
		return toReturn;
	}
	
	sample<hpint> histogram(sample<int> &input)
    {
        int minVal = min(input);
        int maxVal = max(input);
        
		sample<hpint> toReturn;
		
		int* count = (int*)std::calloc(maxVal-minVal+1, sizeof(int));
		
		for (int i = 0; i < input.getN(); i++)
		{
			count[input[i]-minVal]++;
		}
		
		for (int i = 0; i < maxVal-minVal+1; i++)
		{
			hpint toAdd;
			toAdd.count = count[i];
			toAdd.pos = minVal + i;
			toReturn.add(toAdd);
		}
		
		free(count);
		
		return toReturn;
    }
	
	template<>
	void sample<hp>::print()
	{
		int size = getN();
		for (int i = 0; i < size; i++)
		{
			std::cout << dataVector[i].pos << "\t" << dataVector[i].count << std::endl;
		}
	}

	template<>
	void sample<hp>::save(std::string filename)
	{
		int size = getN();
		std::ofstream file(filename.c_str());
		for (int i = 0; i < size; i++)
		{
			file << dataVector[i].pos << "\t" << dataVector[i].count << std::endl;
		}
		file.close();
	}
	
	template<>
	void sample<hpint>::print()
	{
		int size = getN();
		for (int i = 0; i < size; i++)
		{
			std::cout << dataVector[i].pos << "\t" << dataVector[i].count << std::endl;
		}
	}

}
