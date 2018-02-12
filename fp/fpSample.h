#ifndef fpsampleh
#define fpsampleh

#include <vector>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>

namespace fp
{
	template <typename type>
	class sample
	{
		private:
		 std::vector<type> dataVector;
		 type* dataArray;
		 int N;
		 
		public:
			
			sample()
			{
				dataArray = NULL;
			}

			sample(const std::vector<type>& v)
			{
				dataArray = NULL;
				dataVector.insert(dataVector.end(), v.begin(), v.end());
			}
			
			~sample()
			{
				if (dataArray != NULL)
				{
					free(dataArray);
				}
			}
			
			int getN()
			{
				N = dataVector.size();
				return N;
			}
			
			type* getData()
			{
				N = getN();
				
				//update dataArray
				if (dataArray != NULL) {free(dataArray);}
				dataArray = (type*)malloc(N * sizeof(type));
				
				for (int i = 0; i < N; i++)
				{
					dataArray[i] = dataVector[i];
				}
				
				return dataArray;
			}
			
			void add(type newData)
			{
				dataVector.push_back(newData);
			}
			
			void clear()
			{
				if (dataArray != NULL) {free(dataArray); dataArray = NULL;}
				dataVector.clear();
			}
			
			type doStatistics(type (*statFunction)(type*, int))
			{
				type* data = getData();
				int N = getN();
				
				type toReturn = statFunction(data, N);
				
				return toReturn;
			}
			
			type& operator[](const int& pos)
			{
				return dataVector[pos];
			}
			
			void print()
			{
				for (int i = 0; i < dataVector.size(); i++)
				{
					std::cout << dataVector[i] << std::endl;
				}
			}

			void save(std::string a)
			{
				//TODO
				std::cerr << "saving not supported!" << std::endl;
			}
	};
	
	
}

#endif
