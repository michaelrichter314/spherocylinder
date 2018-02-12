#ifndef fpstring
#define fpstring

#include <stringstream>
#include <string>

namespace fp
{
	/*! returns characters after last period */
	char* getExtension(const char* in)
	{
		char current;
		int pos = 0;
		int lastPeriod = -1;
		
		do
		{
			current = in[pos];
			if (current == '.') {lastPeriod = pos;}
			pos++;
		} while (current != '\0');
		
		if (lastPeriod != -1)
		{
			int extensionLength = pos - lastPeriod - 2;
			
			char* out = new char[extensionLength + 1];
			for (int i = 0; i < extensionLength + 1; i++)
			{
				out[i] = in[lastPeriod + 1 + i];
			}
			
			return out;
		}
		else
		{
			char* out = new char[1];
			out[0] = '\0';
			return out;
		}
	}

		
}

#endif
