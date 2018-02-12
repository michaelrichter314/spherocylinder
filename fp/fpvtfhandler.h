#ifndef fpvtfhandlerh
#define fpvtfhandlerh

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <fpdim3.h>
#include <fpStatistics.h>
#include <fpvector2.h>
#include <map>

namespace fp
{
	template <class type>
		class MaybeData
		{
			public:
				bool hasData() const;
				type& getData();
				void setData(type);

			private:
				bool hasdata = false;
				type data;
		};

	template <class type>
		bool MaybeData<type>::hasData() const
		{
			return hasdata;
		}

	template <class type>
		type& MaybeData<type>::getData()
		{
			if (!hasdata)
			{
				std::cerr << "MaybeData::getData(): No data available." << std::endl;
				throw 0;
			}
			return data;
		}
	template <class type>
		void MaybeData<type>::setData(type d)
		{
			data = d;
			hasdata = true;
		}



	struct VTFatom
	{
		int ID;
		MaybeData<std::string> name;
		MaybeData<std::string> type;
		MaybeData<int> resid;
		MaybeData<std::string> resname;
		MaybeData<double> radius;
		MaybeData<std::string> segid;
		MaybeData<std::string> chain;
		MaybeData<double> charge;
		MaybeData<int> atomicnumber;
		MaybeData<std::string> altloc;
		MaybeData<std::string> insertion;
		MaybeData<double> occupancy;
		MaybeData<double> bfactor;
		MaybeData<double> mass;

		std::string getLine();
	};

	struct VTFbond
	{
		int from;
		int to;

		std::string getLine() const;
	};

	struct VTFunitcell
	{
		double3 dsize;
		int3 isize;
		MaybeData<double3> angle;

		std::string getLine(bool);
	};

	struct VTFtimestep
	{
		VTFunitcell PBC;
		std::vector<fp::double3> dpos;
		std::vector<fp::int3> ipos;

		MaybeData<double> E;
		std::vector<VTFbond> dbond;

	};

	struct TwoWords
	{
		std::string A;
		std::string B;
	};

	class VTFhandler
	{
		public:
			VTFhandler(std::string, bool = false);
			VTFhandler();
			//~VTFhandler();
			void open(std::string, bool = false);

			VTFtimestep getTimestep(int = -1);
			std::vector<VTFatom> getAtoms() const;
			std::vector<VTFbond> getBonds() const;
			int getIdFromName(std::string) const;
			int getTimesteps() const;
		[[deprecated]]
		int getIdFromName(MaybeData<std::string>&);

		void writeHeader(std::ostream);


		private:
			std::vector<VTFatom> atom;
			std::vector<VTFbond> bond;
			bool lattice = false;
			int timesteps = 0;

			std::string filename;
			bool opened = false;
			

			//bool atomsValid() const;
			//bool bondsValid() const;
			//std::map<std::string, int> nameToId;
			//void constructNameMap();

			bool fileOpen = false;
			bool indexed = false;
			std::ifstream file;

			bool isNumber(const char&) const;

			void parseLineFirstPass(const char*, bool&, std::streampos);

			// get vector of ID's of atoms, e.g. {0,1,2,5} for 0:2,5
			std::vector<int> parseAidSpecifier(std::string) const;
			void addAtomInfo(VTFatom&, const TwoWords&) const;
			std::vector<VTFatom> parseAtomLine(const char*) const;

			std::vector<VTFbond> parseBondLine(const char*) const;
			void addBonds(std::vector<VTFbond>&, int, int, bool) const;

			VTFunitcell parseUnitcellLine(const char*) const;

			bool isContinuingLine(const char*) const;

			std::vector<std::streampos> timestepPosition;
			std::vector<VTFunitcell> timestepPBC;

			void parseTimestepLine(const char*, std::streampos, bool&);
			void parseTimestepLineNoAdd(const char*);
			int parseCoordinate(const char*, VTFtimestep&, std::vector<bool>&, int);
		
			[[deprecated]]
			void constructNameMap();

			[[deprecated]]
			bool nameMapExists = false;

			[[deprecated]]
			std::map<std::string,int> nameToId;
	};

}

#endif
