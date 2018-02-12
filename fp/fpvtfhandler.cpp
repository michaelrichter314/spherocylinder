#include <fpvtfhandler.h>

namespace fp
{
	std::string VTFatom::getLine()
	{
		std::string toReturn = "atom ";
		toReturn += std::to_string(ID);
		if (name.hasData()) {toReturn += " name "; toReturn += name.getData();}
		if (type.hasData()) {toReturn += " type "; toReturn += type.getData();}
		if (resid.hasData()) {toReturn += " resid "; toReturn += std::to_string(resid.getData());}
		if (resname.hasData()) {toReturn += " resname "; toReturn += resname.getData();}
		if (radius.hasData()) {toReturn += " radius "; toReturn += std::to_string(radius.getData());}
		if (segid.hasData()) {toReturn += " segid "; toReturn += segid.getData();}
		if (chain.hasData()) {toReturn += " chain "; toReturn += chain.getData();}
		if (charge.hasData()) {toReturn += " charge "; toReturn += std::to_string(charge.getData());}
		if (atomicnumber.hasData()) {toReturn += " atomicnumber "; toReturn += std::to_string(atomicnumber.getData());}
		if (altloc.hasData()) {toReturn += " altloc "; toReturn += altloc.getData();}
		if (insertion.hasData()) {toReturn += " insertion "; toReturn += insertion.getData();}
		if (occupancy.hasData()) {toReturn += " occupancy "; toReturn += std::to_string(occupancy.getData());}
		if (bfactor.hasData()) {toReturn += " bfactor "; toReturn += std::to_string(bfactor.getData());}
		if (mass.hasData()) {toReturn += " mass "; toReturn += std::to_string(mass.getData());}

		return toReturn;

	}

	std::string VTFbond::getLine() const
	{
		std::string toReturn = "bond ";
		toReturn += std::to_string(from);
		toReturn += ":";
		toReturn += std::to_string(to);

		return toReturn;
	}

	std::string VTFunitcell::getLine(bool lattice)
	{
		std::string toReturn = "pbc ";
		if (lattice)
		{
			toReturn += std::to_string(isize.x);
			toReturn += " ";
			toReturn += std::to_string(isize.y);
			toReturn += " ";
			toReturn += std::to_string(isize.z);
		}
		else
		{
			toReturn += std::to_string(dsize.x);
			toReturn += " ";
			toReturn += std::to_string(dsize.y);
			toReturn += " ";
			toReturn += std::to_string(dsize.z);
		}

		if (angle.hasData())
		{
			auto a = angle.getData();
			toReturn += " ";
			toReturn += std::to_string(a.x);
			toReturn += " ";
			toReturn += std::to_string(a.y);
			toReturn += " ";
			toReturn += std::to_string(a.z);
		}

		return toReturn;
	}






	VTFhandler::VTFhandler()
	{
	}

	VTFhandler::VTFhandler(std::string fname, bool l)
	{
		open(fname, l);
	}

	/*VTFhandler::~VTFhandler()
	  {
	  file.close();
	  }*/

	void VTFhandler::open(std::string fname, bool lat)
	{
		if (opened) {return;}

		lattice = lat;

		filename = fname;

		std::ifstream file(filename.c_str());

		const int buffersize = 4096;
		char buffer[buffersize];

		bool isContinued = false; // remembers if last line had tailing '\'
		std::stringstream continueBuffer;

		bool parsingHeader = true;

		std::streampos sp = file.tellg();


		while (file.getline(buffer, buffersize))
		{
			if (isContinued)
			{
				int p = 0;
				while (buffer[p] != '\\' && buffer[p] != '\0')
				{
					continueBuffer << buffer[p];
					p++;
				}

				if (!isContinuingLine(buffer))
				{
					// last line of continuing sequence: parse line
					continueBuffer >> buffer;
					continueBuffer.clear();
					isContinued = false;
					parseLineFirstPass(buffer, parsingHeader, sp);
				}
				else
				{
					// this one is another continuing line.
					isContinued = true;
				}
			}
			else
			{
				isContinued = isContinuingLine(buffer);

				if (isContinued)
				{
					// write first line into continue buffer
					int p = 0;
					while (buffer[p] != '\\' && buffer[p] != '\0')
					{
						continueBuffer << buffer[p];
						p++;
					}
				}
				else
				{
					// parse line
					parseLineFirstPass(buffer, parsingHeader, sp);
				}
			}

			sp = file.tellg();
		}


		file.close();
		timesteps = timestepPosition.size();
		opened = true;

		//		constructNameMap();


	}

	void VTFhandler::parseLineFirstPass(const char* s, bool& header, std::streampos sp)
	{
		int p = 0;
		while (s[p] == ' ' || s[p] == '\t') {p++;}

		const char& firstNonWhitespace = s[p];

		if (header)
		{
			switch (firstNonWhitespace)
			{
				case '#':
					break;
				case 'a':
					{
					auto a = parseAtomLine(s);
					atom.insert(atom.end(), a.begin(), a.end());
					}
					break;
				case 'b':
					{
					auto b = parseBondLine(s);
					bond.insert(bond.end(), b.begin(), b.end());
					}
					break;
				case 't':
				case 'c':
					parseTimestepLine(s, sp, header);
					break;
				case 'p':
				case 'u':
					{
					auto u = parseUnitcellLine(s);
					if (timestepPBC.size() == 0)
					{
						timestepPBC.push_back(u);
					}
					else
					{
						timestepPBC[0] = u;
					}
					}
					break;
				default:
					if (isNumber(firstNonWhitespace))
					{
						auto a = parseAtomLine(s);
						atom.insert(atom.end(), a.begin(), a.end());
					}
					else
					{
						std::cerr << "VTFhandler::parseLineFirstPass(...): Cannot parse line. Unknown linetype for header: \"" << s << "\"." << std::endl;
					}
			}
		}
		else
		{
			switch (firstNonWhitespace)
			{
				case '#':
					break;
				case 't':
				case 'c':
					parseTimestepLine(s, sp, header);
					break;
				case 'P':
				case 'p':
				case 'u':
					{
					auto u = parseUnitcellLine(s);
					if (timestepPBC.size() < timestepPosition.size())
					{
						timestepPBC.push_back(u);
					}
					else
					{
						timestepPBC[timestepPBC.size()-1] = u;
					}
					}
					break;
				case 'i':
					indexed = true;
					break;
				default:
					if (!isNumber(firstNonWhitespace))
					{
						std::cerr << "VTFhandler::parseLineFirstPass(...): Cannot parse line. Unknown linetype for body: \"" << s << "\"." << std::endl;
					}
			}
		}
	}


	VTFtimestep VTFhandler::getTimestep(int t)
	{
		VTFtimestep toReturn;
		if (lattice)
		{
			toReturn.ipos.resize(atom.size());
		}
		else
		{
			toReturn.dpos.resize(atom.size());
		}

		if (!opened) {throw 0;}

		// default to last timestep
		if (t < 0) {t = timesteps-1;}
		if (t >= timesteps) {t = timesteps - 1;}

		//std::cout << "getting timestep " << t << std::endl;

		if (!fileOpen)
		{
			file.open(filename.c_str());
			fileOpen = true;
		}

		// go line by line
		const int buffersize = 4096;
		char buffer[buffersize];

		// remember which atoms have been found yet
		std::vector<bool> foundAtom(atom.size(), false);
		bool foundAllAtoms = false;
		int backupSteps = 0;

		while (!foundAllAtoms)
		{
			// go to correct position
			file.clear();
			file.seekg(timestepPosition[t-backupSteps]);
			int lastAtomFound = -1;

			bool doneWithTimestep = false;
			bool startedTimestep = false;

			while (!doneWithTimestep && file.getline(buffer, buffersize))
			{
				int p = 0;
				while (buffer[p] == ' ' || buffer[p] == '\t') {p++;}
				const char& firstNonWhitespace = buffer[p];

				switch (firstNonWhitespace)
				{
					case 't':
					case 'c':
						if (!startedTimestep)
						{
							parseTimestepLineNoAdd(buffer);
							startedTimestep = true;
						}
						else
						{
							doneWithTimestep = true;
						}

						break;
					case 'i':
						indexed = true;
						break;
					case 'P':
					case 'p':
					case 'u':
						break;
					default:
						if (isNumber(firstNonWhitespace))
						{
							// coordinate
							lastAtomFound = parseCoordinate(buffer, toReturn, foundAtom, lastAtomFound);
						}
						else
						{
							std::cerr << "VTFhandler::getTimestep(" << t << "): Cannot read line \"" << buffer << "\"." << std::endl;
						}
						break;
				}
			}

			foundAllAtoms = true;
			// check if all atoms have been foudn
			for (int i = 0; i < foundAtom.size(); i++)
			{
				if (!foundAtom[i])
				{
					foundAllAtoms = false;
					break;
				}
			}

			if (!foundAllAtoms)
			{
				if (backupSteps < t)
				{
					backupSteps++;
				}
				else
				{
					std::cerr << "VTFhandler::getTimestep(" << t << "): Cannot find some positions." << std::endl;
					throw 1;
				}
			}
		}



		toReturn.PBC = timestepPBC[t];
		return toReturn;
	}

	/*bool VTFhandler::atomsValid() const
	  {
	  bool valid = true;
	  std::vector<bool> exists(atom.size(), false);
	  for (const auto& a : atom)
	  {
	  if (a.ID < 0 || a.ID >= atom.size())
	  {
	  valid = false;
	  }
	  else
	  {
	  exists[a.ID] = true;
	  }
	  }

	  for (const auto& e : exists)
	  {
	  if (!e) {valid = false;}
	  }

	  return valid;
	  }

	  bool VTFhandler::bondsValid() const
	  {
	  bool valid = true;

	// check if bonded atoms exist (given atoms are all valid) and if we don't bond an atom to itself
	for (const auto& b : bond)
	{
	if (b.from < 0 || b.from >= atom.size() || b.to < 0 || b.to >= atom.size() || b.from == b.to)
	{
	valid = false;
	}
	}

	// check if we have any double bonds
	for (int i = 0; i < bond.size(); i++)
	{
	for (int j = i+1; j < bond.size(); j++)
	{
	if ((bond[i].from == bond[j].from) && (bond[i].to == bond[j].to)) { valid = false;}
	if ((bond[i].from == bond[j].to) && (bond[i].to == bond[j].from)) { valid = false;}
	}
	}

	return valid;
	}*/

	void VTFhandler::constructNameMap()
	{
		// generate vector with all names
		std::vector<std::string> names;
		for (auto&& a : atom)
		{
			const std::string& n = a.name.getData();

			if (find(names.begin(), names.end(), n) == names.end())
			{
				// not yet in list. add
				names.push_back(n);
			}
		}

		// sort list
		std::sort(names.begin(), names.end());

		// generate map
		int i = 0;
		for (const auto& n : names)
		{
			nameToId[n]=i;
			i++;
		}

		nameMapExists = true;
	}

	//for (const auto& n : nameToId)
	//{
	//cout << n.first << " " << n.second << endl;
	//}

	std::vector<VTFatom> VTFhandler::getAtoms() const
	{
		return atom;
	}

	std::vector<VTFbond> VTFhandler::getBonds() const
	{
		return bond;
	}

	int VTFhandler::getIdFromName(MaybeData<std::string>& d)
	{
		if (!nameMapExists)
		{
			constructNameMap();
		}

		if (d.hasData())
		{
			return nameToId.at(d.getData());
		}
		else
		{
			std::cerr << "VTFhandler::getIdFromName(): No name defined." << std::endl;
			throw 1;
		}
	}

	int VTFhandler::getTimesteps() const
	{
		return timesteps;
	}

	bool VTFhandler::isNumber(const char& c) const
	{
		return (c >= '0' && c <= '9') || c == '-';
	}

	std::vector<int> VTFhandler::parseAidSpecifier(std::string s) const
	{
		std::vector<int> toReturn;
		size_t N = s.size();

		int currentID = 0;
		int lastID = -1;

		for (size_t i = 0; i < N; i++)
		{
			if (isNumber(s[i]))
			{
				currentID *= 10;
				currentID += int(s[i] - '0');
			}
			else if (s[i] == ',')
			{
				// finish last one
				if (lastID != -1) // e.g. 2:13
				{
					int start = std::min(lastID, currentID);
					int end = std::max(lastID, currentID);

					// add last one
					for (int j = start; j <= end; j++)
					{
						// add if it does not exist yet
						if (std::find(toReturn.begin(), toReturn.end(), j) == toReturn.end())
						{
							toReturn.push_back(j);
						}
					}
				}
				else
				{
					// add if it does not exist yet
					if (std::find(toReturn.begin(), toReturn.end(), currentID) == toReturn.end())
					{
						toReturn.push_back(currentID);
					}
				}

				// reset stuff
				currentID = 0;
				lastID = -1;
			}
			else if (s[i] == ':')
			{
				lastID = currentID;
				currentID = 0;
			}
			else if (s[i] == 'd'
					&& s[i+1] == 'e'
					&& s[i+2] == 'f'
					&& s[i+3] == 'a'
					&& s[i+4] == 'u'
					&& s[i+5] == 'l'
					&& s[i+6] == 't')
			{
				toReturn.push_back(-1);
			}
			else
			{
				std::cerr << "VTFhandler::parseAidSpecifier(...): Problems parsing \"" << s << "\"." << std::endl;
			}


		}
		// finish last one
		if (lastID != -1) // e.g. 2:13
		{
			int start = std::min(lastID, currentID);
			int end = std::max(lastID, currentID);

			// add last one
			for (int j = start; j <= end; j++)
			{
				// add if it does not exist yet
				if (std::find(toReturn.begin(), toReturn.end(), j) == toReturn.end())
				{
					toReturn.push_back(j);
				}
			}
		}
		else
		{
			// add if it does not exist yet
			if (std::find(toReturn.begin(), toReturn.end(), currentID) == toReturn.end())
			{
				toReturn.push_back(currentID);
			}
		}

		return toReturn;

	}

	std::vector<VTFatom> VTFhandler::parseAtomLine(const char* c) const
	{

		int p = 0;
		// first: read word "atom", if it's there
		if (c[0] == 'a')
		{
			while (c[p] != '\0' && (c[p] != ' ' && c[p] != '\t')) {p++;}
		}

		// second: read whitespaces after word "atom"
		while (c[p] == ' ' || c[p] == '\t') {p++;}

		// now position is first character of aid specifier. read the whole thing
		std::string aIDs;

		while (c[p] != '\0' && (c[p] != ' ' && c[p] != '\t'))
		{
			aIDs.push_back(c[p]);
			p++;
		}
		std::vector<int> IDs = parseAidSpecifier(aIDs);

		// skip all whitespaces after aid specifier	
		while (c[p] == ' ' || c[p] == '\t') {p++;}

		// parse rest of line
		VTFatom a;
		a.ID = -1;
		TwoWords words;
		bool thisIsSecondWord = false;
		while (c[p] != '\0')
		{
			if (c[p] == ' ' || c[p] == '\t')
			{
				if (thisIsSecondWord)
				{
					addAtomInfo(a, words);
					thisIsSecondWord = false;
					words.A.clear();
					words.B.clear();
				}
				else
				{
					thisIsSecondWord = true;
				}

				// skip all following whitespaces
				while (c[p+1] == ' ' || c[p+1] == '\t') {p++;}
			}
			else
			{
				if (thisIsSecondWord)
				{
					words.B.push_back(c[p]);
				}
				else
				{
					words.A.push_back(c[p]);
				}
			}
			p++;
		}

		if (thisIsSecondWord)
		{
			addAtomInfo(a,words);
		}

		// now, all atom info is saved and we know how many we need. Build them!
		std::vector<VTFatom> toReturn;

		for (const auto& ID : IDs)
		{
			a.ID = ID;
			toReturn.push_back(a);
		}

		return toReturn;

	}

	void VTFhandler::addAtomInfo(VTFatom& a, const TwoWords& w) const
	{
		if (w.A.compare("n") == 0 || w.A.compare("name") == 0)
		{
			a.name.setData(w.B);
		}
		else if (w.A.compare("t") == 0 || w.A.compare("type") == 0)
		{
			a.type.setData(w.B);
		}
		else if (w.A.compare("resid") == 0)
		{
			a.resid.setData(atoi(w.B.c_str()));
		}
		else if (w.A.compare("res") == 0 || w.A.compare("resname") == 0)
		{
			a.resname.setData(w.B);
		}
		else if (w.A.compare("r") == 0 || w.A.compare("radius") == 0)
		{
			a.radius.setData(atof(w.B.c_str()));
		}
		else if (w.A.compare("s") == 0 || w.A.compare("segid") == 0)
		{
			a.segid.setData(w.B);
		}
		else if (w.A.compare("c") == 0 || w.A.compare("chain") == 0)
		{
			a.chain.setData(w.B);
		}
		else if (w.A.compare("q") == 0 || w.A.compare("charge") == 0)
		{
			a.charge.setData(atof(w.B.c_str()));
		}
		else if (w.A.compare("a") == 0 || w.A.compare("atomicnumber") == 0)
		{
			a.atomicnumber.setData(atoi(w.B.c_str()));
		}
		else if (w.A.compare("altloc") == 0)
		{
			a.altloc.setData(w.B);
		}
		else if (w.A.compare("i") == 0 || w.A.compare("insertion") == 0)
		{
			a.insertion.setData(w.B);
		}
		else if (w.A.compare("o") == 0 || w.A.compare("occupancy") == 0)
		{
			a.occupancy.setData(atof(w.B.c_str()));
		}
		else if (w.A.compare("b") == 0 || w.A.compare("bfactor") == 0)
		{
			a.bfactor.setData(atof(w.B.c_str()));
		}
		else if (w.A.compare("m") == 0 || w.A.compare("mass") == 0)
		{
			a.mass.setData(atof(w.B.c_str()));
		}
		else
		{
			std::cerr << "VTFhander::addAtomInfo(...): Unknown atom info: \"" << w.A << "\"." << std::endl;
		}
	}


	std::vector<VTFbond> VTFhandler::parseBondLine(const char* c) const
	{
		std::vector<VTFbond> toReturn;

		int p = 0;
		// first: read word "bond"
		while (c[p] != '\0' && (c[p] != ' ' && c[p] != '\t')) {p++;}

		// second: read whitespaces after word "bond"
		while (c[p] == ' ' || c[p] == '\t') {p++;}

		int currentNumber = 0;
		int lastNumber = -1;
		bool isChain = false;
		while (c[p] != '\0')
		{
			if (isNumber(c[p]))
			{
				currentNumber *= 10;
				currentNumber += (int)(c[p] - '0');
			}
			else if (c[p] == ':')
			{
				// get ready for second number
				lastNumber = currentNumber;
				currentNumber = 0;

				if (c[p+1] == ':')
				{
					isChain = true;
					p++;
				}
			}
			else if (c[p] == ',')
			{
				addBonds(toReturn, lastNumber, currentNumber, isChain);
			}
			else
			{
				std::cerr << "VTFhandler::parseBondLine(): Cannot read bond-specifier in \"" << c << "\" at position " << p << "." << std::endl;
			}

			p++;
		}

		// parse last specifier
		if (lastNumber != -1)
		{
			addBonds(toReturn, lastNumber, currentNumber, isChain);
		}
		else
		{
			std::cerr << "VTFhandler::parseBondLine(): Problem reading last bond-specifier." << std::endl;
		}

		return toReturn;
	}

	void VTFhandler::addBonds(std::vector<VTFbond>& v, int a, int b, bool c) const
	{
		int A = std::min(a,b);
		int B = std::max(a,b);

		if (c)
		{
			for (int i = A; i < B; i++)
			{
				v.push_back({i, i+1});
			}
		}
		else
		{
			v.push_back({A,B});
		}
	}


	VTFunitcell VTFhandler::parseUnitcellLine(const char* s) const
	{
		VTFunitcell toReturn;

		toReturn.dsize = {0.0,0.0,0.0};
		toReturn.isize = {0,0,0};

		std::stringstream stream;
		stream << s;

		std::string sbuf;

		double3 anglebuf;

		int word = 0;
		while (stream >> sbuf)
		{
			if (lattice)
			{
				switch (word)
				{
					case 0:
						break;
					case 1:
						toReturn.isize.x = stoi(sbuf);
						break;
					case 2:
						toReturn.isize.y = stoi(sbuf);
						break;
					case 3:
						toReturn.isize.z = stoi(sbuf);
						break;
				}
			}
			else
			{
				switch (word)
				{
					case 0:
						break;
					case 1:
						toReturn.dsize.x = stod(sbuf);
						break;
					case 2:
						toReturn.dsize.y = stod(sbuf);
						break;
					case 3:
						toReturn.dsize.z = stod(sbuf);
						break;
				}
			}
			switch (word)
			{
				case 0:
				case 1:
				case 2:
				case 3:
					break;
				case 4:
					anglebuf.x = stod(sbuf);
					break;
				case 5:
					anglebuf.y = stod(sbuf);
					break;
				case 6:
					anglebuf.z = stod(sbuf);
					toReturn.angle.setData(anglebuf);
					break;
				default:
					std::cerr << "VTFhander::parseUnitcellLine(): Line too long: \"" << s << "\", there are " << word << " words." << std::endl;
			}
			word++;
		}
		return toReturn;
	}

	bool VTFhandler::isContinuingLine(const char* buffer) const
	{
		int p = 0;

		while (buffer[p] != '\0') {p++;} // skip to end
		do 
		{
			p--;
		} while (buffer[p] == ' ' || buffer[p] == '\t'); // return to last non-whitespace

		return ( buffer[p] == '\\');

	}

	void VTFhandler::parseTimestepLine(const char* s, std::streampos sp, bool& header)
	{
		timestepPosition.push_back(sp);
		parseTimestepLineNoAdd(s);
		header = false;
	}
	void VTFhandler::parseTimestepLineNoAdd(const char* s)
	{
		indexed = false;

		// check if indexed
		int p = 0;
		while (s[p] == ' ' || s[p] == '\t') {p++;} // skip leading spaces
		while (s[p] != ' ' && s[p] != '\t' && s[p] != '\0') {p++;} // skip next word
		while (s[p] != '\0' && (s[p] == ' ' || s[p] == '\t')) {p++;} // skip spaces after first word ("timestep")
		if (s[p] == 'i') {indexed = true;}
	}

	int VTFhandler::parseCoordinate(const char* s, VTFtimestep& timestep, std::vector<bool>& found, int lastFound)
	{
		int thisID;
		std::stringstream stream;
		stream << s;

		if (indexed)
		{
			stream >> thisID;
		}
		else
		{
			thisID = lastFound + 1;
		}
		if (thisID < 0 || thisID >= atom.size())
		{
			std::cerr << "VTFhandler::parseCoordinate(): Invalid atom ID \"" << thisID << "\"." << std::endl;
		}

		// do not override already found atoms
		if (found[thisID]) 
		{
			return thisID;
		}

		bool done = false;

		if (lattice)
		{
			int3 pos;
			double dbuf;
			int ibuf;
			int words = 0;
			while (!done && stream >> dbuf)
			{
				ibuf = std::round(dbuf);
				switch (words)
				{
					case 0:
						pos.x = ibuf;
						break;
					case 1:
						pos.y = ibuf;
						break;
					case 2:
						pos.z = ibuf;
						done = true;
						break;
				}
				words++;
			}

			if (done)
			{
				timestep.ipos[thisID] = pos;
				found[thisID] = true;
				return thisID;
			}
			else
			{
				std::cerr << "VTFhandler::parseCoordinate(): Invalid line \"" << s << "\"." << std::endl;
			}
		}
		else
		{
			double3 pos;
			double dbuf;
			int words = 0;

			while (!done && stream >> dbuf)
			{
				switch (words)
				{
					case 0:
						pos.x = dbuf;
						break;
					case 1:
						pos.y = dbuf;
						break;
					case 2:
						pos.z = dbuf;
						done = true;
						break;
				}
				words++;
			}

			if (done)
			{
				timestep.dpos[thisID] = pos;
				found[thisID] = true;
				return thisID;
			}
			else
			{
				std::cerr << "VTFhandler::parseCoordinate(): Invalid line \"" << s << "\"." << std::endl;
			}
		}

		return -1;
	}

	void VTFhandler::writeHeader(std::ostream stream)
	{
		for (auto&& a : atom)
		{
			stream << a.getLine() << std::endl;
		}

		for (auto&& b : bond)
		{
			stream << b.getLine() << std::endl;
		}
	}
}
