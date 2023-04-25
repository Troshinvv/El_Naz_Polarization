#include <iostream> // std::cout
#include <fstream>  // std::ifstream
#include <map>

#include "MpdAnalysisTask2.h"

using namespace std;

ClassImp(MpdAnalysisTask2);

MpdAnalysisTask2::MpdAnalysisTask2(const char *name, const char *outputName) : MpdAnalysisTask(name, outputName) {}

void MpdAnalysisTask2::readParameters(std::string fname)
{
   // Read from file in format
   // varname  value
   // comments start from #

   std::string fInputFileTtx = fname;
   fInputFileTtx             = fInputFileTtx + ".txt";

   cout << endl << "Read input text file: " << fInputFileTtx << endl;

   std::ifstream ifs(fInputFileTtx);
   if (!ifs) {
      cout << "File " << fInputFileTtx << " can not be opened -> using input parameters from the header file instead"
           << endl;
      return;
   }

   std::string a, b;
   while (ifs.good()) {
      ifs >> a;
      // if comment, skip to the enf of line
      if (a.find_first_of('#') == 0 || a.find_first_of("//") == 0) {
         ifs.ignore(999, '\n');
         continue;
      } else {
	      ifs >> b;
         mMap.insert({a, b});
      }
   }
   ifs.close();
}

void MpdAnalysisTask2::param(std::string name, bool &param, bool defaultvalue)
{
   auto search = mMap.find(name);
   if (search != mMap.end()) {
      if (search->second.compare("true") == 0 || search->second.compare("TRUE") == 0) {
         param = true;
      } else {
         param = false;
      }
      if (printoutparams) cout << "parameter " << name << " set to value: " << param << endl;
   } else {
      param = defaultvalue;
      if (printoutparams) cout << "parameter " << name << " not set. default value used: " << param << endl;
   }
}
void MpdAnalysisTask2::param(std::string name, int &param, int defaultvalue)
{
   auto search = mMap.find(name);
   if (search != mMap.end()) {
      param = atoi(search->second.data());
      if (printoutparams) cout << "parameter " << name << " set to value: " << param << endl;
   } else {
      param = defaultvalue;
      if (printoutparams) cout << "parameter " << name << " not set. default value used: " << param << endl;
   }
}
void MpdAnalysisTask2::param(std::string name, float &param, float defaultvalue)
{
   auto search = mMap.find(name);
   if (search != mMap.end()) {
      param = atof(search->second.data());
      if (printoutparams) cout << "parameter " << name << " set to value: " << param << endl;
   } else {
      param = defaultvalue;
      if (printoutparams) cout << "parameter " << name << " not set. default value used: " << param << endl;
   }
}
void MpdAnalysisTask2::param(std::string name, double &param, double defaultvalue)
{
   auto search = mMap.find(name);
   if (search != mMap.end()) {
      param = atof(search->second.data());
      if (printoutparams) cout << "parameter " << name << " set to value: " << param << endl;
   } else {
      param = defaultvalue;
      if (printoutparams) cout << "parameter " << name << " not set. default value used: " << param << endl;
   }
}
void MpdAnalysisTask2::param(std::string name, std::string &param, std::string defaultvalue)
{
   auto search = mMap.find(name);
   if (search != mMap.end()) {
      param = search->second;
      if (printoutparams) cout << "parameter " << name << " set to value: " << param << endl;
   } else {
      param = defaultvalue;
      if (printoutparams) cout << "parameter " << name << " not set. default value used: " << param << endl;
   }
}