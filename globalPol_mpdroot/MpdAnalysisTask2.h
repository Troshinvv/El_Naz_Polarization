#ifndef MpdAnalysisTask2_H
#define MpdAnalysisTask2_H

#include <map>
#include <string>

#include "Rtypes.h"
#include "TObject.h"

#include "MpdAnalysisTask.h"

class MpdAnalysisTask2 : public MpdAnalysisTask {
   ClassDef(MpdAnalysisTask2, 1);

public:
   MpdAnalysisTask2() {}
   MpdAnalysisTask2(const char *name, const char *outputName);
   ~MpdAnalysisTask2() {}

   /**
    * @brief reads parameters from given file and stores them in map. parameters can be extracted with param(...)
    * function
    *
    * @param fname path to file in txt format. the extension .txt must not be written as it will be attached
    */
   void readParameters(std::string fname);

   void param(std::string name, bool &param, bool defaultvalue);
   void param(std::string name, int &param, int defaultvalue);
   void param(std::string name, float &param, float defaultvalue);
   void param(std::string name, double &param, double defaultvalue);
   void param(std::string name, std::string &param, std::string defaultvalue);

   bool printoutparams = true;

   std::map<std::string, std::string> mMap;

private:
};
#endif
