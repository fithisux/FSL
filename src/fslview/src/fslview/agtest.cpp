#include "atlas.h"

#include <iostream>

using namespace std;

void doLookups(Atlas::Handle& atlas)
{
  cout << "Atlas: " << atlas->inqName() << endl;
  cout << "01: " << atlas->getLabelText(-20.0,   8.0,   8.0) << endl;
  cout << "02: " << atlas->getLabelText(-22.0,   1.0,   9.2) << endl;
  cout << "03: " << atlas->getLabelText(-24.0,  -6.0,  10.4) << endl;
  cout << "04: " << atlas->getLabelText(-26.0, -13.0,  11.6) << endl;
  cout << "05: " << atlas->getLabelText(-28.0, -20.0,  12.8) << endl;
  cout << "06: " << atlas->getLabelText(-30.0, -27.0,  14.0) << endl;
  cout << "07: " << atlas->getLabelText(-32.0, -34.0,  15.2) << endl;
  cout << "08: " << atlas->getLabelText(-34.0, -41.0,  16.4) << endl;
  cout << "09: " << atlas->getLabelText(-36.0, -48.0,  17.6) << endl;
  cout << "10: " << atlas->getLabelText(-38.0, -55.0,  18.8) << endl;
  cout << "11: " << atlas->getLabelText(-40.0, -62.0,  20.0) << endl;
}

int agtest(int argc, char **argv)
{
  AtlasGroup::Handle ag = AtlasGroup::create();

  string atlasdir(string(getenv("FSLDIR")) + "/lib/atlases");

  cout << atlasdir << endl;

  ag->readAtlas(atlasdir, "MNISPMaps.xml");
  ag->readAtlas(atlasdir, "TAL2MNI.xml");
  ag->readAtlas(atlasdir, "JCMaps.xml");
  ag->readAtlas(atlasdir, "JCMax.xml");
  ag->readAtlas(atlasdir, "HOSPM.xml");

  for(AtlasGroup::ConstAtlasIterator it = ag->begin(); it != ag->end(); it++)
    {
      Atlas::Handle b = it->second;
      
      doLookups(b);
    }
     
//   Atlas::Handle b = ag->getAtlasByName("MNI Structural Probability Atlas");
//   doLookups(b);
//   cout << "*******" << endl;
//   b = ag->getAtlasByName("Talairach Daemon Labels");
//   doLookups(b);
//   cout << "*******" << endl;
//   b = ag->getAtlasByName("Juelich Cytoarchitectonic Probability Map");
//   doLookups(b);
//   cout << "*******" << endl;
//   b = ag->getAtlasByName("Juelich Cytoarchitectonic Maximum Probability Map");
//   doLookups(b);
//   cout << "*******" << endl;
//   b = ag->getAtlasByName("Harvard-Oxford Structural Probability Map");
//   doLookups(b);

}
