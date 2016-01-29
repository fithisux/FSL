#include "atlas.h"

#include <iostream>

using namespace std;

int atlastst(int argc, char **argv)
{
  ProbabalisticAtlas::Handle a = 
    ProbabalisticAtlas::create(Image::load("/usr/people/flitney/atlases/mnispm/SPMaps.int.2mm"), "SP Maps");
  
  cout << a->getLabelText(-20.0,   8.0,   8.0) << endl;
  cout << a->getLabelText(-22.0,   1.0,   9.2) << endl;
  cout << a->getLabelText(-24.0,  -6.0,  10.4) << endl;
  cout << a->getLabelText(-26.0, -13.0,  11.6) << endl;
  cout << a->getLabelText(-28.0, -20.0,  12.8) << endl;
  cout << a->getLabelText(-30.0, -27.0,  14.0) << endl;
  cout << a->getLabelText(-32.0, -34.0,  15.2) << endl;
  cout << a->getLabelText(-34.0, -41.0,  16.4) << endl;
  cout << a->getLabelText(-36.0, -48.0,  17.6) << endl;
  cout << a->getLabelText(-38.0, -55.0,  18.8) << endl;
  cout << a->getLabelText(-40.0, -62.0,  20.0) << endl;

  cout << "Adding labels" << endl;

  a->addLabel(0, "Caudate");
  a->addLabel(1, "Cerebellum");
  a->addLabel(2, "Frontal Lobe");
  a->addLabel(3, "Insula");
  a->addLabel(4, "Occipital Lobe");
  a->addLabel(5, "Parietal Lobe");
  a->addLabel(6, "Putamen");
  a->addLabel(7, "Temporal Lobe");
  a->addLabel(8, "Thalamus");

  cout << a->getLabelText(-20.0,   8.0,   8.0) << endl;
  cout << a->getLabelText(-22.0,   1.0,   9.2) << endl;
  cout << a->getLabelText(-24.0,  -6.0,  10.4) << endl;
  cout << a->getLabelText(-26.0, -13.0,  11.6) << endl;
  cout << a->getLabelText(-28.0, -20.0,  12.8) << endl;
  cout << a->getLabelText(-30.0, -27.0,  14.0) << endl;
  cout << a->getLabelText(-32.0, -34.0,  15.2) << endl;
  cout << a->getLabelText(-34.0, -41.0,  16.4) << endl;
  cout << a->getLabelText(-36.0, -48.0,  17.6) << endl;
  cout << a->getLabelText(-38.0, -55.0,  18.8) << endl;
  cout << a->getLabelText(-40.0, -62.0,  20.0) << endl;
}
