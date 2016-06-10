#include "demBaseFoam.H"
#include "pisoControl.H"

class demIcoFoam : public demBaseFoam {
public:
  demIcoFoam();
  virtual ~demIcoFoam();
  virtual void run(double t);

  volVectorField *f_;
  VEC_SETGET(f)

  pisoControl  *piso_;
};
