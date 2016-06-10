#include "demBaseFoam.H"

class demIcoFoam : public demBaseFoam {
public:
  demIcoFoam();
  virtual ~demIcoFoam();

  virtual void run(double t);
  pisoControl  *piso_;
  volVectorField *f_;
  VEC_SETGET(f)

};
