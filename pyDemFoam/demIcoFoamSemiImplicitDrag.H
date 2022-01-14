#include "demBaseFoam.H"
#include "pisoControl.H"

class demIcoFoamSemiImplicitDrag : public demBaseFoam {
public:
  demIcoFoamSemiImplicitDrag();
  virtual ~demIcoFoamSemiImplicitDrag();
  virtual void run(double t);

  volVectorField *ubar_;
  volScalarField *beta_;
  VEC_SETGET(ubar)
  SETGET(beta)

  pisoControl  *piso_;
};
