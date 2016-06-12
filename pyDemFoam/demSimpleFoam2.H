#include "demBaseFoam.H"
#include "simpleControl.H"

class demSimpleFoam2 : public demBaseFoam {
public:
  demSimpleFoam2();
  virtual ~demSimpleFoam2();
  virtual void run(double t);

  volVectorField *ubar_;
  volScalarField *beta_;
  VEC_SETGET(ubar)
  SETGET(beta)

  simpleControl  *simple_;
};
