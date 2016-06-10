#include "demBaseFoam.H"
#include "simpleControl.H"

class demSimpleFoam : public demBaseFoam {
public:
  demSimpleFoam();
  virtual ~demSimpleFoam();
  virtual void run(double t);

  volVectorField *ubar_;
  volScalarField *beta_;
  VEC_SETGET(ubar)
  SETGET(beta)

  simpleControl  *simple_;
};
