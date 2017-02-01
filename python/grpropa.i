%module(directors="1") grpropa
%feature("autodoc", "1"); // automatic docstrings

%{
// workaround for SWIG < 2.0.5 with GCC >= 4.7
#include <cstddef>
using std::ptrdiff_t;
%}

%include stl.i
%include std_set.i
%include std_multiset.i
%include std_map.i
%include std_pair.i
%include std_multimap.i
%include std_vector.i
%include std_string.i
%include std_list.i
%include stdint.i
%include std_container.i
%include "exception.i"

#ifdef grpropa_HAVE_SAGA
    %import (module="saga") saga.i
#endif

%{
#include "grpropa/module/InverseCompton.h"
#include "grpropa/module/PairProduction.h"
#include "grpropa/module/Synchrotron.h"
#include "grpropa/module/Redshift.h"
#include "grpropa/module/BreakCondition.h"
#include "grpropa/module/Boundary.h"
#include "grpropa/module/Observer.h"
#include "grpropa/module/Output.h"
#include "grpropa/module/SimplePropagation.h"
#include "grpropa/module/PropagationCK.h"
#include "grpropa/module/TextOutput.h"
#include "grpropa/module/Tools.h"

#include "grpropa/magneticField/MagneticField.h"
#include "grpropa/magneticField/MagneticFieldGrid.h"
#include "grpropa/magneticField/AMRMagneticField.h"
#include "grpropa/magneticField/JF12Field.h"
#include "grpropa/magneticField/TurbulentMagneticField.h"

#include "grpropa/Referenced.h"
#include "grpropa/Candidate.h"
#include "grpropa/ParticleState.h"
#include "grpropa/Module.h"
#include "grpropa/ModuleList.h"
#include "grpropa/Random.h"
#include "grpropa/Units.h"
#include "grpropa/Vector3.h"
#include "grpropa/Source.h"
#include "grpropa/Common.h"
#include "grpropa/Cosmology.h"
#include "grpropa/PhotonBackground.h"
#include "grpropa/Grid.h"
#include "grpropa/GridTools.h"
%}


%exception
{
 try
 {
   $action
 }
 catch (Swig::DirectorException &e) {
   SWIG_exception(SWIG_RuntimeError, e.getMessage());
 }
 catch (const std::exception& e) {
   SWIG_exception(SWIG_RuntimeError, e.what());
 }
 catch (const char *e) {
   SWIG_exception(SWIG_RuntimeError, e);
 }
}

%ignore operator<<;
%ignore operator>>;
%ignore *::operator=;
%ignore operator grpropa::Source*;
%ignore operator grpropa::SourceFeature*;
%ignore operator grpropa::Candidate*;
%ignore operator grpropa::Module*;
%ignore operator grpropa::ModuleList*;
%ignore operator grpropa::MagneticField*;

%feature("ref")   grpropa::Referenced "$this->addReference();"
%feature("unref") grpropa::Referenced "$this->removeReference();"


%include "grpropa/Vector3.h"
%template(Vector3d) grpropa::Vector3<double>;
%template(Vector3f) grpropa::Vector3<float>;

%include "grpropa/Referenced.h"
%include "grpropa/Units.h"
%include "grpropa/Common.h"
%include "grpropa/Cosmology.h"
%include "grpropa/PhotonBackground.h"
%include "grpropa/Random.h"
%include "grpropa/ParticleState.h"


%extend grpropa::Candidate {
    PyObject * getProperty(PyObject * name){

        std::string value;
        std::string input;

        if (PyString_Check( name )){ 
            input = PyString_AsString( name );
        } else {
            std::cerr << "ERROR: The argument of getProperty() must be a string!" << std::endl;
            return NULL;
        }   
        $self->getProperty( input, value );

        return PyString_FromString( value.c_str() );
    }   
};
%template(CandidateVector) std::vector< grpropa::ref_ptr<grpropa::Candidate> >;
%template(CandidateRefPtr) grpropa::ref_ptr<grpropa::Candidate>;
%include "grpropa/Candidate.h"


%template(ModuleRefPtr) grpropa::ref_ptr<grpropa::Module>;
%template(stdModuleList) std::list< grpropa::ref_ptr<grpropa::Module> >;
%feature("director") grpropa::Module;
%include "grpropa/Module.h"

%implicitconv grpropa::ref_ptr<grpropa::MagneticField>;
%template(MagneticFieldRefPtr) grpropa::ref_ptr<grpropa::MagneticField>;
%include "grpropa/magneticField/MagneticField.h"

%include "grpropa/Grid.h"
%include "grpropa/GridTools.h"

%implicitconv grpropa::ref_ptr<grpropa::Grid<grpropa::Vector3<float> > >;
%template(VectorGridRefPtr) grpropa::ref_ptr<grpropa::Grid<grpropa::Vector3<float> > >;
%template(VectorGrid) grpropa::Grid<grpropa::Vector3<float> >;

%implicitconv grpropa::ref_ptr<grpropa::Grid<float> >;
%template(ScalarGridRefPtr) grpropa::ref_ptr<grpropa::Grid<float> >;
%template(ScalarGrid) grpropa::Grid<float>;

%include "grpropa/magneticField/MagneticFieldGrid.h"
%include "grpropa/magneticField/AMRMagneticField.h"
%include "grpropa/magneticField/JF12Field.h"
%include "grpropa/magneticField/TurbulentMagneticField.h"

%include "grpropa/module/BreakCondition.h"
%include "grpropa/module/Boundary.h"
%include "grpropa/module/Observer.h"
%include "grpropa/module/SimplePropagation.h"
%include "grpropa/module/PropagationCK.h"
%include "grpropa/module/Output.h"
%include "grpropa/module/Synchrotron.h"
%include "grpropa/module/InverseCompton.h"
%include "grpropa/module/PairProduction.h"
%include "grpropa/module/Redshift.h"
%include "grpropa/module/TextOutput.h"
%include "grpropa/module/Tools.h"

%template(SourceRefPtr) grpropa::ref_ptr<grpropa::Source>;
%feature("director") grpropa::Source;
%template(SourceFeatureRefPtr) grpropa::ref_ptr<grpropa::SourceFeature>;
%feature("director") grpropa::SourceFeature;
%include "grpropa/Source.h"

%template(ModuleListRefPtr) grpropa::ref_ptr<grpropa::ModuleList>;
%include "grpropa/ModuleList.h"

%template(IntSet) std::set<int>;
%include "grpropa/module/Tools.h"

__REPR__(grpropa::ParticleState);
__REPR__(grpropa::Candidate);
__REPR__(grpropa::Module);
__REPR__(grpropa::ModuleList);
__REPR__(grpropa::Source);
__REPR__(grpropa::SourceList);
__REPR__(grpropa::SourceFeature);
__REPR__(grpropa::Observer);
__REPR__(grpropa::ObserverFeature);

VECTOR3__REPR__(grpropa::Vector3);

%pythoncode %{
    DeflectionCK = PropagationCK  # legacy name
%}

%pythoncode %{
    def Vector3__repr__(self):
        return "Vector(%.3g, %.3g, %.3g)" % (self.x, self.y, self.z)
    Vector3d.__repr__ = Vector3__repr__
    Vector3f.__repr__ = Vector3__repr__
%}

namespace std {
   %template(Vector1i) vector<int>;
   %template(Vector1d) vector<double>;
   %template(Vector1f) vector<float>;
};
