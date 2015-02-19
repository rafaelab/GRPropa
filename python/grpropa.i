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

#ifdef CRPROPA_HAVE_SAGA
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
#include "grpropa/module/OutputTXT.h"
#include "grpropa/module/OutputShell.h"
#include "grpropa/module/SimplePropagation.h"
#include "grpropa/module/PropagationCK.h"
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
%include "grpropa/module/OutputTXT.h"
%include "grpropa/module/OutputShell.h"
%include "grpropa/module/Synchrotron.h"
%include "grpropa/module/InverseCompton.h"
%include "grpropa/module/PairProduction.h"
%include "grpropa/module/Redshift.h"
%include "grpropa/module/Tools.h"

%template(SourceRefPtr) grpropa::ref_ptr<grpropa::Source>;
%feature("director") grpropa::Source;
%template(SourceFeatureRefPtr) grpropa::ref_ptr<grpropa::SourceFeature>;
%feature("director") grpropa::SourceFeature;
%include "grpropa/Source.h"

%template(ModuleListRefPtr) grpropa::ref_ptr<grpropa::ModuleList>;
%include "grpropa/ModuleList.h"


// pretty print
%pythoncode %{
ParticleState.__repr__ = ParticleState.getDescription
Candidate.__repr__ = Candidate.getDescription

Module.__repr__ = Module.getDescription
ModuleList.__repr__ = ModuleList.getDescription

Source.__repr__ = Source.getDescription
SourceList.__repr__ = SourceList.getDescription
SourceParticleType.__repr__ = SourceParticleType.getDescription
SourceMultipleParticleTypes.__repr__ = SourceMultipleParticleTypes.getDescription
SourceEnergy.__repr__ = SourceEnergy.getDescription
SourcePowerLawSpectrum.__repr__ = SourcePowerLawSpectrum.getDescription
SourceMultiplePositions.__repr__ = SourceMultiplePositions.getDescription
SourcePosition.__repr__ = SourcePosition.getDescription
SourceUniform1D.__repr__ = SourceUniform1D.getDescription
SourceUniformBox.__repr__ = SourceUniformBox.getDescription
SourceUniformShell.__repr__ = SourceUniformShell.getDescription
SourceUniformSphere.__repr__ = SourceUniformSphere.getDescription
SourceDensityGrid.__repr__ = SourceDensityGrid.getDescription
SourceDensityGrid1D.__repr__ = SourceDensityGrid1D.getDescription
SourceDirection.__repr__ = SourceDirection.getDescription
SourceIsotropicEmission.__repr__ = SourceIsotropicEmission.getDescription
SourceEmissionCone.__repr__ = SourceEmissionCone.getDescription
SourceRedshift.__repr__ = SourceRedshift.getDescription
SourceRedshift1D.__repr__ = SourceRedshift1D.getDescription
SourceUniformRedshift.__repr__ = SourceUniformRedshift.getDescription

Observer.__repr__ = Observer.getDescription
ObserverPoint.__repr__ = ObserverPoint.getDescription
ObserverSmallSphere.__repr__ = ObserverSmallSphere.getDescription
ObserverLargeSphere.__repr__ = ObserverLargeSphere.getDescription
ObserverRedshiftWindow.__repr__ = ObserverRedshiftWindow.getDescription
ObserverNeutrinoVeto.__repr__ = ObserverNeutrinoVeto.getDescription
ObserverPhotonVeto.__repr__ = ObserverPhotonVeto.getDescription
ObserverOutput1D.__repr__ = ObserverOutput1D.getDescription
ObserverOutput3D.__repr__ = ObserverOutput3D.getDescription

def Vector3__repr__(self):
    return "Vector(%.3g, %.3g, %.3g)" % (self.x, self.y, self.z)
Vector3d.__repr__ = Vector3__repr__
Vector3f.__repr__ = Vector3__repr__

DeflectionCK = PropagationCK  # legacy name
%}
