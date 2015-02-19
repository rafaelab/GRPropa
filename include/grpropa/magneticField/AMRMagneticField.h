#ifndef GRPROPA_AMRMAGNETICFIELD_H
#define GRPROPA_AMRMAGNETICFIELD_H

#ifdef GRPROPA_HAVE_SAGA

#include <iostream> 
#include <string>
#include <cstdio>

#ifdef _OPENMP
    #include "omp.h"
#endif

#include "grpropa/Units.h"
#include "grpropa/magneticField/MagneticField.h"
#include "grpropa/Vector3.h"

#include "saga/LocalProperties.h"
#include "saga/AMRgrid.h"
#include "saga/MagneticField.h"
#include "saga/Referenced.h"



namespace grpropa {

/**
 @class AMRMagneticField
 @brief Wrapper for saga::MagneticField
 */
class AMRMagneticField: public MagneticField {

private:
	saga::ref_ptr<saga::MagneticField> field;
    double cfLength;
    double cfDensity;
    double cfMagneticField;

public:        
  
    AMRMagneticField(saga::ref_ptr<saga::MagneticField> field_, double convLength, double convDensity, double convMagneticField)            
    {
        field = field_;
        cfLength = convLength;
        cfDensity = convDensity;
        cfMagneticField = convMagneticField;
    }

    Vector3d getField(const Vector3d &position) const {

        double x = position.x/cfLength;
        double y = position.y/cfLength;
        double z = position.z/cfLength;

        std::vector<double> b;
        #ifdef _OPENMP
            #pragma omp critical 
            {    
		        b = field->getField(x, y, z);
            }
        #else 
            b = field->getField(x, y, z);
        #endif

        for(int i=0; i<3; i++)
            b[i]*=cfMagneticField;
        //std::cout << x*cfLength/Mpc << "  " << y*cfLength/Mpc << "  " << z*cfLength/Mpc << "  ||   " << b[0] << " " << b[1] << "  " << b[0] << std::endl;

		return Vector3d(b[0], b[1], b[2]) * tesla;   
    }

};

} // namespace grpropa

#endif // GRPROPA_HAVE_SAGA
#endif // GRPROPA_AMRMAGNETICFIELD_H
