#ifndef SIMPLEPROPAGATION_H
#define SIMPLEPROPAGATION_H

#include "grpropa/Module.h"
#include "grpropa/Units.h"

namespace grpropa {

/**
 @class SimplePropagation
 @brief Simple rectlinear propagation in absence of magnetic fields.

 This module performs a rectlinear propagation.
 The step size is guaranteed to be larger than minStep and smaller than maxStep.
 It always proposes a next step size of maxStep.
 */
class SimplePropagation: public Module {
private:
    double minStep, maxStep;

public:
    SimplePropagation(double minStep = 0, double maxStep = 10 * Mpc);
    void process(Candidate *candidate) const;
    void setMinimumStep(double minStep);
    void setMaximumStep(double maxStep);
    double getMinimumStep() const;
    double getMaximumStep() const;
    std::string getDescription() const;
};

} // namespace grpropa

#endif // SIMPLEPROPAGATION_H

