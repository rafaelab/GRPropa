#ifndef GRPROPA_OBSERVER_H
#define GRPROPA_OBSERVER_H

#include <fstream>
#include <limits>
#include <string>
#include <vector>

#include "../Candidate.h"
#include "../Module.h"
#include "../Referenced.h"
#include "../Vector3.h"

namespace grpropa {

enum DetectionState {
    DETECTED, VETO, NOTHING
};

/**
 @class ObserverFeature
 @brief Abstract base class for features of cosmic ray observers
 */
class ObserverFeature: public Referenced {
protected:
    std::string description;
public:
    virtual DetectionState checkDetection(Candidate *candidate) const;
    virtual void onDetection(Candidate *candidate) const;
    virtual std::string getDescription() const;
};

/**
 @class Observer
 @brief General cosmic ray observer
 */
class Observer: public Module {
private:
    std::vector<ref_ptr<ObserverFeature> > features;
    bool makeInactive;
public:
    Observer(bool makeInactive = true);
    void add(ObserverFeature *property);
    void process(Candidate *candidate) const;
    std::string getDescription() const;
};

/**
 @class ObserverSmallSphere
 @brief Detects particles upon entering a sphere
 */
class ObserverSmallSphere: public ObserverFeature {
private:
    Vector3d center;
    double radius;
public:
    ObserverSmallSphere(Vector3d center = Vector3d(0.), double radius = 0);
    DetectionState checkDetection(Candidate *candidate) const;
    std::string getDescription() const;
};

/**
 @class ObserverLargeSphere
 @brief Detects particles upon exiting a sphere
 */
class ObserverLargeSphere: public ObserverFeature {
private:
    Vector3d center;
    double radius;
public:
    ObserverLargeSphere(Vector3d center = Vector3d(0.), double radius = 0);
    DetectionState checkDetection(Candidate *candidate) const;
    std::string getDescription() const;
};

/**
 @class ObserverRedshiftWindow
 @brief Detects particles in a given redshift window
 */
class ObserverRedshiftWindow: public ObserverFeature {
private:
    double zmin, zmax;
public:
    ObserverRedshiftWindow(double zmin = 0, double zmax = 0.1);
    DetectionState checkDetection(Candidate *candidate) const;
    std::string getDescription() const;
};

/**
 @class ObserverPoint
 @brief Detects particles when reaching x = 0

 Should be renamed to Observer1D
 */
class ObserverPoint: public ObserverFeature {
public:
    DetectionState checkDetection(Candidate *candidate) const;
    std::string getDescription() const;
};

/**
 @class ObserverNeutrinoVeto
 @brief Veto for neutrinos
 */
class ObserverNeutrinoVeto: public ObserverFeature {
public:
    DetectionState checkDetection(Candidate *candidate) const;
    std::string getDescription() const;
};

/**
 @class ObserverChargedLeptonVeto
 @brief Veto for charged leptons
 */
class ObserverChargedLeptonVeto: public ObserverFeature {
public:
    DetectionState checkDetection(Candidate *candidate) const;
    std::string getDescription() const;
};

/**
 @class ObserverPhotonVeto
 @brief Veto for photons
 */
class ObserverPhotonVeto: public ObserverFeature {
public:
    DetectionState checkDetection(Candidate *candidate) const;
    std::string getDescription() const;
};

/**
 @class ObserverOutput3D
 @brief Plain text output of 3D properties
 */
class ObserverOutput3D: public ObserverFeature {
private:
    mutable std::ofstream fout;
public:
    ObserverOutput3D(std::string filename);
    ~ObserverOutput3D();
    void onDetection(Candidate *candidate) const;
};

/**
 @class ObserverOutput1D
 @brief Plain text output of 1D properties
 */
class ObserverOutput1D: public ObserverFeature {
private:
    mutable std::ofstream fout;
public:
    ObserverOutput1D(std::string filename);
    ~ObserverOutput1D();
    void onDetection(Candidate *candidate) const;
};


} // namespace grpropa

#endif // GRPROPA_OBSERVER_H