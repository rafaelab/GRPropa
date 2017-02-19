#include "grpropa/Candidate.h"
#include "grpropa/Cosmology.h"
#include "grpropa/Units.h"

namespace grpropa {


Candidate::Candidate(int id, double E, Vector3d pos, Vector3d dir, double z, double weight) :
        trajectoryLength(0), currentStep(0), nextStep(0), weight(0), active(true) {
    ParticleState state(id, E, pos, dir);
    source = state;
    created = state;
    previous = state;
    current = state;
    setRedshift(z);
    setWeight(weight);
    setCosmicTime(1 / H0() - redshift2LightTravelDistance(z) / c_light);
}

Candidate::Candidate(const ParticleState &state) :
        source(state), created(state), current(state), previous(state), redshift(0), trajectoryLength(0), currentStep(0), nextStep(0), active(true) {
}

bool Candidate::isActive() const {
    return active;
}

void Candidate::setActive(bool b) {
    active = b;
}

double Candidate::getRedshift() const {
    return redshift;
}

double Candidate::getCosmicTime() const {
    return cosmicTime;
}

double Candidate::getTrajectoryLength() const {
    return trajectoryLength;
}

double Candidate::getCurrentStep() const {
    return currentStep;
}

double Candidate::getNextStep() const {
    return nextStep;
}

double Candidate::getWeight() const {
    return weight;
}

void Candidate::setRedshift(double z) {
    redshift = z;
}

void Candidate::setCosmicTime(double T) {
    cosmicTime = T;
}

void Candidate::setTrajectoryLength(double a) {
    trajectoryLength = a;
}

void Candidate::setCurrentStep(double lstep) {
    currentStep = lstep;
    trajectoryLength += lstep;
    cosmicTime += lstep / current.getSpeed();
}

void Candidate::setNextStep(double step) {
    nextStep = step;
}

void Candidate::setWeight(double w) {
    weight = w;
}

void Candidate::limitNextStep(double step) {
    nextStep = std::min(nextStep, step);
}

void Candidate::setProperty(const std::string &name, const std::string &value) {
    properties[name] = value;
}

bool Candidate::getProperty(const std::string &name, std::string &value) const {
    PropertyMap::const_iterator i = properties.find(name);
    if (i == properties.end())
        return false;
    value = i->second;
    return true;
}

bool Candidate::removeProperty(const std::string& name) {
    PropertyMap::iterator i = properties.find(name);
    if (i == properties.end())
        return false;
    properties.erase(i);
    return true;
}

bool Candidate::hasProperty(const std::string &name) const {
    PropertyMap::const_iterator i = properties.find(name);
    if (i == properties.end())
        return false;
    return true;
}

void Candidate::addSecondary(Candidate *c) {
    secondaries.push_back(c);
}

void Candidate::addSecondary(int id, double energy, double weight) {
    ref_ptr<Candidate> secondary = new Candidate;
    secondary->setWeight(weight);
    secondary->setRedshift(redshift);
    secondary->setCosmicTime(cosmicTime);
    secondary->setTrajectoryLength(trajectoryLength);
    secondary->source = source;
    secondary->previous = previous;
    secondary->created = current;
    secondary->current = current;
    secondary->current.setId(id);
    secondary->current.setEnergy(energy);
    secondaries.push_back(secondary);
}

void Candidate::addSecondary(int id, double energy, Vector3d position, double weight) {
    ref_ptr<Candidate> secondary = new Candidate;
    secondary->setRedshift(redshift);
    secondary->setCosmicTime(cosmicTime);
    secondary->setTrajectoryLength(trajectoryLength - (current.getPosition() - position).getR());
    secondary->setWeight(weight);
    secondary->source = source;
    secondary->previous = previous;
    secondary->created = current;
    secondary->current = current;
    secondary->current.setId(id);
    secondary->current.setEnergy(energy);
    secondary->current.setPosition(position);
    secondary->created.setPosition(position);
    secondaries.push_back(secondary);
}

void Candidate::clearSecondaries() {
    secondaries.clear();
}

std::string Candidate::getDescription() const {
    std::stringstream ss;
    ss << "CosmicRay at z = " << getRedshift() << "\n";
    ss << "  source:  " << source.getDescription() << "\n";
    ss << "  current: " << current.getDescription();
    return ss.str();
}

ref_ptr<Candidate> Candidate::clone(bool recursive) const {
    ref_ptr<Candidate> cloned = new Candidate;
    cloned->source = source;
    cloned->created = created;
    cloned->current = current;
    cloned->previous = previous;

    cloned->properties = properties;
    cloned->active = active;
    cloned->weight = weight;
    cloned->redshift = redshift;
    cloned->cosmicTime = cosmicTime;
    cloned->trajectoryLength = trajectoryLength;
    cloned->currentStep = currentStep;
    cloned->nextStep = nextStep;
    if (recursive) {
        cloned->secondaries.reserve(secondaries.size());
        for (size_t i = 0; i < secondaries.size(); i++) {
            cloned->secondaries.push_back(secondaries[i]->clone(recursive));
        }
    }
    return cloned;
}


} // namespace grpropa
