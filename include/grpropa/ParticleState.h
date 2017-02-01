#ifndef GRPROPA_PARTICLE_STATE_H
#define GRPROPA_PARTICLE_STATE_H

#include "grpropa/Vector3.h"

#include <sstream>

namespace grpropa {

/**
 @class ParticleState
 @brief State of the particle: ID, energy, position, direction

 The ParticleState defines the state of an ultra-high energy cosmic ray, which
 is assumed to be travelling at the exact speed of light.
 The cosmic ray state is defined by particle ID, energy and position and
 direction vector.
 For faster lookup mass and charge of the particle are stored as members.
 */
class ParticleState {
private:
    int id; /* particle ID (Particle Data Group numbering scheme) */
    double energy; /* total energy */
    Vector3d position; /* position vector in comoving coordinates */ 
    Vector3d direction; /* unit vector of velocity or momentum */
    double pmass; /* particle rest mass */ 
    double charge; /* particle charge */

public:
    ParticleState(int id = 0, double energy = 0, Vector3d position = Vector3d(0, 0, 0), Vector3d direction = Vector3d(-1, 0, 0));

    std::string getDescription() const;
    void setPosition(const Vector3d &pos); //* Set position in comoving coordinates */
    const Vector3d &getPosition() const; /* Get position in comoving coordinates */
    void setDirection(const Vector3d &dir); /* Set direction unit vector, non unit-vectors are normalized */
    const Vector3d &getDirection() const; /* Get direction unit vector */
    void setEnergy(double newEnergy); /* Set energy in [J] */
    double getEnergy() const; /* Get energy in [J] */
    void setId(int); /* Set particle ID */
    int getId() const; /* Get particle ID */


    double getCharge() const; /* Electric charge of the particle in [C] */
    double getMass() const; /* Mass of the particle in [kg] */
    void setLorentzFactor(double gamma); /* Set Lorentz factor and modify the particle's energy accordingly */
    double getLorentzFactor() const; /* Get Lorentz factor */
    Vector3d getVelocity() const; /* Velocity: direction times the speed of light in [m/s] */
    Vector3d getMomentum() const; /* Momentum: direction times energy divided by the speed of light [kg m/s] */
    double getSpeed() const; /* Returns the speed [m/s] */
};

} // namespace grpropa

#endif // GRPROPA_PARTICLE_STATE_H
