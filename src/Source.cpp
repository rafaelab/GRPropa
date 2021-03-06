#include "grpropa/Source.h"
#include "grpropa/Random.h"
#include "grpropa/Cosmology.h"
#include "grpropa/Common.h"
#include "grpropa/Units.h"

#include <sstream>
#include <stdexcept>

namespace grpropa {

// Source ---------------------------------------------------------------------
void Source::add(SourceFeature* property) {
    features.push_back(property);
}

ref_ptr<Candidate> Source::getCandidate() const {
    ref_ptr<Candidate> candidate = new Candidate();
    for (int i = 0; i < features.size(); i++)
        (*features[i]).prepareCandidate(*candidate);
    return candidate;
}

std::string Source::getDescription() const {
    std::stringstream ss;
    ss << "Cosmic ray source\n";
    for (int i = 0; i < features.size(); i++)
        ss << "    " << features[i]->getDescription();
    return ss.str();
}

// SourceList------------------------------------------------------------------
void SourceList::add(Source* source, double weight) {
    sources.push_back(source);
    if (cdf.size() > 0)
        weight += cdf.back();
    cdf.push_back(weight);
}

ref_ptr<Candidate> SourceList::getCandidate() const {
    if (sources.size() == 0)
        throw std::runtime_error("SourceList: no sources set");
    size_t i = Random::instance().randBin(cdf);
    return (sources[i])->getCandidate();
}

std::string SourceList::getDescription() const {
    std::stringstream ss;
    ss << "List of cosmic ray sources\n";
    for (int i = 0; i < sources.size(); i++)
        ss << "  " << sources[i]->getDescription();
    return ss.str();
}

// SourceFeature---------------------------------------------------------------
void SourceFeature::prepareCandidate(Candidate& candidate) const {
    ParticleState &source = candidate.source;
    prepareParticle(source);
    candidate.created = source;
    candidate.current = source;
    candidate.previous = source;
}

std::string SourceFeature::getDescription() const {
    return description;
}

// ----------------------------------------------------------------------------
SourceParticleType::SourceParticleType(int id) :
        id(id) {
    setDescription();
}

void SourceParticleType::prepareParticle(ParticleState& particle) const {
    particle.setId(id);
}

void SourceParticleType::setDescription() {
    std::stringstream ss;
    ss << "SourceParticleType: " << id << "\n";
    description = ss.str();
}

// ----------------------------------------------------------------------------
SourceEnergy::SourceEnergy(double energy) :
        E(energy) {
    setDescription();
}

void SourceEnergy::prepareParticle(ParticleState& p) const {
    p.setEnergy(E);
}

void SourceEnergy::setDescription() {
    std::stringstream ss;
    ss << "SourceEnergy: " << E / EeV << " EeV\n";
    description = ss.str();
}

// ----------------------------------------------------------------------------
SourcePowerLawSpectrum::SourcePowerLawSpectrum(double Emin, double Emax,
        double index) :
        Emin(Emin), Emax(Emax), index(index) {
    setDescription();
}

void SourcePowerLawSpectrum::prepareParticle(ParticleState& particle) const {
    Random &random = Random::instance();
    double E = random.randPowerLaw(index, Emin, Emax);
    particle.setEnergy(E);
}

void SourcePowerLawSpectrum::setDescription() {
    std::stringstream ss;
    ss << "SourcePowerLawSpectrum: Random energy ";
    ss << "E = " << Emin / EeV << " - " << Emax / EeV << " EeV, ";
    ss << "dN/dE ~ E^" << index  << "\n";
    description = ss.str();
}

// ----------------------------------------------------------------------------
SourcePosition::SourcePosition(Vector3d position) :
        position(position) {
    setDescription();
}

SourcePosition::SourcePosition(double d) :
        position(Vector3d(d, 0, 0)) {
    setDescription();
}

void SourcePosition::prepareParticle(ParticleState& particle) const {
    particle.setPosition(position);
}

void SourcePosition::setDescription() {
    std::stringstream ss;
    ss << "SourcePosition: " << position / Mpc << " Mpc\n";
    description = ss.str();
}

// ----------------------------------------------------------------------------
SourceMultiplePositions::SourceMultiplePositions() {
    setDescription();
}

void SourceMultiplePositions::add(Vector3d pos, double weight) {
    positions.push_back(pos);
    if (cdf.size() > 0)
        weight += cdf.back();
    cdf.push_back(weight);
}

void SourceMultiplePositions::prepareParticle(ParticleState& particle) const {
    if (positions.size() == 0)
        throw std::runtime_error("SourceMultiplePositions: no position set");
    size_t i = Random::instance().randBin(cdf);
    particle.setPosition(positions[i]);
}

void SourceMultiplePositions::setDescription() {
    std::stringstream ss;
    ss << "SourceMultiplePositions: Random position from list\n";
    for (int i = 0; i < positions.size(); i++)
        ss << "  " << positions[i] / Mpc << " Mpc\n";
    description = ss.str();
}

// ----------------------------------------------------------------------------
SourceUniformSphere::SourceUniformSphere(Vector3d center, double radius) :
        center(center), radius(radius) {
    setDescription();
}

void SourceUniformSphere::prepareParticle(ParticleState& particle) const {
    Random &random = Random::instance();
    double r = pow(random.rand(), 1. / 3.) * radius;
    particle.setPosition(center + random.randVector() * r);
}

void SourceUniformSphere::setDescription() {
    std::stringstream ss;
    ss << "SourceUniformSphere: Random position within a sphere at ";
    ss << center / Mpc << " Mpc with";
    ss  << radius / Mpc << " Mpc radius\n";
    description = ss.str();
}

// ----------------------------------------------------------------------------
SourceUniformShell::SourceUniformShell(Vector3d center, double radius) :
        center(center), radius(radius) {
    setDescription();
}

void SourceUniformShell::prepareParticle(ParticleState& particle) const {
    Random &random = Random::instance();
    particle.setPosition(random.randVector() * radius);
}

void SourceUniformShell::setDescription() {
    std::stringstream ss;
    ss << "SourceUniformShell: Random position on a spherical shell at ";
    ss << center / Mpc << " Mpc with ";
    ss << radius / Mpc << " Mpc radius\n";
    description = ss.str();
}

// ----------------------------------------------------------------------------
SourceUniformBox::SourceUniformBox(Vector3d origin, Vector3d size) :
        origin(origin), size(size) {
    setDescription();
}

void SourceUniformBox::prepareParticle(ParticleState& particle) const {
    Random &random = Random::instance();
    Vector3d pos(random.rand(), random.rand(), random.rand());
    particle.setPosition(pos * size + origin);
}

void SourceUniformBox::setDescription() {
    std::stringstream ss;
    ss << "SourceUniformBox: Random uniform position in box with ";
    ss << "origin = " << origin / Mpc << " Mpc and ";
    ss << "size = " << size / Mpc << " Mpc\n";
    description = ss.str();
}

// ----------------------------------------------------------------------------
SourceUniform1D::SourceUniform1D(double minD, double maxD, bool withCosmology) {
    this->withCosmology = withCosmology;
    if (withCosmology) {
        this->minD = comoving2LightTravelDistance(minD);
        this->maxD = comoving2LightTravelDistance(maxD);
    } else {
        this->minD = minD;
        this->maxD = maxD;
    }
    setDescription();
}

void SourceUniform1D::prepareParticle(ParticleState& particle) const {
    Random& random = Random::instance();
    double d = random.rand() * (maxD - minD) + minD;
    if (withCosmology)
        d = lightTravel2ComovingDistance(d);
    particle.setPosition(Vector3d(d, 0, 0));
}

void SourceUniform1D::setDescription() {
    std::stringstream ss;
    ss << "SourceUniform1D: Random uniform position in D = " << minD << " - " << maxD;
    if (withCosmology)
        ss << " (including cosmology)";
    ss << "\n";
    description = ss.str();
}

// ----------------------------------------------------------------------------
SourceDensityGrid::SourceDensityGrid(ref_ptr<ScalarGrid> grid) :
        grid(grid) {
    float sum = 0;
    for (int ix = 0; ix < grid->getNx(); ix++) {
        for (int iy = 0; iy < grid->getNy(); iy++) {
            for (int iz = 0; iz < grid->getNz(); iz++) {
                sum += grid->get(ix, iy, iz);
                grid->get(ix, iy, iz) = sum;
            }
        }
    }
    setDescription();
}

void SourceDensityGrid::prepareParticle(ParticleState& particle) const {
    Random &random = Random::instance();

    // draw random bin
    size_t i = random.randBin(grid->getGrid());
    Vector3d pos = grid->positionFromIndex(i);

    // draw uniform position within bin
    double dx = random.rand() - 0.5;
    double dy = random.rand() - 0.5;
    double dz = random.rand() - 0.5;
    pos += Vector3d(dx, dy, dz) * grid->getSpacing();

    particle.setPosition(pos);
}

void SourceDensityGrid::setDescription() {
    description = "SourceDensityGrid: 3D source distribution according to density grid\n";
}

// ----------------------------------------------------------------------------
SourceDensityGrid1D::SourceDensityGrid1D(ref_ptr<ScalarGrid> grid) :
        grid(grid) {
    if (grid->getNy() != 1)
        throw std::runtime_error("SourceDensityGrid1D: Ny != 1");
    if (grid->getNz() != 1)
        throw std::runtime_error("SourceDensityGrid1D: Nz != 1");

    float sum = 0;
    for (int ix = 0; ix < grid->getNx(); ix++) {
        sum += grid->get(ix, 0, 0);
        grid->get(ix, 0, 0) = sum;
    }
    setDescription();
}

void SourceDensityGrid1D::prepareParticle(ParticleState& particle) const {
    Random &random = Random::instance();

    // draw random bin
    size_t i = random.randBin(grid->getGrid());
    Vector3d pos = grid->positionFromIndex(i);

    // draw uniform position within bin
    double dx = random.rand() - 0.5;
    pos.x += dx * grid->getSpacing();

    particle.setPosition(pos);
}

void SourceDensityGrid1D::setDescription() {
    description = "SourceDensityGrid1D: 1D source distribution according to density grid\n";
}

// ----------------------------------------------------------------------------
SourceIsotropicEmission::SourceIsotropicEmission() {
    setDescription();
}

void SourceIsotropicEmission::prepareParticle(ParticleState& particle) const {
    Random &random = Random::instance();
    particle.setDirection(random.randVector());
}

void SourceIsotropicEmission::setDescription() {
    description = "SourceIsotropicEmission: Random isotropic direction\n";
}

// ----------------------------------------------------------------------------
SourceDirection::SourceDirection(Vector3d direction) :
        direction(direction) {
    setDescription();
}

void SourceDirection::prepareParticle(ParticleState& particle) const {
    particle.setDirection(direction);
}

void SourceDirection::setDescription() {
    std::stringstream ss;
    ss <<  "SourceDirection: Emission direction = " << direction << "\n";
    description = ss.str();
}

// ----------------------------------------------------------------------------
SourceEmissionCone::SourceEmissionCone(Vector3d direction, double aperture) :
        direction(direction), aperture(aperture) {
    setDescription();
}

void SourceEmissionCone::prepareParticle(ParticleState& particle) const {
    Random &random = Random::instance();
    particle.setDirection(random.randConeVector(direction, aperture));
}

void SourceEmissionCone::setDescription() {
    std::stringstream ss;
    ss << "SourceEmissionCone: Jetted emission in ";
    ss << "direction = " << direction << " with ";
    ss << "half-opening angle = " << aperture << " rad\n";
    description = ss.str();
}

// ----------------------------------------------------------------------------
SourceRedshift::SourceRedshift(double z) :
        z(z) {
    setDescription();
}

void SourceRedshift::prepareCandidate(Candidate& candidate) const {
    candidate.setRedshift(z);
}

void SourceRedshift::setDescription() {
    std::stringstream ss;
    ss << "SourceRedshift: Redshift z = " << z << "\n";
    description = ss.str();
}

// ----------------------------------------------------------------------------
SourceUniformRedshift::SourceUniformRedshift(double zmin, double zmax) :
        zmin(zmin), zmax(zmax) {
    setDescription();
}

void SourceUniformRedshift::prepareCandidate(Candidate& candidate) const {
    double z = Random::instance().randUniform(zmin, zmax);
    candidate.setRedshift(z);
}

void SourceUniformRedshift::setDescription() {
    std::stringstream ss;
    ss << "SourceUniformRedshift: Uniform redshift in z = ";
    ss << zmin << " - " << zmax << "\n";
    description = ss.str();
}

// ----------------------------------------------------------------------------
SourceRedshiftEvolution::SourceRedshiftEvolution(double m, double zmin, double zmax) : m(m), zmin(zmin), zmax(zmax) {
    std::stringstream ss;
    ss << "SourceRedshiftEvolution: (1+z)^m, m = " << m;
    ss << ", z = " << zmin << " - " << zmax << "\n";
    description = ss.str();
}

void SourceRedshiftEvolution::prepareCandidate(Candidate& candidate) const {
    double x = Random::instance().randUniform(0, 1);
    double norm, z;

    // special case: m=-1
    if ((std::abs(m+1)) < std::numeric_limits<double>::epsilon()) {
        norm = log(1+zmax) - log(1+zmin);
        z = exp(norm*x) - 1;
    } else {
        norm = ( pow(1+zmax, m+1) - pow(1+zmin, m+1) ) / (m+1);
        z = pow( norm*(m+1)*x + pow(1+zmin, m+1), 1./(m+1)) - 1;
    }
    candidate.setRedshift(z);
}

// ----------------------------------------------------------------------------
SourceRedshift1D::SourceRedshift1D() {
    setDescription();
}

void SourceRedshift1D::prepareCandidate(Candidate& candidate) const {
    double d = candidate.source.getPosition().getR();
    double z = comovingDistance2Redshift(d);
    candidate.setRedshift(z);
}

void SourceRedshift1D::setDescription() {
    description = "SourceRedshift1D: Redshift according to source distance\n";
}

} // namespace grpropa
