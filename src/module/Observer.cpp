#include "grpropa/module/Observer.h"
#include "grpropa/Units.h"
#include "grpropa/Cosmology.h"

namespace grpropa {

DetectionState ObserverFeature::checkDetection(Candidate *candidate) const {
    return NOTHING;
}

void ObserverFeature::onDetection(Candidate *candidate) const {
}

std::string ObserverFeature::getDescription() const {
    return description;
}

Observer::Observer(bool makeInactive) :
        makeInactive(makeInactive) {
}

void Observer::add(ObserverFeature *feature) {
    features.push_back(feature);
}

void Observer::process(Candidate *candidate) const {
    // loop over all features and have them check the particle
    DetectionState state = NOTHING;
    for (int i = 0; i < features.size(); i++) {
        DetectionState s = features[i]->checkDetection(candidate);
        if (s == VETO)
            state = VETO;
        else if ((s == DETECTED) && (state != VETO))
            state = DETECTED;
    }

    if (state == DETECTED) {
        for (int i = 0; i < features.size(); i++) {
            features[i]->onDetection(candidate);
        }

        if (makeInactive)
            candidate->setActive(false);
    }
}

std::string Observer::getDescription() const {
    std::stringstream ss;
    ss << "Observer\n";
    for (int i = 0; i < features.size(); i++)
        ss << "    " << features[i]->getDescription() << "\n";
    return ss.str();
}

ObserverSmallSphere::ObserverSmallSphere(Vector3d center, double radius) :
        center(center), radius(radius) {
}

DetectionState ObserverSmallSphere::checkDetection(Candidate *candidate) const {
    // current distance to observer sphere center
    double d = (candidate->current.getPosition() - center).getR();

    // conservatively limit next step to prevent overshooting
    candidate->limitNextStep(fabs(d - radius));

    // no detection if outside of observer sphere
    if (d > radius)
        return NOTHING;

    // previous distance to observer sphere center
    double dprev = (candidate->previous.getPosition() - center).getR();

    // if particle was inside of sphere in previous step it has already been detected
    if (dprev <= radius)
        return NOTHING;

    // else detection
    return DETECTED;
}

std::string ObserverSmallSphere::getDescription() const {
    std::stringstream ss;
    ss << "ObserverSmallSphere: ";
    ss << "center = " << center / Mpc << " Mpc, ";
    ss << "radius = " << radius / Mpc << " Mpc";
    return ss.str();
}

ObserverLargeSphere::ObserverLargeSphere(Vector3d center, double radius) :
        center(center), radius(radius) {
}

DetectionState ObserverLargeSphere::checkDetection(Candidate *candidate) const {
    // current distance to observer sphere center
    double d = (candidate->current.getPosition() - center).getR();

    // conservatively limit next step size to prevent overshooting
    candidate->limitNextStep(fabs(radius - d));

    // no detection if inside observer sphere
    if (d < radius)
        return NOTHING;

    // previous distance to observer sphere center
    double dprev = (candidate->previous.getPosition() - center).getR();

    // if particle was outside of sphere in previous step it has already been detected
    if (dprev >= radius)
        return NOTHING;

    // else: detection
    return DETECTED;
}

std::string ObserverLargeSphere::getDescription() const {
    std::stringstream ss;
    ss << "ObserverLargeSphere: ";
    ss << "center = " << center / Mpc << " Mpc, ";
    ss << "radius = " << radius / Mpc << " Mpc";
    return ss.str();
}


DetectionState ObserverPoint::checkDetection(Candidate *candidate) const {
    double x = candidate->current.getPosition().x;
    if (x > 0) {
        candidate->limitNextStep(x);
        return NOTHING;
    }
    return DETECTED;
}

std::string ObserverPoint::getDescription() const {
    return "ObserverPoint: observer at x = 0";
}


ObserverRedshiftWindow::ObserverRedshiftWindow(double zmin, double zmax) :
        zmin(zmin), zmax(zmax) {
}

DetectionState ObserverRedshiftWindow::checkDetection(
        Candidate *candidate) const {
    double z = candidate->getRedshift();
    if (z > zmax)
        return VETO;
    if (z < zmin)
        return VETO;
    return NOTHING;
}

std::string ObserverRedshiftWindow::getDescription() const {
    std::stringstream ss;
    ss << "ObserverRedshiftWindow: z = " << zmin << " - " << zmax;
    return ss.str();
}

DetectionState ObserverNeutrinoVeto::checkDetection(Candidate *c) const {
    int id = fabs(c->current.getId());
    if ((id == 12) or (id == 14) or (id == 16))
        return NOTHING;
    return VETO;
}

std::string ObserverNeutrinoVeto::getDescription() const {
    return "ObserverNeutrinoVeto";
}

DetectionState ObserverChargedLeptonVeto::checkDetection(Candidate *c) const {
    int id = fabs(c->current.getId());
    if ((fabs(id) == 11) or (fabs(id) == 13) or (fabs(id == 15)))
        return NOTHING;
    return VETO;
}

std::string ObserverChargedLeptonVeto::getDescription() const {
    return "ObserverChargedLeptonVeto";
}

DetectionState ObserverPhotonVeto::checkDetection(Candidate *c) const {
    if (c->current.getId() == 22)
        return NOTHING;
    return VETO;
}

std::string ObserverPhotonVeto::getDescription() const {
    return "ObserverPhotonVeto";
}

ObserverOutput3D::ObserverOutput3D(std::string fname) {
    description = "ObserverOutput3D: " + fname;
    fout.open(fname.c_str());
    fout
            << "# dT\tD\tID\tID0\tE\tE0\tX\tY\tZ\tX0\tY0\tZ0\tPx\tPy\tPz\tP0x\tP0y\tP0z\tz\n";
    fout << "#\n";
    fout << "# D           Trajectory length [Mpc]\n";
    fout << "# ID          Particle type (PDG MC numbering scheme)\n";
    fout << "# E           Energy [EeV]\n";
    fout << "# X, Y, Z     Position [Mpc]\n";
    fout << "# Px, Py, Pz  Heading (unit vector of momentum)\n";
    fout << "# Initial state: ID0, E0, ...\n";
    fout << "# z           Redshift\n";
    fout << "#\n";
    
}

ObserverOutput3D::~ObserverOutput3D() {
    fout.close();
}

void ObserverOutput3D::onDetection(Candidate *candidate) const {
    char buffer[512];
    size_t p = 0;

    Vector3d xi = candidate->source.getPosition();
    Vector3d xf = candidate->current.getPosition();
    double dt = (xf - xi).getR() / c_light;
    p += sprintf(buffer + p, "%.9e\t", dt);
    p += sprintf(buffer + p, "%10i\t", candidate->current.getId());
    p += sprintf(buffer + p, "%10i\t", candidate->source.getId());
    p += sprintf(buffer + p, "%4.4e\t", candidate->current.getEnergy() / eV);
    p += sprintf(buffer + p, "%4.4e\t", candidate->source.getEnergy() / eV);
    Vector3d pos = candidate->current.getPosition() / Mpc;
    p += sprintf(buffer + p, "%7.6e\t%7.6e\t%7.6e\t", pos.x, pos.y, pos.z);
    Vector3d ipos = candidate->source.getPosition() / Mpc;
    p += sprintf(buffer + p, "%7.6e\t%7.6e\t%7.6e\t", ipos.x, ipos.y, ipos.z);
    Vector3d dir = candidate->current.getDirection();
    p += sprintf(buffer + p, "%.8e\t%.8e\t%.8e\t", dir.x, dir.y, dir.z);
    Vector3d idir = candidate->source.getDirection();
    p += sprintf(buffer + p, "%.8e\t%.8e\t%.8e\t", idir.x, idir.y, idir.z);
    p += sprintf(buffer + p, "%.9e\n", candidate->getRedshift());


#pragma omp critical
    {
        fout.write(buffer, p);
        fout.flush();
    }
}

ObserverOutput1D::ObserverOutput1D(std::string fname) {
    description = "ObserverOutput1D: " + fname;
    fout.open(fname.c_str());
    fout << "#ID\tE\tD\tID0\tE0\n";
    fout << "#\n";
    fout << "# ID  Particle type\n";
    fout << "# E   Energy [EeV]\n";
    fout << "# D   Comoving trajectory length [Mpc]\n";
    fout << "# ID0 Initial particle type\n";
    fout << "# E0  Initial energy [eV]\n";
}

ObserverOutput1D::~ObserverOutput1D() {
    fout.close();
}

void ObserverOutput1D::onDetection(Candidate *candidate) const {
    char buffer[256];
    size_t p = 0;

    p += sprintf(buffer + p, "%10i\t", candidate->current.getId());
    p += sprintf(buffer + p, "%.4e\t", candidate->current.getEnergy() / eV);
    p += sprintf(buffer + p, "%9.4f\t", candidate->getTrajectoryLength() / Mpc);
    p += sprintf(buffer + p, "%10i\t", candidate->source.getId());
    p += sprintf(buffer + p, "%.4e\n", candidate->source.getEnergy() / eV);

#pragma omp critical
    {
        fout.write(buffer, p);
        fout.flush();
    }
}

}// namespace