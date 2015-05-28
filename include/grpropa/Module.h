#ifndef GRPROPA_MODULE_H
#define GRPROPA_MODULE_H

#include "grpropa/Candidate.h"
#include "grpropa/Referenced.h"
#include "grpropa/Common.h"

#include <string>

namespace grpropa {

class Candidate;

/**
 @class Module
 @brief Abstract base class for modules
 */
class Module: public Referenced {
    std::string description;
public:
    Module();
    virtual ~Module() {
    }
    virtual std::string getDescription() const;
    void setDescription(const std::string &description);
    virtual void process(Candidate *candidate) const = 0;
    inline void process(ref_ptr<Candidate> candidate) const {
        process(candidate.get());
    }
};

} // namespace grpropa

#endif /* GRPROPA_MODULE_H */
