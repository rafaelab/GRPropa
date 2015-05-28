#ifndef GRPROPA_OUTPUTSHELL_H
#define GRPROPA_OUTPUTSHELL_H

#include "grpropa/Module.h"
#include "grpropa/AssocVector.h"

namespace grpropa {

/**
 @class ShellOutput
 @brief Show the trajectory in the shell.
 */
class ShellOutput: public Module {
public:
    void process(Candidate *candidate) const;
    std::string getDescription() const;
};

/**
 @class ShellOutput1D
 @brief Show the trajectory in the shell.
 */
class ShellOutput1D: public Module {
public:
    void process(Candidate *candidate) const;
    std::string getDescription() const;
};

/**
 @class ShellPropertyOutput
 @brief Show the candidate properties in the shell.
 */
class ShellPropertyOutput: public Module {
public:
    typedef Loki::AssocVector<std::string, std::string> PropertyMap;
    void process(Candidate *candidate) const;
    std::string getDescription() const;
};

} // namespace cprpropa

#endif // GRPROPA_OUTPUTSHELL_H
