#ifndef GRPROPA_MODULETOOLS_H
#define GRPROPA_MODULETOOLS_H

#include "grpropa/Module.h"

#include <vector>

namespace grpropa {

class PerformanceModule: public Module {
private:
    struct _module_info {
        double time;
        ref_ptr<Module> module;
    };

    mutable std::vector<_module_info> modules;
    mutable size_t calls;

public:
    ~PerformanceModule();
    void add(Module* module);
    void process(Candidate* candidate) const;
    std::string getDescription() const;
};

} // namespace grpropa

#endif // GRPROPA_MODULETOOLS_H
