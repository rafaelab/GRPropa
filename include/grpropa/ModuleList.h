#ifndef GRPROPA_MODULE_LIST_H
#define GRPROPA_MODULE_LIST_H

#include "grpropa/Candidate.h"
#include "grpropa/Module.h"
#include "grpropa/Source.h"

#include <list>
#include <sstream>

namespace grpropa {

/**
 @class ModuleList
 @brief List of modules
 */
class ModuleList: public Referenced {
public:
    typedef std::list<ref_ptr<Module> > module_list_t;
    typedef std::vector<ref_ptr<Candidate> > candidate_vector_t;

    ModuleList();
    virtual ~ModuleList();
    void setShowProgress(bool show);

    void add(Module* module);
    virtual void process(Candidate *candidate);
    void run(Candidate *candidate, bool recursive = true);
    void run(candidate_vector_t &candidates, bool recursive = true);
    void run(Source *source, size_t count, bool recursive = true);

    module_list_t &getModules();
    const module_list_t &getModules() const;

    std::string getDescription() const;
    void showModules() const;

private:
    module_list_t modules;
    bool showProgress;
};

} // namespace grpropa

#endif // GRPROPA_MODULE_LIST_H
