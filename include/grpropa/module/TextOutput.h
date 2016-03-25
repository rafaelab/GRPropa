#ifndef GRPROPA_TEXTOUTPUT_H
#define GRPROPA_TEXTOUTPUT_H

#include "grpropa/module/Output.h"

#include <fstream>

namespace grpropa {

/**
 @class TextOutput
 @brief Configurable plain text output for cosmic ray information.
 */
class TextOutput: public Output {
protected:
	std::ostream *out;
	std::ofstream outfile;
	std::string filename;

	void printHeader() const;

public:
	TextOutput();
	TextOutput(OutputType outputtype);
	TextOutput(std::ostream &out);
	TextOutput(std::ostream &out, OutputType outputtype);
	TextOutput(const std::string &filename);
	TextOutput(const std::string &filename, OutputType outputtype);
	~TextOutput();

	void close();
	void gzip();

	void process(Candidate *candidate) const;
	std::string getDescription() const;
};

} // namespace grpropa

#endif // GRPROPA_TEXTOUTPUT_H