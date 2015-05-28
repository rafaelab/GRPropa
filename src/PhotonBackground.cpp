#include "grpropa/PhotonBackground.h"
#include "grpropa/Common.h"

#include <vector>
#include <iostream>

namespace grpropa {

void listOfPhotonBackgroundModels()
{
    std::cout << "List of available phton background models." << std::endl;
    std::cout << " - CMB" << std::endl;
    std::cout << " - CRB:" << std::endl;
    std::cout << "     ARCADE2" << std::endl;
    std::cout << "     Protheroe96" << std::endl;
    std::cout << " - EBL:" << std::endl;
    std::cout << "      Gilmore+ '12" << std::endl;
    std::cout << "      Dominguez+ '11 (std, upper, lower)" << std::endl;
    std::cout << "      Finke+ '10" << std::endl;
    std::cout << "      Kneiske & Dole '10" << std::endl;
    std::cout << "      Franceschini+ '08" << std::endl;
}

} // namespace grpropa
