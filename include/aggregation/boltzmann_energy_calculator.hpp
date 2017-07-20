#ifndef BOLTZMANN_ENERGY_CALCULATOR_HPP
#define BOLTZMANN_ENERGY_CALCULATOR_HPP

#include <vector>

#include "gromacs/utility/real.h"


/*!
 * Enum for energy units. Boltzmann refers to energy in units of 
 * \f$ k_{\text{B}} T \f$ where \f$ k_{\text{B}} \f$ is the Boltzmann constant 
 * and \f$ T \f$ is thermodynamic temperature in Kelvin. 
 */
enum eEnergyUnit {eEnergyUnitBoltzmann, 
                  eEnergyUnitKiloJoulePerMol,
                  eEnergyUnitKiloCaloriePerMol};


/*!
 * \brief Calculation of energy from Boltzmann statistics.
 *
 * This class is used to calculate an energy from the equilibrium density of
 * solvent particles, i.e. it implements the inversion of
 *
 * \f[
 *      n(s) = \frac{1}{Z} \exp\big( -\beta E(s) \big)
 * \f]
 *
 * where \f$ \beta = 1/k_{\text{B}} T \f$.
 *
 * Note that the left hand side in the above equation will usually be a number
 * density determined by NumberDensityCalculator, but it may also be a 
 * probability density (i.e. normalised to one instead of number of particles)
 * as any constant normalising factor will only enter the energy as an additive
 * constant. For the same reason the partition function is ignored.
 *
 * It is also important to note that the nature of the energy in the above 
 * equation depends on which ensemble was used in the simulation from which the
 * density is calculated. If the NVT ensemble was used, it is Helmholtz free
 * energy \f$ F \f$, if the NPT ensemble was used, it is Gibbs free energy
 * \f$ G \f$.
 *
 * \todo Implement automatic shift of energy so that bulk value is at zero.
 */
class BoltzmannEnergyCalculator
{
    public:

        // constructor:
        BoltzmannEnergyCalculator();

        // public interface for energy calculation:
        std::vector<real> calculate(
                const std::vector<real> &density);


        // set units:
        void setTemperature(const real temperature);
        void setEnergyUnits(const eEnergyUnit &unit);
        

    private:

        // internal variables:
        real energyUnitFactor_;
        real temperature_;
        bool temperatureIsSet_;

        // internal constants:
        const real cGasConstantKiloJoulePerMol_ = 0.0083144598;
        const real cGasConstantKiloCaloriePerMol_ = 0.0019872036;

        // auxiliary functions:
        void mendInfinities(std::vector<real> &energy);
};

#endif

