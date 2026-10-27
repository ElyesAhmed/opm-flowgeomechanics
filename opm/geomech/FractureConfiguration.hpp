#ifndef OPM_FRACTURE_CONFIGURATION_HPP
#define OPM_FRACTURE_CONFIGURATION_HPP

#include <opm/simulators/linalg/PropertyTree.hpp>
#include <string>

namespace Opm {

// Validate the whole parameter file before any solver uses it. Numeric and
// boolean strings remain supported for compatibility with the existing files.
// Does not insert defaults or change the supplied tree. Throws invalid_argument
// with the source and full key paths on failure.
void validateFractureConfiguration(const PropertyTree& settings,
                                   const std::string& source);

// Reports the input tree and the resolved principal mode choices. Solver-local
// defaults are not exhaustively expanded; this is stated in the report itself.
std::string fractureConfigurationReport(const PropertyTree& settings,
                                        const std::string& source);

// Apply and validate the two WSEED overrides on a COPY of the model settings.
// This avoids one seed's overrides leaking into the next seed/model defaults.
PropertyTree fractureSeedConfiguration(const PropertyTree& modelSettings,
                                       double verticalExtent, double flowWidth);

// Exact per-seed tree passed to Fracture::init, with provenance and SI units.
std::string fractureSeedConfigurationReport(const PropertyTree& settings,
                                            const std::string& well,
                                            int parentIndex, int completionIndex,
                                            int gridNumber);

} // namespace Opm
#endif
