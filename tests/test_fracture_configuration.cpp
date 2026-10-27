#include <config.h>
#include <opm/geomech/FractureConfiguration.hpp>
#include <opm/geomech/FractureModel.hpp>

#include <cmath>
#include <iostream>
#include <filesystem>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>

namespace Opm::Detail {
PropertyTree builtinFractureParam(const std::string& filename);
}

namespace {
void require(bool condition, const std::string& message)
{
    if (!condition) throw std::runtime_error(message);
}

template<class Function>
void rejects(Function function, const std::string& diagnostic)
{
    try {
        function();
    } catch (const std::invalid_argument& error) {
        require(std::string(error.what()).find(diagnostic) != std::string::npos,
                "Missing diagnostic '" + diagnostic + "' in: " + error.what());
        return;
    }
    throw std::runtime_error("Expected rejection: " + diagnostic);
}
}

int main(int argc, char** argv)
{
    try {
        require(argc == 2, "Expected source directory argument");
        const std::string source = argv[1];
        const auto defaults = Opm::makeDefaultFractureParam(false);
        Opm::validateFractureConfiguration(defaults, "C++ defaults");
        // This legacy factory variant produces perf_rate, which has no active
        // consumer. Reject it explicitly rather than silently enabling a mode.
        rejects([&] {
            Opm::validateFractureConfiguration(Opm::makeDefaultFractureParam(true), "legacy rate-control defaults");
        }, "perf_rate");
        for (const auto* alias : {"standard", "sequential_implicit"}) {
            Opm::validateFractureConfiguration(Opm::Detail::builtinFractureParam(alias), alias);
        }
        int rootProfiles = 0;
        for (const auto& entry : std::filesystem::directory_iterator(source + "/data")) {
            const auto name = entry.path().filename().string();
            if (entry.path().extension() == ".json" && name.rfind("fracture", 0) == 0) {
                Opm::validateFractureConfiguration(Opm::PropertyTree(entry.path().string()), name);
                ++rootProfiles;
            }
        }
        require(rootProfiles >= 13, "Expected the current root fracture profiles");
        int curatedProfiles = 0;
        for (const auto& entry : std::filesystem::directory_iterator(source + "/data/best_practice")) {
            if (entry.path().extension() != ".json") continue;
            const Opm::PropertyTree profile(entry.path().string());
            if (profile.get_child_optional("fractureparam")) {
                Opm::validateFractureConfiguration(profile, entry.path().string());
                ++curatedProfiles;
            }
        }
        require(curatedProfiles >= 14, "Expected the current curated fracture profiles");
        const Opm::PropertyTree baseline(source + "/data/best_practice/seq_implicit_fb_robust_v2.json");
        Opm::validateFractureConfiguration(baseline, "baseline profile");
        std::ostringstream before;
        baseline.write_json(before, false);
        const auto report = Opm::PropertyTree::fromJsonString(
            Opm::fractureConfigurationReport(baseline, "baseline profile"));
        require(report.get<std::string>("resolved.fracture_flow_mode") == "wi_upscaling", "Flow default");
        require(report.get<std::string>("source") == "baseline profile", "Source provenance");
        require(report.get<std::string>("resolved.vem_face_quadrature") == "kite", "Clean VEM quadrature default");
        require(report.get<std::string>("compatibility_fields.fractureparam.solver.full_intersections").find("no active reader") != std::string::npos,
                "Unused legacy option must be identified, not reported as active");
        for (const auto* quadrature : {"kite", "triangulated"}) {
            auto profile = baseline;
            profile.put("fractureparam.vem_face_quadrature", std::string(quadrature));
            Opm::validateFractureConfiguration(profile, "clean VEM profile");
            const auto modes = Opm::PropertyTree::fromJsonString(
                Opm::fractureConfigurationReport(profile, "clean VEM profile"));
            require(modes.get<std::string>("resolved.vem_face_quadrature") == quadrature,
                    "VEM quadrature must be reported faithfully");
        }
        std::ostringstream after;
        baseline.write_json(after, false);
        require(before.str() == after.str(), "Validation/reporting must not change input");

        const auto invalid = [&](const std::string& key, const auto& value,
                                 const std::string& message) {
            auto settings = baseline;
            if constexpr (std::is_convertible_v<decltype(value), const char*>) {
                settings.put(key, std::string(value));
            } else {
                settings.put(key, value);
            }
            rejects([&] { Opm::validateFractureConfiguration(settings, "invalid-test.json"); }, message);
        };
        invalid("fractureparam.vem_face_quadrature", "triangles", "unsupported value");
        invalid("fractureparam.solver.damping_p", 0.5, "unknown key");
        invalid("fractureparam.solver.damping_w", 0.5, "unknown key");
        invalid("box_from_flow", true, "unknown key"); // Separate mech-sector file.
        invalid("fractureparam.solver.tolerence", 1e-6, "fractureparam.solver.tolerence: unknown key");
        invalid("fractureparam.solver.well_local_solve", "fractured", "unknown key");
        invalid("fractureparam.control.type", "perf_rate", "unsupported value");
        invalid("solver.ctf_change_threshold", 1e-15, "solver.ctf_change_threshold: unknown key");
        invalid("solver.implicit_flow", "perhaps", "solver.implicit_flow: expected boolean");
        invalid("fractureparam.solver.max_iter", 1.5, "expected integer");
        invalid("fractureparam.solver.tolerance", "nan", "finite number");
        invalid("fractureparam.solver.tolerance", "1e9999", "finite number");
        invalid("fractureparam.solver.tolerance", -1.0, "must be positive");
        invalid("fractureparam.config.trires", 0, "must be positive");
        invalid("fractureparam.mech_grid_bc", "rollers", "unsupported value");
        invalid("fractureparam.solver.closed_cell_policy", "stiky", "unsupported value");
        invalid("solver.method", "FullyImplicitMech", "unsupported value");
        invalid("fractureparam.solver.coupling_state_relaxation", 1.5, "must be in (0, 1]");
        invalid("fractureparam.solver.max_closed_cell_toggle_fraction", -0.1, "must be in [0, 1]");
        invalid("fractureparam.solver.min_width", 1.0, "must not exceed");
        invalid("fractureparam.solver.coupling_relax_min", 2.0, "must not exceed");
        invalid("solver.method", "SeqMech", "not implemented");
        invalid("fractureparam.solver.pressure_from_flow", true, "requires fracture_flow_mode=embedded");
        invalid("solver.max_growth_iterations", 0, "requires SeqMechFrac");
        rejects([&] {
            Opm::validateFractureConfiguration(Opm::PropertyTree::fromJsonString(
                R"({"solver":{"method":"PostSolve","method":"SeqMechFrac"}})"), "duplicate.json");
        }, "solver.method: duplicate key");
        rejects([&] {
            Opm::validateFractureConfiguration(Opm::PropertyTree::fromJsonString(
                R"({"solver":{"method":["PostSolve"]}})"), "array.json");
        }, "solver.method: expected a scalar");
        rejects([&] {
            Opm::validateFractureConfiguration(Opm::PropertyTree::fromJsonString(
                R"({"solver.method":"PostSolve"})"), "dots.json");
        }, "without dots");
        rejects([&] {
            Opm::validateFractureConfiguration(Opm::PropertyTree::fromJsonString("{}"), "empty.json");
        }, "empty.json");

        rejects([&] {
            Opm::validateFractureConfiguration(Opm::PropertyTree::fromJsonString(
                R"({"solver":{"method":"PostSolve","implicit_flow":false},"fractureparam":{"verbosity":0}})"), "missing.json");
        }, "fractureparam.solver.method: required key missing");

        auto embedded = defaults;
        embedded.put("fractureparam.solver.fracture_flow_mode", std::string("embedded"));
        embedded.put("fractureparam.solver.pressure_from_flow", true);
        embedded.put("fractureparam.solver.embedded_well_source", std::string("config"));
        embedded.put("fractureparam.solver.leakoff_mobility", std::string("upwind"));
        // FlexibleSolver extensions are deliberately outside this registry.
        embedded.put("fractureparam.solver.linsolver.preconditioner.flow_solver.preconditioner.extension", 42);
        Opm::validateFractureConfiguration(embedded, "embedded profile");

        const auto model = baseline.get_child("fractureparam");
        const auto first = Opm::fractureSeedConfiguration(model, 2.0, 0.0003);
        const auto second = Opm::fractureSeedConfiguration(model, 5.0, 0.001);
        require(model.get<double>("config.min_width") == 0, "Original width was mutated");
        require(!model.get_child_optional("config.axis_scale"), "Original scale was mutated");
        require(first.get<double>("config.axis_scale") == 2.0, "First seed scale");
        require(first.get<double>("config.min_width") == 0.0003, "First seed flow floor");
        require(second.get<double>("config.min_width") == 0.001, "Second seed flow floor");
        require(first.get<double>("solver.min_width") == model.get<double>("solver.min_width"),
                "WSEED must not override the numerical solver width");
        const auto seedReport = Opm::PropertyTree::fromJsonString(
            Opm::fractureSeedConfigurationReport(first, "B-3H", 1875, 6562, 1));
        require(seedReport.get<int>("identity.wseed_parent_index") == 1875, "Parent identity");
        require(seedReport.get<int>("identity.completion_index") == 6562, "Completion identity");
        require(seedReport.get<int>("identity.grid_number") == 1, "Grid identity");
        require(seedReport.get<double>("settings.config.min_width") == 0.0003, "Effective width");
        require(seedReport.get<std::string>("units.length_and_aperture") == "m", "SI units");
        rejects([&] { Opm::fractureSeedConfiguration(model, 0, 0.001); }, "WSEED vertical extent");
        rejects([&] { Opm::fractureSeedConfiguration(model, 1, -0.001); }, "WSEED WIDTH");
        rejects([&] { Opm::fractureSeedConfiguration(model, std::numeric_limits<double>::infinity(), 0); }, "finite");
        rejects([&] { Opm::fractureSeedConfiguration(model, 1, std::numeric_limits<double>::quiet_NaN()); }, "finite");
        std::cout << "Fracture configuration validation and WSEED reports passed: "
                  << rootProfiles << " root profiles, " << curatedProfiles << " curated profiles\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
