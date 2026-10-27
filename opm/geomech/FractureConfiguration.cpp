#include <config.h>
#include <opm/geomech/FractureConfiguration.hpp>

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <cmath>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace Opm {
namespace {
using Tree = boost::property_tree::ptree;
enum class Kind { Boolean, Integer, Number, String };

// Explicit registry of keys read by the fracture/outer-coupling implementation,
// plus legacy fields present in the C++ defaults and built-in profiles. Update
// this registry alongside new parameter reads. FlexibleSolver's flow_solver
// subtree is delegated to that external component rather than frozen here.
const std::map<std::string, Kind>& rules()
{
    static const std::map<std::string, Kind> registry {
        {"add_perfs_to_schedule", Kind::Boolean},
        {"fractureparam.KMax", Kind::Number},
        {"fractureparam.addconnections", Kind::Boolean},
        {"fractureparam.addperfs_to_schedule", Kind::Boolean},
        {"fractureparam.config.axis_scale", Kind::Number},
        {"fractureparam.config.gravity_off", Kind::Boolean},
        {"fractureparam.config.initial_fracture_width", Kind::Number},
        {"fractureparam.config.min_width", Kind::Number},
        {"fractureparam.config.num_exp", Kind::Integer},
        {"fractureparam.config.num_lin", Kind::Integer},
        {"fractureparam.config.scale_filtrate", Kind::Boolean},
        {"fractureparam.config.trires", Kind::Number},
        {"fractureparam.config.type", Kind::String},
        {"fractureparam.config.verbosity", Kind::Integer},
        {"fractureparam.control.WI", Kind::Number},
        {"fractureparam.control.pressure", Kind::Number},
        {"fractureparam.control.rate", Kind::Number},
        {"fractureparam.control.type", Kind::String},
        {"fractureparam.extended_fractures", Kind::Boolean},
        {"fractureparam.fractureWI", Kind::Number},
        {"fractureparam.fracturesolver", Kind::String},
        {"fractureparam.include_fracture_contributions", Kind::Boolean},
        {"fractureparam.max_iter", Kind::Integer},
        {"fractureparam.mech_diagonal_scaling", Kind::Boolean},
        {"fractureparam.mech_grid_bc", Kind::String},
        {"fractureparam.method.iterate", Kind::Boolean},
        {"fractureparam.method.max_it", Kind::Integer},
        {"fractureparam.method.tolerance", Kind::Number},
        {"fractureparam.patch_recovery_stress", Kind::Boolean},
        {"fractureparam.pressuresolver", Kind::String},
        {"fractureparam.reduce_boundary", Kind::Boolean},
        {"fractureparam.require_converged_fracture_for_wi_update", Kind::Boolean},
        {"fractureparam.reservoir.calculate_dist", Kind::Boolean},
        {"fractureparam.reservoir.dist", Kind::Number},
        {"fractureparam.reservoir.interpolate_stress", Kind::Boolean},
        {"fractureparam.reservoir.mobility", Kind::Number},
        {"fractureparam.reservoir.perm", Kind::Number},
        {"fractureparam.smooth_force", Kind::Boolean},
        {"fractureparam.smooth_force_length", Kind::Number},
        {"fractureparam.solver.area_change_fac", Kind::Number},
        {"fractureparam.solver.cellcount_threshold", Kind::Integer},
        {"fractureparam.solver.chatter_flip_limit", Kind::Integer},
        {"fractureparam.solver.check_coupling_fd", Kind::Boolean},
        {"fractureparam.solver.check_coupling_fd_columns", Kind::Integer},
        {"fractureparam.solver.check_coupling_fd_perturbation", Kind::Number},
        {"fractureparam.solver.check_coupling_fd_tolerance", Kind::Number},
        {"fractureparam.solver.check_coupling_fd_verbosity", Kind::Integer},
        {"fractureparam.solver.close_force_tolerance", Kind::Number},
        {"fractureparam.solver.closed_cell_policy", Kind::String},
        {"fractureparam.solver.closing_type", Kind::String},
        {"fractureparam.solver.conservative_propagation", Kind::Boolean},
        {"fractureparam.solver.coupled_step_limit", Kind::Number},
        {"fractureparam.solver.coupling_acceptance_factor", Kind::Number},
        {"fractureparam.solver.coupling_anderson_alpha_max", Kind::Number},
        {"fractureparam.solver.coupling_anderson_alpha_min", Kind::Number},
        {"fractureparam.solver.coupling_convergence_mode", Kind::String},
        {"fractureparam.solver.coupling_dt_area_target", Kind::Number},
        {"fractureparam.solver.coupling_dt_min", Kind::Number},
        {"fractureparam.solver.coupling_dt_pressure_target", Kind::Number},
        {"fractureparam.solver.coupling_dt_quiet_steps", Kind::Integer},
        {"fractureparam.solver.coupling_fallback_legacy", Kind::Boolean},
        {"fractureparam.solver.coupling_legacy_first_damp", Kind::Boolean},
        {"fractureparam.solver.coupling_relax_init", Kind::Number},
        {"fractureparam.solver.coupling_relax_max", Kind::Number},
        {"fractureparam.solver.coupling_relax_min", Kind::Number},
        {"fractureparam.solver.coupling_state_relaxation", Kind::Number},
        {"fractureparam.solver.coupling_tolerance", Kind::Number},
        {"fractureparam.solver.coupling_update_mode", Kind::String},
        {"fractureparam.solver.ctf_change_threshold", Kind::Number},
        {"fractureparam.solver.ctf_change_threshold_rel", Kind::Number},
        {"fractureparam.solver.ctf_relative_scale", Kind::Number},
        {"fractureparam.solver.damping", Kind::Number},
        {"fractureparam.solver.damping_factor_perf", Kind::Number},
        {"fractureparam.solver.damping_factor_wi", Kind::Number},
        {"fractureparam.solver.disable_propagation", Kind::Boolean},
        {"fractureparam.solver.divide_wellidx", Kind::Boolean},
        {"fractureparam.solver.drop_fluid_mech_linearization", Kind::Boolean},
        {"fractureparam.solver.drop_tol_h", Kind::Number},
        {"fractureparam.solver.drop_tol_p", Kind::Number},
        {"fractureparam.solver.dt_limit", Kind::Number},
        {"fractureparam.solver.dump_case_min_nonlin_iter", Kind::Integer},
        {"fractureparam.solver.dump_linear_system_on_failure", Kind::Boolean},
        {"fractureparam.solver.efac", Kind::Number},
        {"fractureparam.solver.embedded_capacity", Kind::Integer},
        {"fractureparam.solver.embedded_cell_dump", Kind::Integer},
        {"fractureparam.solver.embedded_leakoff_report", Kind::Boolean},
        {"fractureparam.solver.embedded_max_dp", Kind::Number},
        {"fractureparam.solver.embedded_max_ds", Kind::Number},
        {"fractureparam.solver.embedded_max_dt", Kind::Number},
        {"fractureparam.solver.embedded_max_temp_change", Kind::Number},
        {"fractureparam.solver.embedded_min_aperture", Kind::Number},
        {"fractureparam.solver.embedded_perf_rw", Kind::Number},
        {"fractureparam.solver.embedded_perf_wi_mode", Kind::String},
        {"fractureparam.solver.embedded_perf_width", Kind::Number},
        {"fractureparam.solver.embedded_remap_state", Kind::Boolean},
        {"fractureparam.solver.embedded_residual_check", Kind::Boolean},
        {"fractureparam.solver.embedded_satnum", Kind::Integer},
        {"fractureparam.solver.embedded_static", Kind::Boolean},
        {"fractureparam.solver.embedded_thermal_distance", Kind::Number},
        {"fractureparam.solver.embedded_wall_conduction", Kind::Boolean},
        {"fractureparam.solver.embedded_well_source", Kind::String},
        {"fractureparam.solver.enable_wi_coupling_update", Kind::Boolean},
        {"fractureparam.solver.failure_on_nonconvergence", Kind::Boolean},
        {"fractureparam.solver.fb_clamp_width", Kind::Boolean},
        {"fractureparam.solver.fb_line_search", Kind::Boolean},
        {"fractureparam.solver.fb_line_search_min", Kind::Number},
        {"fractureparam.solver.fb_mu", Kind::Number},
        {"fractureparam.solver.fb_scaled_convergence", Kind::Boolean},
        {"fractureparam.solver.force_limit", Kind::Number},
        {"fractureparam.solver.fractureWI_mode", Kind::String},
        {"fractureparam.solver.fractureWI_rw", Kind::Number},
        {"fractureparam.solver.fracture_flow_mode", Kind::String},
        {"fractureparam.solver.fracture_storage_term", Kind::Boolean},
        {"fractureparam.solver.freeze_closed_after", Kind::Integer},
        {"fractureparam.solver.full_intersections", Kind::Boolean},
        {"fractureparam.solver.guard_max_abs_pressure", Kind::Number},
        {"fractureparam.solver.guard_max_abs_width", Kind::Number},
        {"fractureparam.solver.guard_state", Kind::Boolean},
        {"fractureparam.solver.iteration_snapshot_cells", Kind::Boolean},
        {"fractureparam.solver.iteration_snapshot_dir", Kind::String},
        {"fractureparam.solver.k1_correlation_cells", Kind::Number},
        {"fractureparam.solver.k1_correlation_length", Kind::Number},
        {"fractureparam.solver.k1_extraction", Kind::String},
        {"fractureparam.solver.leakoff_distance_factor", Kind::Number},
        {"fractureparam.solver.leakoff_mobility", Kind::String},
        {"fractureparam.solver.leakoff_model", Kind::String},
        {"fractureparam.solver.leakoff_sides", Kind::Integer},
        {"fractureparam.solver.legacy_coupling_change_logic", Kind::Boolean},
        {"fractureparam.solver.legacy_parent_setup_iteration", Kind::Boolean},
        {"fractureparam.solver.linear_system_dump_prefix", Kind::String},
        {"fractureparam.solver.linsolver.atol", Kind::Number},
        {"fractureparam.solver.linsolver.debug_checksum", Kind::Boolean},
        {"fractureparam.solver.linsolver.direct_fallback_solver", Kind::String},
        {"fractureparam.solver.linsolver.direct_max_cells", Kind::Integer},
        {"fractureparam.solver.linsolver.failure_policy", Kind::String},
        {"fractureparam.solver.linsolver.fallback_gmres_on_breakdown", Kind::Boolean},
        {"fractureparam.solver.linsolver.max_iter", Kind::Number},
        {"fractureparam.solver.linsolver.preconditioner.diag_flow", Kind::Boolean},
        {"fractureparam.solver.linsolver.preconditioner.diag_mech", Kind::Boolean},
        {"fractureparam.solver.linsolver.preconditioner.fixed_stress", Kind::Boolean},
        {"fractureparam.solver.linsolver.preconditioner.flow_solver.maxiter", Kind::Integer},
        {"fractureparam.solver.linsolver.preconditioner.flow_solver.preconditioner.type", Kind::String},
        {"fractureparam.solver.linsolver.preconditioner.flow_solver.preconditioner.verbosity", Kind::Integer},
        {"fractureparam.solver.linsolver.preconditioner.flow_solver.solver", Kind::String},
        {"fractureparam.solver.linsolver.preconditioner.flow_solver.tol", Kind::Number},
        {"fractureparam.solver.linsolver.preconditioner.flow_solver.verbosity", Kind::Integer},
        {"fractureparam.solver.linsolver.preconditioner.lu_pivoting", Kind::Boolean},
        {"fractureparam.solver.linsolver.preconditioner.mech_first", Kind::Boolean},
        {"fractureparam.solver.linsolver.preconditioner.mech_press_coupling", Kind::Boolean},
        {"fractureparam.solver.linsolver.preconditioner.mode_policy", Kind::String},
        {"fractureparam.solver.linsolver.preconditioner.mode_switch_coupling_threshold", Kind::Number},
        {"fractureparam.solver.linsolver.preconditioner.update_mech_on_reuse", Kind::Boolean},
        {"fractureparam.solver.linsolver.preconditioner.verbosity", Kind::Integer},
        {"fractureparam.solver.linsolver.restart", Kind::Integer},
        {"fractureparam.solver.linsolver.scaling", Kind::String},
        {"fractureparam.solver.linsolver.solver", Kind::String},
        {"fractureparam.solver.linsolver.tol", Kind::Number},
        {"fractureparam.solver.linsolver.verbosity", Kind::Number},
        {"fractureparam.solver.max_area_growth_per_step", Kind::Number},
        {"fractureparam.solver.max_change", Kind::Number},
        {"fractureparam.solver.max_closed_cell_toggle_count", Kind::Integer},
        {"fractureparam.solver.max_closed_cell_toggle_fraction", Kind::Number},
        {"fractureparam.solver.max_dp", Kind::Number},
        {"fractureparam.solver.max_dwidth", Kind::Number},
        {"fractureparam.solver.max_expand_iter", Kind::Integer},
        {"fractureparam.solver.max_iter", Kind::Integer},
        {"fractureparam.solver.max_iter_on_same_level", Kind::Integer},
        {"fractureparam.solver.max_num_coarsening", Kind::Integer},
        {"fractureparam.solver.max_width", Kind::Number},
        {"fractureparam.solver.method", Kind::String},
        {"fractureparam.solver.min_area_growth_guard", Kind::Number},
        {"fractureparam.solver.min_width", Kind::Number},
        {"fractureparam.solver.no_leakof_outercells", Kind::Boolean},
        {"fractureparam.solver.numcell_threshold", Kind::Integer},
        {"fractureparam.solver.onset_dt_limit", Kind::Number},
        {"fractureparam.solver.onset_pressure_rise", Kind::Number},
        {"fractureparam.solver.onset_seed_area", Kind::Number},
        {"fractureparam.solver.perf_pressure_change_threshold", Kind::Number},
        {"fractureparam.solver.perf_pressure_change_threshold_rel", Kind::Number},
        {"fractureparam.solver.perf_pressure_relative_scale", Kind::Number},
        {"fractureparam.solver.pressure_from_flow", Kind::Boolean},
        {"fractureparam.solver.pressure_from_flow_iterations", Kind::Integer},
        {"fractureparam.solver.pressure_from_flow_min_width_factor", Kind::Number},
        {"fractureparam.solver.pressure_from_flow_tolerance", Kind::Number},
        {"fractureparam.solver.pressure_from_flow_verbosity", Kind::Integer},
        {"fractureparam.solver.propagate_converged_cells_only", Kind::Boolean},
        {"fractureparam.solver.propagate_fb_tol", Kind::Number},
        {"fractureparam.solver.propagate_open_front_only", Kind::Boolean},
        {"fractureparam.solver.propagate_stable_iters", Kind::Integer},
        {"fractureparam.solver.rebind_after_growth", Kind::Boolean},
        {"fractureparam.solver.rebind_after_growth_rounds", Kind::Integer},
        {"fractureparam.solver.remap_solution", Kind::Boolean},
        {"fractureparam.solver.reopen_force_tolerance", Kind::Number},
        {"fractureparam.solver.reopen_use_open_pressure", Kind::Boolean},
        {"fractureparam.solver.reopen_width_tolerance", Kind::Number},
        {"fractureparam.solver.retract_established_compressed", Kind::Boolean},
        {"fractureparam.solver.rfac", Kind::Number},
        {"fractureparam.solver.smooth_boundary", Kind::Boolean},
        {"fractureparam.solver.solve_failure_policy", Kind::String},
        {"fractureparam.solver.stress_dt_min", Kind::Number},
        {"fractureparam.solver.stress_dt_target", Kind::Number},
        {"fractureparam.solver.target_cellcount", Kind::Integer},
        {"fractureparam.solver.toggle_guard_mode", Kind::String},
        {"fractureparam.solver.tolerance", Kind::Number},
        {"fractureparam.solver.use_ad_pressure_assembly", Kind::Boolean},
        {"fractureparam.solver.use_fluid_system_water_properties", Kind::Boolean},
        {"fractureparam.solver.use_water_phase_mobility_for_fracture_flow", Kind::Boolean},
        {"fractureparam.solver.use_water_phase_mobility_for_well_fracture_flow", Kind::Boolean},
        {"fractureparam.solver.value_only_wi_update", Kind::Boolean},
        {"fractureparam.solver.verbosity", Kind::Integer},
        {"fractureparam.solver.volume_paced_propagation", Kind::Boolean},
        {"fractureparam.solver.volume_pacing_factor", Kind::Number},
        {"fractureparam.solver.volume_pacing_pressure", Kind::Boolean},
        {"fractureparam.solver.volume_pacing_source", Kind::String},
        {"fractureparam.solver.volume_pacing_verbosity", Kind::Integer},
        {"fractureparam.solver.warn_flux_tol", Kind::Number},
        {"fractureparam.solver.warn_open_area", Kind::Number},
        {"fractureparam.solver.warn_open_not_conducting", Kind::Boolean},
        {"fractureparam.solver.well_consistency_warn_bhp", Kind::Number},
        {"fractureparam.solver.well_source_all_perfs", Kind::Boolean},
        {"fractureparam.solver.well_source_radius", Kind::Number},
        {"fractureparam.solver.wi_flux_normalization", Kind::Boolean},
        {"fractureparam.solver.wi_normalization_max", Kind::Number},
        {"fractureparam.solver.wi_pressure_floor", Kind::Number},
        {"fractureparam.solver.wi_sign_gate", Kind::Boolean},
        {"fractureparam.solver.wi_upscaling", Kind::String},
        {"fractureparam.solver.wi_vector_acceleration", Kind::String},
        {"fractureparam.solver.write_coupled_linear_system", Kind::Boolean},
        {"fractureparam.solver.write_fracture_system", Kind::Boolean},
        {"fractureparam.solver.write_iteration_snapshots", Kind::Boolean},
        {"fractureparam.stab_on_stress", Kind::Boolean},
        {"fractureparam.topology_hysteresis.confirm_iterations", Kind::Integer},
        {"fractureparam.topology_hysteresis.cooldown_iterations", Kind::Integer},
        {"fractureparam.topology_hysteresis.emergency_threshold", Kind::Number},
        {"fractureparam.topology_hysteresis.enable", Kind::Boolean},
        {"fractureparam.topology_hysteresis.off_threshold", Kind::Number},
        {"fractureparam.topology_hysteresis.on_threshold", Kind::Number},
        {"fractureparam.vem_face_quadrature", Kind::String},
        {"fractureparam.vem_force", Kind::Boolean},
        {"fractureparam.vem_source", Kind::Boolean},
        {"fractureparam.vem_stability_choice", Kind::Integer},
        {"fractureparam.vem_stress", Kind::Boolean},
        {"fractureparam.verbose", Kind::Boolean},
        {"fractureparam.verbosity", Kind::Integer},
        {"fractureparam.vtkwritewells", Kind::Boolean},
        {"fractureparam.write_fracture_system", Kind::Boolean},
        {"fractureparam.write_pressure_system", Kind::Boolean},
        {"hasfractures", Kind::Boolean},
        {"solver.area_growth_guard_min_dt", Kind::Number},
        {"solver.growth_driven_switch", Kind::Boolean},
        {"solver.growth_flow_iterations", Kind::Integer},
        {"solver.implicit_flow", Kind::Boolean},
        {"solver.max_frac_it", Kind::Integer},
        {"solver.max_growth_iterations", Kind::Integer},
        {"solver.max_mech_it", Kind::Integer},
        {"solver.max_mech_solves_per_step", Kind::Integer},
        {"solver.mech_skip_coupling_threshold", Kind::Number},
        {"solver.method", Kind::String},
        {"solver.unconverged_fracture_dt_cut", Kind::Boolean},
        {"solver.verbosity", Kind::Integer},
        {"solver.well_local_inner_well_iterations", Kind::Boolean},
        {"solver.well_local_max_iter", Kind::Integer},
        {"solver.well_local_reduction", Kind::Number},
        {"solver.well_local_ring", Kind::Integer},
        {"solver.well_local_solve", Kind::String},
        {"solver.well_local_verbosity", Kind::Integer},
        {"verbosity", Kind::Integer},
    };
    return registry;
}

Tree asTree(const PropertyTree& settings)
{
    std::stringstream stream;
    settings.write_json(stream, false);
    Tree tree;
    boost::property_tree::read_json(stream, tree);
    return tree;
}

std::string json(const Tree& tree)
{
    std::ostringstream stream;
    boost::property_tree::write_json(stream, tree, true);
    return stream.str();
}

const std::string flowSolver = "fractureparam.solver.linsolver.preconditioner.flow_solver";

void inspect(const Tree& node, const std::string& prefix,
             std::vector<std::string>& errors)
{
    std::set<std::string> seen;
    for (const auto& [key, child] : node) {
        const std::string path = prefix.empty() ? key : prefix + "." + key;
        if (!seen.insert(key).second) {
            errors.push_back(path + ": duplicate key");
            continue;
        }
        if (key.empty() || key.find('.') != std::string::npos) {
            errors.push_back(path + ": expected an object key without dots (arrays are unsupported)");
            continue;
        }
        if (key.rfind("__", 0) == 0) {
            if (!child.empty()) errors.push_back(path + ": documentation must be scalar");
            continue;
        }
        const auto rule = rules().find(path);
        if (rule == rules().end()) {
            const auto next = rules().lower_bound(path + ".");
            const bool object = next != rules().end() && next->first.rfind(path + ".", 0) == 0;
            if (!object) {
                errors.push_back(path + ": unknown key (outer solver.* and fractureparam.solver.* are distinct)");
            } else if (child.empty() || !child.data().empty()) {
                errors.push_back(path + ": expected a nonempty object");
            } else if (path != flowSolver) {
                inspect(child, path, errors);
            }
            // flow_solver's keys and values belong to Dune::FlexibleSolver.
            continue;
        }
        if (!child.empty()) {
            errors.push_back(path + ": expected a scalar value");
            continue;
        }
        try {
            switch (rule->second) {
            case Kind::Boolean:
                // Match PropertyTree's bool conversion, including legacy 0/1.
                child.get_value<bool>();
                break;
            case Kind::Integer:
                child.get_value<int>();
                break;
            case Kind::Number:
                if (!std::isfinite(child.get_value<double>())) {
                    throw std::invalid_argument("nonfinite number");
                }
                break;
            case Kind::String:
                if (child.data().empty() || child.data() == "null") {
                    throw std::invalid_argument("empty string");
                }
                break;
            }
        } catch (const std::exception&) {
            const char* expected = rule->second == Kind::Boolean ? "boolean" :
                rule->second == Kind::Integer ? "integer" :
                rule->second == Kind::Number ? "finite number" : "nonempty string";
            errors.push_back(path + ": expected " + expected + ", got '" + child.data() + "'");
        }
    }
}

void choice(const Tree& tree, const std::string& key,
            const std::set<std::string>& values, std::vector<std::string>& errors)
{
    const auto value = tree.get_optional<std::string>(key);
    if (!value || values.count(*value)) return;
    std::string message = key + ": unsupported value '" + *value + "'; expected one of";
    for (const auto& option : values) message += " " + option;
    errors.push_back(message);
}

void positive(const Tree& tree, const std::string& key, bool allowZero,
              std::vector<std::string>& errors)
{
    const auto value = tree.get_optional<double>(key);
    if (value && (allowZero ? *value < 0 : *value <= 0)) {
        errors.push_back(key + (allowZero ? ": must be nonnegative" : ": must be positive"));
    }
}

void ordered(const Tree& tree, const std::string& lower, const std::string& upper,
             std::vector<std::string>& errors)
{
    const auto low = tree.get_optional<double>(lower);
    const auto high = tree.get_optional<double>(upper);
    if (low && high && *low > *high) errors.push_back(lower + ": must not exceed " + upper);
}

void units(Tree& report)
{
    report.put("units.length_and_aperture", "m");
    report.put("units.pressure_and_traction", "Pa unless explicitly listed below");
    report.put("units.rate", "m^3/s (configuration control.rate)");
    report.put("units.permeability", "m^2");
    report.put("units.toughness_KMax", "Pa*sqrt(m)");
    report.put("units.fractureparam_solver_dt_limit", "day");
    report.put("units.fractureparam_solver_onset_dt_limit", "day");
    report.put("units.fractureparam_solver_onset_pressure_rise", "bar");
    report.put("units.fractureparam_solver_pressure_from_flow_tolerance", "bar");
    report.put("units.fractureparam_solver_embedded_max_dp", "bar");
    report.put("units.fractureparam_solver_embedded_max_dt", "K (temperature update)");
}
} // namespace

void validateFractureConfiguration(const PropertyTree& settings, const std::string& source)
{
    const Tree tree = asTree(settings);
    std::vector<std::string> errors;
    inspect(tree, "", errors);
    const auto throwErrors = [&]() {
        if (errors.empty()) return;
        std::string message = "Invalid fracture configuration '" + source + "':";
        for (const auto& error : errors) message += "\n  " + error;
        throw std::invalid_argument(message);
    };
    throwErrors();
    // Required driver settings; branch-specific inner-solver requirements remain
    // with their consumers. Do not merge a different profile into user input.
    for (const auto* key : {"solver.method", "solver.implicit_flow",
         "fractureparam.solver.method", "fractureparam.control.type",
         "fractureparam.addconnections", "fractureparam.reduce_boundary",
         "fractureparam.include_fracture_contributions"}) {
        if (!tree.get_child_optional(key)) errors.push_back(std::string(key) + ": required key missing");
    }
    if (!tree.get_child_optional("fractureparam")) errors.push_back("fractureparam: required object missing");

    choice(tree, "solver.method", {"PostSolve", "SeqMech", "SeqMechFrac"}, errors);
    choice(tree, "solver.well_local_solve", {"none", "fractured", "all"}, errors);
    choice(tree, "fractureparam.vem_face_quadrature", {"kite", "triangulated"}, errors);
    choice(tree, "fractureparam.mech_grid_bc", {"fixed", "roller", "roller_free_top"}, errors);
    choice(tree, "fractureparam.config.type", {"well_seed", "perp_well", "tensile_fracture"}, errors);
    choice(tree, "fractureparam.solver.method", {"nothing", "simple", "only_pressure", "only_width", "iterative", "if", "if_propagate", "if_propagate_trimesh"}, errors);
    choice(tree, "fractureparam.control.type", {"pressure", "rate", "perf_pressure", "rate_well", "bhp_well", "well"}, errors);
    choice(tree, "fractureparam.solver.fracture_flow_mode", {"wi_upscaling", "embedded"}, errors);
    choice(tree, "fractureparam.solver.closed_cell_policy", {"legacy", "sticky", "pdas", "fischer_burmeister", "none"}, errors);
    choice(tree, "fractureparam.solver.closing_type", {"org", "open", "simple"}, errors);
    choice(tree, "fractureparam.solver.coupling_update_mode", {"legacy", "current", "aitken", "anderson"}, errors);
    choice(tree, "fractureparam.solver.k1_extraction", {"aperture_cell", "displacement_correlation"}, errors);
    choice(tree, "fractureparam.solver.toggle_guard_mode", {"count", "chatter"}, errors);
    choice(tree, "fractureparam.solver.wi_vector_acceleration", {"none", "aitken", "anderson"}, errors);
    choice(tree, "fractureparam.solver.coupling_convergence_mode", {"threshold", "residual"}, errors);
    choice(tree, "fractureparam.solver.wi_upscaling", {"legacy", "conductivity"}, errors);
    choice(tree, "fractureparam.solver.fractureWI_mode", {"constant", "estimate"}, errors);
    choice(tree, "fractureparam.solver.leakoff_model", {"legacy", "onesided", "reveal"}, errors);
    choice(tree, "fractureparam.solver.leakoff_mobility", {"reservoir_total", "upwind"}, errors);
    choice(tree, "fractureparam.solver.volume_pacing_source", {"perf", "well"}, errors);
    choice(tree, "fractureparam.solver.solve_failure_policy", {"throw", "warn_keep_state"}, errors);
    choice(tree, "fractureparam.solver.embedded_perf_wi_mode", {"fracture", "estimate"}, errors);
    choice(tree, "fractureparam.solver.embedded_well_source", {"ring", "config"}, errors);
    choice(tree, "fractureparam.solver.linsolver.solver", {"bicgstab", "gmres", "fgmres", "direct"}, errors);
    choice(tree, "fractureparam.solver.linsolver.direct_fallback_solver", {"bicgstab", "gmres", "fgmres"}, errors);
    choice(tree, "fractureparam.solver.linsolver.scaling", {"none", "blockmax"}, errors);
    choice(tree, "fractureparam.solver.linsolver.failure_policy", {"throw", "ladder"}, errors);
    choice(tree, "fractureparam.solver.linsolver.preconditioner.mode_policy", {"manual", "coupling_auto"}, errors);

    for (const auto* key : {"fractureparam.config.axis_scale", "fractureparam.config.trires",
         "fractureparam.solver.tolerance", "fractureparam.solver.max_width",
         "fractureparam.solver.max_dp", "fractureparam.solver.max_dwidth",
         "fractureparam.solver.max_iter", "fractureparam.solver.target_cellcount",
         "fractureparam.solver.cellcount_threshold", "fractureparam.solver.numcell_threshold",
         "fractureparam.solver.linsolver.tol", "fractureparam.solver.linsolver.max_iter",
         "fractureparam.solver.linsolver.restart", "fractureparam.solver.embedded_capacity",
         "fractureparam.solver.embedded_min_aperture", "fractureparam.solver.fractureWI_rw"}) {
        positive(tree, key, false, errors);
    }
    for (const auto* key : {"fractureparam.config.min_width", "fractureparam.config.initial_fracture_width",
         "fractureparam.solver.min_width", "fractureparam.solver.linsolver.atol",
         "fractureparam.solver.dt_limit", "fractureparam.solver.onset_dt_limit",
         "fractureparam.solver.max_num_coarsening", "fractureparam.reservoir.perm",
         "fractureparam.reservoir.mobility", "fractureparam.fractureWI",
         "solver.max_growth_iterations", "solver.max_mech_solves_per_step"}) {
        positive(tree, key, true, errors);
    }
    ordered(tree, "fractureparam.solver.min_width", "fractureparam.solver.max_width", errors);
    ordered(tree, "fractureparam.solver.coupling_relax_min", "fractureparam.solver.coupling_relax_max", errors);
    ordered(tree, "fractureparam.solver.coupling_anderson_alpha_min", "fractureparam.solver.coupling_anderson_alpha_max", errors);
    ordered(tree, "fractureparam.topology_hysteresis.off_threshold", "fractureparam.topology_hysteresis.on_threshold", errors);
    const auto fraction = tree.get_optional<double>("fractureparam.solver.max_closed_cell_toggle_fraction");
    if (fraction && (*fraction < 0 || *fraction > 1)) errors.push_back("fractureparam.solver.max_closed_cell_toggle_fraction: must be in [0, 1]");
    const auto relaxation = tree.get_optional<double>("fractureparam.solver.coupling_state_relaxation");
    if (relaxation && (*relaxation <= 0 || *relaxation > 1)) errors.push_back("fractureparam.solver.coupling_state_relaxation: must be in (0, 1]");

    const auto method = tree.get<std::string>("solver.method", "");
    if (method == "SeqMech" && tree.get<bool>("solver.implicit_flow", false)) {
        errors.push_back("solver.implicit_flow: true is not implemented with solver.method=SeqMech");
    }
    if (tree.get<bool>("fractureparam.solver.pressure_from_flow", false)
        && tree.get<std::string>("fractureparam.solver.fracture_flow_mode", "wi_upscaling") != "embedded") {
        errors.push_back("fractureparam.solver.pressure_from_flow: requires fracture_flow_mode=embedded");
    }
    if (tree.get<bool>("solver.growth_driven_switch", false)
        && (method != "SeqMechFrac" || tree.get<int>("solver.max_growth_iterations", 0) <= 0)) {
        errors.push_back("solver.growth_driven_switch: requires SeqMechFrac and positive max_growth_iterations");
    }
    throwErrors();
}

std::string fractureConfigurationReport(const PropertyTree& settings, const std::string& source)
{
    Tree report;
    report.put("schema_version", 1);
    report.put("source", source);
    report.put("scope", "Validated input and principal resolved modes; solver-local defaults are not exhaustively expanded. WSEED overrides are reported per seed.");
    for (const auto* key : {"add_perfs_to_schedule", "hasfractures",
         "fractureparam.control.WI", "fractureparam.solver.full_intersections",
         "fractureparam.method.iterate",
         "fractureparam.method.max_it", "fractureparam.method.tolerance",
         "fractureparam.config.verbosity", "fractureparam.fracturesolver",
         "fractureparam.solver.write_fracture_system"}) {
        if (settings.get_child_optional(key)) {
            report.put(std::string("compatibility_fields.") + key,
                "Retained legacy field; no active reader in the current simulator path");
        }
    }
    report.put("activation_source", "Deck FRAC, not the legacy hasfractures JSON field");
    report.add_child("settings", asTree(settings));
    report.put("resolved.fracture_flow_mode", settings.get<std::string>("fractureparam.solver.fracture_flow_mode", "wi_upscaling"));
    report.put("resolved.vem_face_quadrature", settings.get<std::string>("fractureparam.vem_face_quadrature", "kite"));
    report.put("resolved.closed_cell_policy", settings.get<std::string>("fractureparam.solver.closed_cell_policy", "sticky"));
    report.put("delegated_subtree", flowSolver + " (Dune::FlexibleSolver validates its own options)");
    units(report);
    return json(report);
}

PropertyTree fractureSeedConfiguration(const PropertyTree& modelSettings,
                                       double verticalExtent, double flowWidth)
{
    if (!std::isfinite(verticalExtent) || verticalExtent <= 0) {
        throw std::invalid_argument("WSEED vertical extent: must be positive and finite (m)");
    }
    if (!std::isfinite(flowWidth) || flowWidth < 0) {
        throw std::invalid_argument("WSEED WIDTH: must be nonnegative and finite (m)");
    }
    PropertyTree seed = modelSettings;
    seed.put("config.axis_scale", verticalExtent);
    seed.put("config.min_width", flowWidth);
    return seed;
}

std::string fractureSeedConfigurationReport(const PropertyTree& settings, const std::string& well,
                                            int parentIndex, int completionIndex, int gridNumber)
{
    Tree report;
    report.put("schema_version", 1);
    report.put("scope", "Exact per-seed parameter tree passed to Fracture::init; solver-local defaults are not exhaustively expanded");
    report.put("well", well);
    report.put("identity.wseed_parent_index", parentIndex);
    report.put("identity.completion_index", completionIndex);
    report.put("identity.grid_number", gridNumber);
    report.put("overrides.config_axis_scale.source", "WSEED vertical extent, converted to SI by deck parser");
    report.put("overrides.config_axis_scale.value_m", settings.get<double>("config.axis_scale"));
    report.put("overrides.config_min_width.source", "WSEED WIDTH, converted to SI by deck parser");
    report.put("overrides.config_min_width.value_m", settings.get<double>("config.min_width"));
    if (settings.get_child_optional("solver.min_width")) {
        report.put("configured_solver_min_width_m", settings.get<double>("solver.min_width"));
    }
    report.add_child("settings", asTree(settings));
    units(report);
    return json(report);
}
} // namespace Opm
