#include "trm_trans.h"

using namespace std;
using namespace tinyxml2;

transport_module::transport_module()
{
    volume.set_metadata("", "", "denseVector", "double", "3");
    volume.set_default_value(1);
    porosity.set_metadata("", "", "denseVector", "double", "0");
    porosity.set_default_value(0.2);
    conductivity_long.set_metadata("", "", "denseVector", "double", "3");
    conductivity_long.set_default_value(1);
    conductivity_trans.set_metadata("", "", "denseVector", "double", "3");
    conductivity_trans.set_default_value(1);
    position.set_metadata("", "", "sparseMatrix", "int", "0");

    save_csv_debug = false;
    numerical_tolerance = numeric_limits<double>::quiet_NaN();
}

/**
 * - sets name, version and other module metadata
 * @param name module name
 * @param version module version
 * @param type module type
 * @param next_config next config file
 */
void transport_module::set_metadata(const string& name, const string& version, const string& type, double num_tolerance)
{
    trm_module::set_metadata(name, version);
    this->type = type;
    this->numerical_tolerance = num_tolerance;
}

/**
 * - checks if version in input file matches with the TRM version (version of the transport module represents version of whole TRM)
 */
void transport_module::check_version()
{
    if (this->version != TRM_VERSION) {
        throw trm_error("Input file version '" + version + "' does not match with the software version '" + TRM_VERSION + "'.");
    }
}

/**
 * - gets component for given variable
 * @param var_name variable name in camel case, e.g. "Volume" or "HydraulicHead"
 * @param bc_time boundary conditions change time (negative if not relevant for this variable)
 * @return pointer to the component, nullptr if component for given time is not initialized
 */
component* transport_module::get_variable(const string& var_name, const int bc_time)
{
    if (var_name == "InOutFlow" || var_name == "HydraulicHead" || var_name == "GridFlow") {
        if (bc_time < 0) {
            throw trm_error("Invalid boundary conditions change time for variable " + var_name + ".");
        }
        unsigned bc_time_unsigned = static_cast<unsigned>(bc_time);
        if (var_name == "HydraulicHead") {
            if (hydraulic_head.find(bc_time_unsigned) == hydraulic_head.end()) {
                return nullptr;
            }
            return &hydraulic_head[bc_time_unsigned];
        }
        else if (var_name == "InOutFlow") {
            return &in_out_flow[bc_time_unsigned];
        }
        else if (var_name == "GridFlow") {
            if (grid_flow.find(bc_time_unsigned) == grid_flow.end()) {
                return nullptr;
            }
            return &grid_flow[bc_time_unsigned];
        }
    }
    else {
        if (var_name == "Volume") {
            return &volume;
        }
        else if (var_name == "Porosity") {
            return &porosity;
        }
        else if (var_name == "Position") {
            return &position;
        }
        else if (var_name == "HydraulicConductivityLongitudinal") {
            return &conductivity_long;
        }
        else if (var_name == "HydraulicConductivityTransverse") {
            return &conductivity_trans;
        }
        else {
            throw trm_error("Invalid variable name '" + var_name + "'.");
        }
    }
    return nullptr;
}

/**
 * - gets variable with all componets of given category
 * @param category category name (one of "InComponents", "TransportComponents")
 * @return pointer to the time-based map of ID-based maps of components
 */
map<unsigned, map<string, component>>* transport_module::get_components(const string& category)
{
    if (category == "InComponents") {
        return &in_components;
    }
    else if (category == "TransportComponents") {
        return &transport_components;
    }
    else {
        throw trm_error("Invalid component category name '" + category + "'.");
    }
}

bool transport_module::set_grid_size_component(unsigned nrows, unsigned ncols, const string& var_name, int bc_time)
{
    component *curr_component = get_variable(var_name, bc_time);
    bool init_variables = false;
    if (curr_component != nullptr && !curr_component->is_empty()) {
        if (var_name == "InOutFlow" || var_name == "Volume" || var_name == "Porosity" ||
            var_name == "HydraulicConductivityLongitudinal" || var_name == "HydraulicConductivityTransverse") {
            // these are vectors for filling the grid -> the old grid, not the vector needs to be resized
            curr_component->resize_vector(this->nrows, this->ncols, nrows, ncols);
        }
        // resizing makes no sense, clean the component (to be initialized during next calculation)
        else if (var_name == "Position" || var_name == "GridFlow" || var_name == "HydraulicHead") {
            curr_component->set_size(0, 0);
            init_variables = true;
        }
    }
    return init_variables;
}

/**
 * - sets number of cells, resizes corresponding matrices and removes results
 * @param nrows number of rows
 * @param ncols number of columns
 */
void transport_module::set_grid_size(unsigned nrows, unsigned ncols)
{
    if (this->nrows != nrows || this->ncols != ncols) {
        vector<string> var_names = {
            "HydraulicHead", "Volume", "Porosity", "HydraulicConductivityLongitudinal", "HydraulicConductivityTransverse",
            "InOutFlow", "Position", "GridFlow"
        };
        bool init_variables = false;
        for (auto& var_name : var_names) {
            if (var_name == "InOutFlow" || var_name == "HydraulicHead" || var_name == "GridFlow") {
                for (auto& bc_time_step : bc_times_steps) {
                    set_grid_size_component(nrows, ncols, var_name, static_cast<int>(bc_time_step.first));
                }
            }
            else {
                set_grid_size_component(nrows, ncols, var_name, -1);
            }
        }
        remove_results();
        for (auto& in_comp_for_time : in_components) {
            for (auto& in_comp : in_comp_for_time.second) {
                in_comp.second.resize(nrows * ncols, 1);
            }
        }
        for (auto& comp_for_time : transport_components) {
            for (auto& comp : comp_for_time.second) {
                comp.second.resize(nrows * ncols, 1);
            }
        }

        this->nrows = nrows;
        this->ncols = ncols;
        ncells = nrows * ncols;

        // initialize only if there were any values before
        // otherwise this can be setting when loading from file when output during initialization is not ready
        if (init_variables) {
            init_position();
            init_hydraulic_head(0);
            init_grid_flow(0);
            modify_grid_flow(0);
        }
    }
}

/**
 * - sets one characteristics value for given component
 * @param category component category ("TransportComponents" only)
 * @param component_name name of component
 * @param char_id ID of characteristics
 * @param value characteristics value to be set
 */
void transport_module::set_component_char(const std::string& category, const std::string& component_name, const std::string& char_id, double value)
{
    if (category == "TransportComponents") {
        string component_id = find_id_for_component_name(transport_components.begin()->second, component_name);
        transport_components_chars[component_id][char_id] = value;
    }
    else {
        throw trm_error("Invalid name of component '" + category + "' for characteristics.");
    }
}

/**
 * - gets one characteristics value for given component
 * @param category component category ("TransportComponents" only)
 * @param component_name name of component
 * @param char_id ID of characteristics
 * @return characteristics value
 */
double transport_module::get_component_char(const std::string& category, const std::string& component_name, const std::string& char_id)
{
    map<string, map<string, double>> *chars;
    string component_id;
    if (category == "TransportComponents") {
        chars = &transport_components_chars;
        component_id = find_id_for_component_name(transport_components.begin()->second, component_name);
    }
    else {
        throw trm_error("Invalid name of component '" + category + "' for characteristics.");
    }
    if (chars->find(component_id) != chars->end() && (*chars)[component_id].find(char_id) != (*chars)[component_id].end()) {
        return (*chars)[component_id][char_id];
    }
    else {
        throw trm_error("Characteristics '" + char_id + "' for component '" + component_id + "' and category '" + category + "' does not exist.");
    }
}

/**
 * - adds component
 * @param comp component instance
 * @param category component category ("InComponents", "TransportComponents")
 * @param time_step time step for the instance (boundary conditions change time for "InComponents")
 */
void transport_module::add_component(component& comp, const string& category, unsigned time_step)
{
    map<unsigned, map<string, component>> *components = get_components(category);
    (*components)[time_step][comp.get_id()] = comp;
}

/**
 * - gets names of all components from given category belonging to the first time (if relevant)
 * @param comp component_category category of transport module ("InComponents", "TransportComponents")
 * @return component names
 */
vector<string> transport_module::get_component_names(const string& component_category)
{
    vector<string> comp_names;
    map<unsigned, map<string, component>> *components = get_components(component_category);
    if (!components->empty()) {
        for (auto& comp : components->begin()->second) {
            comp_names.push_back(comp.second.get_name());
        }
    }
    return comp_names;
}

/**
 * - gets times for given component category
 * @param comp component_category category of transport module ("InComponents", "TransportComponents")
 * @return time IDs
 */
vector<unsigned> transport_module::get_component_times(const string& component_category)
{
    vector<unsigned> time_ids;
    map<unsigned, map<string, component>> *components = get_components(component_category);
    if (!components->empty()) {
        for (auto& comp : (*components)) {
            time_ids.push_back(comp.first);
        }
    }
    return time_ids;
}

/**
 * - gets data for given component
 * @param component_category category of transport module ("InComponents", "TransportComponents")
 *     or variable ("Volume", "Porosity", "InOutFlow" etc.)
 * @param component_name name of particular component
 * @param time_id ID of required time, i.e. boundary conditions change time for for InComponents or InOutFlow
 * @return required data as a simple vector
 */
vector<double> transport_module::get_component_data(const string& component_category, const string& component_name, unsigned time_id)
{
    arma::mat tmp_data;
    if (component_category == "Volume" || component_category == "Porosity" ||
        component_category == "HydraulicHead" || component_category == "InOutFlow" ||
        component_category == "HydraulicConductivityLongitudinal" || component_category == "HydraulicConductivityTransverse") {
        component *tmp_component = get_variable(component_category, time_id);
        if (tmp_component == nullptr) {
            ostringstream oss;
            oss << time_id;
            throw trm_error("Component for variable '" + component_category + "' and time " + oss.str() + " is not available.");
        }
        else if (component_category == "InOutFlow") {
            tmp_data = get_variable(component_category, time_id)->get_data_sparse();
        }
        else {
            tmp_data = get_variable(component_category, time_id)->get_data();
        }
    }
    else {
        map<unsigned, map<string, component>> *components = get_components(component_category);
        string comp_id = trm_module::find_id_for_component_name((*components)[time_id], component_name);
        tmp_data = get_component_data(component_category, time_id, comp_id);
    }
    return arma::conv_to<vector<double>>::from(tmp_data);
}

/**
 * - sets data for given component
 * @param data data to be set as a simple vector
 * @param component_category category of module ("InComponents", "TransportComponents")
 *     or variable ("Volume", "Porosity", "InOutFlow" etc.)
 * @param component_name name of particular component (not relevant for variables)
 * @param time_id ID of required time
 */
void transport_module::set_component_data(const vector<double>& data, const string& component_category, const string& component_name, unsigned time_id)
{
    if (data.size() != nrows * ncols) {
        throw trm_error("Provided data size does not match with the component size.");
    }
    if (component_category == "InOutFlow" || component_category == "HydraulicHead" ||
        component_category == "Volume" || component_category == "Porosity" ||
        component_category == "HydraulicConductivityLongitudinal" || component_category == "HydraulicConductivityTransverse") {
        component *tmp_component = get_variable(component_category, time_id);
        if (tmp_component == nullptr) {
            ostringstream oss;
            oss << time_id;
            throw trm_error("Component for variable '" + component_category + "' and time " + oss.str() + " is not available.");
        }
        arma::mat tmp_data(data);
        if (tmp_component->is_sparse_type()) {
            arma::sp_mat new_data = arma::conv_to<arma::sp_mat>::from(tmp_data);
            tmp_component->set_data_sparse(new_data, false);
        }
        else {
            tmp_component->set_data(tmp_data, false);
        }
    }
    else if (component_category == "InComponents" || component_category == "TransportComponents") {
        string comp_id;

        map<unsigned, map<string, component>> *components = get_components(component_category);
        comp_id = trm_module::find_id_for_component_name((*components)[time_id], component_name);
        if (components->find(time_id) == components->end()) {
            ostringstream oss;
            oss << "Invalid time '" << time_id << "' for component '" << component_name << "' of transport module.";
            throw trm_error(oss.str());
        }
        else if ((*components)[time_id].find(comp_id) == (*components)[time_id].end()) {
            throw trm_error("Invalid name '" + component_name + "' for component of transport module.");
        }

        arma::mat tmp_data(data);
        if ((*components)[time_id][comp_id].is_sparse_type()) {
            arma::sp_mat new_data = arma::conv_to<arma::sp_mat>::from(tmp_data);
            (*components)[time_id][comp_id].set_data_sparse(new_data, false);
        }
        else {
            (*components)[time_id][comp_id].set_data(tmp_data, false);
        }
    }
    else {
        throw trm_error("Invalid transport component category '" + component_category + "'.");
    }
}

/**
 * - creates a new component from scratch
 * @param category component category ("InComponents", "TransportComponents")
 * @param time_id time ID (for in component boundary conditions change time)
 * @param time actual value of time
 * @param id component ID
 * @param name component name
 * @param type matrix type
 * @param data_type data type
 * @param quantity_id not used
 * @param data component data
 */
void transport_module::create_component(const string& category, unsigned time_id, double time, const string& id, const string& name, const string& type, const string& data_type, const string& quantity_id, arma::mat data)
{
    map<unsigned, map<std::string, component>> *components = get_components(category);
    (*components)[time_id][id].set_metadata(id, name, type, data_type, quantity_id);
    if (category == "InComponents") {
        arma::sp_mat sp_data(data);
        (*components)[time_id][id].set_data_sparse(sp_data, true);
    }
    else {
        (*components)[time_id][id].set_data(data, true);
    }
    (*components)[time_id][id].set_time(time);
}

/**
 * - calculates temporary position matrices needed for initialization of variables
 * @param side_size size of side of the squared matrices
 */
void transport_module::calc_positions(unsigned side_size)
{
    arma::sp_mat tmp_empty(side_size, side_size);
    position1 = tmp_empty;
    position4 = position3 = position2 = position1;

    unsigned i, j, ii, jj;
    arma::mat tmp_conductivity_long = conductivity_long.get_data();
    arma::mat tmp_conductivity_trans = conductivity_trans.get_data();
    for (ii = 0; ii < nrows; ii++) {
        for (jj = 0; jj < ncols; jj++) {
            i = ii * ncols + jj;
            for (j = 0; j < side_size; j++) {
                double conductivity_long_value = (tmp_conductivity_long(i) + tmp_conductivity_long(j)) / 2;
                double conductivity_trans_value = (tmp_conductivity_trans(i) + tmp_conductivity_trans(j)) / 2;
                if (j == i - ncols && ii > 0) {
                    position3(i, j) = conductivity_long_value;
                    position3(i, i) = position3(i, i) - conductivity_long_value;
                }
                if (j == i - 1 && jj > 0) {
                    position1(i, j) = conductivity_long_value;
                    position1(i, i) = position1(i, i) - conductivity_long_value;
                }
                if (j == i + 1 && jj < ncols - 1) {
                    position2(i, j) = conductivity_trans_value;
                    position2(i, i) = position2(i, i) - conductivity_trans_value;
                }
                if (j == i + ncols && ii < nrows - 1) {
                    position4(i, j) = conductivity_trans_value;
                    position4(i, i) = position4(i, i) - conductivity_trans_value;
                }
            }
        }
    }
}

/**
 * - checks if conductivity is available, if it is not, initializes it with default value 1
 */
void transport_module::check_conductivity()
{
    if (conductivity_long.get_data().n_elem != get_cell_count("HydraulicConductivityLongitudinal")) {
        conductivity_long.resize(ncells, 1);
    }
    if (conductivity_trans.get_data().n_elem != get_cell_count("HydraulicConductivityTransverse")) {
        conductivity_trans.resize(ncells, 1);
    }
}

/**
 * - calculates matrix of positions if it is not already defined with corresponding size
 */
void transport_module::init_position()
{
    check_conductivity();
    if (position.get_data_sparse().n_elem != get_cell_count("Position")) {
        arma::sp_mat tmp_position(ncells, ncells);

        calc_positions(ncells);
        tmp_position = position1 + position2 + position3 + position4;
        position.set_data_sparse(tmp_position, true);
        log(BASIC) << "Position values initialized.";
    }
}

/**
 * - calculates matrix of diagonal inverse volume to be used in transport calculation
 * - before that, volume of liquid part is calculated from initial volume and porosity
 */
void transport_module::init_volume()
{
    if (porosity.get_data().n_elem != get_cell_count("Porosity")) {
        porosity.resize(ncells, 1);
    }
    volume_diag_inverse = volume.get_data() % porosity.get_data();
    volume_diag_inverse = diagmat(arma::ones<arma::vec>(ncells) / volume_diag_inverse, 0);
}

/**
 * - calculates vector of hydraulic head if it is not already defined with corresponding size
 * - initializes positions if they are not ready
 * @param bc_time_id boundary conditions change time for which the head will be calculated
 */
void transport_module::init_hydraulic_head(unsigned bc_time)
{
    if (hydraulic_head.find(bc_time) == hydraulic_head.end()) {
        component tmp_component;
        tmp_component.set_metadata("", "", "denseVector", "double", "2");
        hydraulic_head[bc_time] = tmp_component;
    }
    if (position.get_data_sparse().n_elem != get_cell_count("Position")) {
        init_position();
    }
    unsigned cells_head = get_cell_count("HydraulicHead");
    if (hydraulic_head[bc_time].get_data().n_elem != cells_head) {
        arma::sp_mat position_red(position.get_data_sparse());
        position_red.shed_row(0);
        position_red.shed_col(0);
        arma::mat in_out_flow_red(in_out_flow[bc_time].get_data_sparse());
        in_out_flow_red.shed_row(0);
        arma::vec head_red(cells_head - 1);
        head_red = arma::spsolve(position_red, in_out_flow_red, "lapack"); //TDD use superlu
        head_red.insert_rows(0, 1);
        hydraulic_head[bc_time].set_data(head_red, false, true);

        log(BASIC) << "Hydraulic head values for boundary conditions change time " << bc_time << " initialized.";
    }
}

/**
 * - converts sp_imat to sp_mat (direct conversion is not provided by Armadillo)
 * @param sp_imatrix sparse matrix of integers
 * @return sparse matrix of doubles
 */
arma::sp_mat transport_module::conv_sp_imat_to_sp_mat(arma::sp_imat sp_imatrix)
{
    arma::imat imatrix(sp_imatrix);
    arma::mat matrix = arma::conv_to<arma::mat>::from(imatrix);
    arma::sp_mat sp_matrix(matrix);
    return sp_matrix;
}

/**
 * - calculates matrix of grid flow if it is not already defined with corresponding size
 * @param bc_time grid flow will be initialized for time of this boundary condition change
 */
void transport_module::init_grid_flow(unsigned bc_time)
{
    if (grid_flow.find(bc_time) == grid_flow.end()) {
        component tmp_component;
        tmp_component.set_metadata("", "", "sparseMatrix", "double", "4");
        grid_flow[bc_time] = tmp_component;
    }
    arma::sp_mat tmp_grid_flow = grid_flow[bc_time].get_data_sparse();
    if (tmp_grid_flow.n_elem != get_cell_count("GridFlow")) {
        check_conductivity();
        if (position1.n_elem == 0) {
            calc_positions(ncells);
        }
        if (hydraulic_head[bc_time].get_data().n_elem == 0) {
            init_hydraulic_head(bc_time);
        }
        tmp_grid_flow.set_size(ncells, ncells);
        tmp_grid_flow.zeros();
        arma::sp_mat head_as_sp = arma::conv_to<arma::sp_mat>::from(hydraulic_head[bc_time].get_data());
        arma::sp_mat flow1 = position1 * head_as_sp;
        arma::sp_mat flow2 = position2 * head_as_sp;
        arma::sp_mat flow3 = position3 * head_as_sp;
        arma::sp_mat flow4 = position4 * head_as_sp;

        // for 1D grid, two conditions in this loop can be fulfilled at once
        // for one of them, zero is assigned -> adding works here - if non-zero
        // value is the result of the first condition, it is not replaced by zero in the second one
        for (unsigned i = 0; i < ncells; i++) {
            for (unsigned j = 0; j < ncells; j++) {
                if (j == i - ncols && i > ncols - 1) {
                    tmp_grid_flow(i, j) += -flow3(j + ncols, 0);
                }
                if (j == i - 1 && i > 0) {
                    tmp_grid_flow(i, j) += -flow1(j + 1, 0);
                }
                if (j == i + 1) {
                    tmp_grid_flow(i, j) += -flow2(j - 1, 0);
                }
                if (j == i + ncols) {
                    tmp_grid_flow(i, j) += -flow4(j - ncols, 0);
                }
            }
        }
        grid_flow[bc_time].set_data_sparse(tmp_grid_flow, false, true);
        log(BASIC) << "Grid flow values for boundary conditions change time " << bc_time << " initialized.";
    }
}

/**
 * - checks if inflow/outflow values sum up to 0
 * - value of numerical tolerance or default value based on numeric limits is used
 * @param bc_time boundary condition time for which inflow/outflow is checked
 */
void transport_module::check_in_out_flow(unsigned bc_time)
{
    double mat_sum = arma::accu(in_out_flow[bc_time].get_data_sparse());
    double tolerance = isnan(numerical_tolerance) ? 10 * numeric_limits<double>::epsilon() : numerical_tolerance;
    if (!(mat_sum > -tolerance && mat_sum < tolerance)) {
        ostringstream oss;
        oss << tolerance;
        throw trm_error("Values of inflow/outflow have to sum up to 0 (with tolerance " + oss.str() + ").");
    }
}

/**
 * - modifies matrix of the grid flow to be prepared for run
 * @param bc_time_id boundary conditions change time whose grid flow will be used
 */
void transport_module::modify_grid_flow(unsigned bc_time)
{
    check_in_out_flow(bc_time);
    arma::sp_mat tmp_in_out_flow = in_out_flow[bc_time].get_data_sparse();
    // use only negative values of in/outflow
    tmp_in_out_flow.transform([](double value) { return (value > 0 ? 0 : value); });

    grid_flow_modif = grid_flow[bc_time].get_data_sparse();
    grid_flow_modif.transform([](double value) { return (value > 0 ? 0 : value); });

    arma::vec row_sums(ncells);
    for (unsigned i = 0; i < ncells; i++) {
        row_sums(i) = arma::accu(grid_flow_modif.row(i));
    }
    grid_flow_modif.diag() = -row_sums - tmp_in_out_flow;
    grid_flow_modif = -grid_flow_modif;

    in_out_flow_modif.set_size(ncells, ncells);
    in_out_flow_modif.diag() = tmp_in_out_flow.col(0);
    in_out_flow_modif = -in_out_flow_modif;

    if (save_csv_debug) {
        ostringstream oss;
        oss << "GridFlowModified" << bc_time << ".csv";
        vector<string> empty;
        output->serialize_variable(grid_flow_modif, "TransportModule/GridElements2GridElements/GridFlow", empty, empty, nullptr, oss.str(), "");
    }

    log(BASIC) << "Grid flow values for boundary conditions change time " << bc_time << " modified for run.";
}

/**
 * - initializes transport components by values obtained from the activation of reaction
 * - called only for the first time step
 * - components are stored (and outputted then) always, even if the initial time is not chosen in output settings
 * @param concentrations concentrations of transport components from the reaction part
 * @param component_names names of transport components
 * @param initial_time_index index of the first time step (for saving components)
 */
void transport_module::init(vector<double>& concentrations, vector<string>& component_names, unsigned initial_time_index)
{
    vector<string>::size_type ncomps = component_names.size();
    for (vector<string>::size_type comp = 0; comp < ncomps; comp++) {
        arma::mat conc_mat(ncells, 1);
        for (unsigned c = 0; c < ncells; c++) {
            conc_mat(c) = concentrations[comp * ncells + c];
        }
        component trans_comp;
        trans_comp.copy_metadata(get_first_component());
        ostringstream oss;
        oss << comp + 1;
        trans_comp.set_id(oss.str());
        trans_comp.set_name(component_names[comp]);
        trans_comp.set_data(conc_mat, true);
        add_component(trans_comp, "TransportComponents", initial_time_index);
    }
    set_time_end_id();

    log(BASIC) << "Initialization of transport module ended.";
}

/**
 * - activates transport module, i.e. calculates positions, hydraulic head and grid flow
 */
void transport_module::activate()
{
    if (has_state(DISABLED)) {
        return;
    }
    log(BASIC) << "Activation of transport module started.";

    init_position();
    init_volume();
    for (auto bc_time_step : bc_times_steps) {
        unsigned bc_time = bc_time_step.first;
        init_hydraulic_head(bc_time);
        init_grid_flow(bc_time);
    }

    log(BASIC) << "Transport module activated.";
}

/**
 * - sets boundary conditions for transport module, i.e. adds InComponents
 * @param component_names names of transport components
 * @param boundary_concs boundary concentrations by cell index (unsigned map key) and component (position in vector)
 * @param bc_time boundary conditions change time used for initialization of grid components
 */
void transport_module::set_boundary_conditions(vector<string>& component_names, map<unsigned, vector<double>>& boundary_concs, unsigned bc_time)
{
    if (has_state(DISABLED)) {
        return;
    }
    log(BASIC) << "Initialization of transport module started.";

    vector<string>::size_type ncomps = component_names.size();
    for (vector<string>::size_type comp = 0; comp < ncomps; comp++) {
        // in components setting is compulsory
        arma::sp_mat bound_conc_mat(ncells, 1);
        component in_comp;
        for (auto& bcond : boundary_concs) {
            bound_conc_mat(bcond.first) = bcond.second[comp];
        }
        ostringstream oss;
        oss << comp + 1;
        in_comp.set_id(oss.str());
        in_comp.set_name(component_names[comp]);
        in_comp.set_data_type("double");
        in_comp.set_quantity_id("101");
        in_comp.set_data_sparse(bound_conc_mat, true);
        add_component(in_comp, "InComponents", bc_time);

    }

    unsigned t_index = get_time_init_id();
    if (transport_components.find(t_index) == transport_components.end()) {
        ostringstream oss;
        oss << "No transport component available for initial time step with ID '" << t_index << "'.";
        throw trm_error(oss.str());
    }
    transport_components_mat.set_size(ncells, ncomps + 1);

    log(BASIC) << "Boundary conditions set for time step " << bc_time << ".";
}

/**
 * - runs the model for given time
 * @param t_index index of time to be calculated
 * @param time corresponding time
 * @param concentrations concentrations of transport components
 * @param component_names names of transport components
 * @param component_names names of transport components
 * @param curr_bc_time current boundary conditions change time, values of in components for it will be used
 * @param mapping values of mapping component
 */
void transport_module::calculate(unsigned t_index, double time, vector<double>& concentrations, vector<string>& component_names, unsigned curr_bc_time, vector<double>& mapping)
{
    if (has_state(DISABLED) || has_state(DISABLED_CALCULATION)) {
        return;
    }

    unsigned col, row;
    // initialization for t_index 1 is enough, in next steps content of the first column of transport_component_mat is kept
    if (t_index == 1) {
        for (row = 0; row < ncells; row++) {
            transport_components_mat(row, 0) = get_component_data("TransportComponents", t_index - 1, "0")(row, 0);
        }
    }
    unsigned c = 0;
    for (col = 1; col < transport_components_mat.n_cols; col++) {
        for (row = 0; row < ncells; row++) {
            transport_components_mat(row, col) = concentrations[c];
            c++;
        }
    }

    // concentrations of in components taken from boundary conditions
    arma::sp_mat in_comps(ncells, transport_components_mat.n_cols);
    for (col = 0; col < transport_components_mat.n_cols; col++) {
        ostringstream oss;
        oss << col;
        string col_str = oss.str();
        arma::sp_mat comp_data = in_components[curr_bc_time][col_str].get_data_sparse();
        for (unsigned row = 0; row < ncells; row++) {
            in_comps(row, col) = comp_data(row, 0);
        }
    }
    transport_components_mat += volume_diag_inverse * (in_out_flow_modif * in_comps * bc_times_steps[curr_bc_time] + grid_flow_modif * transport_components_mat * bc_times_steps[curr_bc_time]);

    // change mapping sent to the reaction module
    mapping = arma::conv_to<vector<double>>::from(transport_components_mat.col(0));

    // change concentration sent to the reaction module
    c = 0;
    for (col = 1; col < transport_components_mat.n_cols; col++) {
        for (row = 0; row < ncells; row++) {
            concentrations[c] = transport_components_mat(row, col);
            c++;
        }
    }

    if (output->is_step_for_output(t_index, time_end_id) && output->get_cell_output_type() != "none") {
        vector<string> tmp_component_names = component_names;
        tmp_component_names.insert(tmp_component_names.begin(), 1, get_first_component().get_name());

        vector<string>::size_type ncomps = tmp_component_names.size();
        for (vector<string>::size_type cn = 0; cn < ncomps; cn++) {
            arma::mat tmp_comp(ncells, 1);
            for (row = 0; row < ncells; row++) {
                tmp_comp(row, 0) = transport_components_mat(row, cn);
            }

            component trans_comp;
            trans_comp.copy_metadata(get_first_component());
            ostringstream oss;
            oss << cn;
            trans_comp.set_id(oss.str());
            trans_comp.set_name(tmp_component_names[cn]);
            trans_comp.set_data(tmp_comp, true);
            trans_comp.set_time(time);
            add_component(trans_comp, "TransportComponents", t_index);
        }
    }
    log(BASIC) << "Transport for time step ID '" << t_index << "' (time " << time << ") calculated.";
}

/**
 * - gets first transport component (the first one for the first time)
 */
component& transport_module::get_first_component()
{
    if (transport_components.size() == 0) {
        throw trm_error("The first component is required but no one is available.");
    }
    return transport_components.begin()->second.begin()->second;
}

/**
 * - gets given transport component for the first time
 * @param id ID of required transport component
 */
component& transport_module::get_first_component_of_id(string id)
{
    if (transport_components.size() == 0) {
        throw trm_error("The first component of '" + id + "'is required but no transport component is available.");
    }
    if (transport_components.begin()->second.find(id) == transport_components.begin()->second.end()) {
        throw trm_error("The first component of '" + id + "' is required but this ID is not available.");
    }
    return transport_components.begin()->second[id];
}

/**
 * - removes calculated transport components
 */
void transport_module::remove_results()
{
    // last read from input file with time_init_id has to be kept
    for (unsigned t_index = time_init_id + 1; t_index <= time_end_id; t_index++) {
        transport_components.erase(t_index);
    }
}

/**
 * - resets variables initialized during activation (position, hydraulic head, grid flow)
 */
void transport_module::remove_activation()
{
    position.set_size(0, 0);
    for (auto bc_time_step : bc_times_steps) {
        unsigned bc_time = bc_time_step.first;
        hydraulic_head[bc_time].set_size(0, 0);
        grid_flow[bc_time].set_size(0, 0);
    }
}

/**
 * - gets data for component of transport module
 * @param component_category one of "InComponents", "TransportComponents"
 * @param time_index ID of time, i.e. boundary conditions change time for InComponents
 * @param component_id ID of component
 * @return matrix of component values
 */
arma::mat transport_module::get_component_data(const string& component_category, unsigned time_index, string component_id)
{
    map<unsigned, map<string, component>> *components = get_components(component_category);
    if (components->find(time_index) == components->end()) {
        ostringstream oss;
        oss << "Invalid time '" << time_index << "' for component of transport module.";
        throw trm_error(oss.str());
    }
    else if ((*components)[time_index].find(component_id) == (*components)[time_index].end()) {
        throw trm_error("Invalid ID '" + component_id + "' for transport component.");
    }
    if (component_category == "InComponents") {
        arma::mat tmp_mat((*components)[time_index][component_id].get_data_sparse());
        return tmp_mat;
    }
    else {
        return (*components)[time_index][component_id].get_data();
    }
}

/**
 * - serializes all transport components to XML tree
 * @param xmldoc XML document needed for creating nodes
 * @param before_elem serialized components will be added after this element
 */
void transport_module::serialize_transport_components(XMLDocument *xmldoc, XMLElement *before_elem)
{
    serialize_components(xmldoc, before_elem, transport_components, "TransportComponents", "transport");
}

/**
 * - serializes all in components to XML tree
 * @param before_elem serialized components will be added after this element
 */
void transport_module::serialize_in_components(XMLElement *before_elem)
{
    for (auto& in_comp_for_bc_time : in_components) {
        XMLElement *in_comps_elem = output->get_xmldoc()->NewElement("InComponents");
        in_comps_elem->SetAttribute("idTimeStep", in_comp_for_bc_time.first);
        before_elem->Parent()->InsertAfterChild(before_elem, in_comps_elem);
        for (auto& in_comp : in_comp_for_bc_time.second) {
            ostringstream oss;
            oss << "InComponent" << in_comp_for_bc_time.first << ".csv";
            in_comp.second.change_to_vector();
            in_comp.second.serialize(oss.str(), in_comps_elem, output, "transport", "both");
        }
        before_elem = before_elem->NextSiblingElement();
    }
}

void transport_module::serialize_grid_component(const string& var_name, int bc_time, const string& xml_path)
{
    component *curr_component;
    curr_component = get_variable(var_name, bc_time);
    if (curr_component != nullptr && !curr_component->is_empty()) {
        string var_name_bc_time = var_name;
        if (bc_time > -1) {
            ostringstream oss;
            oss << var_name << bc_time;
            var_name_bc_time = oss.str();
        }

        curr_component->serialize(
            var_name_bc_time + ".csv", nullptr, output, "transport", "both", xml_path + var_name, bc_time, false);
    }
}

/**
 * - serializes grid components (volume, in/out flow etc.) to XML tree
 * @param parent_elem serialized components will be added as children of this element
 * @param replace_elems whether the existing elements in document (if they exist) will be replaced (this argument will not be needed when save_inputs is complete)
 */
void transport_module::serialize_grid_components(XMLElement *parent_elem, bool replace_elems)
{
    XMLElement *grid_elements_elem;
    if (replace_elems) {
        grid_elements_elem = output->get_child_elem(parent_elem, "GridElements");
    }
    else {
        grid_elements_elem = output->create_element(parent_elem, "GridElements", true);
    }
    vector<string> var_names = { "Volume", "Porosity", "HydraulicHead", "InOutFlow", "HydraulicConductivityLongitudinal", "HydraulicConductivityTransverse" };
    for (auto &var_name : var_names) {
        if (replace_elems) {
            output->delete_child_elems(grid_elements_elem, var_name);
        }
        if (var_name == "HydraulicHead" || var_name == "InOutFlow") {
            for (auto &bc_time_step : bc_times_steps) {
                serialize_grid_component(var_name, static_cast<int>(bc_time_step.first), "TransportModule/GridElements/");
            }
        }
        else {
            serialize_grid_component(var_name, -1, "TransportModule/GridElements/");
        }
    }
    XMLElement *grid2grid_elements_elem;
    if (replace_elems) {
        grid2grid_elements_elem = output->get_child_elem(parent_elem, "GridElements2GridElements");
    }
    else {
        grid2grid_elements_elem = output->create_element(parent_elem, "GridElements2GridElements", true);
    }
    var_names = { "Position", "GridFlow" };
    for (auto &var_name : var_names) {
        if (replace_elems) {
            output->delete_child_elems(grid2grid_elements_elem, var_name);
        }
        if (var_name == "GridFlow") {
            for (auto &bc_time_step : bc_times_steps) {
                serialize_grid_component(var_name, static_cast<int>(bc_time_step.first), "TransportModule/GridElements2GridElements/");
            }
        }
        else {
            serialize_grid_component(var_name, -1, "TransportModule/GridElements2GridElements/");
        }
    }
}

/**
 * - serializes whole transport module XML tree to the document associated with the transport part (xmldoc attribute)
 */
void transport_module::serialize()
{
    XMLDocument *xmldoc = output->get_xmldoc();
    XMLElement *trans_module_elem = output->create_element(xmldoc->FirstChildElement("TRM"), "TransportModule", true);
    trans_module_elem->SetAttribute("name", name.c_str());
    trans_module_elem->SetAttribute("version", version.c_str());
    trans_module_elem->SetAttribute("type", type.c_str());
    if (!isnan(numerical_tolerance)) {
        trans_module_elem->SetAttribute("numTolerance", numerical_tolerance);
    }
    XMLElement *time_elem = output->create_element(trans_module_elem, "Time", true);
    serialize_time(xmldoc, time_elem);

    XMLElement *grid_size_elem = output->create_element(trans_module_elem, "GridSize", true);
    serialize_grid_size(grid_size_elem);
    serialize_grid_components(trans_module_elem, false);
    serialize_in_components(trans_module_elem->LastChildElement());

    XMLElement *in_comps_elem = trans_module_elem->LastChildElement("InComponents");
    serialize_transport_components(xmldoc, in_comps_elem != nullptr ? in_comps_elem : grid_size_elem);
    XMLElement *child = nullptr, *parent = xmldoc->FirstChildElement("TRM");
    output->find_in_xml_tree("TransportModule/TransportComponents", &parent, &child);
    if (child == nullptr) {
        output->create_element(trans_module_elem, "TransportComponents", true);
    }
    else {
        output->get_child_elem(trans_module_elem, "TransportComponents", false, true);
    }
}
