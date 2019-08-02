#include <sstream>
#include <iostream>
#include <stdexcept>
#include <cuda_runtime.h>
#include <multem/multem_ext.h>
#include <types.cuh>
#include <input_multislice.cuh>
#include <device_functions.cuh>
#include <multem.cu>

bool is_gpu_available_wrapper() {
  return mt::is_gpu_available();
}

int number_of_gpu_available_wrapper() {
  return mt::number_of_gpu_available();
}

mt::eDevice string_to_device_enum(const std::string &device) {
  if (device == "host") {
    return mt::e_host;
  } else if (device == "device") {
    return mt::e_device;
  } else if (device == "host_device") {
    return mt::e_host_device;
  }
  throw std::runtime_error("Unknown device");
}

mt::ePrecision string_to_precision_enum(const std::string &precision) {
  if (precision == "float") {
    return mt::eP_float;
  } else if (precision == "double") {
    return mt::eP_double;
  }
  throw std::runtime_error("Unknown precision");
}

mt::eElec_Spec_Int_Model string_to_interation_model_enum(const std::string &interaction_model) {
  if (interaction_model == "Multislice") {
    return mt::eESIM_Multislice;
  } else if (interaction_model == "Phase_Object") {
    return mt::eESIM_Phase_Object;
  } else if (interaction_model == "Weak_Phase_Object") {
    return mt::eESIM_Weak_Phase_Object;
  }
  throw std::runtime_error("Invalid interaction model");
}

mt::ePotential_Type string_to_potential_type_enum(const std::string &potential_type) {
  if (potential_type == "Doyle_0_4") {
    return mt::ePT_Doyle_0_4;
  } else if (potential_type == "Peng_0_4") {
    return mt::ePT_Peng_0_4;
  } else if (potential_type == "Peng_0_12") {
    return mt::ePT_Peng_0_12;
  } else if (potential_type == "Kirkland_0_12") {
    return mt::ePT_Kirkland_0_12;
  } else if (potential_type == "Weickenmeier_0_12") {
    return mt::ePT_Weickenmeier_0_12;
  } else if (potential_type == "Lobato_0_12") {
    return mt::ePT_Lobato_0_12;
  } else if (potential_type == "none") {
    return mt::ePT_none;
  }
  throw std::runtime_error("Invalid potential type");
}

mt::ePhonon_Model string_to_phonon_model_enum(const std::string &pn_model) {
  if (pn_model == "Still_Atom") {
    return mt::ePM_Still_Atom;
  } else if (pn_model == "Absorptive_Model") {
    return mt::ePM_Absorptive_Model;
  } else if (pn_model == "Frozen_Phonon") {
    return mt::ePM_Frozen_Phonon;
  }
  throw std::runtime_error("Invalid phonon model");
}

mt::eRot_Point_Type string_to_rot_point_type_enum(const std::string &spec_rot_center_type) {
  if (spec_rot_center_type == "geometric_center") {
    return mt::eRPT_geometric_center;
  } else if (spec_rot_center_type == "User_Define") {
    return mt::eRPT_User_Define;
  }
  throw std::runtime_error("Invalid spec rot center type");
}

mt::eThick_Type string_to_thick_type_enum(const std::string &thick_type) {
  if (thick_type == "Whole_Spec") {
    return mt::eTT_Whole_Spec;
  } else if (thick_type == "Through_Thick") {
    return mt::eTT_Through_Thick;
  } else if (thick_type == "Through_Slices") {
    return mt::eTT_Through_Slices;
  }
  throw std::runtime_error("Invalid thickness type");
}

mt::ePotential_Slicing string_to_potential_slicing_enum(const std::string &potential_slicing) {
  if (potential_slicing == "Planes") {
    return mt::ePS_Planes;
  } else if (potential_slicing == "dz_Proj") {
    return mt::ePS_dz_Proj;
  } else if (potential_slicing == "dz_Sub") {
    return mt::ePS_dz_Sub;
  } else if (potential_slicing == "Auto") {
    return mt::ePS_Auto;
  }
  throw std::runtime_error("Invalid potential slicing");
}

mt::eTEM_Sim_Type string_to_tem_sim_type_enum(const std::string &simulation_type) {
  if (simulation_type == "STEM") {
    return mt::eTEMST_STEM ;
  } else if (simulation_type == "ISTEM") {
    return mt::eTEMST_ISTEM ;
  } else if (simulation_type == "CBED") {
    return mt::eTEMST_CBED ;
  } else if (simulation_type == "CBEI") {
    return mt::eTEMST_CBEI ;
  } else if (simulation_type == "ED") {
    return mt::eTEMST_ED ;
  } else if (simulation_type == "HRTEM") {
    return mt::eTEMST_HRTEM ;
  } else if (simulation_type == "PED") {
    return mt::eTEMST_PED ;
  } else if (simulation_type == "HCTEM") {
    return mt::eTEMST_HCTEM ;
  } else if (simulation_type == "EWFS") {
    return mt::eTEMST_EWFS ;
  } else if (simulation_type == "EWRS") {
    return mt::eTEMST_EWRS ;
  } else if (simulation_type == "EELS") {
    return mt::eTEMST_EELS ;
  } else if (simulation_type == "EFTEM") {
    return mt::eTEMST_EFTEM ;
  } else if (simulation_type == "IWFS") {
    return mt::eTEMST_IWFS ;
  } else if (simulation_type == "IWRS") {
    return mt::eTEMST_IWRS ;
  } else if (simulation_type == "PPFS") {
    return mt::eTEMST_PPFS ;
  } else if (simulation_type == "PPRS") {
    return mt::eTEMST_PPRS ;
  } else if (simulation_type == "TFFS") {
    return mt::eTEMST_TFFS ;
  } else if (simulation_type == "TFRS") {
    return mt::eTEMST_TFRS ;
  } else if (simulation_type == "PropFS") {
    return mt::eTEMST_PropFS ;
  } else if (simulation_type == "PropRS") {
    return mt::eTEMST_PropRS ;
  }
  throw std::runtime_error("Invalid simulation type");
}

mt::eIncident_Wave_Type string_to_incident_wave_type_enum(const std::string &iw_type) {
  if (iw_type == "eIWT_Plane_Wave") {
    return mt::eIWT_Plane_Wave;
  } else if (iw_type == "eIWT_Convergent_Wave") { 
    return mt::eIWT_Convergent_Wave;
  } else if (iw_type == "eIWT_User_Define_Wave") { 
    return mt::eIWT_User_Define_Wave;
  } else if (iw_type == "eIWT_Auto") { 
    return mt::eIWT_Auto;
  }
  throw std::runtime_error("Invalid iw type");
}

mt::eIllumination_Model string_to_illumination_model_enum(const std::string &illumination_model) {
  if (illumination_model == "Coherent") {
    return mt::eIM_Coherent;
  } else if (illumination_model == "Partial_Coherent") {
    return mt::eIM_Partial_Coherent;
  } else if (illumination_model == "Trans_Cross_Coef") {
    return mt::eIM_Trans_Cross_Coef;
  } else if (illumination_model == "Full_Integration") {
    return mt::eIM_Full_Integration;
  } else if (illumination_model == "none") {
    return mt::eIM_none;
  }
  throw std::runtime_error("Invalid illumination model");
}
  
mt::eOperation_Mode string_to_operation_model_enum(const std::string &operation_mode) {
  if (operation_mode == "Normal") {
    return mt::eOM_Normal;
  } else if (operation_mode == "Advanced") {
    return mt::eOM_Advanced;
  }
  throw std::runtime_error("Invalid operation mode");
}

mt::eLens_Var_Type string_to_lens_var_type_enum(const std::string &cdl_var_type) {
  if (cdl_var_type == "off") {
    return mt::eLVT_off;
  } else if (cdl_var_type == "m") {
    return mt::eLVT_m;
  } else if (cdl_var_type == "f") {
    return mt::eLVT_f;
  } else if (cdl_var_type == "Cs3") {
    return mt::eLVT_Cs3;
  } else if (cdl_var_type == "Cs5") {
    return mt::eLVT_Cs5;
  } else if (cdl_var_type == "mfa2") {
    return mt::eLVT_mfa2;
  } else if (cdl_var_type == "afa2") {
    return mt::eLVT_afa2;
  } else if (cdl_var_type == "mfa3") {
    return mt::eLVT_mfa3;
  } else if (cdl_var_type == "afa3") {
    return mt::eLVT_afa3;
  } else if (cdl_var_type == "inner_aper_ang") {
    return mt::eLVT_inner_aper_ang;
  } else if (cdl_var_type == "outer_aper_ang") {
    return mt::eLVT_outer_aper_ang;
  }
  throw std::runtime_error("Invalid cdl_var_type");
}


SystemConfigurationWrapper::SystemConfigurationWrapper()
    : handle_(new mt::System_Configuration()) {
  handle_->validate_parameters();
}

SystemConfigurationWrapper::SystemConfigurationWrapper(
    const char *device,
    const char *precision,
    int cpu_ncores,
    int cpu_nthread,
    int gpu_device,
    int gpu_stream)
    : handle_(new mt::System_Configuration()) {
  handle_->device = string_to_device_enum(device);
  handle_->precision = string_to_precision_enum(precision);
  handle_->cpu_ncores = cpu_ncores;
  handle_->cpu_nthread = cpu_nthread;
  handle_->gpu_device = gpu_device;
  handle_->gpu_nstream = gpu_stream;
  handle_->active = true;
  handle_->validate_parameters();
  handle_->set_device();
}

int SystemConfigurationWrapper::get_device() {
  return handle_->get_device();
}

bool SystemConfigurationWrapper::is_host() const {
  return handle_->is_host();
}

bool SystemConfigurationWrapper::is_device() const {
  return handle_->is_device();
}

bool SystemConfigurationWrapper::is_float() const {
  return handle_->is_float();
}

bool SystemConfigurationWrapper::is_double() const {
  return handle_->is_double();
}

bool SystemConfigurationWrapper::is_float_host() const {
  return handle_->is_float_host();
}

bool SystemConfigurationWrapper::is_double_host() const {
  return handle_->is_double_host();
}

bool SystemConfigurationWrapper::is_float_device() const {
  return handle_->is_float_device();
}

bool SystemConfigurationWrapper::is_double_device() const {
  return handle_->is_double_device();
}






InputMultisliceWrapper::InputMultisliceWrapper()
  : handle_(new InputMultisliceType()) {}


void InputMultisliceWrapper::set_system_conf(SystemConfigurationWrapper system_conf) {
  handle_->system_conf = *(system_conf.handle());
}

void InputMultisliceWrapper::set_interaction_model(const char *interaction_model) {
  handle_->interaction_model = string_to_interation_model_enum(interaction_model);
}

void InputMultisliceWrapper::set_potential_type(const char *potential_type) {
  handle_->potential_type = string_to_potential_type_enum(potential_type);
}

void InputMultisliceWrapper::set_pn_model(const char *pn_model) {
  handle_->pn_model = string_to_phonon_model_enum(pn_model);
}

void InputMultisliceWrapper::set_pn_coh_contrib(bool pn_coh_contrib) {
  handle_->pn_coh_contrib = pn_coh_contrib;
}

void InputMultisliceWrapper::set_pn_single_conf(bool pn_single_conf) {
  handle_->pn_single_conf = pn_single_conf;
}

void InputMultisliceWrapper::set_fp_dist(int fp_dist) {
  handle_->fp_dist = fp_dist;
}

void InputMultisliceWrapper::set_pn_seed(int pn_seed) {
  handle_->pn_seed = pn_seed;
}

void InputMultisliceWrapper::set_pn_nconf(int pn_nconf) {
  handle_->pn_nconf = pn_nconf;
}

void InputMultisliceWrapper::set_fp_iconf_0(int fp_iconf_0) {
  handle_->fp_iconf_0 = fp_iconf_0;
}

void InputMultisliceWrapper::set_is_crystal(bool is_crystal) {
  handle_->is_crystal = is_crystal;
}

void InputMultisliceWrapper::set_spec_rot_theta(FloatType spec_rot_theta) {
  handle_->spec_rot_theta = spec_rot_theta;
}
  
void InputMultisliceWrapper::set_spec_rot_center_type(const char *spec_rot_center_type) {
  handle_->spec_rot_center_type = string_to_rot_point_type_enum(spec_rot_center_type);
}
  
void InputMultisliceWrapper::set_thick_type(const char *thick_type) {
  handle_->thick_type = string_to_thick_type_enum(thick_type);
}
  
void InputMultisliceWrapper::set_potential_slicing(const char *potential_slicing) {
  handle_->potential_slicing = string_to_potential_slicing_enum(potential_slicing);
}

void InputMultisliceWrapper::set_simulation_type(const char *simulation_type) {
  handle_->simulation_type = string_to_tem_sim_type_enum(simulation_type);
}

void InputMultisliceWrapper::set_iw_type(const char *iw_type) {
  handle_->iw_type = string_to_incident_wave_type_enum(iw_type);
}
  
void InputMultisliceWrapper::set_illumination_model(const char *illumination_model) {
  handle_->illumination_model = string_to_illumination_model_enum(illumination_model);
}

void InputMultisliceWrapper::set_operation_mode(const char *operation_mode) {
  handle_->operation_mode = string_to_operation_model_enum(operation_mode);
}

void InputMultisliceWrapper::set_slice_storage(bool slice_storage) {
  handle_->slice_storage = slice_storage;
}

void InputMultisliceWrapper::set_mul_sign(int mul_sign) {
  handle_->mul_sign = mul_sign;
}

void InputMultisliceWrapper::set_Vrl(FloatType Vrl) {
  handle_->Vrl = Vrl;
}

void InputMultisliceWrapper::set_nR(int nR) {
  handle_->nR = nR;
}

void InputMultisliceWrapper::set_nrot(int nrot) {
  handle_->nrot = nrot;
}

void InputMultisliceWrapper::set_cdl_var_type(const char *cdl_var_type) {
  handle_->cdl_var_type = string_to_lens_var_type_enum(cdl_var_type);
}
  
void InputMultisliceWrapper::set_islice(int islice) {
  handle_->islice = islice;
}

void InputMultisliceWrapper::set_dp_Shift(bool dp_Shift) {
  handle_->dp_Shift = dp_Shift;
}

void InputMultisliceWrapper::set_E_0(FloatType E_0) {
  handle_->E_0 = E_0;
}

void InputMultisliceWrapper::set_lambda(FloatType lambda) {
  handle_->lambda = lambda;
}

void InputMultisliceWrapper::set_theta(FloatType theta) {
  handle_->theta = theta;
}

void InputMultisliceWrapper::set_phi(FloatType phi) {
  handle_->phi = phi;
}


void InputMultisliceWrapper::assign(InputMultisliceWrapper &input_multislice) {
  handle_->assign(*(input_multislice.handle_));
}
  
void InputMultisliceWrapper::validate_parameters() {
  handle_->validate_parameters();
}
  
void InputMultisliceWrapper::validate_output_area() {
  handle_->validate_output_area();
}

void InputMultisliceWrapper::set_iscan_beam_position() {
  handle_->set_iscan_beam_position();
}

void InputMultisliceWrapper::set_reverse_multislice(bool rm) {
  handle_->set_reverse_multislice(rm);
}

/* InputMultisliceWrapper::FloatType InputMultisliceWrapper::Rx_exp_factor() { */
/*   return handle_->Rx_exp_factor(); */
/* } */

/* InputMultisliceWrapper::FloatType InputMultisliceWrapper::Ry_exp_factor() { */
/*   return handle_->Ry_exp_factor(); */
/* } */

InputMultisliceWrapper::FloatType InputMultisliceWrapper::set_incident_angle(const
    InputMultisliceWrapper::FloatType &theta) const {
  return handle_->set_incident_angle(theta);
}

InputMultisliceWrapper::FloatType InputMultisliceWrapper::get_phonon_rot_weight() const {
  return handle_->get_phonon_rot_weight();
}

void InputMultisliceWrapper::set_phi(const int &irot) {
  handle_->set_phi(irot);
}

InputMultisliceWrapper::FloatType InputMultisliceWrapper::get_propagator_factor(const
    InputMultisliceWrapper::FloatType &z) const {
  return handle_->get_propagator_factor(z);
}

InputMultisliceWrapper::FloatType InputMultisliceWrapper::Vr_factor() const {
  return handle_->Vr_factor();
}

InputMultisliceWrapper::FloatType InputMultisliceWrapper::gx_0() const {
  return handle_->gx_0();
}

InputMultisliceWrapper::FloatType InputMultisliceWrapper::gy_0() const {
  return handle_->gy_0();
}
  
int InputMultisliceWrapper::number_conf() {
  return handle_->number_conf();
}

int InputMultisliceWrapper::number_of_beams() {
  return handle_->number_of_beams();
}

bool InputMultisliceWrapper::is_multi_beam() {
  return handle_->is_multi_beam();
}
  
bool InputMultisliceWrapper::is_spec_rot_active() const {
  return handle_->is_spec_rot_active();
}

bool InputMultisliceWrapper::is_multislice() const {
  return handle_->is_multislice();
}

bool InputMultisliceWrapper::is_phase_object() const {
  return handle_->is_phase_object();
}

bool InputMultisliceWrapper::is_weak_phase_object() const {
  return handle_->is_weak_phase_object();
}

bool InputMultisliceWrapper::is_still_atom() const {
  return handle_->is_still_atom();
}

bool InputMultisliceWrapper::is_absorptive_model() const {
  return handle_->is_absorptive_model();
}

bool InputMultisliceWrapper::is_frozen_phonon() const {
  return handle_->is_frozen_phonon();
}

bool InputMultisliceWrapper::is_frozen_phonon_single_conf() const {
  return handle_->is_frozen_phonon_single_conf();
}

bool InputMultisliceWrapper::is_whole_spec() const {
  return handle_->is_whole_spec();
}

bool InputMultisliceWrapper::is_through_slices() const {
  return handle_->is_through_slices();
}

bool InputMultisliceWrapper::is_through_thick() const {
  return handle_->is_through_thick();
}

bool InputMultisliceWrapper::is_slicing_by_planes() const {
  return handle_->is_slicing_by_planes();
}

bool InputMultisliceWrapper::is_slicing_by_dz() const {
  return handle_->is_slicing_by_dz();
}

bool InputMultisliceWrapper::is_subslicing() const {
  return handle_->is_subslicing();
}

bool InputMultisliceWrapper::is_subslicing_whole_spec() const {
  return handle_->is_subslicing_whole_spec();
}

bool InputMultisliceWrapper::is_plane_wave() const {
  return handle_->is_plane_wave();
}

bool InputMultisliceWrapper::is_convergent_wave() const {
  return handle_->is_convergent_wave();
}

bool InputMultisliceWrapper::is_user_define_wave() const {
  return handle_->is_user_define_wave();
}

bool InputMultisliceWrapper::is_STEM() const {
  return handle_->is_STEM();
}

bool InputMultisliceWrapper::is_ISTEM() const {
  return handle_->is_ISTEM();
}

bool InputMultisliceWrapper::is_CBED() const {
  return handle_->is_CBED();
}

bool InputMultisliceWrapper::is_CBEI() const {
  return handle_->is_CBEI();
}

bool InputMultisliceWrapper::is_ED() const {
  return handle_->is_ED();
}

bool InputMultisliceWrapper::is_HRTEM() const {
  return handle_->is_HRTEM();
}

bool InputMultisliceWrapper::is_PED() const {
  return handle_->is_PED();
}

bool InputMultisliceWrapper::is_HCTEM() const {
  return handle_->is_HCTEM();
}

bool InputMultisliceWrapper::is_EWFS() const {
  return handle_->is_EWFS();
}

bool InputMultisliceWrapper::is_EWRS() const {
  return handle_->is_EWRS();
}

bool InputMultisliceWrapper::is_EWFS_SC() const {
  return handle_->is_EWFS_SC();
}

bool InputMultisliceWrapper::is_EWRS_SC() const {
  return handle_->is_EWRS_SC();
}

bool InputMultisliceWrapper::is_EELS() const {
  return handle_->is_EELS();
}

bool InputMultisliceWrapper::is_EFTEM() const {
  return handle_->is_EFTEM();
}

bool InputMultisliceWrapper::is_IWFS() const {
  return handle_->is_IWFS();
}

bool InputMultisliceWrapper::is_IWRS() const {
  return handle_->is_IWRS();
}

bool InputMultisliceWrapper::is_PPFS() const {
  return handle_->is_PPFS();
}

bool InputMultisliceWrapper::is_PPRS() const {
  return handle_->is_PPRS();
}

bool InputMultisliceWrapper::is_TFFS() const {
  return handle_->is_TFFS();
}

bool InputMultisliceWrapper::is_TFRS() const {
  return handle_->is_TFRS();
}

bool InputMultisliceWrapper::is_PropFS() const {
  return handle_->is_PropFS();
}

bool InputMultisliceWrapper::is_PropRS() const {
  return handle_->is_PropRS();
}

bool InputMultisliceWrapper::is_STEM_ISTEM() const {
  return handle_->is_STEM_ISTEM();
}

bool InputMultisliceWrapper::is_CBED_CBEI() const {
  return handle_->is_CBED_CBEI();
}

bool InputMultisliceWrapper::is_ED_HRTEM() const {
  return handle_->is_ED_HRTEM();
}

bool InputMultisliceWrapper::is_PED_HCTEM() const {
  return handle_->is_PED_HCTEM();
}

bool InputMultisliceWrapper::is_EWFS_EWRS() const {
  return handle_->is_EWFS_EWRS();
}

bool InputMultisliceWrapper::is_EWFS_EWRS_SC() const {
  return handle_->is_EWFS_EWRS_SC();
}

bool InputMultisliceWrapper::is_EELS_EFTEM() const {
  return handle_->is_EELS_EFTEM();
}

bool InputMultisliceWrapper::is_IWFS_IWRS() const {
  return handle_->is_IWFS_IWRS();
}

bool InputMultisliceWrapper::is_PPFS_PPRS() const {
  return handle_->is_PPFS_PPRS();
}

bool InputMultisliceWrapper::is_TFFS_TFRS() const {
  return handle_->is_TFFS_TFRS();
}

bool InputMultisliceWrapper::is_PropFS_PropRS() const {
  return handle_->is_PropFS_PropRS();
}

bool InputMultisliceWrapper::is_grid_FS() const {
  return handle_->is_grid_FS();
}

bool InputMultisliceWrapper::is_grid_RS() const {
  return handle_->is_grid_RS();
}

bool InputMultisliceWrapper::is_simulation_type_FS() const {
  return handle_->is_simulation_type_FS();
}

bool InputMultisliceWrapper::is_simulation_type_RS() const {
  return handle_->is_simulation_type_RS();
}

bool InputMultisliceWrapper::is_specimen_required() const {
  return handle_->is_specimen_required();
}

bool InputMultisliceWrapper::is_ISTEM_CBEI_HRTEM_HCTEM_EFTEM() const {
  return handle_->is_ISTEM_CBEI_HRTEM_HCTEM_EFTEM();
}

bool InputMultisliceWrapper::is_CBED_ED_EWFS_PED() const {
  return handle_->is_CBED_ED_EWFS_PED();
}

bool InputMultisliceWrapper::is_obj_lens_temp_spat() const {
  return handle_->is_obj_lens_temp_spat();
}

bool InputMultisliceWrapper::is_cond_lens_temp_spat() const {
  return handle_->is_cond_lens_temp_spat();
}

bool InputMultisliceWrapper::is_scanning() const {
  return handle_->is_scanning();
}

bool InputMultisliceWrapper::is_illu_mod_coherent() const {
  return handle_->is_illu_mod_coherent();
}

bool InputMultisliceWrapper::is_illu_mod_partial_coherent() const {
  return handle_->is_illu_mod_partial_coherent();
}

bool InputMultisliceWrapper::is_illu_mod_trans_cross_coef() const {
  return handle_->is_illu_mod_trans_cross_coef();
}

bool InputMultisliceWrapper::is_illu_mod_full_integration() const {
  return handle_->is_illu_mod_full_integration();
}

bool InputMultisliceWrapper::is_incoh_temporal_spatial() const {
  return handle_->is_incoh_temporal_spatial();
}

bool InputMultisliceWrapper::is_incoh_temporal() const {
  return handle_->is_incoh_temporal();
}

bool InputMultisliceWrapper::is_incoh_spatial() const {
  return handle_->is_incoh_spatial();
}

bool InputMultisliceWrapper::is_detector_circular() const {
  return handle_->is_detector_circular();
}

bool InputMultisliceWrapper::is_detector_radial() const {
  return handle_->is_detector_radial();
}

bool InputMultisliceWrapper::is_detector_matrix() const {
  return handle_->is_detector_matrix();
}

bool InputMultisliceWrapper::is_slice_storage() const {
  return handle_->is_slice_storage();
}

bool InputMultisliceWrapper::is_operation_mode_normal() const {
  return handle_->is_operation_mode_normal();
}

bool InputMultisliceWrapper::is_operation_mode_advanced() const {
  return handle_->is_operation_mode_advanced();
}

bool InputMultisliceWrapper::is_lvt_off() const {
  return handle_->is_lvt_off();
}

bool InputMultisliceWrapper::is_lvt_m() const {
  return handle_->is_lvt_m();
}

bool InputMultisliceWrapper::is_lvt_Cs3() const {
  return handle_->is_lvt_Cs3();
}

bool InputMultisliceWrapper::is_lvt_Cs5() const {
  return handle_->is_lvt_Cs5();
}

bool InputMultisliceWrapper::is_lvt_mfa2() const {
  return handle_->is_lvt_mfa2();
}

bool InputMultisliceWrapper::is_lvt_afa2() const {
  return handle_->is_lvt_afa2();
}

bool InputMultisliceWrapper::is_lvt_mfa3() const {
  return handle_->is_lvt_mfa3();
}

bool InputMultisliceWrapper::is_lvt_afa3() const {
  return handle_->is_lvt_afa3();
}

bool InputMultisliceWrapper::is_lvt_inner_aper_ang() const {
  return handle_->is_lvt_inner_aper_ang();
}

bool InputMultisliceWrapper::is_lvt_outer_aper_ang() const {
  return handle_->is_lvt_outer_aper_ang();
}






/* typedef float FloatType; */

/* void simulate( */
/*     mt::System_Configuration &system_conf, */
/*     mt::Input_Multislice<FloatType> &input_multislice, */ 
/*     mt::Output_Multislice<FloatType> &output_multislice) { */

/*   mt_run_multislice<FloatType, mt::e_host>( */
/*     system_conf, */
/*     input_multislice, */ 
/*     output_multislice); */
/* } */
