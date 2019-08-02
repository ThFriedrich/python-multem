
#ifndef MULTEM_EXT_H
#define MULTEM_EXT_H

#include <memory>

namespace mt {
  class System_Configuration;

  template<typename T>
  class Input_Multislice;
}


bool is_gpu_available_wrapper();
int number_of_gpu_available_wrapper();

class SystemConfigurationWrapper {
public:
  
  SystemConfigurationWrapper();
  SystemConfigurationWrapper(
    const char *device,
    const char *precision,
    int cpu_ncores,
    int cpu_nthread,
    int gpu_device,
    int gpu_stream);
  int get_device();
  bool is_host() const;
  bool is_device() const;
  bool is_float() const;
  bool is_double() const;
  bool is_float_host() const;
  bool is_double_host() const;
  bool is_float_device() const;
  bool is_double_device() const;

  std::shared_ptr<mt::System_Configuration> handle() {
    return handle_;
  }

private:
    
  std::shared_ptr<mt::System_Configuration> handle_; 
 
};



class InputMultisliceWrapper {
public:

  typedef float FloatType;
  typedef mt::Input_Multislice<FloatType> InputMultisliceType;

  void set_system_conf(SystemConfigurationWrapper system_conf);
  void set_interaction_model(const char *interaction_model);
  void set_potential_type(const char *potential_type);
  void set_pn_model(const char *pn_model);
  void set_pn_coh_contrib(bool pn_coh_contrib);
  void set_pn_single_conf(bool pn_single_conf);
  /* FP_Dim pn_dim;                    // Phonon dimensions */
  void set_fp_dist(int fp_dist);
  void set_pn_seed(int pn_seed);
  void set_pn_nconf(int pn_nconf);
  void set_fp_iconf_0(int fp_iconf_0);

  /* Atom_Data<T> atoms;                 // atoms */
  void set_is_crystal(bool is_crystal);
  void set_spec_rot_theta(FloatType spec_rot_theta);
  
  /* r3d<T> spec_rot_u0;                 // unitary vector */     
  void set_spec_rot_center_type(const char *spec_rot_center_type);
  /* r3d<T> spec_rot_center_p;               // rotation point */

  void set_thick_type(const char *thick_type);
  /* host_vector<T> thick;                 // Array of thickes */

  void set_potential_slicing(const char *potential_slicing);

  /* Grid_2d<T> grid_2d;                 // grid information */

  /* Range_2d output_area;               // Output region information */

  void set_simulation_type(const char *simulation_type);
  void set_iw_type(const char *iw_type);
  /* host_vector<complex<T>> iw_psi;           // user define incident wave */
  /* host_vector<T> iw_x;                // x position */
  /* host_vector<T> iw_y;                // y position */

  void set_E_0(FloatType E_0);
  void set_lambda(FloatType lambda);
  void set_theta(FloatType theta);
  void set_phi(FloatType phi);

  void set_illumination_model(const char *illumination_model);

  /* Lens<T> cond_lens;                  // Condenser lens */
  /* Lens<T> obj_lens;                   // Objective lens */

  /* Scanning<T> scanning;                 // Scanning */

  /* Detector<T, e_host> detector;             // STEM Detectors */

  /* EELS<T> eels_fr;                  // EELS */

  void set_operation_mode(const char *operation_mode);
  void set_slice_storage(bool slice_storage);
  void set_mul_sign(int mul_sign);
  void set_Vrl(FloatType Vrl);
  void set_nR(int nR);
  void set_nrot(int nrot);
  void set_cdl_var_type(const char *cdl_var_type);
  /* host_vector<T> cdl_var;               // Array of thickes */
  /* host_vector<int> iscan; */
  /* host_vector<T> beam_x;                // temporal variables for */
  /* host_vector<T> beam_y; */

  void set_islice(int islice);
  void set_dp_Shift(bool dp_Shift);                    // Shift diffraction pattern

  InputMultisliceWrapper();
  void assign(InputMultisliceWrapper &input_multislice);
  void validate_parameters();
  void set_reverse_multislice(bool rm);
  void validate_output_area();
  void set_iscan_beam_position();
  /* template<class TVector>; */
  /* void set_beam_position(TVector &x, TVector &y); */
  /* eSpace get_simulation_space() const; */
  /* void set_incident_wave_type(eIncident_Wave_Type iw_type_i); */
  /* void set_eels_fr_atom(const int &iatoms, const Atom_Data<T> &atoms); */
  
  /* FloatType Rx_exp_factor(); */
  /* FloatType Ry_exp_factor(); */
  FloatType set_incident_angle(const FloatType &theta) const;
  FloatType get_phonon_rot_weight() const;
  void set_phi(const int &irot);
  FloatType get_propagator_factor(const FloatType &z) const;
  FloatType Vr_factor() const;
  FloatType gx_0() const;
  FloatType gy_0() const;
  int number_conf();
  int number_of_beams();
  bool is_multi_beam();
  bool is_spec_rot_active() const;
  bool is_multislice() const;
  bool is_phase_object() const;
  bool is_weak_phase_object() const;
  bool is_still_atom() const;
  bool is_absorptive_model() const;
  bool is_frozen_phonon() const;
  bool is_frozen_phonon_single_conf() const;
  bool is_whole_spec() const;
  bool is_through_slices() const;
  bool is_through_thick() const;
  bool is_slicing_by_planes() const;
  bool is_slicing_by_dz() const;
  bool is_subslicing() const;
  bool is_subslicing_whole_spec() const;
  bool is_plane_wave() const;
  bool is_convergent_wave() const;
  bool is_user_define_wave() const;
  bool is_STEM() const;
  bool is_ISTEM() const;
  bool is_CBED() const;
  bool is_CBEI() const;
  bool is_ED() const;
  bool is_HRTEM() const;
  bool is_PED() const;
  bool is_HCTEM() const;
  bool is_EWFS() const;
  bool is_EWRS() const;
  bool is_EWFS_SC() const;
  bool is_EWRS_SC() const;
  bool is_EELS() const;
  bool is_EFTEM() const;
  bool is_IWFS() const;
  bool is_IWRS() const;
  bool is_PPFS() const;
  bool is_PPRS() const;
  bool is_TFFS() const;
  bool is_TFRS() const;
  bool is_PropFS() const;
  bool is_PropRS() const;
  bool is_STEM_ISTEM() const;
  bool is_CBED_CBEI() const;
  bool is_ED_HRTEM() const;
  bool is_PED_HCTEM() const;
  bool is_EWFS_EWRS() const;
  bool is_EWFS_EWRS_SC() const;
  bool is_EELS_EFTEM() const;
  bool is_IWFS_IWRS() const;
  bool is_PPFS_PPRS() const;
  bool is_TFFS_TFRS() const;
  bool is_PropFS_PropRS() const;
  bool is_grid_FS() const;
  bool is_grid_RS() const;
  bool is_simulation_type_FS() const;
  bool is_simulation_type_RS() const;
  bool is_specimen_required() const;
  bool is_ISTEM_CBEI_HRTEM_HCTEM_EFTEM() const;
  bool is_CBED_ED_EWFS_PED() const;
  bool is_obj_lens_temp_spat() const;
  bool is_cond_lens_temp_spat() const;
  bool is_scanning() const;
  bool is_illu_mod_coherent() const;
  bool is_illu_mod_partial_coherent() const;
  bool is_illu_mod_trans_cross_coef() const;
  bool is_illu_mod_full_integration() const;
  bool is_incoh_temporal_spatial() const;
  bool is_incoh_temporal() const;
  bool is_incoh_spatial() const;
  bool is_detector_circular() const;
  bool is_detector_radial() const;
  bool is_detector_matrix() const;
  bool is_slice_storage() const;
  bool is_operation_mode_normal() const;
  bool is_operation_mode_advanced() const;
  bool is_lvt_off() const;
  bool is_lvt_m() const;
  bool is_lvt_Cs3() const;
  bool is_lvt_Cs5() const;
  bool is_lvt_mfa2() const;
  bool is_lvt_afa2() const;
  bool is_lvt_mfa3() const;
  bool is_lvt_afa3() const;
  bool is_lvt_inner_aper_ang() const;
  bool is_lvt_outer_aper_ang() const;

private:

  std::shared_ptr<InputMultisliceType> handle_; 
};



#endif
