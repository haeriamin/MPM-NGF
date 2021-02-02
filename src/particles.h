/*******************************************************************************
    Copyright (c) The Taichi MPM Authors (2018- ). All Rights Reserved.
    The use of this software is governed by the LICENSE file.
*******************************************************************************/

#pragma once

#include <iostream>
#include <taichi/util.h>
#include <taichi/math/svd.h>
#include <taichi/math/array.h>
#include <taichi/math/levelset.h>
#include "mpm_fwd.h"

TC_NAMESPACE_BEGIN

template <int dim>
class MPMParticle : public Unit {
 public:
  using Vector = VectorND<dim, real>;
  using VectorP = VectorND<dim + 1, real>;
  using Matrix = MatrixND<dim, real>;
  using Region = RegionND<dim>;

 private:
  VectorP v_and_m;

 public:
  Vector pos;
  Matrix dg_e;  // Elastic deformation gradient
  Matrix apic_b;  // Affine momemtum (APIC), c234
  Vector boundary_normal;
  real boundary_distance;
  real vol;
  int dt_limit;
  int stiffness_limit;
  int cfl_limit;
  bool near_boundary_;
  bool sticky;
  bool is_rigid_;
  uint32 states;
  int32 id;
  real debug;
  Matrix dg_t;  // added: total deformation gradient
  Matrix dg_p;  // added: plastic deformation gradient
  Matrix T;  // added: stress tensor
  real p;  // added: normal stress (pressure)
  real tau;  // added: shear stress
  real gf;  // added: granular fluidity
  // real mu_visual;  // added: friction coeff (FOR VISUALIZATION ONLY)
  // bool is_free;  // added (FOR VISUALIZATION ONLY)
  Vector rigid_impulse;  // added: Impulses on rigid boundary particles
  // Matrix apic_c;  // added: c567
  // real dg_p_det;
  // real dg_p_inv_det;

  // keep the order
  // TC_IO_DEF_VIRT(v_and_m,
  //                pos,
  //                dg_e,
  //                apic_b,
  //                boundary_normal,
  //                boundary_distance,
  //                vol,
  //                dt_limit,
  //                stiffness_limit,
  //                cfl_limit,
  //                sticky,
  //                near_boundary_,
  //                states,
  //                id,
  //                is_rigid_,
  //                dg_t,
  //                dg_p,
  //                T,
  //                p,
  //                tau,
  //                gf,
  //                mu_visual,
  //                is_free,
  //                rigid_impulse,
  //                dg_p_det,
  //                dg_p_inv_det,
  //               );
  TC_IO_DEF_VIRT(v_and_m,
                 pos,
                 dg_e,
                 apic_b,
                 boundary_normal,
                 boundary_distance,
                 vol,
                 dt_limit,
                 stiffness_limit,
                 cfl_limit,
                 sticky,
                 near_boundary_,
                 states,
                 id,
                 is_rigid_,
                 dg_t,
                 dg_p,
                 T,
                 p,
                 tau,
                 gf,
                 rigid_impulse
                );

  TC_FORCE_INLINE bool is_rigid() const {
    return is_rigid_;
  }

  TC_FORCE_INLINE real get_mass() const {
    return v_and_m[dim];
  }

  TC_FORCE_INLINE void set_mass(real mass) {
    v_and_m[dim] = mass;
  }

  TC_FORCE_INLINE Vector get_velocity() const {
    return Vector(v_and_m);
  }

  TC_FORCE_INLINE void set_velocity(const Vector &v) {
    v_and_m = VectorP(v, v_and_m[dim]);
  }

  TC_FORCE_INLINE void set_velocity(const __m128 &v) {
    v_and_m = _mm_blend_ps(v, v_and_m, 0x7);
  }

  // TC_FORCE_INLINE real get_gf() const {
  //   return gf;
  // }

  MPMParticle() {
    dg_e = Matrix(1.0_f);
    apic_b = Matrix(0);
    v_and_m = VectorP(0.0f);
    vol = 1.0f;
    states = 0;
    boundary_normal = Vector(0.0_f);
    boundary_distance = 0.0_f;
    dt_limit = 1;
    stiffness_limit = 1;
    cfl_limit = 1;
    near_boundary_ = false;
    id = 0;
    is_rigid_ = false;
    dg_t = Matrix(1.0_f);
    dg_p = Matrix(1.0_f);
    T = Matrix(0.0_f);
    p = 1.0_f;
    tau = 0.0_f;
    gf = 1.0_f;  // it should be set correctly, was 1 or 10
    // mu_visual = 0.0_f;
    // is_free = false;
    rigid_impulse = Vector(0.0_f);
    // apic_c = Matrix(0);
    // dg_p_det = 1.0_f;
    // dg_p_inv_det = 1.0_f;
  }

 public:
  bool near_boundary() const {
    return near_boundary_;
  }

  uint32 get_state(int i) const {
    return (states >> (i * 2)) % 4;
  }

  void set_state(int i, int state) {
    states = states ^ ((state ^ get_state(i)) << (i * 2));
  }

  virtual real get_allowed_dt(const real &dx) const = 0;

  virtual void initialize(const Config &config) {
    if (config.has_key("compressibility")) {
      TC_ERROR("'compressibility' is deprecated. Use 'initial_dg' instead");
    }
    dg_e = Matrix(config.get("initial_dg", 1.0_f));
  }

  virtual real potential_energy() const {
    TC_NOT_IMPLEMENTED
    return 0;
  }

  virtual Matrix first_piola_kirchhoff() {
    TC_NOT_IMPLEMENTED
    return Matrix(0.0f);
  }

  virtual Matrix calculate_force() {
    TC_NOT_IMPLEMENTED
    return Matrix(0.0f);
  }

  virtual int plasticity(const Matrix &cdg, const real &laplacian_gf) {
    return 0;
  }

  virtual Matrix get_first_piola_kirchoff_differential(const Matrix &dF) {
    return Matrix(0.0f);
  }

  virtual ~MPMParticle() {
  }

  virtual real get_stiffness() const {
    TC_NOT_IMPLEMENTED
    return 0.0f;
  }

  virtual Vector3 get_debug_info() const {
    return Vector3(0);
  }

  virtual void set_mu_to_zero() {
  }

  virtual void set_lambda_and_mu_to_zero() {
  }

};

using MPMParticle2D = MPMParticle<2>;
using MPMParticle3D = MPMParticle<3>;
TC_INTERFACE(MPMParticle2D);
TC_INTERFACE(MPMParticle3D);

#define TC_REGISTER_MPM_PARTICLE(name)                                \
  using name##Particle2D = name##Particle<2>;                         \
  using name##Particle3D = name##Particle<3>;                         \
  static_assert(                                                      \
      sizeof(name##Particle2D) <= get_particle_size_upper_bound<2>(), \
      "2D MPM particle (" #name ") cannot exceed 192B");              \
  static_assert(                                                      \
      sizeof(name##Particle3D) <= get_particle_size_upper_bound<3>(), \
      "3D MPM particle (" #name ") cannot exceed 256B");              \
  TC_IMPLEMENTATION(MPMParticle2D, name##Particle2D,                  \
                    name##Particle2D().get_name());                   \
  TC_IMPLEMENTATION(MPMParticle3D, name##Particle3D,                  \
                    name##Particle2D().get_name());

TC_NAMESPACE_END
