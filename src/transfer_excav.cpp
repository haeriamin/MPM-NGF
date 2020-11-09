/*******************************************************************************
    Copyright (c) The Taichi MPM Authors (2018- ). All Rights Reserved.
    The use of this software is governed by the LICENSE file.
*******************************************************************************/

#include <taichi/system/profiler.h>
#include <taichi/math/sparse.h>
#include "mpm.h"
#include "kernel.h"
#include "taichi/dynamics/rigid_body.h"
#include "boundary_particle.h"
#include <taichi/common/testing.h>

#ifndef MPM_TRANSFER_OPT

TC_NAMESPACE_BEGIN

// #define TC_MPM_USE_LOCKS
#ifdef TC_MPM_USE_LOCKS
#define LOCK_GRID grid_locks[ind].lock();
#define UNLOCK_GRID grid_locks[ind].unlock();
#else
#define LOCK_GRID
#define UNLOCK_GRID
#endif

#define MLSMPM

#if defined(MLSMPM)
constexpr bool use_mls_mpm = true;
#else
constexpr bool use_mls_mpm = false;
#endif

// x64 Intrinsics
// Do not move them to taichi main lib since MSVC forces
// CPU intrinsics to have compile-time const parameters
TC_FORCE_INLINE float32 extract_float32(const __m128 &s, int i) {
  int ret = _mm_extract_ps(s, i);
  return reinterpret_cast<float32 *>(&ret)[0];
}
/*
TC_FORCE_INLINE __m128 broadcast(const __m128 &s, int i) {
  return _mm_shuffle_ps(s, s, 0x55 * i);
}
*/
#define broadcast(s, i) _mm_shuffle_ps((s), (s), 0x55 * (i))

// grid cache ------------------------------------------------------------------
template <typename MPM, bool v_and_m_only = false>
struct GridCache {
  static_assert(mpm_kernel_order == 2, "Only supports quadratic kernel");

  using SparseMask = typename MPM::SparseMask;

  static constexpr int dim = 3;
  static constexpr int scratch_x_size = (1 << SparseMask::block_xbits) + 2;
  static constexpr int scratch_y_size = (1 << SparseMask::block_ybits) + 2;
  static constexpr int scratch_z_size = (1 << SparseMask::block_zbits) + 2;
  static constexpr int scratch_size =
      scratch_x_size * scratch_y_size * scratch_z_size;

  static constexpr int num_nodes = pow<dim>(mpm_kernel_order + 1); // quad3D: 27

  using ElementType =
      std::conditional_t<v_and_m_only, Vector4f, GridState<dim>>;

  using GridCacheType =
      ElementType[scratch_x_size][scratch_y_size][scratch_z_size];
  using GridCacheLinearizedType = ElementType[scratch_size];

  using SparseGrid = typename MPM::SparseGrid;

  SparseGrid &grid;
  uint64 block_offset;
  bool write_back;

  TC_ALIGNED(64) GridCacheType blocked;
  GridCacheLinearizedType &linear =
      *reinterpret_cast<GridCacheLinearizedType *>(&blocked[0][0][0]);

  static constexpr int kernel_linearized(int x) {
    return ((x / 9) * scratch_y_size * scratch_z_size +
            (x / 3 % 3) * scratch_z_size + x % 3);
  }

  // constructor
  TC_FORCE_INLINE GridCache(SparseGrid &grid,
                            const uint64 &block_offset,
                            bool write_back)
      : grid(grid), block_offset(block_offset), write_back(write_back) {
    Vector3i block_base_coord(MPM::SparseMask::LinearToCoord(block_offset));
    auto grid_array = grid.Get_Array();
    for (int i = 0; i < scratch_x_size; i++) {
      for (int j = 0; j < scratch_y_size; j++) {
        for (int k = 0; k < scratch_z_size; k++) {
          TC_STATIC_IF(v_and_m_only) {
            id(blocked[i][j][k]) =
                grid_array(to_std_array(block_base_coord + Vector3i(i, j, k)))
                    .velocity_and_mass;
          }
          TC_STATIC_ELSE {
            id(blocked[i][j][k]) =
                grid_array(to_std_array(block_base_coord + Vector3i(i, j, k)));
          }
          TC_STATIC_END_IF
        }
      }
    }
  }

  ~GridCache() {
    if (!write_back) {
      return;
    }
    Vector3i block_base_coord(MPM::SparseMask::LinearToCoord(block_offset));
    auto grid_array = grid.Get_Array();
    for (int i = 0; i < scratch_x_size; i++) {
      for (int j = 0; j < scratch_y_size; j++) {
        for (int k = 0; k < scratch_z_size; k++) {
          TC_STATIC_IF(v_and_m_only) {
            id(grid_array(to_std_array(block_base_coord + Vector3i(i, j, k)))
                .velocity_and_mass) = blocked[i][j][k];
          }
          TC_STATIC_ELSE {
            id(grid_array(to_std_array(block_base_coord + Vector3i(i, j, k)))) =
                blocked[i][j][k];
          }
          TC_STATIC_END_IF
        }
      }
    }
  }

  TC_FORCE_INLINE static constexpr int linearized_offset(int x, int y, int z) {
    return x * scratch_y_size * scratch_z_size + y * scratch_z_size + z;
  }

  TC_FORCE_INLINE static Vector3i spgrid_block_linear_to_vector(int elem) {
    int elem_x = (elem >> (SparseMask::block_zbits + SparseMask::block_ybits));
    int elem_y = ((elem >> SparseMask::block_zbits) &
                  ((1 << SparseMask::block_ybits) - 1));
    int elem_z = elem & ((1 << SparseMask::block_zbits) - 1);
    return Vector3i(elem_x, elem_y, elem_z);
  }

  TC_FORCE_INLINE static constexpr int spgrid_block_to_grid_cache_block(
      int elem) {
    int elem_x = (elem >> (SparseMask::block_zbits + SparseMask::block_ybits));
    int elem_y = ((elem >> SparseMask::block_zbits) &
                  ((1 << SparseMask::block_ybits) - 1));
    int elem_z = elem & ((1 << SparseMask::block_zbits) - 1);
    return linearized_offset(elem_x, elem_y, elem_z);
  }
};

TC_FORCE_INLINE __m128 make_float4(float32 a, float32 b, float32 c, float32 d) {
  return _mm_set_ps(d, c, b, a);
}

// MLS-MPM fast kernel ---------------------------------------------------------
struct MLSMPMFastKernel32 {
  static constexpr int dim = 3;
  using Base = MPMKernelBase<dim, 2>;
  using Vector = typename Base::Vector;
  TC_ALIGNED(64) __m128 kernels[3][3];

  TC_FORCE_INLINE MLSMPMFastKernel32(const __m128 &pos, real inv_delta_x) {
    __m128 p_fract = _mm_sub_ps(pos, _mm_set1_ps(0.5f));
    TC_ALIGNED(64) __m128 w_cache[dim];
    for (int k = 0; k < dim; k++) {
      __m128 t = _mm_sub_ps(_mm_set1_ps(p_fract[k]),
                            make_float4(-0.5f, 0.5f, 1.5f, 0.0f));
      __m128 tt = _mm_mul_ps(t, t);
      w_cache[k] =
          _mm_fmadd_ps(make_float4(0.5f, -1.0f, 0.5f, 0.0f), tt,
                       _mm_fmadd_ps(make_float4(-1.5f, 0.0f, 1.5f, 0.0f), t,
                                    make_float4(1.125f, 0.75f, 1.125f, 0.0f)));
    }
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        kernels[i][j] =
            _mm_mul_ps(_mm_set1_ps(w_cache[0][i] * w_cache[1][j]), w_cache[2]);
      }
    }
  }

  TC_FORCE_INLINE static int get_stencil_start(real x) {
    return int(x - 0.5f);
  }
};

// rasterize ------------------------------------------------------------- : OFF
template <int dim>
void MPM<dim>::rasterize(real delta_t, bool with_force) {
  for (auto &r : this->rigids) {
    r->reset_tmp_velocity();
  }
  parallel_for_each_particle([&](Particle &p) {
    if (p.is_rigid()) {
      return;
    }
    if (particle_gravity) {
      p.set_velocity(p.get_velocity() + gravity * delta_t);
    }
    // Note, pos is magnified grid pos
    const Vector pos = p.pos * inv_delta_x;
    const Vector v = p.get_velocity();
    const real mass = p.get_mass();
    // Note, apic_b has delta_x issue
    const Matrix apic_b_inv_d_mass = p.apic_b * (Kernel::inv_D() * mass);
    const Vector mass_v = mass * v;
    Matrix delta_t_tmp_force(0.0f);
    if (with_force) {
      delta_t_tmp_force = delta_t * p.calculate_force();
    }
    RegionND<dim> region(VectorI(0), VectorI(Kernel::kernel_size));

    Vectori grid_base_pos = get_grid_base_pos(pos);
    Kernel kernel(pos, inv_delta_x);

    for (auto &ind : region) {
      auto i = ind.get_ipos() + grid_base_pos;
      Vector dpos = pos - i.template cast<real>();
      VectorP dw_w = kernel.get_dw_w(ind.get_ipos());
      GridState<dim> &g = get_grid(i);

      // Coloring
      uint64 grid_state = g.get_states(), particle_state = p.states;
      uint64 mask = (grid_state & particle_state & state_mask) >> 1;

      if ((grid_state & mask) != (particle_state & mask)) {
        // Different color **********
        // Directly project instead of writing to the grid

        RigidBody<dim> *r = get_rigid_body_ptr(get_grid(i).get_rigid_body_id());
        if (r == nullptr)
          continue;

        Vector grid_pos = i.template cast<real>() * delta_x;

        Vector rigid_v = r->get_velocity_at(grid_pos);

        Vector v_acc = v;  // + apic_b_inv_d_mass * dpos / Vector(mass);

        Vector velocity_change =
            v_acc -
            friction_project(v_acc, rigid_v, p.boundary_normal,
                             r->frictions[(particle_state >> (2 * r->id)) % 2]);

        Vector impulse = mass * dw_w[dim] * velocity_change +
                         delta_t_tmp_force * Vector(dw_w);
        r->apply_tmp_impulse(impulse, grid_pos);
        continue;
      }

      VectorP delta;
      if (mpm_kernel_order == 1) {
        delta = dw_w[dim] * VectorP(mass_v, mass) +
                VectorP(mass * p.apic_b * Vector(dw_w) * delta_x) +
                VectorP(delta_t_tmp_force * Vector(dw_w));
      } else {
#ifdef MLSMPM
        delta = dw_w[dim] *
                (VectorP(mass_v + apic_b_inv_d_mass * dpos, mass) +
                 VectorP(-delta_t_tmp_force * dpos * 4.0_f * inv_delta_x));
#else
        delta = dw_w[dim] * VectorP(mass_v + apic_b_inv_d_mass * dpos, mass) +
                VectorP(delta_t_tmp_force * Vector(dw_w));
#endif
      }
      g.velocity_and_mass += delta;
    }
  });
  for (auto &r : rigids) {
    r->apply_tmp_velocity();
  }
  // TC_P(E);
}

template <>
void MPM<2>::rasterize_optimized(real delta_t) {
  rasterize(delta_t, true);
}
TC_FORCE_INLINE __m128 make_float3(float a, float b, float c) {
  return make_float4(a, b, c, 0);
}

// clang-format off
TC_ALIGNED(64) const static __m128 grid_pos_offset_[27] = {
    make_float3(0, 0, 0),  // 0
    make_float3(0, 0, 1),  // 1
    make_float3(0, 0, 2),  // 2
    make_float3(0, 1, 0),  // 3
    make_float3(0, 1, 1),  // 4, i-1
    make_float3(0, 1, 2),  // 5
    make_float3(0, 2, 0),  // 6
    make_float3(0, 2, 1),  // 7
    make_float3(0, 2, 2),  // 8
    make_float3(1, 0, 0),  // 9
    make_float3(1, 0, 1),  // 10, j-1
    make_float3(1, 0, 2),  // 11
    make_float3(1, 1, 0),  // 12, k-1
    make_float3(1, 1, 1),  // 13, i,j,k
    make_float3(1, 1, 2),  // 14, k+1
    make_float3(1, 2, 0),  // 15
    make_float3(1, 2, 1),  // 16, j+1
    make_float3(1, 2, 2),  // 17
    make_float3(2, 0, 0),  // 18
    make_float3(2, 0, 1),  // 19
    make_float3(2, 0, 2),  // 20
    make_float3(2, 1, 0),  // 21
    make_float3(2, 1, 1),  // 22, i+1
    make_float3(2, 1, 2),  // 23
    make_float3(2, 2, 0),  // 24
    make_float3(2, 2, 1),  // 25
    make_float3(2, 2, 2),  // 26
};
// clang-format on

// clang-format off
TC_ALIGNED(64) const static Vector3f grid_pos_offset[27] = {
    Vector3f(0, 0, 0),
    Vector3f(0, 0, 1),
    Vector3f(0, 0, 2),
    Vector3f(0, 1, 0),
    Vector3f(0, 1, 1),
    Vector3f(0, 1, 2),
    Vector3f(0, 2, 0),
    Vector3f(0, 2, 1),
    Vector3f(0, 2, 2),
    Vector3f(1, 0, 0),
    Vector3f(1, 0, 1),
    Vector3f(1, 0, 2),
    Vector3f(1, 1, 0),
    Vector3f(1, 1, 1),
    Vector3f(1, 1, 2),
    Vector3f(1, 2, 0),
    Vector3f(1, 2, 1),
    Vector3f(1, 2, 2),
    Vector3f(2, 0, 0),
    Vector3f(2, 0, 1),
    Vector3f(2, 0, 2),
    Vector3f(2, 1, 0),
    Vector3f(2, 1, 1),
    Vector3f(2, 1, 2),
    Vector3f(2, 2, 0),
    Vector3f(2, 2, 1),
    Vector3f(2, 2, 2),
};
// clang-format on

TC_TEST("grid_pos_offset") {
  for (int i = 0; i < 27; i++) {
    CHECK(grid_pos_offset[i].x == i / 9);
    CHECK(grid_pos_offset[i].y == i / 3 % 3);
    CHECK(grid_pos_offset[i].z == i % 3);
  }
}

// optimized rasterization function --------------------------------------- : ON
template <>
void MPM<3>::rasterize_optimized(real delta_t) {
  constexpr int dim = 3;
  for (auto &r : this->rigids) {
    r->reset_tmp_velocity();
  }

  // block_op_rigid, called from block_op_switch -------------------------------
  auto block_op_rigid = [&](uint32 b, uint64 block_offset, GridState<dim> *g_) {
    using Cache = GridCache<MPM<dim>>;
    Cache grid_cache(*grid, block_offset, true);
    int particle_begin;
    int particle_end = block_meta[b].particle_offset;

    for (uint32 t = 0; t < SparseMask::elements_per_block; t++) {
      particle_begin = particle_end;
      particle_end += g_[t].particle_count;
      int grid_cache_offset = grid_cache.spgrid_block_to_grid_cache_block(t);

      Vectori grid_base_pos = Vectori(SparseMask::LinearToCoord(block_offset)) +
                              grid_cache.spgrid_block_linear_to_vector(t);
      Vector grid_base_pos_f = Vector(grid_base_pos);

      Vector grid_pos[27];
      for (int i = 0; i < 27; i++) {
        grid_pos[i] = grid_pos_offset[i] + grid_base_pos_f;
      }

      // added: Reset forces on rigid body boundary particles
      if (config_backup.get("visualize_particle_impulses", false))
      {
        for (int r_p_i = particle_begin; r_p_i < particle_end; r_p_i++) {
          Particle &r_p = *allocator[particles[r_p_i]];
          if (r_p.is_rigid())
            r_p.rigid_impulse = Vector(0.0_f);
        }
      }

      for (int p_i = particle_begin; p_i < particle_end; p_i++) {
        Particle &p = *allocator[particles[p_i]];
        if (p.is_rigid()) {
          continue;
        }
        // add particle gravity ------------------------------------------------
        if (particle_gravity) {
          p.set_velocity(p.get_velocity() + gravity * delta_t);
        }
        // Note, pos is magnified (0-res) grid pos ------------- PARTICLE POS ??
        const Vector pos = p.pos * inv_delta_x;

        // particle kernel
        Kernel kernel(pos, inv_delta_x);
        const VectorP(&kernels)[3][3][3] = kernel.kernels;
        using KernelLinearized = VectorP[27];
        const KernelLinearized &kernels_linearized =
            *reinterpret_cast<const KernelLinearized *>(&kernels[0][0][0]);

        const Vector v = p.get_velocity();
        const real mass = p.get_mass();

        // added: Disconnection handling
        real gf = 0.0_f;
        if (p.p > 0.0_f)  // pressure @ n
          gf = p.gf;
        Matrix delta_t_tmp_force(0.0_f);
        delta_t_tmp_force = (delta_t * p.calculate_force());

        // Note, apic_b has delta_x issue
        const Matrix apic_b_inv_d_mass = p.apic_b * (Kernel::inv_D() * mass);
        const Vector mass_v = mass * v;

        for (int node_id = 0; node_id < Cache::num_nodes; node_id++) {
          Vector dpos = pos - grid_pos[node_id];

          GridState<dim> &g =
              grid_cache.linear[grid_cache.kernel_linearized(node_id) +
                                grid_cache_offset];

          const VectorP &dw_w = kernels_linearized[node_id];

          // Coloring
          uint64 grid_state = g.get_states(), particle_state = p.states;
          uint64 mask = (grid_state & particle_state & state_mask) >> 1;

          // incompatible grid and particle ------------------------------------
          if ((grid_state & mask) != (particle_state & mask)) {
            // Different color **********
            // Directly project instead of writing to the grid

            // calculate impulse on rigid bodies -------------------------------
            if (config_backup.get("compute_particle_impulses", false))
            {
              RigidBody<dim> *r = get_rigid_body_ptr(g.get_rigid_body_id());
              if (r == nullptr)
                continue;
              Vector impulse(0.0_f);
              Vector force_tmp(0.0_f);
              if (p.boundary_distance <= 0.05_f * delta_x){
                Vector rigid_v = r->get_velocity_at(delta_x * grid_pos[node_id]);
                Vector v_acc = v;  // + apic_b_inv_d_mass * dpos / Vector(mass);
                Vector velocity_change =
                  v_acc -
                  friction_project(
                    v_acc,
                    rigid_v,
                    p.boundary_normal,
                    r->frictions[(particle_state >> (2 * r->id)) % 2]);
                impulse = mass * dw_w[dim] * velocity_change
                        + delta_t_tmp_force * Vector(dw_w);
                // added: Force and torque on rigid bodies' center of mass
                force_tmp = impulse / delta_t;
                r->rigid_force_tmp += force_tmp;
                r->rigid_torque_tmp += cross(delta_x * grid_pos[node_id] - r->position, force_tmp);
              }

              // Impulse on rigid body boundary particles
              // TODO: Make it more realistic for visualization only
              if (config_backup.get("visualize_particle_impulses", false))
              {
                for (int r_p_i = particle_begin; r_p_i < particle_end; r_p_i++)
                {
                  Particle &r_p = *allocator[particles[r_p_i]];
                  if (r_p.is_rigid()){
                    // Vector dpos2 = pos - r_p.pos*inv_delta_x;
                    // real w2 = 1.0_f - sqrt(pow(dpos2[0],2)+pow(dpos2[1],2)+pow(dpos2[2],2)) / sqrt(2);
                    r_p.rigid_impulse = force_tmp;
                    // r_p.rigid_impulse += force_tmp * w2;
                  }
                }
              }

              // Apply impulses on rigid body
              if (config_backup.get("affect_particle_impulses", false))
                r->apply_tmp_impulse(impulse, delta_x * grid_pos[13]);
            }

            continue;
          }

          VectorP delta;
          real delta_gf;  // added
#ifdef MLSMPM
          delta = dw_w[dim]*(VectorP(mass_v + apic_b_inv_d_mass * dpos, mass) +
                  VectorP(-delta_t_tmp_force * dpos * 4.0_f * inv_delta_x));
          delta_gf = dw_w[dim] * gf;
#else
          delta = dw_w[dim]* VectorP(mass_v + apic_b_inv_d_mass * dpos, mass) +
                  VectorP(delta_t_tmp_force * Vector(dw_w));
#endif
          g.velocity_and_mass += delta;
          g.granular_fluidity += delta_gf;  // added
        }
      }
    }
  };

  __m128 S = _mm_set1_ps(-4.0_f * inv_delta_x * delta_t);

  // block_op_normal, called from block_op_switch ------------------------------
  auto block_op_normal = [&](uint32 b, uint64 block_offset,
                             GridState<dim> *g_) {
    // using Cache = GridCache<MPM<dim>, true>;
    using Cache = GridCache<MPM<dim>>;  // added
    Cache grid_cache(*grid, block_offset, true);
    int particle_begin;
    int particle_end = block_meta[b].particle_offset;

    // grid loop
    for (uint32 t = 0; t < SparseMask::elements_per_block; t++) {
      particle_begin = particle_end;
      particle_end += g_[t].particle_count;
      int grid_cache_offset = grid_cache.spgrid_block_to_grid_cache_block(t);

      Vectori grid_base_pos = Vectori(SparseMask::LinearToCoord(block_offset)) +
                              grid_cache.spgrid_block_linear_to_vector(t);
      Vector grid_base_pos_f = Vector(grid_base_pos);

      // particle loop
      for (int p_i = particle_begin; p_i < particle_end; p_i++) {
        Particle &p = *allocator[particles[p_i]];
        if (particle_gravity) {
          p.set_velocity(p.get_velocity() + gravity * delta_t);
        }

        // Note, pos is magnified grid pos
        __m128 pos_ = _mm_mul_ps(p.pos.v, _mm_set1_ps(inv_delta_x));

#if defined(MLSMPM)
        MLSMPMFastKernel32 kernel(_mm_sub_ps(pos_, grid_base_pos_f),
                                  inv_delta_x);
        const __m128(&kernels)[3][3] = kernel.kernels;
        using KernelLinearized = real[3 * 3 * 4];
        const KernelLinearized &kernels_linearized =
            *reinterpret_cast<const KernelLinearized *>(&kernels[0][0][0]);
#else
        Kernel kernel(Vector(pos_), inv_delta_x);
        const VectorP(&kernels)[3][3][3] = kernel.kernels;
        using KernelLinearized = VectorP[27];
        const KernelLinearized &kernels_linearized =
            *reinterpret_cast<const KernelLinearized *>(&kernels[0][0][0]);
#endif
        const __m128 v = p.get_velocity().v;
        const real mass = p.get_mass();
        __m128 mass_ = _mm_set1_ps(p.get_mass());
        // Note, apic_b has delta_x issue
        const Matrix apic_b_inv_d_mass = p.apic_b * (Kernel::inv_D() * mass);
        const __m128 mass_v = _mm_mul_ps(_mm_set1_ps(mass), v);

        // added: Disconnection handling
        __m128 gf = _mm_setzero_ps();
        if (p.p > 0.0_f)  // pressure @ n
          gf = _mm_set_ss(p.gf);
        Matrix stress(0.0_f);
        stress = p.calculate_force();

        __m128 delta_t_tmp_force_[3];
        Matrix &delta_t_tmp_force =
            reinterpret_cast<Matrix &>(delta_t_tmp_force_);
        for (int i = 0; i < 3; i++) {
          delta_t_tmp_force_[i] = _mm_mul_ps(_mm_set1_ps(delta_t), stress[i]);
        }

        __m128 rela_pos = _mm_sub_ps(pos_, grid_base_pos_f);
        __m128 affine[3];

        for (int i = 0; i < 3; i++)
          affine[i] = _mm_fmadd_ps(stress[i], S, apic_b_inv_d_mass[i]);

// Loop start
#ifdef MLSMPM
#define LOOP(node_id)                                                          \
  {                                                                            \
    __m128 dpos = _mm_sub_ps(rela_pos, grid_pos_offset_[node_id]);             \
    __m128 g =                                                                 \
        grid_cache                                                             \
            .linear[grid_cache.kernel_linearized(node_id) + grid_cache_offset] \
            .velocity_and_mass;                                                \
    __m128 weight =                                                            \
        _mm_set1_ps(kernels[node_id / 9][node_id / 3 % 3][node_id % 3]);       \
    __m128 affine_prod = _mm_fmadd_ps(                                         \
        affine[2], broadcast(dpos, 2),                                         \
        _mm_fmadd_ps(affine[1], broadcast(dpos, 1),                            \
                     _mm_fmadd_ps(affine[0], broadcast(dpos, 0), mass_v)));    \
    __m128 contrib = _mm_blend_ps(mass_, affine_prod, 0x7);                    \
    __m128 delta = _mm_mul_ps(weight, contrib);                                \
    g = _mm_add_ps(g, delta);                                                  \
    grid_cache                                                                 \
        .linear[grid_cache.kernel_linearized(node_id) + grid_cache_offset]     \
        .velocity_and_mass = g;                                                \
    __m128 delta_gf = _mm_mul_ss(weight, gf);                                  \
    __m128 gg =                                                                \
        _mm_set_ss(grid_cache                                                  \
            .linear[grid_cache.kernel_linearized(node_id) + grid_cache_offset] \
            .granular_fluidity);                                               \
    gg = _mm_add_ss(gg, delta_gf);                                             \
    _mm_store_ss((float *)&grid_cache                                          \
        .linear[grid_cache.kernel_linearized(node_id) + grid_cache_offset]     \
        .granular_fluidity, gg);                                               \
  }
#else
#define LOOP(node_id)                                                          \
  {                                                                            \
    __m128 dpos = _mm_sub_ps(rela_pos, grid_pos_offset_[node_id]);             \
    __m128 g =                                                                 \
        grid_cache                                                             \
            .linear[grid_cache.kernel_linearized(node_id) + grid_cache_offset] \
            .velocity_and_mass;                                                \
    const VectorP &dw_w = kernels_linearized[node_id];                         \
    __m128 delta =                                                             \
        dw_w[dim] *                                                            \
            VectorP(Vector(mass_v) + apic_b_inv_d_mass * Vector(dpos), mass) + \
        VectorP(delta_t_tmp_force * Vector(dw_w));                             \
    g = _mm_add_ps(g, delta);                                                  \
    grid_cache                                                                 \
        .linear[grid_cache.kernel_linearized(node_id) + grid_cache_offset]     \
        .velocity_and_mass = g;                                                \
  }
#endif
        TC_REPEAT27(LOOP);
#undef LOOP
      }
    }
  };

  // block_op_switch -----------------------------------------------------------
  auto block_op_switch = [&](uint32 b, uint64 block_offset, GridState<dim> *g) {
    if (rigid_page_map->Test_Page(block_offset)) {
      block_op_rigid(b, block_offset, g);
    } else {
      block_op_normal(b, block_offset, g);
    }
  };

  // calls block_op_switch
  parallel_for_each_block_with_index(block_op_switch, false, true);
  // apply impulses from particles on rigid bodies
  for (auto &r : rigids) {
    r->apply_tmp_velocity();

    // added: Reset force and torque on rigid bodies' center of mass
    r->rigid_force = r->rigid_force_tmp;
    r->rigid_torque = r->rigid_torque_tmp;
    r->rigid_force_tmp = Vector(0.0_f);
    r->rigid_torque_tmp = Vector(0.0_f);
    
  }
}
// end -------------------------------------------------------------------------

// resample -------------------------------------------------------------- : OFF
static_assert(mpm_kernel_order == 2, "");
template <int dim>
void MPM<dim>::resample() {
  for (auto &r : rigids) {
    r->reset_tmp_velocity();
  }

  // particle
  parallel_for_each_particle([&](Particle &p) {
    if (p.is_rigid()) {
      return;
    }
    real delta_t = base_delta_t;
    Vector v(0.0f), bv(0.0f);
    Matrix b(0.0f);
    Matrix cdg(0.0f);
    Vector pos = p.pos * inv_delta_x;

    RegionND<dim> region(VectorI(0), VectorI(Kernel::kernel_size));
    Vectori grid_base_pos = get_grid_base_pos(pos);
    Kernel kernel(pos, inv_delta_x);

    int rigid_id = -1;

    // grid
    for (auto &ind : region) {
      auto i = ind.get_ipos() + grid_base_pos;

      Vector grid_pos = i.template cast<real>() * delta_x;

      GridState<dim> &g = get_grid(i);

      auto grid_vel = grid_velocity(i);
      Vector dpos = pos - i.template cast<real>();
      VectorP dw_w = kernel.get_dw_w(ind.get_ipos());

      // Coloring
      uint64 grid_state = g.get_states();
      uint64 particle_state = p.states;
      uint64 mask = (grid_state & particle_state & state_mask) >> 1;
      if ((grid_state & mask) != (particle_state & mask)) {
        // different color
        Vector fake_v = p.get_velocity();
        // Vector fake_v(0.0_f);
        RigidBody<dim> *r = get_rigid_body_ptr(g.get_rigid_body_id());
        Vector v_g(0.0_f);
        real friction = 0.0_f;
        if (r != nullptr) {
          v_g = r->get_velocity_at(grid_pos);
          rigid_id = g.get_rigid_body_id();
          friction = r->frictions[(particle_state >> (2 * r->id)) % 2];
        }
        if (p.near_boundary()) {
          if (p.sticky) {
            friction = -1;
          }
          fake_v = friction_project(
                       p.get_velocity(),  // + p.apic_b * dpos * kernel.inv_D(),
                       v_g, p.boundary_normal, friction) +
                   p.boundary_normal * (delta_t * delta_x * pushing_force);
        }
        grid_vel = fake_v;
      }
      v += dw_w[dim] * grid_vel;
      b += Matrix::outer_product(dw_w[dim] * grid_vel, dpos);
      cdg += Matrix::outer_product(grid_vel, Vector(dw_w));
    }
    // end grid

    if (p.near_boundary()) {
      p.apic_b = Matrix(0);
      if (p.sticky) {
        // v *= 0;
      }
    } else {
      p.apic_b = damp_affine_momemtum(b);
    }
    p.set_velocity(v);

#ifdef MLSMPM
    cdg = b * (-4 * inv_delta_x);
#endif
    cdg = Matrix(1.0f) + delta_t * cdg;
    plasticity_counter += p.plasticity(cdg, 0.0f);

    p.pos += delta_t * p.get_velocity();

    // Position correction
    p.pos = (p.pos * inv_delta_x)
                .clamp(Vector(0.0_f), res.template cast<real>() - Vector(eps)) *
            delta_x;
    if (p.near_boundary()) {
      if (p.boundary_distance < -0.05 * delta_x &&
          p.boundary_distance > -delta_x * 0.3) {
        Vector delta_velocity =
            p.boundary_distance * p.boundary_normal * penalty;
        p.set_velocity(p.get_velocity() - delta_velocity);
        if (rigid_id != -1) {
          RigidBody<dim> *r = get_rigid_body_ptr(rigid_id);
          r->apply_tmp_impulse(delta_velocity * p.get_mass(), p.pos);
        }
      }
    }
  });
  for (auto &r : rigids) {
    r->apply_tmp_velocity();
  }
}

template void MPM<2>::rasterize(real delta_t, bool);
template void MPM<2>::resample();
template void MPM<3>::rasterize(real delta_t, bool);
template void MPM<3>::resample();
template <>
void MPM<2>::resample_optimized() {resample();}

// optimized resampling --------------------------------------------------- : ON
template <>
void MPM<3>::resample_optimized() {
  constexpr int dim = 3;

  // block_op_rigid ------------------------------------------------------------
  auto block_op_rigid = [&](uint32 b, uint64 block_offset, GridState<dim> *g) {
    using Cache = GridCache<MPM<dim>>;
    Cache grid_cache(*grid, block_offset, false);
    int particle_begin;
    int particle_end = block_meta[b].particle_offset;

    // element loop
    for (uint32 t = 0; t < SparseMask::elements_per_block; t++) {
      particle_begin = particle_end;

      // number of particles in grid
      particle_end += g[t].particle_count;

      int grid_cache_offset = grid_cache.spgrid_block_to_grid_cache_block(t);

      Vectori grid_base_pos = Vectori(SparseMask::LinearToCoord(block_offset)) +
                              grid_cache.spgrid_block_linear_to_vector(t);
      Vector grid_base_pos_f = Vector(grid_base_pos);

      Vector grid_pos[27];
      for (int i = 0; i < 27; i++) {
        grid_pos[i] = grid_pos_offset[i] + grid_base_pos_f;
      }

      // for each (non-rigid) particle in grid
      for (int k = particle_begin; k < particle_end; k++) {
        Particle &p = *allocator[particles[k]];
        if (p.is_rigid()) {
          continue;
        }
        real delta_t = base_delta_t;
        Vector v(0.0f);
        Matrix b(0.0f);
        Matrix cdg(0.0f);
        Vector pos = p.pos * this->inv_delta_x;

        Vector v_g(0.0_f);
        real friction = 0;
        
        RegionND<dim> region(VectorI(0), VectorI(Kernel::kernel_size));

        Kernel kernel(pos, inv_delta_x);
        const VectorP(&kernels)[3][3][3] = kernel.kernels;
        using KernelLinearized = VectorP[27];
        const KernelLinearized &kernels_linearized =
            *reinterpret_cast<const KernelLinearized *>(&kernels[0][0][0]);

        int rigid_id = -1;

        // added:
        real del_lap_gf = 0.0_f;
        real del_cen_lap_gf = 0.0_f;
        int c = 0;

        // for each node
        for (int node_id = 0; node_id < Cache::num_nodes; node_id++) {
          Vector dpos = pos - grid_pos[node_id];

          GridState<dim> &g =
              grid_cache.linear[grid_cache.kernel_linearized(node_id) +
                                grid_cache_offset];

          auto grid_vel = Vector(g.velocity_and_mass);
          const VectorP &dw_w = kernels_linearized[node_id];

          // Coloring
          uint64 grid_state = g.get_states();
          uint64 particle_state = p.states;
          uint64 mask = (grid_state & particle_state & state_mask) >> 1;

          // added: Modified rb interaction approach
          if (node_id == 13)
          {
            RigidBody<dim> *r = get_rigid_body_ptr(g.get_rigid_body_id());
            if (r != nullptr)
            {
              v_g = r->get_velocity_at(grid_pos[node_id] * delta_x);
              rigid_id = g.get_rigid_body_id();
              friction = r->frictions[(particle_state >> (2 * r->id)) % 2];
            }
          }

          if ((grid_state & mask) != (particle_state & mask)){
             grid_vel = p.get_velocity();
          }else{
            if (node_id == 13)
            {
              del_cen_lap_gf = g.granular_fluidity;
            }else if (node_id == 22 || node_id == 4 || node_id == 16 ||
                      node_id == 10 || node_id == 14 || node_id == 12)
            { 
              del_lap_gf += g.granular_fluidity;
              c += 1;
            }
          }

          v = fused_mul_add(grid_vel, Vector(dw_w[dim]), v); // a*b+c where c=v=0

          // *** b += Matrix::outer_product(dw_w[dim] * grid_vel, dpos);
          //     cdg += Matrix::outer_product(grid_vel, Vector(dw_w));
          Vector w_grid_vel = dw_w[dim] * grid_vel;
          for (int r = 0; r < dim; r++) {
            b[r]   = fused_mul_add(w_grid_vel, Vector(dpos[r]), b[r]);
            cdg[r] = fused_mul_add(grid_vel, Vector(dw_w[r]), cdg[r]);
          }

        }

        if (p.near_boundary()) {
          p.apic_b = Matrix(0);
        } else {
          p.apic_b = damp_affine_momemtum(b);
        }

        // set non-rigid particle velocity (v) 
        p.set_velocity(v);
        Vector delta_t_vec(delta_t);

        // added: Modified rb interaction approach
        real dist = 0.05_f * delta_x;
        if (p.near_boundary() && p.boundary_distance <= dist)
        {
          Vector v_collision(0.0_f);
          v_collision =
            friction_project(
              p.get_velocity(),  // + p.apic_b * dpos * kernel.inv_D(),
              v_g,  // grid node (in rigidbody) vel
              p.boundary_normal, friction) - v_g
            + (dist-p.boundary_distance)/dist * v_g
            + p.boundary_normal * (delta_t * delta_x * pushing_force);  // pushing force
          p.set_velocity(v_collision);
        }

        // *** cdg = Matrix(1.0f) + delta_t * cdg;
#ifdef MLSMPM
        cdg = b * (-4 * inv_delta_x);
#endif
        for (int i = 0; i < dim; i++){
          cdg[i] = fused_mul_add(delta_t_vec, cdg[i], Vector::axis(i));
        }

        // added
        // laplacian of granular fluidity using central FD scheme
        real laplacian_gf = 0.0_f;
        if (c == 6){
          // laplacian_gf = this->inv_delta_x * this->inv_delta_x * (del_lap_gf - del_cen_lap_gf*c);
          laplacian_gf = inv_delta_x * inv_delta_x * (
            + grid_cache.linear[grid_cache.kernel_linearized(22) + grid_cache_offset].granular_fluidity
            + grid_cache.linear[grid_cache.kernel_linearized(4)  + grid_cache_offset].granular_fluidity
            + grid_cache.linear[grid_cache.kernel_linearized(16) + grid_cache_offset].granular_fluidity
            + grid_cache.linear[grid_cache.kernel_linearized(10) + grid_cache_offset].granular_fluidity
            + grid_cache.linear[grid_cache.kernel_linearized(14) + grid_cache_offset].granular_fluidity
            + grid_cache.linear[grid_cache.kernel_linearized(12) + grid_cache_offset].granular_fluidity
            -(grid_cache.linear[grid_cache.kernel_linearized(13) + grid_cache_offset].granular_fluidity * 6.0_f));
        }

        // added: Update granular fluidity and deformation gradient
        plasticity_counter += p.plasticity(cdg, laplacian_gf);

        p.pos = fused_mul_add(p.get_velocity(), delta_t_vec, p.pos);

        // Position correction: apply penalty on boundary particle
        if (p.near_boundary()) {
          if (p.boundary_distance < -0.05 * delta_x &&
              p.boundary_distance > -delta_x * 0.3) {
            Vector delta_velocity =
                p.boundary_distance * p.boundary_normal * penalty;
            p.set_velocity(p.get_velocity() - delta_velocity);
            if (rigid_id != -1) {
              RigidBody<dim> *r = get_rigid_body_ptr(rigid_id);
              r->apply_tmp_impulse(delta_velocity * p.get_mass(), p.pos);
            }
          }
        } 

      }  // particle loop end
    }
  };

  auto block_op_normal = [&](uint32 b, uint64 block_offset, GridState<dim> *g) {
    // using Cache = GridCache<MPM<dim>, true>;
    using Cache = GridCache<MPM<dim>>;  // added
    Cache grid_cache(*grid, block_offset, false);
    int particle_begin;
    int particle_end = block_meta[b].particle_offset;

    real inv_delta_x = this->inv_delta_x;

    // grid loop
    // elements_per_block = 8 x 8 x 4 = 256
    for (uint32 t = 0; t < SparseMask::elements_per_block; t++) {
      particle_begin = particle_end;
      particle_end += g[t].particle_count;

      int grid_cache_offset = grid_cache.spgrid_block_to_grid_cache_block(t);

      Vectori grid_base_pos = Vectori(SparseMask::LinearToCoord(block_offset)) +
                              grid_cache.spgrid_block_linear_to_vector(t);
      Vector grid_base_pos_f = Vector(grid_base_pos);

      // added
      // it's the same for each element
      // laplacian of granular fluidity using central FD scheme
      real laplacian_gf = inv_delta_x * inv_delta_x * (
          + grid_cache.linear[grid_cache.kernel_linearized(22) + grid_cache_offset].granular_fluidity
          + grid_cache.linear[grid_cache.kernel_linearized(4)  + grid_cache_offset].granular_fluidity
          + grid_cache.linear[grid_cache.kernel_linearized(16) + grid_cache_offset].granular_fluidity
          + grid_cache.linear[grid_cache.kernel_linearized(10) + grid_cache_offset].granular_fluidity
          + grid_cache.linear[grid_cache.kernel_linearized(14) + grid_cache_offset].granular_fluidity
          + grid_cache.linear[grid_cache.kernel_linearized(12) + grid_cache_offset].granular_fluidity
          -(grid_cache.linear[grid_cache.kernel_linearized(13) + grid_cache_offset].granular_fluidity * 6.0_f)
        );

      // particle loop
      for (int k = particle_begin; k < particle_end; k++) {
        Particle &p = *allocator[particles[k]];
        real delta_t = base_delta_t;
        Vector pos = p.pos * inv_delta_x;

#if defined(MLSMPM)
        MLSMPMFastKernel32 kernel(pos - grid_base_pos_f, inv_delta_x);
        const __m128(&kernels)[3][3] = kernel.kernels;
        using KernelLinearized = real[3 * 3 * 4];
        const KernelLinearized &kernels_linearized =
            *reinterpret_cast<const KernelLinearized *>(&kernels[0][0][0]);
#else
        Kernel kernel(pos, inv_delta_x);
        const VectorP(&kernels)[3][3][3] = kernel.kernels;
        using KernelLinearized = VectorP[27];
        const KernelLinearized &kernels_linearized =
            *reinterpret_cast<const KernelLinearized *>(&kernels[0][0][0]);
#endif

        /// SIMD Def Region
        __m128 b_[3];
        __m128 cdg_[3];
        for (int i = 0; i < dim; i++) {
          b_[i] = _mm_setzero_ps();
          cdg_[i] = _mm_setzero_ps();
        }
        __m128 v_ = _mm_setzero_ps();
        __m128 pos_ = pos.v;
        __m128 rela_pos = _mm_sub_ps(pos_, grid_base_pos_f);
///

#if defined(MLSMPM)
#define LOOP(node_id)                                                          \
  {                                                                            \
    __m128 dpos = _mm_sub_ps(rela_pos, grid_pos_offset_[node_id]);             \
    float *addr =                                                              \
        (float *)&grid_cache                                                   \
            .linear[grid_cache.kernel_linearized(node_id) + grid_cache_offset] \
            .velocity_and_mass;                                                \
    __m128 grid_vel = _mm_load_ps(addr);                                       \
    __m128 w =                                                                 \
        _mm_set1_ps(kernels[node_id / 9][node_id / 3 % 3][node_id % 3]);       \
    v_ = _mm_fmadd_ps(grid_vel, w, v_);                                        \
    __m128 w_grid_vel = _mm_mul_ps(w, grid_vel);                               \
    for (int r = 0; r < dim; r++) {                                            \
      b_[r] = _mm_fmadd_ps(w_grid_vel, broadcast(dpos, r), b_[r]);             \
    }                                                                          \
  }
#else
#define LOOP(node_id)                                                          \
  {                                                                            \
    __m128 dpos = _mm_sub_ps(rela_pos, grid_pos_offset_[node_id]);             \
    float *addr =                                                              \
        (float *)&grid_cache                                                   \
            .linear[grid_cache.kernel_linearized(node_id) + grid_cache_offset] \
            .velocity_and_mass;                                                \
    __m128 grid_vel = _mm_load_ps(addr);                                       \
    __m128 dw_w = kernels_linearized[node_id].v;                               \
    v_ = _mm_fmadd_ps(grid_vel, broadcast(dw_w, dim), v_);                     \
    __m128 w_grid_vel = _mm_mul_ps(broadcast(dw_w, dim), grid_vel);            \
    for (int r = 0; r < dim; r++) {                                            \
      b_[r] = _mm_fmadd_ps(w_grid_vel, broadcast(dpos, r), b_[r]);             \
      cdg_[r] = _mm_fmadd_ps(grid_vel, broadcast(dw_w, r), cdg_[r]);           \
    }                                                                          \
  }
#endif
        TC_REPEAT27(LOOP);
#undef LOOP
        if (rpic_damping != 0 && apic_damping != 0) {
          p.apic_b = damp_affine_momemtum(b);
        } else {
          for (int i = 0; i < dim; i++) {
            p.apic_b[i].v = b_[i];
          }
        }

        p.set_velocity(Vector(v_));

        __m128 delta_t_vec = _mm_set1_ps(delta_t);
#ifdef MLSMPM
        // cdg = -b * 4 * inv_delta_x;
        __m128 scale = _mm_set1_ps(-4.0_f * inv_delta_x * delta_t);

        cdg_[0] = _mm_fmadd_ps(scale, b_[0], _mm_set_ps(0, 0, 0, 1));
        cdg_[1] = _mm_fmadd_ps(scale, b_[1], _mm_set_ps(0, 0, 1, 0));
        cdg_[2] = _mm_fmadd_ps(scale, b_[2], _mm_set_ps(0, 1, 0, 0));
#else
        // cdg = Matrix(1.0f) + delta_t * cdg;
        cdg_[0] = _mm_fmadd_ps(delta_t_vec, cdg_[0], _mm_set_ps(0, 0, 0, 1));
        cdg_[1] = _mm_fmadd_ps(delta_t_vec, cdg_[1], _mm_set_ps(0, 0, 1, 0));
        cdg_[2] = _mm_fmadd_ps(delta_t_vec, cdg_[2], _mm_set_ps(0, 1, 0, 0));
#endif
        Matrix &cdg = reinterpret_cast<Matrix &>(cdg_[0]);

        // added: Update granular fluidity and deformation gradient
        p.plasticity(cdg, laplacian_gf);

        // advect particles
        p.pos.v = _mm_fmadd_ps(v_, delta_t_vec, p.pos.v);

      }
    }
  };

  for (auto &r : rigids) {
    r->reset_tmp_velocity();
  }

  auto block_op_switch = [&](uint32 b, uint64 block_offset, GridState<dim> *g) {
    if (rigid_page_map->Test_Page(block_offset)) {
      block_op_rigid(b, block_offset, g);
    } else {
      block_op_normal(b, block_offset, g);
    }
  };

  parallel_for_each_block_with_index(block_op_switch, false, false);

  for (auto &r : rigids) {
    r->apply_tmp_velocity();
  }
}
// optimized resampling end ----------------------------------------------------

template void MPM<2>::resample_optimized();
template void MPM<3>::resample_optimized();
TC_TEST("mls_kernel") {
  for (int t = 0; t < 10000; t++) {
    Vector3 pos = Vector3::rand() + Vector3(0.5_f);
    MPMFastKernel32 gt(pos, 1);
    MLSMPMFastKernel32 fast(pos, 1);
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
          TC_CHECK_EQUAL(gt.get_w(Vector3i(i, j, k)), fast.kernels[i][j][k],
                         1e-6_f);
        }
      }
    }
  }
}

TC_NAMESPACE_END
#endif
