/*******************************************************************************
    Copyright (c) The Taichi MPM Authors (2018- ). All Rights Reserved.
    The use of this software is governed by the LICENSE file.
*******************************************************************************/

#include <taichi/system/profiler.h>
#include <taichi/dynamics/rigid_body.h>
#include "mpm.h"
#include "kernel.h"
#include "boundary_particle.h"

#ifndef MPM_TRANSFER_OPT
// Non-optimized Transfer

TC_NAMESPACE_BEGIN

// in this section (rP2G), grid state has issue with corners *******************
// rP2G ------------------------------------------------------------------------
// writes to grid_state, grid_cdf, grid_sdf_last_rigid_particle, rigid_markers
template <int dim>
void MPM<dim>::rasterize_rigid_boundary() {
  
  // find nearest rigid particle to grid node ----------------------------------
  // Note: we assume particles from the same rigid body are contiguous (adjacent)
  parallel_for_each_particle([&](Particle &p_) {
    // ONLY rigid particles
    auto *p = dynamic_cast<RigidBoundaryParticle<dim> *>(&p_);
    if (!p)
      return;

    // cdf_kernel_order_rasterize <- 2 (from mpm_fwd.h)
    constexpr int kernel_size = cdf_kernel_order_rasterize + 1;  // was 1
    // from 0 to 3 (grid) in each direction
    RegionND<dim> region(VectorI(0), VectorI(kernel_size));
    // grid base pos from rigid particle pos
    Vectori grid_base_pos =
        get_grid_base_pos_with<kernel_size>(p->pos * inv_delta_x);

    // triangular (3D) element where particle is in
    auto elem = p->get_world_space_element();
    Matrix world_to_elem = world_to_element(elem);

    // first vertex of element (counter-clockwise mesh)
    Vector p0 = elem.v[0];

    // grid node loop (ind from 0 to 3 in each direction)
    for (auto &ind : region) {

      // grid node pos
      // if ind=[1-0-1] then ind.get_ipos()=[1,0,1]
      auto i = ind.get_ipos() + grid_base_pos;
      auto grid_pos = i.template cast<real>() * delta_x;;

      // transformed dist from [grid node] to [first element's vertex]
      Vector coord = world_to_elem * (grid_pos - p0);

      // if grid node is inside the rigid body then negative 
      bool negative = coord[dim - 1] < 0;

      // abs of this dist in z direction
      real dist_triangle = abs(coord[dim - 1]);

      if (config_backup.get("cdf_3d_modified", false))
      {
        // added: CPIC modification for corners, less penetration (leakage)
        bool valid = false;
        TC_STATIC_IF(dim == 3)
        {    
          // Option 1 for generic configuration with O(1):
          Vector grid_pos_elem = world_to_elem * grid_pos;
          Vector e0 = world_to_elem * elem.v[0];
          Vector e1 = world_to_elem * elem.v[1];
          Vector e2 = world_to_elem * elem.v[2];
          real A = abs(e0[0]*(e1[1]-e2[1]) +
                        e1[0]*(e2[1]-e0[1]) +
                        e2[0]*(e0[1]-e1[1]))/2.0_f;
          real A0 = abs(grid_pos_elem[0]*(e1[1]-e2[1]) +
                        e1[0]*(e2[1]-grid_pos_elem[1]) +
                        e2[0]*(grid_pos_elem[1]-e1[1]))/2.0_f;
          real A1 = abs(e0[0]*(grid_pos_elem[1]-e2[1]) +
                        grid_pos_elem[0]*(e2[1]-e0[1]) +
                        e2[0]*(e0[1]-grid_pos_elem[1]))/2.0_f;
          real A2 = abs(e0[0]*(e1[1]-grid_pos_elem[1]) +
                        e1[0]*(grid_pos_elem[1]-e0[1]) +
                        grid_pos_elem[0]*(e0[1]-e1[1]))/2.0_f;
            
          Vector u_world = elem.v[1] - elem.v[0];
          Vector v_world = elem.v[2] - elem.v[0];
          Vector n_world = cross(u_world, v_world); 
          Vector dist_world = grid_pos - elem.v[0];
          real dot_prod_world = 0.0_f;
          for (int i = 0; i < 2; i++)
            dot_prod_world += dist_world[i] * n_world[i];

          if (abs(A-(A0+A1+A2)) <= 1e-5_f && abs(dot_prod_world) <= 1e-5_f)
            valid = true;

          // Option 2 for box only (works well):
          // if (negative && dist_triangle <= delta_x)
          // {
          //   Vector3 u = elem.v[1] - elem.v[0];
          //   Vector3 v = elem.v[2] - elem.v[0];
          //   Vector3 n = cross(u, v);
          //   if (n[1] != 0.0_f || n[2] != 0.0_f)
          //     valid = false;
          // }
        }
        TC_STATIC_END_IF
        if (!valid)
          continue;
      }
      else
      {
        bool in_range = false;
        // in-range grid node
        TC_STATIC_IF(dim == 2) {
          if (-0.02 <= coord[0] && coord[0] <= 1.02)
            in_range = true;
        }
        TC_STATIC_ELSE {
          if (0 <= coord[0] && 0 <= coord[1] && coord[0]+coord[1] <= 1)
            in_range = true;
        }
        TC_STATIC_END_IF
        if (!in_range)
          continue;
      }

      dist_triangle *= inv_delta_x;

      GridState<dim> &g = get_grid(i);
      g.get_lock().lock();

      // set minimum dist to the grid in rigid particle's kernel
      if (g.get_rigid_body_id() == -1 || dist_triangle < get_grid(i).get_distance()) {
        g.set_distance(dist_triangle);
        // g.set_distance(5.0_f);
        g.set_rigid_body_id(p->rigid->id);
      }

      // TODO: what happens if the relative position of the grid point to a
      // rigid body is ambiguous???
      // Should be fine for most meshes, though?
      g.set_states(g.get_states() | (2 + (int)(negative)) << (p->rigid->id * 2));

      g.get_lock().unlock();
    }
  });

  // scale the distance
  parallel_for_each_active_grid([&](GridState<dim> &g){
    g.set_distance(g.get_distance() * delta_x);
  });

  // extra steps in 2D 
  TC_STATIC_IF(dim == 2) {
    ArrayND<2, uint32> grid_states_tmp(id(this->res + Vectori(1)), 0);
    for (int e = 0; e < config_backup.get<int>("cdf_expand", 0); e++) {
      for (int k = 0; k < dim; k++) {
        Region region = Region(Vectori::axis(k), res - Vectori::axis(k));
        grid_states_tmp.reset_zero();
        for (auto &ind : region) {
          id(grid_states_tmp)[ind] = this->get_grid(ind.get_ipos()).states;
        }
        for (auto &ind : region) {
          auto update = [&](Vectori offset) {
            auto nei_state = this->get_grid(ind + offset).states;
            auto state = id(grid_states_tmp)[ind];
            auto states_to_add =
                ((nei_state & ~state) & state_mask) & GridState<dim>::tag_mask;
            id(grid_states_tmp)[ind] =
                state | (nei_state & (states_to_add | (states_to_add >> 1)));
          };
          update(Vectori::axis(k));
          update(-Vectori::axis(k));
        }
        for (auto &ind : region) {
          this->get_grid(ind.get_ipos()).states = id(grid_states_tmp)[ind];
        }
      }
      int count = 0;
      Region region = Region(Vectori(0), this->res);
      for (auto &ind : region) {
        count += (0 != this->get_grid(ind.get_ipos()).states);
      }
      // TC_P(count);
    }
  }  

  // added: ???
  TC_STATIC_ELSE {
    if (config_backup.get("cdf_modified_cube", false)){
      // Option 1:
      // ArrayND<3, uint32> grid_states_tmp(id(this->res + Vectori(1)), 0);
      // for (int k = 0; k < dim; k++) {
      //   Region region = Region(
      //     Vectori(int(minXVerPos*inv_delta_x)-1,int(minYVerPos*inv_delta_x)-1,int(minZVerPos*inv_delta_x)-1),
      //     Vectori(int(maxXVerPos*inv_delta_x)+1,int(maxYVerPos*inv_delta_x)+1,int(maxZVerPos*inv_delta_x)+1)
      //   );
      //   grid_states_tmp.reset_zero();
      //   for (auto &ind : region) {
      //     id(grid_states_tmp)[ind] = this->get_grid(ind.get_ipos()).states;
      //   }
      //   for (auto &ind : region) {
      //     auto update = [&](Vectori offset) {
      //       auto nei_state = this->get_grid(ind + offset).states;
      //       auto state = id(grid_states_tmp)[ind];
      //       auto states_to_add =
      //           ((nei_state & ~state) & state_mask) & GridState<dim>::tag_mask;
      //       id(grid_states_tmp)[ind] =
      //           state | (nei_state & (states_to_add | (states_to_add >> 1)));
      //     };
      //     update(Vectori::axis(k));
      //     update(-Vectori::axis(k));
      //   }
      //   for (auto &ind : region) {
      //     this->get_grid(ind.get_ipos()).states = id(grid_states_tmp)[ind];
      //   }
      // }

      // Option 2:
      // Region region = Region(Vectori(int(minXVerPos*inv_delta_x)-1,int(minYVerPos*inv_delta_x)-1,int(minZVerPos*inv_delta_x)-1),
      //                        Vectori(int(maxXVerPos*inv_delta_x)+1,int(maxYVerPos*inv_delta_x)+1,int(maxZVerPos*inv_delta_x)+1));
      // constexpr int kernel_size = cdf_kernel_order_rasterize + 1;
      // Vectori min_grid_base_pos = get_grid_base_pos_with<kernel_size>(minYParPos);
      // VectorI dx1(1/this->res);
      // GridState<dim> &g_aux = get_grid(min_grid_base_pos+Vectori(0,3,0));
      // // TC_P(int(minYVerPos*inv_delta_x)-1);
      // // TC_P(minYParPos);
      // for (auto &ind : region){
      //   auto i = ind.get_ipos();
      //   if ( (i[1] >= min_grid_base_pos[1] && i[1] <= min_grid_base_pos[1]+dx1[1] ) &&
      //       this->get_grid(i).states == 0){
      //     // (minXVerPos <= i[0] && i[0] <= minXVerPos+) &&   
      //     // TC_P(i);
      //     if (g_aux.get_distance() != 0){
      //       for (int c = 0; c <= 2; c++){
      //         GridState<dim> &g = get_grid(i+Vectori(0,c,0));
      //         g.set_distance(g_aux.get_distance());
      //         g.set_rigid_body_id(rbID);
      //         g.get_lock().lock();
      //         g.set_states(g.get_states() | (2 + (int)(false)) << (rbID*2));
      //         g.get_lock().unlock();
      //       }
      //     }
      //   }
      // }
    }

  }

  TC_STATIC_END_IF
}

template void MPM<2>::rasterize_rigid_boundary();

template void MPM<3>::rasterize_rigid_boundary();

// G2P --------------------------------------------------------------------------
// Construct particle states
template <int dim>
void MPM<dim>::gather_cdf() {
  // Gain new color
  // do not modify existing colors
  /*
    p.states = (particle_state | (new_color_mask << 1) |
                //((grid_state & new_color_mask) * bit));
                (new_color_mask * bit));
                */

  // non-rigid particle loop
  parallel_for_each_particle([&](Particle &p) {
    p.sticky = false;

    // ignore rigid particles
    if (p.is_rigid()) {
      return;
    }

    using CDFPrecision = float64;  // was 32

    using VectorPd = VectorND<dim + 1, CDFPrecision>;
    using MatrixPd = MatrixND<dim + 1, CDFPrecision>;
    using Vectord = VectorND<dim, CDFPrecision>;

    p.boundary_distance = 0;
    p.boundary_normal = Vector(0.0_f);
    p.near_boundary_ = false;

    // particle pos
    auto pos = p.pos * inv_delta_x;
    auto linearized_offset =
        SparseMask::Linear_Offset(pos.template cast<int>());
    if (dim == 3 && !rigid_page_map->Test_Page(linearized_offset)) {
      return;
    }

    // particle state
    AffinityType &particle_states = p.states;

    // 3
    constexpr int kernel_size = cdf_kernel_order_gather + 1;

    // 0 to 3 in each direction
    RegionND<dim> region(VectorI(0), VectorI(kernel_size));

    // grid base pos
    Vectori grid_base_pos = get_grid_base_pos_with<kernel_size>(pos);

    // 2
    MPMKernel<dim, kernel_size - 1> kernel(pos, inv_delta_x);

    // if there is only one one-state grid node in kernel
    AffinityType all_boundaries = 0;
    for (auto &ind : region) {
      // grid node pos (in kernel)
      auto i = ind.get_ipos() + grid_base_pos;
      all_boundaries |= (get_grid(i).get_states() & state_mask);
    }

    // Unset states that the particle does not touch
    particle_states &= (all_boundaries + (all_boundaries >> 1));

    // change zero-state particles
    AffinityType all_states_to_add = all_boundaries & (~particle_states);

    // State 1 ---------------------
    while (all_states_to_add != 0) {
      //
      AffinityType state_to_add = (all_states_to_add & -all_states_to_add);

      //
      all_states_to_add ^= state_to_add;

      //
      CDFPrecision weighted_distances[2] = {0, 0};
      for (auto &ind : region) {
        // grid node pos
        auto i = ind.get_ipos() + grid_base_pos;

        // distance from particle to grid node
        Vectord dpos =
            (pos - i.template cast<real>()).template cast<CDFPrecision>();

        // weight func
        VectorP dw_w = kernel.get_dw_w(ind.get_ipos());

        // Coloring
        GridState<dim> &g = get_grid(i);
        uint64 grid_state = g.get_states();

        if (g.get_rigid_body_id() == -1) {
          continue;
        }

        // distance from grid node to the nearest rigid body particle
        CDFPrecision d = g.get_distance() * inv_delta_x;
        VectorPd xp(-dpos, 1);

        CDFPrecision weight;
        if (mpm_use_weighted_reconstruction) {
          weight = dw_w[dim];
        } else {
          weight = 1.0f;
        }

        // assert((state_to_add & (state_to_add - 1)) == 0);
        if (grid_state & state_to_add) {
          // determine if the distance is from left or right
          int sign = (int)((grid_state & (state_to_add >> 1)) != 0);
          weighted_distances[sign] += d * weight;
        }
      }

      //
      if (weighted_distances[0] + weighted_distances[1] > 1e-7_f) {
        p.states |=
            (state_to_add | ((state_to_add >> 1) * int(weighted_distances[0] <
                                                       weighted_distances[1])));
        cutting_counter++;
      }
    }

    // State 2 ----------------
    if (particle_states != 0) {
      MatrixPd XtX(0);
      VectorPd XtY(0);
      for (auto &ind : region) {
        auto i = ind.get_ipos() + grid_base_pos;
        Vectord dpos =
            (pos - i.template cast<real>()).template cast<CDFPrecision>();

        VectorP dw_w = kernel.get_dw_w(ind.get_ipos());

        GridState<dim> &g = get_grid(i);

        if (g.get_rigid_body_id() == -1) {
          continue;
        }

        // Coloring
        // calculated from rP2G
        uint64 grid_state = g.get_states();

        //
        uint64 mask = (grid_state & particle_states & state_mask) >> 1;

        // calculated from rP2G
        // distance from grid node to the nearest rigid body particle
        CDFPrecision d = g.get_distance() * inv_delta_x;

        VectorPd xp(-dpos, 1);
        CDFPrecision weight;
        if (mpm_use_weighted_reconstruction) {
          weight = dw_w[dim];
        } else {
          weight = 1.0f;
        }

        // ?
        if (grid_state != 0) {
          // if compatible
          if ((grid_state & mask) == (particle_states & mask) &&
              (grid_state != 0)) {
            // same color
            // M in equ (4)
            XtX += MatrixPd::outer_product(xp, xp) * weight;
            // the rest in equ (4)
            XtY += VectorPd(-d * dpos, d) * weight;
            // if (i[1]<=31)
            //   TC_P(i);
          // true 
          } else if (cdf_use_negative) {
            // Only one color different, use negative
            uint64 diff = ((grid_state & mask) ^ (particle_states & mask));
            if (diff > 0 && 0 == (diff & (diff - 1))) {
              XtX += MatrixPd::outer_product(xp, xp) * weight;
              XtY += VectorPd(d * dpos, -d) * weight;
            }
          }
        }
      }
      if (std::abs(determinant(XtX)) > mpm_reconstruction_guard<dim>()) {
        // gradient of d
        VectorP r = (inversed(XtX) * XtY).template cast<real>();
        p.near_boundary_ = true;
        // r[dim] -= 1 * (1 - length2(Vector(r)));
        p.boundary_distance = r[dim] * delta_x;
        p.debug = length2(Vector(r));
        if (length2(Vector(r)) > 1e-4_f) {
          p.boundary_normal = normalized(Vector(r));
        } else {
          p.boundary_normal = Vector(0);

          // p.sticky = true;
          
          // TC_WARN("0 boundary normal detected");
        }
      } else {
        // p.states = 0;
        p.boundary_distance = 0;
        p.boundary_normal = Vector(0);
        // p.near_boundary_ = false;
      }
    }

  });
}

template void MPM<2>::gather_cdf();

template void MPM<3>::gather_cdf();

TC_NAMESPACE_END

#endif
