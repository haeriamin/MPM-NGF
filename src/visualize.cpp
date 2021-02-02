/*******************************************************************************
    Copyright (c) The Taichi MPM Authors (2018- ). All Rights Reserved.
    The use of this software is governed by the LICENSE file.
*******************************************************************************/

#include <taichi/system/threading.h>
#include <taichi/visual/texture.h>
#include <taichi/system/profiler.h>
#include <taichi/visual/scene.h>
#include <Partio.h>

#include "mpm.h"
#include "kernel.h"  // added

TC_NAMESPACE_BEGIN

// write_partio
template <int dim>
void MPM<dim>::write_partio(const std::string &file_name) const {
  Partio::ParticlesDataMutable *parts = Partio::create();
  Partio::ParticleAttribute posH, vH, mH, typeH, normH, statH, boundH, distH,
      debugH, indexH, limitH, apicH,
      gfH, prH, muH, rigid_impulseH, isFreeH; // added

  bool verbose = config_backup.get("verbose_bgeo", false);

  posH = parts->addAttribute("position", Partio::VECTOR, 3);
  typeH = parts->addAttribute("type", Partio::INT, 1);
  indexH = parts->addAttribute("index", Partio::INT, 1);
  limitH = parts->addAttribute("limit", Partio::INT, 3);
  vH = parts->addAttribute("v", Partio::VECTOR, 3);

  // added:
  gfH = parts->addAttribute("gf", Partio::FLOAT, 1); // granular fluidity
  prH = parts->addAttribute("pr", Partio::FLOAT, 1); // pressure
  muH = parts->addAttribute("mu", Partio::FLOAT, 1); // friction coeff
  rigid_impulseH = parts->addAttribute("rigid_impulse", Partio::VECTOR, 3);
  isFreeH = parts->addAttribute("free", Partio::INT, 1);
  // distH = parts->addAttribute("boundary_distance", Partio::FLOAT, 1);

  if (verbose) {
    mH = parts->addAttribute("m", Partio::VECTOR, 1);
    normH = parts->addAttribute("boundary_normal", Partio::VECTOR, 3);
    debugH = parts->addAttribute("debug", Partio::VECTOR, 3);
    statH = parts->addAttribute("states", Partio::INT, 1);
    distH = parts->addAttribute("boundary_distance", Partio::FLOAT, 1);
    boundH = parts->addAttribute("near_boundary", Partio::INT, 1);
    apicH = parts->addAttribute("apic_frobenius_norm", Partio::FLOAT, 1);
  }
  auto particles_sorted = particles;
  std::sort(particles_sorted.begin(), particles_sorted.end(),
            [&](ParticlePtr a, ParticlePtr b) {
              return allocator.get_const(a)->id < allocator.get_const(b)->id;
            });
  for (auto p_i : particles_sorted) {
    const Particle *p = allocator.get_const(p_i);
    int idx = parts->addParticle();
    if (verbose) {
      float32 *m_p = parts->dataWrite<float32>(mH, idx);
      float32 *bn_p = parts->dataWrite<float32>(normH, idx);
      float32 *debug_p = parts->dataWrite<float32>(debugH, idx);
      int *st_p = parts->dataWrite<int>(statH, idx);
      int *nb_p = parts->dataWrite<int>(boundH, idx);
      float32 *dist_p = parts->dataWrite<float32>(distH, idx);
      float32 *apic_p = parts->dataWrite<float32>(apicH, idx);
      real mass = p->get_mass();
      Vector norm = p->boundary_normal;
      m_p[0] = mass;
      for (int k = 0; k < 3; k++)
        bn_p[k] = 0.f;
      for (int k = 0; k < 3; k++)
        debug_p[k] = 0.f;
      for (int k = 0; k < dim; k++)
        bn_p[k] = norm[k];
      for (int k = 0; k < 3; k++)
        debug_p[k] = p->get_debug_info()[k];

      st_p[0] = p->states;
      nb_p[0] = p->near_boundary();
      dist_p[0] = p->boundary_distance * inv_delta_x;
      apic_p[0] =
          (0.5_f * (p->apic_b - p->apic_b.transposed())).frobenius_norm();
    }

    Vector vel = p->get_velocity();
    float32 *v_p = parts->dataWrite<float32>(vH, idx);
    for (int k = 0; k < 3; k++)
      v_p[k] = 0.f;
    for (int k = 0; k < dim; k++)
      v_p[k] = vel[k];
    int *type_p = parts->dataWrite<int>(typeH, idx);
    int *index_p = parts->dataWrite<int>(indexH, idx);
    int *limit_p = parts->dataWrite<int>(limitH, idx);
    float32 *p_p = parts->dataWrite<float32>(posH, idx);

    // added:
    float32 *gf_p = parts->dataWrite<float32>(gfH, idx);  // granular fluidity
    float32 *pr_p = parts->dataWrite<float32>(prH, idx);  // pressure
    float32 *mu_p = parts->dataWrite<float32>(muH, idx);  // friction coeff
    float32 *rigid_impulse_p = parts->dataWrite<float32>(rigid_impulseH, idx);
    int *isFree_p = parts->dataWrite<int>(isFreeH, idx);
    // float32 *dist_p = parts->dataWrite<float32>(distH, idx);

    Vector pos = p->pos;

    for (int k = 0; k < 3; k++)
      p_p[k] = 0.f;

    for (int k = 0; k < dim; k++)
      p_p[k] = pos[k];
    type_p[0] = int(p->is_rigid());
    index_p[0] = p->id;
    limit_p[0] = p->dt_limit;
    limit_p[1] = p->stiffness_limit;
    limit_p[2] = p->cfl_limit;

    // added:
    gf_p[0] = p->gf;  // granular fluidity
    pr_p[0] = p->p;  // pressure
    // mu_p[0] = p->mu_visual;  // friction coeff
    Vector rf = p->rigid_impulse;
    for (int k = 0; k < dim; k++)
      rigid_impulse_p[k] = rf[k];
    // isFree_p[0] = int(p->is_free);
    // dist_p[0] = p->boundary_distance * inv_delta_x;
  }
  Partio::write(file_name.c_str(), *parts);
  parts->release();
}

// write_rigid_body
template <int dim>
void MPM<dim>::write_rigid_body(RigidBody<dim> const *rigid,
                                const std::string &file_name) const {
  TC_STATIC_IF(dim == 2) {
    auto &mesh = rigid->mesh;
    auto const &trans = rigid->get_mesh_to_world();

    std::string fn = file_name + std::string(".poly");
    FILE *f = fopen(fn.c_str(), "w");

    fmt::print(f, "POINTS\n");

    int index = 0;
    for (auto &elem : mesh->elements) {
      for (int i = 0; i < dim; i++) {
        Vector v = transform(trans, elem.v[i]);
        fmt::print(f, "{}: {} {} 0.0\n", ++index, v.x, v.y);
      }
    }

    fmt::print(f, "POLYS\n");

    for (int i = 1; i <= index / 2; ++i) {
      fmt::print(f, "{}: {} {}\n", i, i * 2 - 1, i * 2);
    }

    fmt::print(f, "END\n");

    fclose(f);
  }
  // 3D
  TC_STATIC_ELSE {
    auto &mesh = rigid->mesh;
    auto const &trans = rigid->get_mesh_to_world();

    std::string fn = file_name + std::string(".obj");

    FILE *f = std::fopen(fn.c_str(), "w");

    for (auto &elem : mesh->elements) {
      for (int i = 0; i < dim; i++) {
        Vector v = transform(trans, elem.v[i]);
        fmt::print(f, "v {} {} {}\n", v[0], v[1], v[2]);
      }
    }
    int counter = 0;
    for (auto &_ : mesh->elements) {
      fmt::print(f, "f {} {} {}\n", dim * counter + 1, dim * counter + 2,
                 dim * counter + 3);
      counter++;
      trash(_);
    }
    std::fclose(f);

    // added
    std::string gn = file_name + std::string(".txt");
    FILE *g = std::fopen(gn.c_str(), "w");
    Vector force = rigid->rigid_force;
    Vector torque = rigid->rigid_torque;
    Vector position = rigid->position;
    fmt::print(g, "{} {} {} {} {} {} {} {} {}\n",
      force[0], force[1], force[2],
      torque[0], torque[1], torque[2],
      position[0], position[1], position[2]
    );
    std::fclose(g);  // end
  }
  TC_STATIC_END_IF
}

// added: write_particle
template <int dim>
void MPM<dim>::write_particle(const std::string &file_name) const {
    auto particles_sorted = particles;
    std::sort(particles_sorted.begin(), particles_sorted.end(),
              [&](ParticlePtr a, ParticlePtr b) {
                return allocator.get_const(a)->id < allocator.get_const(b)->id;
              });
    std::string fn = file_name + std::string(".txt");
    FILE *f = std::fopen(fn.c_str(), "w");
    for (auto p_i : particles_sorted) {
      const Particle *p = allocator.get_const(p_i);
      // int index_p  = p->id;
      Vector pos = p->pos;
      Vector vel = p->get_velocity();
      int type_p = int(p->is_rigid());
      int bound_p = p->near_boundary();
      float dist_p = p->boundary_distance;
      // real dg_p_det = p->dg_p_det;
      // real dg_p_inv_det = p->dg_p_inv_det;
      real tau = p->tau;
      real pres = p->p;
      // real mu_visual = p->mu_visual;
      real gf = p->gf;

      // fmt::print(f, "{} {} {} {} {} {} {} {} {} {} {} {} {} {} {}\n",
      fmt::print(f, "{} {} {} {} {} {} {} {} {} {} {} {}\n",
        type_p,
        pos[0], pos[1], pos[2],
        vel[0], vel[1], vel[2],
        bound_p, dist_p,
        // dg_p_det, dg_p_inv_det,
        tau, pres,
        // mu_visual,
        gf
        );
    }
    std::fclose(f);
}  // end

// added: write_dataset
template <int dim>
void MPM<dim>::write_dataset(RigidBody<dim> const *rigid,
                             const std::string &file_name) const {
    std::string fn = file_name + std::string(".csv");
    FILE *f = std::fopen(fn.c_str(), "w");

    Vector force = rigid->rigid_force;
    real rigid_vel_mag = length(rigid->velocity);
    fmt::print(f, "{}, {}, {}, {}\n",
      force[0], force[1], force[2], rigid_vel_mag);

    auto particles_sorted = particles;
    std::sort(particles_sorted.begin(), particles_sorted.end(),
              [&](ParticlePtr a, ParticlePtr b) {
                return allocator.get_const(a)->id < allocator.get_const(b)->id;
              });

    for (auto p_i : particles_sorted) {
      const Particle *p = allocator.get_const(p_i);
      // int index_p  = p->id;
      int type_p = int(p->is_rigid());
      Vector pos = p->pos;
      real vel_mag = length(p->get_velocity());

      if (type_p == 0) {
        fmt::print(f, "{}, {}, {}, {}, {}\n",
          type_p, pos[0], pos[1], pos[2], vel_mag);
      } else {
        fmt::print(f, "{}, {}, {}, {}\n",
          type_p, pos[0], pos[1], pos[2]);
      }
    }

    std::fclose(f);
}  // end

template <>
void MPM<3>::visualize() const {
  write_bgeo();
}

template <>
void MPM<2>::visualize() const {
  write_bgeo();
}

template void MPM<2>::write_partio(const std::string &file_name) const;
template void MPM<3>::write_partio(const std::string &file_name) const;

template void MPM<2>::write_rigid_body(RigidBody<2> const *rigid,
                                       const std::string &file_name) const;
template void MPM<3>::write_rigid_body(RigidBody<3> const *rigid,
                                       const std::string &file_name) const;

// added
template void MPM<2>::write_particle(const std::string &file_name) const;
template void MPM<3>::write_particle(const std::string &file_name) const;

// added
template void MPM<2>::write_dataset(RigidBody<2> const *rigid,
                                    const std::string &file_name) const;
template void MPM<3>::write_dataset(RigidBody<3> const *rigid,
                                    const std::string &file_name) const;

TC_NAMESPACE_END
