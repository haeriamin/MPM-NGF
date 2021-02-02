/*******************************************************************************
    Copyright (c) The Taichi MPM Authors (2018- ). All Rights Reserved.
    The use of this software is governed by the LICENSE file.
*******************************************************************************/

#ifdef TC_USE_MPI
#include <mpi.h>
#endif

#include <taichi/system/threading.h>
#include <taichi/visual/texture.h>
#include <taichi/math/svd.h>
#include <taichi/math.h>
#include <taichi/common/asset_manager.h>
#include <taichi/common/testing.h>
#include <taichi/system/profiler.h>

#include "articulation.h"
#include "mpm.h"
#include "poisson_disk_sampler.h"
#include "particle_allocator.h"
#include "boundary_particle.h"

TC_NAMESPACE_BEGIN

// initialize ------------------------------------------------------------------
template <int dim>
void MPM<dim>::initialize(const Config &config) {
  TC_P(grid_block_size());
  TC_TRACE("BaseParticle size: {} B", sizeof(Particle));
  Simulation<dim>::initialize(config);
  config_backup = config;
  res = config.get<Vectori>("res");
  apic_damping = config.get("apic_damping", 0.0f);
  rpic_damping = config.get("rpic_damping", 0.0f);
  penalty = config.get("penalty", 0.0f);
  delta_x = config.get("delta_x", delta_x);
  inv_delta_x = 1.0f / delta_x;
  gravity = config.get("gravity", Vector::axis(1) * (-10.0_f));
  apic = config.get("apic", true);
  pushing_force = config.get("pushing_force", 20000.0_f);
  TC_ASSERT_INFO(!config.has_key("delta_t"),
                 "Please use 'base_delta_t' instead of 'delta_t'");
  base_delta_t = config.get("base_delta_t", 1e-4_f);
  base_delta_t *= config.get("dt-multiplier", 1.0_f);
  reorder_interval = config.get<int>("reorder_interval", 1000);
  // cfl number
  cfl = config.get("cfl", 1.0f);
  particle_gravity = config.get<bool>("particle_gravity", true);
  // affine damping
  TC_LOAD_CONFIG(affine_damping, 0.0f);

  // spgrid --------------------------------------------------------------------
  spgrid_size = 4096;
  while (spgrid_size / 2 > (res.max() + 1)) {
    spgrid_size /= 2;
  }
  TC_INFO("Created SPGrid of size {}", spgrid_size);
  // 2D
  TC_STATIC_IF(dim == 2) {
    grid = std::make_unique<SparseGrid>(spgrid_size, spgrid_size);
  }
  // 3D
  TC_STATIC_ELSE {
    grid = std::make_unique<SparseGrid>(spgrid_size, spgrid_size, spgrid_size);
  }
  TC_STATIC_END_IF
  page_map = std::make_unique<SPGrid_Page_Map<log2_size>>(*grid);
  rigid_page_map = std::make_unique<SPGrid_Page_Map<log2_size>>(*grid);
  fat_page_map = std::make_unique<SPGrid_Page_Map<log2_size>>(*grid);
  grid_region = Region(Vectori(0), res + VectorI(1), Vector(0)); // start, end, offset

  /*
  // Restart?
  pakua = create_instance<Pakua>("json");
  Config config_;
  config_.set("frame_directory", config.get_string("frame_directory"));
  pakua->initialize(config_);
  */

  if (config.get("energy_experiment", false)) {
    TC_ASSERT(config.get("optimized", true) == false);
  }

  // Add a background rigid body
  rigids.emplace_back(std::make_unique<RigidBody<dim>>());
  rigids.back()->set_as_background();
}

// add particles ---------------------------------------------------------------
template <int dim>
std::string MPM<dim>::add_particles(const Config &config) {
  auto region = RegionND<dim>(Vectori(0), res);

  if (config.get_string("type") == "rigid") {
    add_rigid_particle(config);
    return std::to_string((int)rigids.size() - 1);
  } else {

    // check if the particle is inside levelset --------------------------------
    auto detect_particles_inside_levelset = [&](const Vector &pos) {
      // check if there is a levelset
      if (this->levelset.levelset0 == nullptr)
        return;
      real phi = this->levelset.sample(pos*inv_delta_x, this->current_t);
      if (phi < 0) {
        TC_WARN("Particles inside levelset generate.\n");
      }
    };

    // global ------------------------------------------------------------------
    auto create_particle = [&](const Vector &coord, real maximum,
                               const Config &config) {
      ParticlePtr p_i;
      Particle *p;
      std::string type = config.get<std::string>("type");
      std::tie(p_i, p) = allocator.allocate_particle(type);

      Config config_new = config;
      std::string param_string[3] = {"cohesion_tex", "theta_c_tex",
                                     "theta_s_tex"};
      for (auto &p : param_string)
        if (config.has_key(p)) {
          std::shared_ptr<Texture> texture =
              AssetManager::get_asset<Texture>(config.get<int>(p));
          real value = texture->sample(coord).x;
          config_new.set(p.substr(0, p.length() - 4), value);
        }

      p->initialize(config_new);
      p->pos = coord;

      if (config_backup.get("sand_climb", false)) {
        std::shared_ptr<Texture> texture = AssetManager::get_asset<Texture>(
            config_backup.get<int>("sand_texture"));
        real speed = config_backup.get("sand_speed", 0.0_f);
        real radius = 15.0_f / 180 * (real)M_PI;
        real x = p->pos[0] + this->current_t * speed * cos(radius);
        real y = p->pos[1] + this->current_t * speed * sin(radius);
        real z = p->pos[2];
        y -= 0.1 / cos(radius);
        y -= x * tan(radius);
        Vector coord(0.5_f);
        coord.x = x;
        coord.y = z;
        if (texture->sample(coord).x < y)
          return;
      }

      // check if particle is near boundary in "mpm.h" -------------------------
      if (this->near_boundary(*p)) {
        TC_WARN("particle out of box or near boundary. Ignored.");
        return;
      }

      // check if the particle is inside levelset ------------------------------
      detect_particles_inside_levelset(coord);

      // global ----------------------------------------------------------------
      p->vol = pow<dim>(delta_x) / maximum;
      // p->vol = (4.0_f / 3.0_f * (real)M_PI * pow<3>(delta_x/2)) / maximum;

      p->set_mass(p->vol * config.get("density", 400.0f) * config.get("packing_fraction", 1.0_f));
      p->set_velocity(config.get("initial_velocity", Vector(0.0f)));

      // stork nod -------------------------------------------------------------
      if (config.get("stork_nod", 0.0_f) > 0.0_f) {
        real x = coord[0];
        real y = coord[1];
        real d = std::max(-0.5_f * x + y - 0.25_f, 0_f);
        real linear = config.get("stork_nod", 0.0_f);
        Vector v(0.0_f);
        v[0] = -d * linear;
        v[1] = -0.5_f * d * linear;
        p->set_velocity(v);
      }

      particles.push_back(p_i); // p_i: particle index
    };

    // global ------------------------------------------------------------------
    // benchmark
    int benchmark = config.get("benchmark", 0);
    if (benchmark) {
      TC_INFO("Generating particles for benchmarking");
      TC_STATIC_IF(dim == 3) {
        // Make sure it is a cube
        TC_ASSERT(res == Vectori(res[0]));
        real s;
        if (benchmark == 125) {
          s = 0.1;
        } else if (benchmark == 8000) {
          s = 0.4;
        } else {
          s = 0;
          TC_ERROR("s must be 125 or 8000");
        }
        int lower = (int)std::round(res[0] * (0.5_f - s));
        int higher = lower + (int)std::round(res[0] * 2 * s);
        //TC_P(lower);
        //TC_P(higher);
        Vector offset(0.25_f * this->delta_x);
        //----------------------------------------------------------------------
        // create 8 particle per cell
        for (auto ind : Region(Vector3i(res[0]*.2,res[0]*.125,res[0]*.2),
                               Vector3i(res[0]*.4,res[0]*.175,res[0]*.4))) {
          for (int i = 0; i < 8; i++) {
            Vector sign(1);
            if (i % 2 == 0)
              sign[0] = -1;
            if (i / 2 % 2 == 0)
              sign[1] = -1;
            if (i / 4 % 2 == 0)
              sign[2] = -1;
            // the followig func is called 8M times
            create_particle(ind.get_pos() * delta_x + offset * sign, 1, config);
          }
        }
      }
      TC_STATIC_END_IF;
      TC_ASSERT(dim == 3);
      TC_P(particles.size());
      return "";
    }

    // point cloud
    if (config.get("point_cloud", false)) {
      TC_STATIC_IF(dim == 3) {
        TC_INFO("Reading point cloud, fn={}",
                config.get<std::string>("mesh_fn"));
        std::string mesh_fn = config.get<std::string>("mesh_fn");
        std::string full_fn = absolute_path(mesh_fn);
        std::shared_ptr<Mesh> mesh = std::make_shared<Mesh>();
        Config mesh_config;
        mesh_config.set("filename", full_fn);
        mesh->initialize(mesh_config);
        Vector scale = config.get("scale", Vector(1.0f));
        Vector translate = config.get("translate", Vector(0.0f));
        for (auto &coord : mesh->vertices) {
          create_particle(id(coord) * id(scale) + translate, 8, config);
        }
      }
      TC_STATIC_ELSE{TC_NOT_IMPLEMENTED} TC_STATIC_END_IF
    } else {

    // else (global) -----------------------------------------------------------
      std::shared_ptr<Texture> density_texture =
          AssetManager::get_asset<Texture>(config.get<int>("density_tex"));
      real maximum = 0.0f;
      // region
      for (auto &ind : region) {
        Vector coord = (ind.get_ipos().template cast<real>() + Vector(0.5f)) *
                       this->delta_x;
        real sample = density_texture->sample(coord).x;
        maximum = std::max(sample, maximum);
      }
      // pd sampler
      if (config.get("pd", true)) {
        PoissonDiskSampler<dim> sampler;
        std::vector<Vector> samples;
        // write periodic pd
        if (config.get("pd_write_periodic_data", false)) {
          sampler.write_periodic_data();
        }
        // periodic pd
        if (config.get("pd_periodic", true)) {
          // source pd
          if (config.get("pd_source", false)) {
            Vector sample_offset =
                config.get("initial_velocity", Vector(0.0f)) * this->current_t;
            real delta_t = config.get("delta_t", 1e-3f);
            Vector sample_advection =
                config.get("initial_velocity", Vector(0.0f)) * delta_t +
                0.5_f * gravity * (delta_t + base_delta_t) * delta_t;

            sampler.sample_from_source(density_texture, region, this->delta_x,
                                       sample_offset, sample_advection,
                                       samples);
          // packed pd
          } else if (config.get("pd_packed", false)) {
            std::shared_ptr<Texture> local_texture =
                AssetManager::get_asset<Texture>(config.get<int>("local_tex"));
            real radius = config.get("packed_radius", 0.01_f);
            real gap = config.get("packed_gap", 0.002_f);
            sampler.sample_packed(density_texture, local_texture, region,
                                  this->delta_x, radius, gap, samples);
          } else {
            // non-source & non-packed pd
            sampler.sample_from_periodic_data(density_texture, region,
                                              this->delta_x, samples);
          }
        } else {
          // non-periodic pd (from paper)
          real minDistAH = config.get("minDistAH", -1.0f);
          sampler.sample(density_texture, region,
             this->delta_x, samples, minDistAH);
        }
        for (auto &coord : samples) {
          // create particles
          create_particle(coord, maximum, config);
          if (config.get<bool>("only_one", false)) {
            break;
          }
        }
      } else {
      // random sampler start
        TC_P(maximum);
        for (auto &ind : region) {
          for (int l = 0; l < maximum; l++) {
            Vector coord =
                (ind.get_ipos().template cast<real>() + Vector::rand()) *
                this->delta_x;
            if (rand() > density_texture->sample(coord).x / maximum) {
              continue;
            }
            create_particle(coord, maximum, config);
          }
        }
      }
      // end of random sampler
    }
  }
  TC_P(particles.size());
  return "";
}

// render particles -------------------------------------------------------- OFF
template <int dim>
std::vector<RenderParticle> MPM<dim>::get_render_particles() const {
  return std::vector<RenderParticle>();
  /*
  std::vector<RenderParticle> render_particles;
  render_particles.reserve(particles.size());
  Vector center(res.template cast<real>() * 0.5f);
  for (auto p_p : particles) {
    Particle &p = *p_p;
    // at least synchronize the position
    Vector pos = p.pos * inv_delta_x - center;
    render_particles.push_back(
        RenderParticle(Vector3(pos), Vector4(0.8f, 0.9f, 1.0f, 0.5f)));
  }
  return render_particles;
  */
}

// added: Reset grid granular fluidity
template <int dim>
void MPM<dim>::reset_grid_granular_fluidity() {
  parallel_for_each_active_grid([&](GridState<dim> &g) {
    g.granular_fluidity = 0.0_f;
  });
}

// normalize grid & apply external force ---------------------------------------
template <int dim>
void MPM<dim>::normalize_grid_and_apply_external_force(
    Vector velocity_increment_) {
  VectorP velocity_increment(velocity_increment_, 0);
  parallel_for_each_active_grid([&](GridState<dim> &g) {
    // dim'th elements of "velocity_and_mass" store mass
    real mass = g.velocity_and_mass[dim];
    if (mass > 0) {
      real inv_mass = 1.0_f / mass;
      VectorP alpha(Vector(inv_mass), 1);
      // Original:
      // g.velocity_and_mass *= alpha;
      // g.velocity_and_mass += velocity_increment;
      g.velocity_and_mass =
          fused_mul_add(g.velocity_and_mass,
            alpha, velocity_increment); // a*b + c from "vector.h"
    }
  });
}

// apply grid boundary conditions ----------------------------------------------
template <int dim>
void MPM<dim>::apply_grid_boundary_conditions(
    const DynamicLevelSet<dim> &levelset,
    real t) {
  int expr_leaky_levelset = config_backup.get<int>("expr_leaky_levelset", 0);
  real hack_velocity = config_backup.get<real>("hack_velocity", 0.0_f);

  int grid_block_size_max = grid_block_size().max();

  auto block_op = [&](uint32 b, uint64 block_offset, GridState<dim> *g) {
    Vectori block_base_coord(SparseMask::LinearToCoord(block_offset));
    Vector center = Vector(block_base_coord + grid_block_size() / Vectori(2));

    //
    if (!expr_leaky_levelset && levelset.inside(center) &&
        std::abs(levelset.sample(center, t)) >= (real)grid_block_size_max) {
      return;
    }

    //
    Region region(Vectori(0), grid_block_size());
    for (auto &ind_ : region) {
      Vectori ind = block_base_coord + ind_.get_ipos();

      if (grid_mass(ind) == 0.0f) {
        continue;
      }

      // GRID position:
      Vector pos = Vector(ind);
      real phi;
      Vector n;
      Vector boundary_velocity;
      real mu;

      //------------------------------------------------------------------------
      // if grid node's -3 <= phi <= 0 (boundary grid) -------------------------
      // and if not leaky levelset
      if (expr_leaky_levelset == 0) {
        phi = levelset.sample(pos, t);
        if (phi < -3 || 0 < phi)  // was 0 
          continue;
        // normall to the levelset which its phi<0
        n = levelset.get_spatial_gradient(pos, t);

        // if hack velocity is ON
        if (hack_velocity != 0.0_f) {
          if (0.5_f < pos.y * delta_x && pos.y * delta_x < 0.7_f &&
              t <= config_backup.get("hack_time", 0.0_f)) {
            boundary_velocity = hack_velocity * Vector::axis(0);
          } else {
            boundary_velocity = Vector(0);
          }
        // main ----------------------------------------------------------------
        } else {
          // for non-dynamic levelset, d(phi)/dt=0
          boundary_velocity =
              -levelset.get_temporal_derivative(pos, t) * n * delta_x;

          // added: Grid granular fluidity boundary condition
          get_grid(ind).granular_fluidity = 0.0_f;
        }

        // boundary friction is the same as levelset friction ------------------
        mu = levelset.levelset0->friction;

        // sand speed ----------------------------------------------------------
        if (config_backup.get("sand_speed", 0.0_f) > 0) {
          real speed = config_backup.get("sand_speed", 0.0_f);
          real radius = 15.0_f / 180 * (real)M_PI;
          boundary_velocity = Vector(0);
          boundary_velocity.x = speed * (-cos(radius));
          boundary_velocity.y = speed * (-sin(radius));
        }

        // gravity cutting -----------------------------------------------------
        if (config_backup.get("gravity_cutting", false)) {
          if (real(ind.y) > 0.7_f * res[1])
            mu = -1;
        }

        // sand crawler --------------------------------------------------------
        if (config_backup.get("sand_crawler", false)) {
          if (real(ind.y) < 0.535_f * res[1])
            mu = -1;
        }

      // leaky levelset (?) ----------------------------------------------------
      } else {
        int y = ind.y;
        if (res[1] / 2 - expr_leaky_levelset <= y && y < res[1] / 2) {
          n = Vector::axis(1);
        } else {
          continue;
        }
        boundary_velocity = Vector(0);
        mu = -1;
      }

      // friction project returns: projected_relative_vel + boundary_velocity
      // in "mpm_fwd.h"
      Vector v = friction_project(grid_velocity(ind), boundary_velocity, n, mu);

      // VectorP has dim+1 number of elements (vector plus!)
      VectorP &v_and_m = get_grid(ind).velocity_and_mass;

      // set v as boundary grid velocity
      // dim'th elements of "v_and_m" stores mass
      v_and_m = VectorP(v, v_and_m[dim]);
    }
  };
  parallel_for_each_block_with_index(block_op, true);
}

// apply dirichlet boundary conditions (like sticky bc) ------------------------
// 2D
template <>
void MPM<2>::apply_dirichlet_boundary_conditions() {
  real distance   = config_backup.get("dirichlet_boundary_radius", 0.0_f);
  real distance_l = config_backup.get("dirichlet_distance_left", distance);
  real distance_r = config_backup.get("dirichlet_distance_right", distance);

  real velocity   = config_backup.get("dirichlet_boundary_velocity", 0.0_f);
  real velocity_l = config_backup.get("dirichlet_boundary_left", velocity);
  real velocity_r = config_backup.get("dirichlet_boundary_right", velocity);

  Vector vl(0.0_f), vr(0.0_f);
  vl[0] = velocity_l;
  vr[0] = velocity_r;

  for (auto &ind_ : grid_region) {
    Vectori ind = ind_.get_ipos();
    Vector pos = Vector(ind) * delta_x;
    if (pos[0] < distance_l) {
      get_grid(ind).velocity_and_mass =
          VectorP(vl, get_grid(ind).velocity_and_mass[3]);
    } else if (pos[0] > 1.0_f - distance_r) {
      get_grid(ind).velocity_and_mass =
          VectorP(vr, get_grid(ind).velocity_and_mass[3]);
    }
  }
}
// 3D
template <>
void MPM<3>::apply_dirichlet_boundary_conditions() {
  for (auto &ind_ : grid_region) {
    Vectori ind = ind_.get_ipos();
    Vector pos = Vector(ind) * delta_x;
    Vector v(0);
    if (pos.y > 0.525_f) {
      get_grid(ind).velocity_and_mass =
          VectorP(v, get_grid(ind).velocity_and_mass[3]);
    }
  }
  //  real radius = config_backup.get("dirichlet_boundary_radius", 0.0_f);
  //  real linear_v =
  //      config_backup.get("dirichlet_boundary_linear_velocity", 0.0_f);
  //
  //  TC_ASSERT_INFO(!config_backup.has_key("dirichlet_boundary_angular_velocity"),
  //                 "Use 'dirichlet_boundary_rotation_cycle' instead of "
  //                 "'dirichlet_boundary_angular_velocity'.")
  //  real rotation_cycle =
  //      config_backup.get("dirichlet_boundary_rotation_cycle", 0.0_f);
  //  for (auto &ind_ : grid_region) {
  //    Vectori ind = ind_.get_ipos();
  //    Vector pos = Vector(ind) * delta_x;
  //    real dist_to_x = sqrt(sqr(pos[1] - 0.5_f) + sqr(pos[2] - 0.5_f));
  //    if (sqrt(sqr(pos[0] - 0.5f) + sqr(pos[1] - 0.5_f) + sqr(pos[2] - 0.5_f))
  //    <
  //        radius) {
  //      Vector v;
  //      // avoid atan2(0, 0)
  //      if (dist_to_x < eps) {
  //        v = Vector(linear_v, 0.0_f, 0.0_f);
  //      } else {
  //        real rotation_angle = atan2(pos[1] - 0.5_f, pos[2] - 0.5_f);
  //        real rotation_speed;
  //        if (rotation_cycle > 0.0_f) {
  //          rotation_speed = dist_to_x * 2.0_f * M_PI / rotation_cycle;
  //        } else {
  //          rotation_speed = 0.0_f;
  //        }
  //        v = Vector(linear_v, rotation_speed * cos(rotation_angle),
  //                   -rotation_speed * sin(rotation_angle));
  //      }
  //      get_grid(ind).velocity_and_mass =
  //          VectorP(v, get_grid(ind).velocity_and_mass[3]);
  //    }
  //  }
}

// particle-levelset interaction -------------------------------------- : ON/OFF
template <int dim>
void MPM<dim>::particle_collision_resolution(real t) {
  parallel_for_each_particle([&](Particle &p) {
    Vector pos = p.pos * inv_delta_x;
    real phi = this->levelset.sample(pos, t);
    // if there is collision (phi<0)
    if (phi <= 0.25) {
      Vector gradient = this->levelset.get_spatial_gradient(pos, t);
      p.pos -= gradient * phi * delta_x;
      p.set_velocity(p.get_velocity()-dot(gradient, p.get_velocity())*gradient);
    }
  });
}

// added: particle_bc_at_levelset ------------------------------------- : ON/OFF
template <int dim>
void MPM<dim>::particle_bc_at_levelset(real t) {
  parallel_for_each_particle([&](Particle &p) {
    Vector pos = p.pos * inv_delta_x;
    real phi = this->levelset.sample(pos, t);
    if (phi < 0.25)
      p.gf = 0.0_f;
  });
}

// step ------------------------------------------------------------------------
template <int dim>
void MPM<dim>::step(real dt) {
  if (dt < 0) {
    substep();
    request_t = this->current_t;
  } else {
    request_t += dt;
    while (this->current_t + base_delta_t < request_t) {
      update_counter += this->particles.size();
      substep();
    }
  }
  TC_TRACE("#Particles {}", particles.size());
  auto average = std::accumulate(rigid_block_fractions.begin(),
                                 rigid_block_fractions.end(), 0.0) /
                 rigid_block_fractions.size();
  TC_TRACE("Average rigid block fraction: {:.2f}%", 100 * average);
  step_counter += 1;
  if (config_backup.get("print_energy", false)) {
    TC_P(calculate_energy());
  }
  TC_WARN("Times of particle updating : {}", update_counter);
}

//------------------------------------------------------------------------------
// MPM subset (MAIN LOOP) ------------------------------------------------------
template <int dim>
void MPM<dim>::substep() {
  Profiler _p("mpm_substep");
  cutting_counter    = 0;
  plasticity_counter = 0;
  const real delta_t = base_delta_t;
  if (particles.empty()) {
    // TC_DEBUG("dt = {}, No particles", base_delta_t);
  } else {
    // TC_DEBUG("dt = {}, Do substep for {} particles, ", base_delta_t,
    // particles.size());
  }

  TC_PROFILE("sort_particles_and_populate_grid",
             sort_particles_and_populate_grid());

  // added: Reset grid granular fluidity
  TC_PROFILE("reset_grid_granular_fluidity",
              reset_grid_granular_fluidity());

  // articulate ----------------------------------------------------------------
  if (has_rigid_body()) {
    for (int i = 0; i < config_backup.get("coupling_iterations", 1); i++) {
      // check rigidBody collision --------------------------------------- : OFF
      TC_PROFILE("rigidify", rigidify(delta_t));
      // rigid body articulation in "mpm.h" ------------------------------------
      TC_PROFILE("articulate", articulate(delta_t));
      // modified (CDF)
      TC_PROFILE("rasterize_rigid_boundary", rasterize_rigid_boundary());
    }
  }

  // visualize CDF ------------------------------------------------------- : OFF
  if (config_backup.get("visualize_cdf", false)) {
    int counter = 0;
    for (auto &ind : grid_region) {
      Particle *p = allocator[counter];
      GridState<dim> &g = get_grid(ind);
      p->pos = ind.get_pos() * delta_x;
      p->boundary_distance = g.distance;
      p->states = g.states;
      counter++;
    }
    TC_PROFILE("advect_rigid_bodies", advect_rigid_bodies(delta_t));
    this->current_t += delta_t;
    substep_counter += 1;
    return;
  }

  // visualize particle CDF ---------------------------------------------- : OFF
  if (config_backup.get("visualize_particle_cdf", false)) {
    int counter = 0;
    for (auto &ind : grid_region) {
      Region region(Vectori(0), Vectori(4));
      for (auto &sub_ind : region) {
        Particle *p = allocator[counter];
        TC_ASSERT((std::size_t)counter < allocator.pool.size());
        p->states = 0;
        p->pos = (ind.get_pos() + sub_ind.get_pos() * 0.25_f) * delta_x;
        counter++;
      }
    }
    TC_PROFILE("gather_cdf", gather_cdf());
    TC_PROFILE("advect_rigid_bodies", advect_rigid_bodies(delta_t));
    this->current_t += delta_t;
    substep_counter += 1;
    return;
  }

  // gather CDF ----------------------------------------------------------------
  if (has_rigid_body()) {
    TC_PROFILE("gather_cdf", gather_cdf());
  }

  // added: particle bc near levelsets -------------------------------- : On/OFF
  if (config_backup.get("particle_bc_at_levelset", false)) {
    TC_PROFILE("particle_bc_at_levelset",
      particle_bc_at_levelset(this->current_t));
  }

  // rasterize (particle to grid) ----------------------------------------------
  if (!config_backup.get("benchmark_rasterize", false)) {
    // optimized : ON
    if (config_backup.get("optimized", true)) {
      TC_PROFILE_TPE("P2G optimized",
      rasterize_optimized(delta_t), particles.size());
    // else : OFF
    } else {
      TC_PROFILE_TPE("P2G", rasterize(delta_t), particles.size());
    }
  } else {
    while (true) {
      Time::Timer _("Rasterize x 20");
      base_delta_t = 0;
      for (int i = 0; i < 20; i++) {
        rasterize_optimized(delta_t);
      }
    }
  }

  // add gravity to grid instead of particles ----------------------------------
  Vector gravity_velocity_increment = gravity * delta_t;
  if (particle_gravity) {
    gravity_velocity_increment = Vector(0);
  }
  TC_PROFILE(
      "normalize_grid_and_apply_external_force",
      normalize_grid_and_apply_external_force(gravity_velocity_increment));

  // rigidBody-levelset collision ---------------------------------------- : OFF
  if (config_backup.get("rigid_body_levelset_collision", false)) {
    TC_PROFILE("rigid_body_levelset_collision",
      rigid_body_levelset_collision(this->current_t, delta_t));
  }

  // boundary condition --------------------------------------------------------
  TC_PROFILE("boundary_condition",
    apply_grid_boundary_conditions(this->levelset, this->current_t));

  // ---------------------------------------------------------------------------
  if (config_backup.get("dirichlet_boundary_radius", 0.0_f) > 0.0_f) {
    TC_PROFILE("apply_dirichlet_boundary_conditions",
      apply_dirichlet_boundary_conditions());
  }

  // resample (grid to particle) -----------------------------------------------
  if (!config_backup.get("benchmark_resample", false)) {
    // optimized : ON
    if (config_backup.get("optimized", true)) {
      TC_PROFILE_TPE("G2P optimized", resample_optimized(), particles.size());
    // else : OFF
    } else {
      TC_PROFILE_TPE("G2P", resample(), particles.size());
    }
  } else {
    // NOTE: debugging here
    while (true) {
      Time::Timer _("Resample x 20");
      base_delta_t = 0;
      for (int i = 0; i < 20; i++) {
        resample_optimized();
      }
    }
  }

  // clean boundary particles --------------------------------------------------
  if (config_backup.get("clean_boundary", true)) {
    TC_PROFILE("clean boundary", clear_boundary_particles());
  }

  // particle collision ------------------------------------------------ : On/OFF
  if (config_backup.get("particle_collision", false)) {
    TC_PROFILE("particle_collision",
      particle_collision_resolution(this->current_t));
  }

  // advect rigid body ---------------------------------------------------------
  if (has_rigid_body()) {
    TC_PROFILE("advect_rigid_bodies", advect_rigid_bodies(delta_t));
  }

  this->current_t += delta_t;
  substep_counter += 1;
}
// subset end ------------------------------------------------------------------
//------------------------------------------------------------------------------

template <int dim>
bool MPM<dim>::test() const {
  return true;
}

// clear boundary particles ----------------------------------------------------
template <int dim>
void MPM<dim>::clear_boundary_particles() {
  std::vector<ParticlePtr> particles_new;
  static bool has_deleted = false;
  int remove_particles = config_backup.get("remove_particles", 0);
  real remove_height = config_backup.get("remove_height", 0.02_f);

  // Do not use bool here for concurrency
  std::vector<uint8> particle_remaining(particles.size(), 0);

  tbb::parallel_for(0, (int)particles.size(), [&](int i) {
    auto p_i = particles[i];
    auto p = allocator[p_i];
    if (near_boundary(*p) || p->pos.abnormal() ||
        p->get_velocity().abnormal()) {
      return;
    }
    if (remove_particles) {
      int kind = remove_particles;
      real h = remove_height;
      if (!has_deleted && this->current_t >= 0.1) {
        if (p->pos.x >= 0.45 && p->pos.y >= 0.6 - 0.5 * h &&
            p->pos.y <= 0.6 + 0.5 * h) {
          if (kind == 1) {
            return;
          } else if (kind == 2) {
            p->set_mu_to_zero();
          } else if (kind == 3) {
            p->set_lambda_and_mu_to_zero();
          }
        }
      }
    }
    particle_remaining[i] = 1;
  });

  particles_new.reserve(particles.size());
  for (int i = 0; i < (int)particles.size(); i++) {
    if (particle_remaining[i])
      particles_new.push_back(particles[i]);
  }

  if (!has_deleted && this->current_t >= 0.1)
    has_deleted = true;
  int deleted = (int)particles.size() - (int)particles_new.size();
  if (deleted != 0 && config_backup.get("warn_particle_deletion", true)) {
    TC_WARN(
        "{} boundary (or abnormal) particles deleted.\n{} Particles remained\n",
        deleted, particles_new.size());
  }
  particles = particles_new;
}

// get debug information -------------------------------------------------------
template <int dim>
std::string MPM<dim>::get_debug_information() {
  // return std::to_string(rigids[1].velocity.x);
  return "";
}

template <int dim>
MPM<dim>::~MPM() {
}

template <typename T>
Array2D<T> horizontal_stack(const Array2D<T> &a, const Array2D<T> &b) {
  TC_ASSERT(a.get_res()[1] == b.get_res()[1]);
  Array2D<T> merged(Vector2i(a.get_res()[0] + b.get_res()[0], a.get_res()[1]));
  for (auto &ind : a.get_region()) {
    merged[ind] = a[ind];
  }
  for (auto &ind : b.get_region()) {
    merged[ind + Vector2i(a.get_res()[0], 0)] = b[ind];
  }
  return merged;
}

template <typename T>
Array2D<T> vertical_stack(const Array2D<T> &a, const Array2D<T> &b) {
  TC_ASSERT(a.get_res()[0] == b.get_res()[0]);
  Array2D<T> merged(Vector2i(a.get_res()[0], a.get_res()[1] + b.get_res()[1]));
  for (auto &ind : a.get_region()) {
    merged[ind] = a[ind];
  }
  for (auto &ind : b.get_region()) {
    merged[ind + Vector2i(0, a.get_res()[1])] = b[ind];
  }
  return merged;
}

// draw cdf --------------------------------------------------------------------
// 2D
template <>
void MPM<2>::draw_cdf(const Config &config) {
  // Step 1: sample particles
  int scale = 5;

  for (int i = 7; i < (res[0] - 7) * scale; i++) {
    for (int j = 7; j < (res[1] - 7) * scale; j++) {
      auto alloc = allocator.allocate_particle("snow");
      Particle &p = *alloc.second;
      p.pos =
          Vector(delta_x * (i + 0.5_f) / scale, delta_x * (j + 0.5_f) / scale);
      p.vol = pow<dim>(delta_x);
      p.set_mass(p.vol * 400);
      p.set_velocity(config.get("initial_velocity", Vector(0.0f)));
      particles.push_back(alloc.first);
    }
  }
  base_delta_t = 0.001;
  TC_P("seeded");
  for (int i = 0; i < 5; i++) {
    substep();
  }
  TC_P("stepped");
  Array2D<Vector3> img_dist, img_affinity, img_normal;
  img_dist.initialize(res * scale);
  img_affinity.initialize(res * scale);
  img_normal.initialize(res * scale);

  // Step 2: compute sample color and gather from grid
  for (auto &p_i : particles) {
    auto &p = *allocator[p_i];
    if (p.near_boundary()) {
      Vectori ipos = (p.pos * real(scale) * inv_delta_x).template cast<int>();
      TC_ASSERT(img_affinity.inside(ipos));
      // img[ipos] = Vector3(p->near_boundary);
      img_dist[ipos] = Vector3(p.boundary_distance * inv_delta_x * 0.2f);
      Vector3 affinity_color;
      affinity_color.x = p.states % 4 * 0.33_f;
      affinity_color.y = p.states / 4 % 4 * 0.33_f;
      affinity_color.z = p.states / 16 % 4 * 0.33_f;
      img_affinity[ipos] = Vector3(affinity_color);
      img_normal[ipos] = Vector3(p.boundary_normal * 0.5_f + Vector(0.5_f));
    }
  }
  auto merged_lower =
      horizontal_stack(horizontal_stack(img_dist, img_affinity), img_normal);

  // Step 3: Grid CDF
  Region region(Vectori(0), res * scale);
  for (auto &ind : region) {
    Vectori ipos = ind.get_ipos() / Vectori(scale);
    img_dist[ind] = Vector3(get_grid(ipos).get_distance() * inv_delta_x * 0.2f);
    Vector3 affinity_color;
    auto states = get_grid(ipos).get_states();
    affinity_color.x = states % 4 * 0.33_f;
    affinity_color.y = states / 4 % 4 * 0.33_f;
    affinity_color.z = states / 16 % 4 * 0.33_f;
    img_affinity[ind] = Vector3(affinity_color);
  }
  auto merged_upper =
      horizontal_stack(horizontal_stack(img_dist, img_affinity), img_normal);
  auto merged = vertical_stack(merged_lower, merged_upper);
  merged.write_as_image("/tmp/cdf.png");
}

// 3D
template <>
void MPM<3>::draw_cdf(const Config &config) {
  TC_NOT_IMPLEMENTED
}
template <int dim>
void MPM<dim>::sort_allocator() {
  TC_TRACE("Reordering particles");
  Time::Timer _("reorder");

  std::swap(allocator.pool_, allocator.pool);
  if (allocator.pool.size() < allocator.pool_.size()) {
    allocator.pool.resize(allocator.pool_.size());
  }

  for (uint32 i = 0; i < particles.size(); i++) {
    int j = particles[i];
    allocator.pool[i] = allocator.pool_[j];
    particles[i] = i;
  }
  allocator.gc(particles.size());
}

// sort particles & populate grid ----------------------------------------------
template <int dim>
void MPM<dim>::sort_particles_and_populate_grid() {
  // Profiler::disable();
  constexpr int index_bits = (32 - SparseMask::block_bits);

  TC_ASSERT(particles.size() < (1 << index_bits));

  if (particles.size() > particle_sorter.size()) {
    particle_sorter.resize(particles.size());
  }

  auto grid_array = grid->Get_Array();

  {
    Profiler _("prepare array to sort");
    tbb::parallel_for(0, (int)particles.size(), [&](int i) {
      uint64 offset = SparseMask::Linear_Offset(to_std_array(
          get_grid_base_pos(allocator[particles[i]]->pos * inv_delta_x)));
      particle_sorter[i] =
          ((offset >> SparseMask::data_bits) << index_bits) + i;
    });
  }

  TC_PROFILE("parallel_sort",
             tbb::parallel_sort(particle_sorter.begin(),
                                particle_sorter.begin() + particles.size()));

  {
    Profiler _("reorder particle pointers");
    // Reorder particles
    std::swap(particles, particles_);
    // if (particles.size() < particles_.size()) {
    // TODO: make it more efficient
    particles.resize(particles_.size());
    //}
    for (int i = 0; i < (int)particles.size(); i++) {
      particles[i] = particles_[particle_sorter[i] & ((1ll << index_bits) - 1)];
    }
  }

  // Reorder particles
  if (reorder_interval > 0 && substep_counter % reorder_interval == 0) {
    sort_allocator();
  }

  // Reset page_map
  {
    Profiler _("reset page map");
    page_map->Clear();
    for (int i = 0; i < (int)particles.size(); i++) {
      uint64 offset = (particle_sorter[i] >> index_bits)
                      << SparseMask::data_bits;
      page_map->Set_Page(offset);
    }
    page_map->Update_Block_Offsets();
  }

  // Update particle offset
  auto blocks = page_map->Get_Blocks();

  {
    Profiler _("fat page map");
    // Reset fat_page_map
    fat_page_map->Clear();
    for (int b = 0; b < (int)blocks.second; b++) {
      auto base_offset = blocks.first[b];
      TC_STATIC_IF(dim == 2) {
        auto x = 1 << SparseMask::block_xbits;
        auto y = 1 << SparseMask::block_ybits;
        auto c = SparseMask::LinearToCoord(base_offset);
        for (int i = -1 + (c[0] == 0); i < 2; i++) {
          for (int j = -1 + (c[1] == 0); j < 2; j++) {
            fat_page_map->Set_Page(SparseMask::Packed_Add(
                base_offset, SparseMask::Linear_Offset(x * i, y * j)));
          }
        }
      }
      TC_STATIC_ELSE {
        auto x = 1 << SparseMask::block_xbits;
        auto y = 1 << SparseMask::block_ybits;
        auto z = 1 << id(SparseMask::block_zbits);
        auto c = SparseMask::LinearToCoord(base_offset);
        for (int i = -1 + (c[0] == 0); i < 2; i++) {
          for (int j = -1 + (c[1] == 0); j < 2; j++) {
            for (int k = -1 + (c[2] == 0); k < 2; k++) {
              fat_page_map->Set_Page(SparseMask::Packed_Add(
                  base_offset, SparseMask::Linear_Offset(x * i, y * j, z * k)));
            }
          }
        }
      }
      TC_STATIC_END_IF
    }
    fat_page_map->Update_Block_Offsets();
  }

  auto fat_blocks = fat_page_map->Get_Blocks();
  {
    Profiler _("reset grid");
    for (int i = 0; i < (int)fat_blocks.second; i++) {
      auto offset = fat_blocks.first[i];
      std::memset(&grid_array(offset), 0, 1 << log2_size);
    }
  }

  {
    Profiler _("block particle offset");
    block_meta.clear();
    uint64 last_offset = -1;
    for (uint32 i = 0; i < particles.size(); i++) {
      if (last_offset != (particle_sorter[i] >> 32)) {
        block_meta.push_back({i, 0});
      }
      last_offset = particle_sorter[i] >> 32;
    }

    block_meta.push_back({(uint32)particles.size(), 0});
    TC_ASSERT(block_meta.size() == blocks.second + 1);
  }

  {
    Profiler _("grid particle offset");
    parallel_for_each_block_with_index(
        [&](uint32 b, uint64 base_offset, GridState<dim> *g) {
          auto particle_begin = block_meta[b].particle_offset;
          auto particle_end = block_meta[b + 1].particle_offset;

          int count = 0;
          for (uint32 i = particle_begin; i < particle_end; i++) {
            auto base_pos =
                get_grid_base_pos(allocator[particles[i]]->pos * inv_delta_x);
            uint64 offset = SparseMask::Linear_Offset(to_std_array(base_pos));
            auto index = (offset >> SparseMask::data_bits) &
                         ((1 << SparseMask::block_bits) - 1);
            g[index].particle_count += 1;
            count++;
          }
          // TC_ASSERT(count == particle_end - particle_begin);
        },
        false);
  }
  // Profiler::enable();
  TC_STATIC_IF(dim == 3) {
    Profiler _("update rigid page map");
    this->update_rigid_page_map();
  }
  TC_STATIC_END_IF
}

// general actions -------------------------------------------------------------
template <int dim>
std::string MPM<dim>::general_action(const Config &config) {
  std::string action = config.get<std::string>("action");

  // add articulation ----------------------------------------------------------
  if (action == "add_articulation") {
    Config new_config = config;
    new_config.set("obj0", rigids[config.get<int>("obj0")].get());
    if (config.has_key("obj1")) {
      new_config.set("obj1", rigids[config.get<int>("obj1")].get());
    } else {
      new_config.set("obj1", rigids[0].get());
    }
    auto art = create_instance_unique<Articulation<dim>>(
        config.get<std::string>("type"), new_config);
    articulations.push_back(std::move(art));

  // draw cdf ------------------------------------------------------------------
  } else if (action == "cdf") {
    draw_cdf(config);

  // save ----------------------------------------------------------------------
  } else if (action == "save") {
    TC_P(this->get_name());
    write_to_binary_file_dynamic(this, config.get<std::string>("file_name"));

  // calculate energy ----------------------------------------------------------
  } else if (action == "calculate_energy") {
    return fmt::format("{}", calculate_energy());

  // load rigid body from binary file ------------------------------------------
  } else if (action == "load") {
    read_from_binary_file_dynamic(this, config.get<std::string>("file_name"));
    for (auto &r : rigids) {
      if (r->pos_func_id != -1) {
        typename RigidBody<dim>::PositionFunctionType *f =
            (typename RigidBody<dim>::PositionFunctionType *)config.get<uint64>(
                fmt::format("script{:05d}", r->pos_func_id));
        r->pos_func = *f;
        TC_INFO("scripted position loaded");
      }
      if (r->rot_func_id != -1) {
        typename RigidBody<dim>::RotationFunctionType *f =
            (typename RigidBody<dim>::RotationFunctionType *)config.get<uint64>(
                fmt::format("script{:05d}", r->rot_func_id));
        r->rot_func = *f;
        TC_INFO("scripted rotation loaded");
      }
    }

  // delete particles inside level set -----------------------------------------
  } else if (action == "delete_particles_inside_level_set") {
    std::vector<ParticlePtr> particles_new;
    int deleted = 0;
    for (auto p_i : particles) {
      auto p = allocator[p_i];
      if (this->levelset.sample(p->pos * inv_delta_x, this->current_t) < 0) {
        deleted += 1;
        continue;
      }
      particles_new.push_back(p_i);
    }
    TC_INFO(
        "{} boundary (or abnormal) particles deleted.\n{} Particles remained\n",
        deleted, particles_new.size());
    particles = particles_new;
  } else {
    TC_ERROR("Unknown action: {}", action);
  }
  return "";
}
template std::string MPM<2>::general_action(const Config &config);
template std::string MPM<3>::general_action(const Config &config);
using MPM2D = MPM<2>;
using MPM3D = MPM<3>;
TC_IMPLEMENTATION(Simulation2D, MPM2D, "mpm");
TC_IMPLEMENTATION(Simulation3D, MPM3D, "mpm");

// mpm serialization -----------------------------------------------------------
TC_TEST("mpm_serialization") {
  std::string fn = "/tmp/snapshot.tcb", fn2;
  MPM<3> mpm, mpm_loaded;
  mpm.penalty = 2013011374;
  mpm.res = Vector3i(42, 95, 63);
  mpm.delta_x = 1024;
  mpm.inv_delta_x = 1025;
  write_to_binary_file(mpm, fn);
  read_from_binary_file(mpm_loaded, fn);
  if (remove(fn.c_str()) != 0) {
    TC_WARN("Error occur when deleting binary file.");
  }
#define CHECK_LOADED(x) CHECK(mpm_loaded.x == mpm.x)
  CHECK_LOADED(res);
  CHECK_LOADED(penalty);
  CHECK_LOADED(delta_x);
  CHECK_LOADED(inv_delta_x);
  // CHECK_LOADED(grid_region);
  std::vector<Element<2>> elems, elems2;
  elems.resize(1000);
  write_to_binary_file(elems, fn);
  read_from_binary_file(elems2, fn);
  if (remove(fn.c_str()) != 0) {
    TC_WARN("Error occur when deleting binary file.");
  }
}

// update rigid page map -------------------------------------------------------
// 2D
template <>
void MPM<2>::update_rigid_page_map() {
  TC_NOT_IMPLEMENTED
}
// 3D
template <>
void MPM<3>::update_rigid_page_map() {
  constexpr int dim = 3;
  auto block_size = grid_block_size();
  rigid_page_map->Clear();
  std::vector<int> block_has_rigid(page_map->Get_Blocks().second, 0);
  auto block_op = [&](uint32 b, uint64 block_offset, GridState<dim> *g) {
    int particle_begin;
    int particle_end = block_meta[b].particle_offset;

    bool has_rigid = false;
    for (uint32 t = 0; t < SparseMask::elements_per_block; t++) {
      particle_begin = particle_end;
      particle_end += g[t].particle_count;
      for (int k = particle_begin; k < particle_end; k++) {
        Particle &p = *allocator[particles[k]];
        has_rigid = has_rigid || p.is_rigid();
      }
    }
    if (has_rigid) {
      block_has_rigid[b] = true;
    }
  };
  parallel_for_each_block_with_index(block_op, false, false);
  auto blocks = page_map->Get_Blocks();
  // We do not have thread-safe Set_Page...
  for (uint i = 0; i < blocks.second; i++) {
    auto offset = blocks.first[i];
    if (!block_has_rigid[i]) {
      continue;
    }
    Vectori base_pos(SparseMask::LinearToCoord(offset));
    for (auto &ind : Region(Vectori(-1), Vectori(2))) {
      if (VectorI(0) <= VectorI(ind) && VectorI(ind) < VectorI(spgrid_size)) {
        auto nei_offset =
            SparseMask::Linear_Offset(base_pos + ind.get_ipos() * block_size);
        TC_ASSERT(fat_page_map->Test_Page(offset));
        rigid_page_map->Set_Page(nei_offset);
      }
    }
  }
  rigid_page_map->Update_Block_Offsets();
  rigid_block_fractions.push_back(1.0_f * rigid_page_map->Get_Blocks().second /
                                  fat_page_map->Get_Blocks().second);
}

// calculate energy ------------------------------------------------------------
template <int dim>
real MPM<dim>::calculate_energy() {
  sort_particles_and_populate_grid();
  rasterize(0, false);
  float64 kinetic = 0;
  float64 potential = 0;
  auto blocks = fat_page_map->Get_Blocks();
  auto grid_array = grid->Get_Array();
  for (int b = 0; b < (int)blocks.second; b++) {
    GridState<dim> *g =
        reinterpret_cast<GridState<dim> *>(&grid_array(blocks.first[b]));
    for (int i = 0; i < (int)SparseMask::elements_per_block; i++) {
      real mass = g[i].velocity_and_mass[dim];
      if (mass == 0)
        continue;
      Vector velocity(g[i].velocity_and_mass / mass);
      kinetic += mass * velocity.length2() * 0.5_f;
    }
  }
  std::mutex mut;
  parallel_for_each_particle([&](Particle &p) {
    // We do not care about performance here.
    mut.lock();
    potential += p.potential_energy();
    mut.unlock();
  });
  TC_WARN("t={} Energy: kinetic={}, potential={}", this->current_t, kinetic,
          potential);
  auto mechanical_energy = kinetic + potential;
  TC_INFO("Visualize energy: t={}, energy={}", this->current_t,
          mechanical_energy);
  return mechanical_energy;
}
TC_NAMESPACE_END
