## Excavation model prepared for generating training dataset
## Amin Haeri [ahaeri92@gmail.com]
## Early 2021

#!/usr/bin/env python

import math as m
import sys

import taichi


if __name__ == '__main__':

    # Dimensional analysis
    scale_dim = 0.25
    scale_time = scale_dim * 2
    scale_vel = scale_dim * 2

    # Get inputs
    c = int(sys.argv[1])
    depth = float(sys.argv[2])
    speed = float(sys.argv[3])
    angle = float(sys.argv[4])
    motion = float(sys.argv[5])
    print(c, '-', depth/scale_dim, '-', speed/scale_vel, '-', angle, '-', motion)

    # Discretization
    r = 301  # Spatial resolution
    dx = 1/r
    dt = 4e-5
    ppc = 1

    # Gravity
    G = 9.81

    # Data acquisition rate
    frame_rate = 60 / scale_time  # [Hz]

    # Bin
    friction_ls = -1
    offset = 0.05/2
    x_bin = 2
    y_bin = 0.375
    z_bin = 1.5

    # Soil
    modulus_elastic = 15e6 * scale_dim
    ratio_poisson = 0.3
    modulus_shear = modulus_elastic / 2 / (1 + ratio_poisson)
    modulus_bulk = modulus_elastic / 3 / (1 - 2 * ratio_poisson)
    friction_static = 0.7600
    friction_2 = 0.9616
    grain_diameter = 0.0003 * scale_dim
    grain_density = 2583
    packing_fraction = 0.67
    critical_density = packing_fraction * grain_density
    
    # Bucket
    x_buck = 0.0070
    y_buck = 0.0763
    z_buck = 0.1142
    friction_rb = 0.3
    dist_horizontal = 0.4 * scale_dim  # [m]
    speed_vertical = 0.025 * scale_vel  # [m/s]
    time_rest = 3 * scale_time  # [s]

    # Start MPM
    time_vertical = depth / speed_vertical
    time_horizontal = round(dist_horizontal / speed, 1)

    if motion == 1:
        time_total = time_vertical + time_horizontal + time_rest
    elif motion == 2:
        time_total =  time_horizontal + time_rest

    mpm = taichi.dynamics.MPM(
        task_id=str(c)+'_D'+str(depth/scale_dim)+'m_S'+str(speed/scale_vel)+
            'ms-1_A'+str(angle)+'deg_M'+str(motion),
        res=(r, r, r),
        base_delta_t=dt,
        frame_dt=1/frame_rate,
        num_frames=time_total*frame_rate,
        num_threads=-1,
        gravity=(0, -G, 0),
        particle_gravity=True,
        rigidBody_gravity=False,
        rigid_body_collision=False,
        rigid_body_levelset_collision=False,
        particle_collision=False,
        pushing_force=0,
        print_rigid_body_state=False,
        clean_boundary=False,
        warn_particle_deletion=False,
        cdf_3d_modified=True,
        compute_particle_impulses=True,
        visualize_particle_impulses=False,
        affect_particle_impulses=False,
        particle_bc_at_levelset=False, 

        snapshots=False,
        verbose_bgeo=False,
        write_particle=False,
        write_rigid_body=False,
        write_partio=False,
        write_dataset=True,
    )

    # level-set --
    levelset = mpm.create_levelset()
    levelset.add_plane(taichi.Vector(1, 0, 0), -offset)
    levelset.add_plane(taichi.Vector(0, 1, 0), -offset)
    levelset.add_plane(taichi.Vector(0, 0, 1), -offset)
    levelset.add_plane(taichi.Vector(-1, 0, 0), offset+x_bin/10)
    levelset.add_plane(taichi.Vector(0, 0, -1), offset+z_bin/10)
    levelset.set_friction(friction_ls)
    mpm.set_levelset(levelset, False)

    # soil --
    tex = taichi.Texture(
        'mesh',
        scale=taichi.Vector(x_bin, y_bin, z_bin),
        translate=(offset+x_bin/20, offset+y_bin/20, offset+z_bin/20),
        resolution=(2*r, 2*r, 2*r),
        mesh_accuracy=3,
        filename='projects/mpm/data/cube_smooth.obj',
    ) * ppc

    mpm.add_particles(
        type='nonlocal',
        pd=True,
        density_tex=tex.id,
        density=grain_density,
        critical_density=critical_density,
        packing_fraction=packing_fraction,
        S_mod=modulus_shear,
        B_mod=modulus_bulk,
        dia=grain_diameter,
        mu_s=friction_static,
        mu_2=friction_2,
        A_mat=0.48,
        I_0=0.278,
        t_0=1e-4,
    )

    # bucket --
    def position_function(t):
        x_percentage = 0.9
        z_percentage = 0.5

        if motion == 1:
            if (t <= time_vertical):
                return taichi.Vector(
                    offset + x_percentage*x_bin/10 + m.sin(m.radians(angle))*y_buck/2
                        - speed_vertical * t,
                    offset + y_bin/10 + m.cos(m.radians(angle))*y_buck/2 + dx/2
                        - speed_vertical * t,
                    offset + z_percentage * z_bin/10
                )
            elif (t > time_vertical and t <= time_vertical + time_horizontal):
                return taichi.Vector(
                    offset + x_percentage*x_bin/10 + m.sin(m.radians(angle))*y_buck/2
                        - speed_vertical * time_vertical
                        - speed * (t - time_vertical),
                    offset + y_bin/10 + m.cos(m.radians(angle))*y_buck/2 + dx/2
                        - speed_vertical * time_vertical,
                    offset + z_percentage*z_bin/10
                )
            elif (t > time_vertical + time_horizontal and t <= time_vertical*2 + time_horizontal):
                return taichi.Vector(
                    offset + x_percentage*x_bin/10 + m.sin(m.radians(angle))*y_buck/2
                        - speed_vertical * time_vertical
                        - speed * time_horizontal
                        - speed_vertical * (t - (time_vertical + time_horizontal)),
                    offset + y_bin/10 + m.cos(m.radians(angle))*y_buck/2 + dx/2
                        - speed_vertical * time_vertical
                        + speed_vertical * (t - (time_vertical + time_horizontal)),
                    offset + z_percentage*z_bin/10
                )
            else:
                return taichi.Vector(
                    offset + x_percentage*x_bin/10 + m.sin(m.radians(angle))*y_buck/2
                        - speed_vertical * time_vertical
                        - speed * time_horizontal
                        - speed_vertical * time_vertical,
                    offset + y_bin/10 + m.cos(m.radians(angle))*y_buck/2 + dx/2
                        - speed_vertical * time_vertical
                        + speed_vertical * (t - (time_vertical + time_horizontal)),
                    offset + z_percentage*z_bin/10
                )

        elif motion == 2:
            # speed *= m.exp(-100 * dt)
            return taichi.Vector(
                offset + x_percentage*x_bin/10 + m.sin(m.radians(angle))*y_buck/2
                    - speed * min(max(t, 0), time_horizontal*0.95),
                offset + y_bin/10 + m.cos(m.radians(angle))*y_buck/2 + dx/2
                    - (2 * depth / time_horizontal) * m.sqrt(-m.pow(min(max(t, 0), time_horizontal*0.95), 2)
                        + min(max(t, 0), time_horizontal*0.95) * time_horizontal),
                offset + z_percentage * z_bin/10
            )

    mpm.add_particles(
        type='rigid',
        density=1e5,
        packing_fraction=1,  
        rotation_axis=(0, 0, 1),
        friction=friction_rb,
        scripted_position=taichi.function13(position_function),
        scripted_rotation=taichi.constant_function13(taichi.Vector(0, 0, 90-angle)),
        scale=(y_buck, x_buck, z_buck),
        codimensional=False,
        mesh_fn='projects/mpm/data/plate_houdini.obj',
    )

    # --
    mpm.simulate(
        clear_output_directory=True,
        print_profile_info=False,
    )