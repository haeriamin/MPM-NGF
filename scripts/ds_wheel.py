## Wheel-soil model prepared for generating training dataset
## Amin Haeri [ahaeri92@gmail.com]
## Mid 2021

# MPM          EXP
# 1 u      ->  4 m
# 0.03125  <-  0.125 m (wheel's width)

#!/usr/bin/env python

import sys
import math

import taichi


if __name__ == '__main__':

    # Dimensional analysis
    scale_dim = 0.25
    scale_time = scale_dim * 2
    scale_lin_vel = scale_dim * 2
    scale_ang_vel = scale_lin_vel / scale_dim
    scale_force = scale_dim**3

    # Get inputs
    c = int(sys.argv[1])
    gravity = float(sys.argv[2])
    slip = float(sys.argv[3])
    wload = float(sys.argv[4])
    wdia = float(sys.argv[5])/100
    sfangle = float(sys.argv[6])
    print(c, '-', gravity, '-', slip, '-', wload, '-', wdia, '-', sfangle)

    # Discretization
    r = 301  # Spatial resolution
    dx = 1/r
    dt = 1e-5
    ppc = 4

    # Data acquisition rate
    frame_rate = 60 / scale_time  # [Hz]

    # Bin
    friction_ls = -1
    offset = 0.05/2
    x_bin = 0.400 * scale_dim
    y_bin = 0.350/6 * scale_dim
    z_bin = 0.190 * scale_dim

    # Soil
    modulus_elastic = 15e6 * scale_dim
    ratio_poisson = 0.3
    modulus_shear = modulus_elastic / 2 / (1 + ratio_poisson)
    modulus_bulk = modulus_elastic / 3 / (1 - 2 * ratio_poisson)
    friction_static = math.tan(math.pi/180 * sfangle)
    friction_2 = 0.9616
    grain_diameter = 0.0003 * scale_dim
    grain_density = 2583
    packing_fraction = 0.67
    critical_density = packing_fraction * grain_density
    
    # Wheel
    wwidth = 0.125  # [m]
    angular_velocity = 0.13

    if wdia == 0.15:
        wwidth = wwidth * 2/3
        angular_velocity = angular_velocity * 2
    elif wdia == 0.05:
        wwidth = wwidth * 2/3
        wload = wload / 9
        angular_velocity = angular_velocity * 6

    wheel_volume = math.pi * (wdia/2)**2 * wwidth  # [m3]
    wheel_density = wload / wheel_volume / gravity  # = F/V/G
    friction_rb = -5

    if slip == 20:
        horizontal_velocity = 0.01560 * scale_lin_vel  # [m/s]
    elif slip == 40:
        horizontal_velocity = 0.01200 * scale_lin_vel  # [m/s]
    elif slip == 70:
        horizontal_velocity = 0.00585 * scale_lin_vel  # [m/s]

    angular_velocity = 180/math.pi * angular_velocity * scale_ang_vel  # [deg/s]

    wheel_diameter = wdia * scale_dim 
    x_wheel = offset + x_bin * 0.7
    y_wheel = offset + y_bin + wheel_diameter/2 + dx
    z_wheel = offset + z_bin * 0.5
    scale_x_wheel = scale_dim * (wdia/0.30)
    scale_y_wheel = scale_dim * (wdia/0.30)
    scale_z_wheel = scale_dim * (wwidth/0.125)
    time_rest = 4 * scale_time  # 2 sec for OneG and 12 sec for LunarG
    time_forward = 10 * scale_time
    time_total =  time_rest + time_forward

    # Start MPM
    mpm = taichi.dynamics.MPM(
        task_id=str(c)+'_G'+str(gravity)+'ms-2_S'+str(int(slip))+
            'perc_L'+str(int(wload))+'N_D'+str(int(wdia*100))+'cm_A'+str(int(sfangle))+'deg',
        res=(r, r, r),
        base_delta_t=dt,
        frame_dt=1/frame_rate,
        num_frames=time_total*frame_rate,
        num_threads=-1,
        gravity=(0, -gravity, 0),
        particle_gravity=True,
        rigidBody_gravity=True,
        free_axis_in_position=2,  # 1: x, 2: y, 3: z
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
        affect_particle_impulses=True,
        particle_bc_at_levelset=False, 

        snapshots=False,
        verbose_bgeo=False,
        write_particle=False,
        write_rigid_body=False,  # True for test
        write_partio=False, # True for test
        write_dataset=True,
    )

    # level-set --
    levelset = mpm.create_levelset()
    levelset.add_plane(taichi.Vector(1, 0, 0), -offset)
    levelset.add_plane(taichi.Vector(0, 1, 0), -offset)
    levelset.add_plane(taichi.Vector(0, 0, 1), -offset)
    levelset.add_plane(taichi.Vector(-1, 0, 0), offset+x_bin)
    levelset.add_plane(taichi.Vector(0, 0, -1), offset+z_bin)
    levelset.set_friction(friction_ls)
    mpm.set_levelset(levelset, False)

    # soil --
    tex = taichi.Texture(
        'mesh',
        scale=taichi.Vector(x_bin*10, y_bin*10, z_bin*10),
        translate=(offset+x_bin/2, offset+y_bin/2, offset+z_bin/2),
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

    # Wheel --
    def position_function(t):
        if t <= time_rest:
            return taichi.Vector(x_wheel, y_wheel, z_wheel)
        else:
            return taichi.Vector(x_wheel-horizontal_velocity*(t-time_rest), y_wheel, z_wheel)
    
    def rotation_function(t):
        if t <= time_rest:
            return taichi.Vector(0, 0, 0)
        else:
            return taichi.Vector(0, 0, angular_velocity*(t-time_rest))

    mpm.add_particles(
        type='rigid',
        linear_damping=1e3, #2e3, 1e-1/dt
        density=wheel_density,
        packing_fraction=1,  
        friction=friction_rb,
        scripted_position=taichi.function13(position_function),
        scripted_rotation=taichi.function13(rotation_function),
        scale=(scale_x_wheel, scale_y_wheel, scale_z_wheel),
        codimensional=False,
        mesh_fn='projects/mpm/data/wheel_houdini.obj',
        # mesh_fn='projects/mpm/data/wheel_houdini_closed.obj',
    )

    # --
    mpm.simulate(
        clear_output_directory=True,
        print_profile_info=False,
    )