# Wheel-soil model, by Amin Haeri [ahaeri92@gmail.com]

# Bin
    # Dimensions = 0.90 m long x 0.193 m wide x 0.35 m tall (2/3 soil)
# Wheel
    # Diameter = 0.300 m
    # Width = 0.125 m 
    # [Vertical load = 164 N] & [Volume = 0.00883572933 m3] => Density = F/V/g = 1890 kg/m3 in 1G
    # Commanded velocity = 0.02 m/s => Angular velocity = 0.13 rad/s 
    # Testing 20% slip and 70% slip (so horizontal velocities of 0.0156 m/s and 0.00585 m/s, respectively) 
# Soil 
    # Bulk density = 1730  [kg/m3] => Volume fraction = 0.67
    # Internal friction angle = 37.2
    # 71 kg of soil in the bin 

# MPM          EXP
# 1 u      ->  4 m
# 0.22500  <-  0.900   m (bin's length)  => (k = 4)
# 0.08750  <-  0.350   m (soil's height)
# 0.04750  <-  0.190   m (bin's width)
# 0.07500  <-  0.300   m (wheel's diameter)
# 0.03125  <-  0.125   m (wheel's width)
# 13.5     <-  27      s (forward time)
# ?  <-  4*(grid size)


import math as m
import taichi as tc


if __name__ == '__main__':
    r = 301 #251, 301, 401
    dx = 1/r

    dim_scale = 0.25
    time_scale = 0.5
    lin_velocity_scale = 0.5
    ang_velocity_scale = 2

    restTime = 2 * time_scale
    finalTime = 27 * time_scale + restTime
    frameRate = 60

    FrictionRB = 0.4
    FrictionLS = -1

    G = 9.81
    wheel_density = 164 / 0.00883572933 / G

    offset = 0.05

    x_bin = 0.900 * dim_scale
    y_bin = 0.350/2 * dim_scale
    z_bin = 0.190 * dim_scale #0.190 #0.35

    # Testing 20%, 40%, and 70% slip so horizontal velocities of 0.0156, 0.012, and 0.00585 m/s, respectively
    horizontal_velocity = 0.012 * lin_velocity_scale  # [m/s]
    angular_velocity = 180/3.14 * 0.1333 * ang_velocity_scale  # [deg/s] 
    wheel_diameter = 0.300 * dim_scale 
    wheel_width = 0.125 * dim_scale
    x_wheel = offset + x_bin * 0.8
    y_wheel = offset + y_bin + wheel_diameter/2 + dx
    z_wheel = offset + z_bin - wheel_width/2 + dx


    # --------------------------------------------------------------------------
    mpm = tc.dynamics.MPM(
        res=(r, r, r),
        base_delta_t=1e-4,#3.5e-4,#1e-4,  # can't be larger
        frame_dt=1/frameRate,
        num_frames=finalTime*frameRate,
        num_threads=-1,
        gravity=(0, -G, 0),
        particle_gravity=True,
        rigidBody_gravity=True,
        free_axis_in_position=2,  # 1: x, 2: y, 3: z
        rigid_body_collision=False,
        rigid_body_levelset_collision=False,  # apply impulse on rB if its particle's phi<0
        particle_collision=False,  # update particle pos and vel if its phi<0
        # if things do not separate
        # G2P, add normal-to-boundary velocity to boundary particle, def: 20000
        pushing_force=0,
        # due to numerical advection error (penetration penalty):
        # G2P, add normal-to-boundary velocity to boundary particle, typical: 1e3, def: 0
        # penalty=1e3,
        # sand_crawler=True,
        # rigid_penalty=0,  # def: 1e3, for rigidBody collision
        print_rigid_body_state=False,  # print pos, vel and angVel on terminal
        clean_boundary=True,  # def: true, clear boundary particles
        warn_particle_deletion=False,
        snapshots=False,
        verbose_bgeo=True,  # write other particles attributes, "visualize.cpp"
        write_particle_info=True,

        write_Houdini_and_rigidbody_info=True,
        cdf_3d_modified=True,
        compute_particle_impulses=True,

        visualize_particle_impulses=True,
        affect_particle_impulses=True,
        particle_bc_at_levelset=False, 
    )


    # Level-set ----------------------------------------------------------------
    levelset = mpm.create_levelset()
    levelset.add_plane(tc.Vector(1, 0, 0), -offset)
    levelset.add_plane(tc.Vector(0, 1, 0), -offset)
    levelset.add_plane(tc.Vector(0, 0, 1), -offset)
    levelset.add_plane(tc.Vector(-1, 0, 0), offset+x_bin)
    levelset.add_plane(tc.Vector(0, 0, -1), offset+z_bin)
    levelset.set_friction(FrictionLS)
    mpm.set_levelset(levelset, False)


    # Wheel -------------------------------------------------------------------
    def position_function(t):
        if t < restTime:
            return tc.Vector(x_wheel, y_wheel, z_wheel)
        else:
            return tc.Vector(x_wheel-horizontal_velocity*(t-restTime), y_wheel, z_wheel)
    
    def rotation_function(t):
        if t < restTime:
            return tc.Vector(0, 0, 0)
        else:
            return tc.Vector(0, 0, angular_velocity*(t-restTime))

    wheel = mpm.add_particles(
        type='rigid',
        linear_damping=1e3, #1e2, 1e3
        density=wheel_density,
        packing_fraction=1,  
        friction=FrictionRB,
        scripted_position=tc.function13(position_function),
        scripted_rotation=tc.function13(rotation_function),
        scale=(0.25, 0.25, 0.25),
        codimensional=False,
        mesh_fn='projects/mpm/data/wheel_houdini_closed.obj',
    )

    # Soil ---------------------------------------------------------------------
    tex = tc.Texture(
        'mesh',
        scale=tc.Vector(x_bin*10, y_bin*10, z_bin*10),
        translate=(offset+x_bin/2, offset+y_bin/2, offset+z_bin/2),
        resolution=(2*r, 2*r, 2*r),
        mesh_accuracy=3,
        filename='projects/mpm/data/cube_smooth.obj',
    ) * 4  # 8, 4

    grain_density = 2583
    packing_fraction = 0.67
    critical_density = packing_fraction * grain_density
    mpm.add_particles(
        type='nonlocal',
        pd=True,
        density_tex=tex.id,
        density=grain_density,  # verified
        critical_density=critical_density,
        packing_fraction=packing_fraction, ## verified

        S_mod=6.0e4,
        B_mod=1.0e5, #1.0e5
        A_mat=0.48,
        dia=0.0003*dim_scale,  # ?

        # These 3 parameters are corrolated [local, mu_2-mu_s=0.2616] (but not from [PhD nonlocal])
        # mu_s should be larger than sqrt(3)*(1-(2*nu))/(1+nu)
        mu_s=.7600, #.7000 ~ 35deg
        mu_2=.9616, #.9616 ~ 44deg
        I_0=0.278,

        # Larger this larger oscillations & denser sand
        t_0=1e-4, #1e-4, #5e-5, 
    )


    # --------------------------------------------------------------------------
    mpm.simulate(
        clear_output_directory=True,
        print_profile_info=True,
    )