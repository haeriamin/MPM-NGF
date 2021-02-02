## Industrial excavator model

# MPM         EXP
# 1 u     ->  4 m
# 0.5000  <-  2.00  m (bin's length)  => (k = 4)
# 0.0375  <-  0.15  m (soil's height)
# 0.2450  <-  1.00  m (bin's width)
# 0.0100  <-  0.04  m/s (excavator's forward velocity)
# 0.0033  <-  4*(grid size)

import math as m
import taichi as tc


if __name__ == '__main__':
    r = 301  # can't be smaller | scale the problem instead
    dx = 1/r
    Scale = 1
    time_scale = 0.5 * m.sqrt(Scale)
    velocity_scale = 0.5 * m.sqrt(Scale)
    verticalTime = 2.5*time_scale
    forwardTime = 8.5*time_scale
    stopTime = 0*time_scale
    finalTime = forwardTime + 2*verticalTime + 2*stopTime
    frameRate = 60

    FrictionRB = 0.1 #-1
    FrictionLS = -1

    G = 9.81

    offset = 0.05/2
    x_bin = 2.5 * Scale  #5.
    y_bin = 0.5 * Scale  #0.414
    z_bin = 1.0 * Scale  #2.5
    bin_size = tc.Vector(x_bin, y_bin, z_bin)

    # --------------------------------------------------------------------------
    mpm = tc.dynamics.MPM(
        res=(r, r, r),
        base_delta_t=1e-4,  # can't be larger
        frame_dt=1/frameRate,
        num_frames=finalTime*frameRate,
        num_threads=-1,
        gravity=(0, -G, 0),
        particle_gravity=True,
        rigidBody_gravity=False,
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
        verbose_bgeo=True,  # write other particles attributes, "visualize.cpp"
        Houdini=True,
        snapshots=False,

        particle_bc_at_levelset=False, 
        cdf_3d_modified=True,
        compute_particle_impulses=True,
        visualize_particle_impulses=True,
        affect_particle_impulses=False,
    )

    # level-set ----------------------------------------------------------------
    levelset = mpm.create_levelset()
    levelset.add_plane(tc.Vector(1, 0, 0), -offset)
    levelset.add_plane(tc.Vector(0, 1, 0), -offset)
    levelset.add_plane(tc.Vector(0, 0, 1), -offset)
    levelset.add_plane(tc.Vector(-1, 0, 0), offset+x_bin/10*Scale)
    levelset.add_plane(tc.Vector(0, 0, -1), offset+z_bin/10*Scale)
    levelset.set_friction(FrictionLS)
    mpm.set_levelset(levelset, False)

    # soil ---------------------------------------------------------------------
    tex = tc.Texture(
        'mesh',
        scale=bin_size,
        translate=(offset+x_bin/20, offset+y_bin/20, offset+z_bin/20),
        resolution=(2*r, 2*r, 2*r),
        mesh_accuracy=3,
        filename='projects/mpm/data/cube_smooth.obj',
    ) * 8

    mpm.add_particles(
        type='nonlocal',
        pd=True,
        density_tex=tex.id,
        density=2583,  # verified
        packing_fraction=0.6, #0.67  # verified

        S_mod=6.0e4,
        B_mod=1.0e5,
        A_mat=0.48,
        dia=0.0003,  # verified

        # These 3 parameters are corrolated [local, mu_2-mu_s=0.2616] (but not from [PhD nonlocal])
        # mu_s should be larger than sqrt(3)*(1-(2*nu))/(1+nu)
        mu_s=.7000, #35deg
        mu_2=.9616, #44deg
        I_0=0.278,

        # Larger this larger oscillations & denser sand
        t_0=1e-4, #1e-4, #5e-5, 
    )

    # bucket -------------------------------------------------------------------
    def position_function(t):
        speed = 0.04*velocity_scale  # m/sec
        if (t <= verticalTime):
            return tc.Vector(offset+0.975*x_bin/10-speed*(t),
                             offset+1.9*y_bin/10-speed*(t),
                             offset+0.5*z_bin/10)
        elif (t > finalTime-(verticalTime+1)):
            return tc.Vector(offset+0.975*x_bin/10-speed*(t),
                             offset+1.9*y_bin/10-speed*(verticalTime)+speed*(t-(finalTime-(verticalTime+1))),
                             offset+0.5*z_bin/10)
        else:
            return tc.Vector(offset+0.975*x_bin/10-speed*(t),
                             offset+1.9*y_bin/10-speed*(verticalTime),
                             offset+0.5*z_bin/10)

    def rotation_function(t):
        return tc.Vector(0, 0, (30-(-60)) / (-finalTime) * t + 30)

    excav = mpm.add_particles(
        type='rigid',
        # scale=(0.25*Scale, 0.25*Scale, 0.25*Scale),
        scale=(0.2*Scale, 0.2*Scale, 0.2*Scale),
        rotation_axis=(0, 0, 1),
        scripted_position=tc.function13(position_function),
        scripted_rotation=tc.function13(rotation_function),
        density=1e7,
        codimensional=False,  # not thin shell
        angular_damping=0,
        friction=FrictionRB,
        # mesh_fn='projects/mpm/data/bucket_fine.obj',
        mesh_fn='projects/mpm/data/bucket.obj',
    )

    # # bucket -------------------------------------------------------------------
    # def position_function(t):
    #     offset_coeff = 1.07 #1.06 #0.9
    #     depth = 4.7625*dx*Scale #=.0125/dx+1 #0
    #     if (t <= verticalTime):
    #         speed = (0.0250*velocity_scale) # m/sec
    #         return tc.Vector(offset+offset_coeff*x_bin/10
    #                             # -speed*(t),
    #                             -speed*(0),
    #                          offset+(1.925*y_bin/10-depth)  #2.04 for 3.8 deg
    #                             -speed*(t),
    #                          offset+0.50*z_bin/10)
    #     elif (t > verticalTime and t <= verticalTime+stopTime):
    #         return tc.Vector(offset+offset_coeff*x_bin/10
    #                             -((0.0250*velocity_scale))*(verticalTime),
    #                          offset+(1.925*y_bin/10-depth)
    #                             -((0.0250*velocity_scale))*(verticalTime),
    #                          offset+0.50*z_bin/10)
    #     elif (t > verticalTime+stopTime and t <= verticalTime+stopTime+forwardTime):
    #         speed = 0.04*velocity_scale  # m/sec
    #         return tc.Vector(offset+offset_coeff*x_bin/10
    #                             # -((0.0250*velocity_scale))*(verticalTime)
    #                             -((0*velocity_scale))*(verticalTime)
    #                             -speed*(t-(verticalTime+stopTime)),
    #                         offset+(1.925*y_bin/10-depth)
    #                             -((0.0250*velocity_scale))*(verticalTime),
    #                         offset+0.50*z_bin/10)
    #     elif (t > verticalTime+stopTime+forwardTime and t <= verticalTime+stopTime+forwardTime+stopTime):
    #         return tc.Vector(offset+offset_coeff*x_bin/10
    #                             # -((0.0250*velocity_scale))*(verticalTime)
    #                             -((0*velocity_scale))*(verticalTime)
    #                             -(0.04*velocity_scale)*(forwardTime),
    #                          offset+(1.925*y_bin/10-depth)
    #                             -((0.0250*velocity_scale))*(verticalTime),
    #                          offset+0.50*z_bin/10)
    #     elif (t > verticalTime+stopTime+forwardTime+stopTime):
    #         speed = (0.0250*velocity_scale)  # m/sec
    #         return tc.Vector(offset+offset_coeff*x_bin/10
    #                             # -((0.0250*velocity_scale))*(verticalTime)
    #                             -((0*velocity_scale))*(verticalTime)
    #                             -(0.04*velocity_scale)*(forwardTime)
    #                             -speed*(t-(verticalTime+stopTime+forwardTime+stopTime)),
    #                          offset+(1.925*y_bin/10-depth)
    #                             -((0.0250*velocity_scale))*(verticalTime)
    #                             +speed*(t-(verticalTime+stopTime+forwardTime+stopTime)),
    #                          offset+0.50*z_bin/10)

    # mpm.add_particles(
    #     type='rigid',
    #     density=1e5,
    #     packing_fraction=1,  
    #     rotation_axis=(0, 0, 1),
    #     friction=FrictionRB,
    #     scripted_position=tc.function13(position_function),
    #     scripted_rotation=tc.function13(rotation_function),
    #     scale=(0.0763*Scale, 0.0070, 0.1142*Scale),
    #     codimensional=False,
    #     mesh_fn='projects/mpm/data/plate_houdini.obj',
    # )

    # --------------------------------------------------------------------------
    mpm.simulate(
        clear_output_directory=True,
        print_profile_info=True,
    )
