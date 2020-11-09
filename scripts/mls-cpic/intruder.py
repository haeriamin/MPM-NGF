import taichi as tc
from math import *
from random import *


if __name__ == '__main__':
    r = 301 #301 #251
    dx = 1/r
    finalTime = 0.5  # sec
    frameRate = 40000  # Hz
    Scale = 1  # problem scale
    oneG = 9.81
    micG = .02*oneG

    FrictionLS = -1

    offset = 0.2
    x_bin = 1 * Scale
    y_bin = 1 * Scale
    z_bin = 1 * Scale
    bin_size = tc.Vector(x_bin, y_bin, z_bin)

    mpm = tc.dynamics.MPM(
        res=(r, r, r),
        # delta_x=.003,  # def: 1/res
        base_delta_t=1e-6,
        frame_dt=1/frameRate,  # inverse frame rate (sec or 1/Hz)
        num_frames=finalTime*frameRate,
        num_threads=-1,
        gravity=(0, -oneG, 0),
        particle_gravity=True,
        rigidBody_gravity=True,
        rigid_body_collision=True,
        rigid_body_levelset_collision=False,  # apply impulse on rB if its particle's phi<0
        particle_collision=True,  # update particle pos and vel if its phi<0
        pushing_force=0,  # G2P, add normal-to-boundary velocity to boundary particle, "transfer.cpp"
        # due to numerical advection error (DOES NOT AFFECT NONLOCAL):
        penalty=0,  # 1e4,  # G2P, add normal-to-boundary velocity to boundary particle, typical: 1e3, "transfer.cpp"
        # rigid_penalty=0,     # def: 1e3, for rigidBody collision
        print_rigid_body_state=False,  # print pos, vel and angVel on terminal
        clean_boundary=True,  # def: true, clear boundary particles
        warn_particle_deletion=True,
        verbose_bgeo=True,  # write other particles attributes, "visualize.cpp"
        snapshots=False,
        write_Houdini_and_rigidbody_info=True,
        write_particle_info=True,

        particle_bc_at_levelset=True, 
        cdf_3d_modified=True,
        compute_particle_impulses=True,
        visualize_particle_impulses=True,
        affect_particle_impulses=True,
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

    # Lower this looser sand
    # Up to some value (1e5Pa) due to large Me_tr!
    Y_mod = 15e4 #15e4 #1e6
    nu = 0.3
    mpm.add_particles(
        type='nonlocal',
        pd=True,
        density_tex=tex.id,
        density=2583, #2450
        critical_density=2583,
        packing_fraction=1,
        initial_velocity=(0, 0, 0),

        S_mod=Y_mod/2/(1+nu),
        B_mod=Y_mod/3/(1-2*nu),
        # S_mod=5.0e4,#5.0e4,  # in FEM was 5.0e6
        # B_mod=2.0e4,#2.0e4,  # in FEM was 2.0e7

        A_mat=.48,#0.48,#0.58,#0.51
        dia=0.0004*Scale,

        # These 3 parameters are corrolated [local, mu_2-mu_s=0.2616] (but not from [PhD nonlocal])
        # mu_s should be larger than sqrt(3)*(1-(2*nu))/(1+nu)
        mu_s=.3819,  # .7000 = 35deg
        mu_2=1.5, #.6435,  # .9616 = 44deg
        I_0=0.278,

        # Larger this larger oscillations & denser sand
        t_0=1e-4, #1e-1 #1e-3 #1e-4
    )


    # intruder -------------------------------------------------------------------
    mpm.add_particles(
        type='rigid',
        density=8910,
        packing_fraction=1,  
        rotation_axis=(1, 1, 1),
        friction=0,
        initial_position=(offset+x_bin/10/2, offset+y_bin/10/2+0.1, offset+z_bin/10/2),
        initial_velocity=(0, -5, 0),
        scale=(0.01*Scale, 0.01*Scale, 0.01*Scale),
        codimensional=False,
        mesh_fn='projects/mpm/data/sphere_small.obj',
    )

    # --------------------------------------------------------------------------
    mpm.simulate(
        clear_output_directory=True,
        print_profile_info=True,
    )
