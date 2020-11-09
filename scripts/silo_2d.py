import taichi as tc
from math import *
from random import *


if __name__ == '__main__':
    r = 251#301
    dx = 1/r
    finalTime = 20  # sec
    frameRate = 60  # Hz
    oneG = 9.81
    micG = .02*oneG

    FrictionLS = -1

    offset = 0.2
    length_bin = 0.25
    height_bin = 0.25
    orifice = 0.025  #0.01
    length_soil = 0.985*length_bin
    height_soil = 0.125

    mpm = tc.dynamics.MPM(
        res=(r, r),
        # delta_x=.003,  # def: 1/res
        base_delta_t=.0001,
        frame_dt=1/frameRate,  # inverse frame rate (sec or 1/Hz)
        num_frames=finalTime*frameRate,
        num_threads=-1,
        gravity=(0, -oneG),
        particle_gravity=True,
        rigidBody_gravity=False,
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
        # write_particle_info=True,

        particle_bc_at_levelset=True, 
        # cdf_3d_modified=True,
        # compute_particle_impulses=True,
        # visualize_particle_impulses=True,
        # affect_particle_impulses=False,
    )

    # --------------------------------------------------------------------------
    # add level-set
    # problem region
    levelset = mpm.create_levelset()  # creates levelset with the dim of [0 to res] (dummy levelset)
    # bottom floor:
    levelset.add_plane(tc.Vector(0, +1), -offset)
    # level-set friction
    levelset.set_friction(FrictionLS)
    # dynamic level-set: False
    mpm.set_levelset(levelset, False)

    # --------------------------------------------------------------------------
    # Top floor
    mesh1 = tc.SegmentMesh()
    segment = ((0, 1),
               (length_bin, 1),)
    mesh1.add_segment(segment)
    mpm.add_particles(
        type='rigid',
        density=1e5,
        friction=-1,
        codimensional=True,
        scripted_position=tc.constant_function((offset+length_bin/2+orifice, offset+height_bin)),
        scripted_rotation=tc.constant_function((0, 0)),
        segment_mesh=mesh1)

    # Right wall
    mesh2 = tc.SegmentMesh()
    segment = ((1, 0),
               (1, .15),)
    mesh2.add_segment(segment)
    mpm.add_particles(
        type='rigid',
        density=1e5,
        friction=-1,
        codimensional=True,
        scripted_position=tc.constant_function((offset+length_bin, offset+height_bin+(height_soil)/2)),
        scripted_rotation=tc.constant_function((0, 0)),
        segment_mesh=mesh2)

    # Left wall
    mesh3 = tc.SegmentMesh()
    segment = ((1, 0),
               (1, .45),)
    mesh3.add_segment(segment)
    mpm.add_particles(
        type='rigid',
        density=1e5,
        friction=-2,
        codimensional=True,
        scripted_position=tc.constant_function((offset, offset + (height_bin+(height_soil))/2)),
        scripted_rotation=tc.constant_function((0, 0)),
        segment_mesh=mesh3)

    # --------------------------------------------------------------------------
    # add soil
    os1 = 0.2025
    os2 = 0.454
    tex = tc.Texture(
        'rect_bound',
        lower_bounds=(os1, os2, 0),
        upper_bounds=(os1+length_soil, os2+height_soil, 1),
        # translate=(offset, offset, 1),
    ) * 8

    # Lower this looser sand
    # Up to some value (1e5Pa) due to large Me_tr and stability issue!
    # Has issue from state 3 to state 4
    Y_mod = 1e6 #1e6
    nu = 0.45
    mpm.add_particles(
        type='nonlocal',
        pd=True,
        density_tex=tex.id,
        density=2450, #2583,
        critical_density=1700,
        packing_fraction=1,
        initial_velocity=(0, 0),

        S_mod=Y_mod/2/(1+nu),
        B_mod=Y_mod/3/(1-2*nu),
        # S_mod=5.0e4,#5.0e4,  # in FEM was 5.0e6
        # B_mod=2.0e4,#2.0e4,  # in FEM was 2.0e7

        A_mat=.48,#0.48,#0.58,#0.51
        dia=0.0003,

        # These 3 parameters are corrolated [local, mu_2-mu_s=0.2616] (but not from [PhD nonlocal])
        # mu_s should be larger than sqrt(3)*(1-(2*nu))/(1+nu)
        mu_s=.3819,  # .7000 = 35deg
        mu_2=1.5, #.6435,  # .9616 = 44deg
        I_0=0.278,

        # Larger this larger oscillations & denser sand
        t_0=1e-1, #1e-3 makes g go to infinity because of soil position
    )

    # --------------------------------------------------------------------------
    mpm.simulate(
        clear_output_directory=True,
        print_profile_info=True)
