# Nonlocal test, by Amin Haeri [ahaeri92@gmail.com]

import taichi as tc

if __name__ == '__main__':
    r = 301
    dx = 1/r
    finalTime = 1  # sec
    frameRate = 60  # Hz
    Scale = 1  # problem scale
    oneG = 9.81
    micG = .02*oneG

    FrictionRB = 0.3
    FrictionLS = 0.3

    Y_mod = 100e6  # lower this looser sand; up to some value due to large Me_tr! (has issue from state 3 to state 4)
    nu = 0.45

    # --------------------------------------------------------------------------
    mpm = tc.dynamics.MPM(
        res=(r, r, r),
        # delta_x=.003,  # def: 1/res
        base_delta_t=.0001,
        frame_dt=1/frameRate,  # inverse frame rate (sec or 1/Hz)
        num_frames=finalTime*frameRate,
        num_threads=-1,
        gravity=(0, -oneG, 0),
        particle_gravity=True,
        rigidBody_gravity=False,
        rigid_body_collision=True,
        rigid_body_levelset_collision=False,  # apply impulse on rB if its particle's phi<0
        particle_collision=False,  # update particle pos and vel if its phi<0
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
    # floor:
    levelset.add_plane(tc.Vector(0, +1, 0), -(.1499-(.05*Scale)))
    # level-set friction
    levelset.set_friction(FrictionLS)
    # dynamic level-set: False
    mpm.set_levelset(levelset, False)

    # --------------------------------------------------------------------------
    # def position_function(t):
    #     speed = 10*dx  # m/sec
    #     return tc.Vector(0.5, 0.1+speed*t, 0.5)

    # mpm.add_particles(
    #     type='rigid',
    #     density=400,#1e5,
    #     packing_fraction=1,  
    #     rotation_axis=(0, 0, 1),
    #     friction=FrictionRB,
    #     scripted_position=tc.function13(position_function),#tc.constant_function13(tc.Vector(0.5, 0.1, 0.5)),

    #     scale=(0.2, 0.0070, 0.2),  # thickness should be at least equal to 2dx
    #     scripted_rotation=tc.constant_function13(tc.Vector(0, 0, 0)),
    #     codimensional=False,
    #     mesh_fn='projects/mpm/data/plate_houdini.obj',
    # )

    # --------------------------------------------------------------------------
    # add soil
    tex = tc.Texture(
        'mesh',
        scale=tc.Vector(0.4, 0.8, 0.4),
        translate=(0.5, 0.20, 0.5), #0.15
        resolution=(2*r, 2*r, 2*r),
        mesh_accuracy=3,
        filename='projects/mpm/data/cube_smooth.obj',
    ) * 8

    mpm.add_particles(
        type='nonlocal',
        pd=True,
        density_tex=tex.id,
        density=2450,
        packing_fraction=1,#.79,
        initial_velocity=(0, 0, 0),

        S_mod=Y_mod/2/(1+nu),
        B_mod=Y_mod/3/(1-2*nu),
        # S_mod=5.0e4,#5.0e4,  # in FEM was 5.0e6
        # B_mod=2.0e4,#2.0e4,  # in FEM was 2.0e7

        A_mat=.48,#0.48,#0.58,#0.51
        dia=0.005,

        # These 3 parameters are corrolated [local, mu_2-mu_s=0.2616] (but not from [PhD nonlocal])
        # mu_s should be larger than sqrt(3)*(1-(2*nu))/(1+nu)
        mu_s=.3819,  # .7000 = 35deg
        mu_2=.6435,  # .9616 = 44deg
        I_0=0.278,

        # Larger this larger oscillations & denser sand
        t_0=1e-1, #1e-3 makes g go to infinity because of soil position
    )

    # --------------------------------------------------------------------------
    mpm.simulate(
        clear_output_directory=True,
        print_profile_info=True,
    )
