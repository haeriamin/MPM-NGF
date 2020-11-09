import taichi as tc


if __name__ == '__main__':
    r = 301  # should be large enough (due to accuracy of rigidbody-particle interaction)
    dx = 1/r
    frameRate = 60
    finalTime = 3

    FrictionRB = 0.75
    FrictionLS = 1.5
    oneG = 9.81

    offset = 30*dx
    x_bin = 5.0
    y_bin = 1.25
    z_bin = 2.5
    bin_size = tc.Vector(x_bin, y_bin, z_bin)

    # --------------------------------------------------------------------------
    mpm = tc.dynamics.MPM(
        res=(r, r, r),
        base_delta_t=1e-4,
        frame_dt=1/frameRate,
        num_frames=finalTime*frameRate,
        num_threads=-1,
        gravity=(0, -oneG, 0),
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
        # penalty=1e1,
        # sand_crawler=True,
        # rigid_penalty=0,  # def: 1e3, for rigidBody collision
        print_rigid_body_state=False,  # print pos, vel and angVel on terminal
        clean_boundary=True,  # def: true, clear boundary particles
        warn_particle_deletion=False,
        verbose_bgeo=True,  # write other particles attributes, "visualize.cpp"
        Houdini=True,
        snapshots=False,
        cdf_3d_modified=True,
        visualize_particle_impulses=True,
    )

    # level-set ----------------------------------------------------------------
    levelset = mpm.create_levelset()
    levelset.add_plane(tc.Vector(0, 1, 0), -offset)
    levelset.set_friction(FrictionLS)
    mpm.set_levelset(levelset, False)

    # fixed soil ---------------------------------------------------------------
    tex = tc.Texture(
        'mesh',
        scale=tc.Vector(x_bin/10, 50*dx, z_bin/8),
        translate=(150*dx, offset+3*dx, 150*dx),
        resolution=(2*r, 2*r, 2*r),
        mesh_accuracy=3,
        filename='projects/mpm/data/cube_smooth.obj',
    ) * 8

    mpm.add_particles(
       type='nonlocal',
        pd=True,
        density_tex=tex.id,
        density=2583,
        packing_fraction=0.67,

        S_mod=6.0e4, # in FEM was 5.0e6
        B_mod=1.0e5, # in FEM was 2.0e7
        A_mat=0.48,
        dia=0.0003,

        # These 3 parameters are corrolated [local, mu_2-mu_s=0.2616] (but not from [PhD nonlocal])
        # mu_s should be larger than sqrt(3)*(1-(2*nu))/(1+nu)
        mu_s=.7000, #35deg  # 0.3819 is 20 deg
        mu_2=.9616, #44deg  # 0.6435 is 32 deg
        I_0=0.278,

        # Larger this larger oscillations & denser sand
        t_0=1e-3,#5e-5,#1e-3, 
    )

    # plate -----------------------------------------------------------------
    def position_function(t):
        speed = 10*dx  # m/sec
        vert_time = 0.4
        if (t <= vert_time):
            return tc.Vector(160*dx-speed*t, offset+11*dx-speed*t, 150*dx)
        elif (t >= 3-vert_time):
            return tc.Vector(160*dx-speed*t, offset+11*dx-speed*vert_time+speed*(t-(3-vert_time)), 150*dx)
        else:
            return tc.Vector(160*dx-speed*t, offset+11*dx-speed*vert_time, 150*dx)

    main=mpm.add_particles(
        type='rigid',
        density=1e5,
        packing_fraction=1,  
        scale=(10*dx, 0.0050, 30*dx),  # thickness should be larger than dx
        rotation_axis=(0, 0, 1),
        friction=FrictionRB,
        scripted_position=tc.function13(position_function),
        scripted_rotation=tc.constant_function13(tc.Vector(0, 0, 90-0)),#3.8
        codimensional=False,
        mesh_fn='projects/mpm/data/plate_houdini.obj',
        # mesh_fn='projects/mpm/data/cube.obj',
        # mesh_fn='projects/mpm/data/plate_round.obj',
        # reverse_vertices=True,
    )

    # --------------------------------------------------------------------------
    mpm.simulate(
        clear_output_directory=True,
        print_profile_info=True,
    )
