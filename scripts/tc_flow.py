## Taylor-Couette Flow

import taichi as tc


if __name__ == '__main__':
    r = 251
    dx = 1/r
    Omega = 0.025  # rad/sec
    finalTime = 22  # sec
    frameRate = 50  # Hz
    Scale = 1 #0.75  # problem scale, has issue with 0.5!
    FrictionRB = -1  # rigidBody friction
    FrictionFloor = 0 # smooth floor
    FrictionLS = -1  # frictional contact makes motion!
    oneG = 9.81
    micG = .02*oneG

    # --------------------------------------------------------------------------
    mpm = tc.dynamics.MPM(
        res=(r, r, r),
        # delta_x=.003,  # def: 1/res
        base_delta_t=1e-4, #1e-4,#5e-5,
        frame_dt=1/frameRate,  # inverse frame rate (sec or 1/Hz)
        num_frames=finalTime*frameRate,
        num_threads=-1,
        gravity=(0., -oneG, 0.),
        particle_gravity=True,
        rigidBody_gravity=False,
        rigid_body_collision=False,
        rigid_body_levelset_collision=False,  # apply impulse on rB if its particle's phi<0
        particle_collision=False,  # update particle pos and vel if its phi<0
        pushing_force=0,  # G2P, add normal-to-boundary velocity to boundary particle, "transfer.cpp"
        # due to numerical advection error (if zero causes leakage):
        penalty=0,  # 1e4,  # G2P, add normal-to-boundary velocity to boundary particle, typical: 1e3, "transfer.cpp"
        # rigid_penalty=0,     # def: 1e3, for rigidBody collision
        print_rigid_body_state=False,  # print pos, vel and angVel on terminal
        dirichlet_boundary_radius=0,     # 1/0 means dirichletBoundary is ON/OFF (apply velocity!)
        clean_boundary=False,  # def: true, clear boundary particles
        warn_particle_deletion=True,

        particle_bc_at_levelset=True, 
        cdf_3d_modified=True,
        # compute_particle_impulses=True,
        # visualize_particle_impulses=True,
        # affect_particle_impulses=False,

        snapshots=False,
        verbose_bgeo=True,  # write other particles attributes, "visualize.cpp"
        write_particle=True,
        write_rigid_body=True,
        write_partio=True,
        write_dataset=False,
    )

    # --------------------------------------------------------------------------
    # add level-set
    # problem region
    levelset = mpm.create_levelset()  # creates levelset with the dim of [0 to res] (dummy levelset)
    # inner cylinder:
    # levelset.add_cylinder(tc.Vector(.3, 0, .3), .0998*Scale, False)  # 0 means it's not important 8
    # outer cylinder:
    levelset.add_cylinder(tc.Vector(.3, 0, .3), .2001*Scale, True)
    # levelset.add_cylinder(tc.Vector(.3, 0, .3), .22*Scale, True)
    # floor:
    # levelset.add_plane(tc.Vector(0, +1, 0), -(.15-(.051*Scale)))
    # top:
    # levelset.add_plane(tc.Vector(0, -1, 0), +(.15+(.050*Scale)))
    # level-set friction
    levelset.set_friction(FrictionLS)
    # dynamic level-set: False
    mpm.set_levelset(levelset, False)

    # --------------------------------------------------------------------------
    # add soil
    tex = tc.Texture(
        'mesh',
        resolution=(2*r, 2*r, 2*r),
        mesh_accuracy=3,
        translate=(.3, .15, .3),
        # translate=(.3, .17, .3),
        scale=(.2*Scale, .1*Scale, .2*Scale),
        filename='projects/mpm/data/tcSoilFine.obj'
    ) * 4 #8  # number of particles per cell (ppc, max sample density)

    ## Option 1:
    E_mod = 1e6 #1e6, 15e4
    nu = 0.3
    ## Option 2:
    # Here E is 3e6 Pa and nu is 0.45
    S_mod=1e6
    B_mod=10e6

    grain_density = 2550
    packing_fraction = 0.645 #0.645, 1
    critical_density = packing_fraction * grain_density

    ## Check 1: dt < (1/res^2)*t0/2/A^2/d^2
    ## Check 2: dt < (1/res)*sqrt(rho/E)
    ## Check 3: mu_s > sqrt(3)*(1-(2*nu))/(1+nu)

    mpm.add_particles(
        type='nonlocal',
        pd=True,
        density_tex=tex.id,
        density=grain_density,
        critical_density=critical_density,
        packing_fraction=packing_fraction,
        initial_velocity=(0, 0, 0),

        # S_mod=E_mod/2/(1+nu),
        # B_mod=E_mod/3/(1-2*nu),
        S_mod=S_mod,  # in FEM was 5.0e6
        B_mod=B_mod,  # in FEM was 2.0e7

        A_mat=0.48,  # 0.48 | 0.58 | 0.51
        dia=0.004,  # 0.004, 0.0003

        # These 3 parameters are corrolated [local, mu_2-mu_s=0.2616] (but not from [PhD nonlocal])
        # mu_s should be larger than sqrt(3)*(1-(2*nu))/(1+nu)
        # mu_s=0.7,
        # mu_2=0.9616,
        mu_s=0.3819,
        mu_2=0.6435,
        I_0=0.278,

        # Larger this larger oscillations & denser sand
        t_0=1e-4, #1e-4, 1e-3
    )

    # --------------------------------------------------------------------------
    # Two rigid bodies can't be close to each other
    # inner cylinder
    # def rotation_function(t):
    #     return tc.Vector(0, Omega*t, 0)

    inner = mpm.add_particles(
        type='rigid',
        density=400,
        scale=(.0999*2*Scale, .15*Scale, .0999*2*Scale),  # ri=coef/2
        # scale=(.08*2*Scale, .15*Scale, .08*2*Scale),  # ri=coef/2
        # scripted_rotation=tc.function13(rotation_function),
        scripted_position=tc.constant_function13(tc.Vector(.3, .15, .3)),
        rotation_axis=(0, 1, 0),
        initial_angular_velocity=Omega,
        friction=FrictionRB,
        codimensional=False,  # Thin shell has leakage issue
        mesh_fn='projects/mpm/data/cylinder_houdini.obj',
    )
    innerAux = mpm.add_particles(
        type='rigid',
        scale=(.0019*Scale, .001*Scale, .0019*Scale),
        scripted_position=tc.constant_function13(tc.Vector(.3, .15, .3)),
        scripted_rotation=tc.constant_function13(tc.Vector(0, 0, 0)),
        rotation_axis=(0, 1, 0),
        friction=0.0001,
        codimensional=True,
        mesh_fn='projects/mpm/data/cylinder_houdini.obj',
    )
    mpm.add_articulation(
        type='stepper',
        obj0=inner,
        obj1=innerAux,
        angular_velocity=Omega,
        axis=(0, 1, 0),
    )

    # --------------------------------------------------------------------------
    # Confining plate in Micro-G only
    # def position_function(t):
    #     time = 0.1
    #     if (t <= time):
    #         return tc.Vector(0.3, (.048*Scale-.051*Scale)/time*t+(.15+(.051*Scale)), 0.3)
    #     else:
    #         return tc.Vector(0.3, .15+(.048*Scale), 0.3)
    # confPlate = mpm.add_particles(
    #     type='rigid',
    #     density=400,
    #     scale=(1*Scale, 1*Scale, 1*Scale),
    #     scripted_position=tc.function13(position_function),
    #     scripted_rotation=tc.constant_function13(tc.Vector(0, 0, 0)),
    #     friction=-2,
    #     codimensional=True,
    #     mesh_fn='projects/mpm/data/conf_plate.obj'
    # )

    # --------------------------------------------------------------------------
    # change in rotation direction if required
    # def frame_update(t, frame_dt):
    #     if t < (finalTime/2)-dt or t > (finalTime/2)+dt:
    #         return
    #     mpm.add_articulation(
    #         type='stepper',
    #         obj0=inner,
    #         obj1=innerAux,
    #         angular_velocity=-Omega,
    #         axis=(0, 1, 0),
    #     )

    # --------------------------------------------------------------------------
    # ## outer cylinder
    # outer = mpm.add_particles(
    #         type              = 'rigid',
    #         density           = 100000,
    #         scale             = (.4000*Scale, .1020*Scale, .4000*Scale),
    #         scripted_rotation = tc.constant_function13(tc.Vector(0, 0, 0)),
    #         scripted_position = tc.constant_function13(tc.Vector(.3, .15, .3)),
    #         friction          = Friction,
    #         codimensional     = False,
    #         mesh_fn           = 'projects/mpm/data/cylinderInner.obj'
    # )
    # --------------------------------------------------------------------------
    ## floor
    floor = mpm.add_particles(
            type              = 'rigid',
            density           = 100000,
            scripted_rotation = tc.constant_function13(tc.Vector(0, 0, 0)),
            scale             = (2.15*Scale, .1*Scale, 2.15*Scale),
            scripted_position = tc.constant_function13(tc.Vector(.3, (.15-(.052*Scale)), .3)),
            friction          = FrictionFloor,
            codimensional     = False,
            mesh_fn           = 'projects/mpm/data/cylinder_jet.obj'
    )
    # --------------------------------------------------------------------------
    # ## ceil
    # ceil = mpm.add_particles(
    #         type              = 'rigid',
    #         density           = 100000,
    #         scripted_rotation = tc.constant_function13(tc.Vector(180, 0, 0)),
    #         scale             = (2.15*Scale, .1*Scale, 2.15*Scale),
    #         scripted_position = tc.constant_function13(tc.Vector(.3, (.15+(.05*Scale)), .3)),
    #         friction          = Friction,
    #         codimensional     = False,
    #         mesh_fn           = 'projects/mpm/data/cylinder_jet.obj'
    # )

    # --------------------------------------------------------------------------
    mpm.simulate(
        clear_output_directory=True,
        print_profile_info=True,
        # frame_update=frame_update,
    )
