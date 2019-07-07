# Taylor-Couette Flow, by Amin Haeri [email](ahaeri92@gmail.com)
import taichi as tc
from math import *

if __name__ == '__main__':
    r           = 301
    dx          = 1/r
    Omega       = .025 # rad/sec
    finalTime   = 22   # sec
    frameRate   = 50   # Hz
    Scale       = 1    # problem scale
    # if not -1, causes segmentation fault:
    FrictionRB  = -1  # rigidBody friction (not particle's)
    FrictionLS  = 0.4  # levelset friction (not particle's)
    E           = 35e4  # Young's modulus
    nu          = .3    # Poisson ratio
    oneG        = 9.81
    micG        = .02*oneG

    #---------------------------------------------------------------------------
    mpm = tc.dynamics.MPM(
    res                       = (r, r, r),
  # delta_x                   = .003,        # def: 1/res
    base_delta_t              = .0001,
    frame_dt                  = 1/frameRate, # inverse frame rate (sec or 1/Hz)
    num_frames                = finalTime*frameRate,
    gravity                   = (0, -oneG, 0),
    particle_gravity          = True,
    rigidBody_gravity         = False,
    rigid_body_collision      = False,
    rigid_body_levelset_collision = False, # apply impulse on rB if its particle's phi<0
    particle_collision        = True,  # update particle pos and vel if its phi<0
    pushing_force             = 0,     # G2P, add normal-to-boundary velocity to boundary particle, "transfer.cpp"
    # due to numerical advection error:
    penalty                   = 1e4,   # G2P, add normal-to-boundary velocity to boundary particle, typical: 1e3, "transfer.cpp"
  # rigid_penalty             = 0,     # def: 1e3, for rigidBody collision
    print_rigid_body_state    = False, # print pos, vel and angVel on terminal
    dirichlet_boundary_radius = 0,     # 1/0 means dirichletBoundary is ON/OFF (apply velocity!)
    clean_boundary            = False, # def: true, clear boundary particles
    warn_particle_deletion    = True,
    verbose_bgeo              = False,  # write other particles attributes, "visualize.cpp"
    Houdini                   = False,
    snapshots                 = False,
    )

    #---------------------------------------------------------------------------
    ## add level-set
    # problem region
    levelset = mpm.create_levelset() # creates levelset with the dim of [0 to res] (dummy levelset)
    # inner cylinder:
    levelset.add_cylinder(tc.Vector(.3, 0, .3), .0998*Scale, False) # 0 means it's not important
    # outer cylinder:
    levelset.add_cylinder(tc.Vector(.3, 0, .3), .2001*Scale, True)
    # floor:
    levelset.add_plane(tc.Vector(0, +1, 0), -(.1499-(.05*Scale)) )
    # top:
    levelset.add_plane(tc.Vector(0, -1, 0), +(.1501+(.05*Scale)) )
    # level-set friction
    levelset.set_friction(FrictionLS)
    # dynamic level-set: False
    mpm.set_levelset(levelset, False)

    #---------------------------------------------------------------------------
    ## add soil
    tex = tc.Texture(
        'mesh',
        resolution    = (2*r, 2*r, 2*r),
        mesh_accuracy = 3,
        translate     = (.3, .15, .3),
        scale         = (.2*Scale, .1*Scale, .2*Scale),
        filename      = 'projects/mpm/data/tcSoilFine.obj'
    ) * 8 # number of particles per cell (ppc, max sample density)

    mpm.add_particles(
        type             = 'sand',
        minDistAH        = -1,    # .325*dx, -1 means use default
        pd_periodic      = False, # periodic domain
        density_tex      = tex.id,
        initial_velocity = (0, 0, 0),
        lambda_0         = E*nu/((1+nu)*(1-2*nu)), # first Lame parameter
        mu_0             = E/(2*(1+nu)), # shear modulus (G)
        alpha            = 1,
        beta             = 1,
        cohesion         = 0,
        logJp            = 0,
        friction_angle   = 40,
        density          = 2550,
    )

    #---------------------------------------------------------------------------
    ## inner cylinder
    # def rotation_function(t):
    #   return tc.Vector(0, Omega*t, 0)
    inner = mpm.add_particles(
        type              = 'rigid',
        density           = 1e6,
        scale             = (.0999*2*Scale, .1*Scale, .0999*2*Scale), # ri=coef/2
      # scripted_rotation = tc.function13(rotation_function),
        scripted_position = tc.constant_function13(tc.Vector(.3, .15, .3)),
        rotation_axis     = (0, 1, 0),
        initial_angular_velocity = Omega,
        friction          = FrictionRB,
        codimensional     = False, # not thin shell
        mesh_fn           = 'projects/mpm/data/cylinderInner.obj'
    )
    innerAux = mpm.add_particles(
        type              = 'rigid',
        scale             = (.19*Scale, .1*Scale, .19*Scale),
        scripted_position = tc.constant_function13(tc.Vector(.3, .15, .3)),
        scripted_rotation = tc.constant_function13(tc.Vector(0, 0, 0)),
        rotation_axis     = (0, 1, 0),
        friction          = 0.0001,
        codimensional     = True,
        mesh_fn           = 'projects/mpm/data/cylinderInner.obj'
    )
    mpm.add_articulation(
        type             = 'stepper',
        obj0             = inner,
        obj1             = innerAux,
        angular_velocity = Omega,
        axis             = (0, 1, 0),
    )

    #---------------------------------------------------------------------------
    ## change in rotation direction if required
    def frame_update(t, frame_dt):
        if t<(finalTime/2)-dt or t >(finalTime/2)+dt:
            return
        mpm.add_articulation(
            type             = 'stepper',
            obj0             = inner,
            obj1             = innerAux,
            angular_velocity = -Omega,
            axis             = (0, 1, 0),
        )

    #---------------------------------------------------------------------------
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
    #---------------------------------------------------------------------------
    # ## floor
    # floor = mpm.add_particles(
    #         type              = 'rigid',
    #         density           = 100000,
    #         scripted_rotation = tc.constant_function13(tc.Vector(0, 0, 0)),
    #         scale             = (2.15*Scale, .1*Scale, 2.15*Scale),
    #         scripted_position = tc.constant_function13(tc.Vector(.3, (.15-(.05*Scale)), .3)),
    #         friction          = Friction,
    #         codimensional     = False,
    #         mesh_fn           = 'projects/mpm/data/cylinder_jet.obj'
    # )
    #---------------------------------------------------------------------------
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
    #---------------------------------------------------------------------------
    # ## confining plate
    # def position_function(t):
    #     if (t < 1):
    #         return tc.Vector(0.3, 0.2001-(t/100), 0.3)
    #     else:
    #         return tc.Vector(0.3, 0.2001-(1/100), 0.3)
    # confPlate = mpm.add_particles(
    #         type              = 'rigid',
    #         density           = 40,
    #         scale             = (0.201, 0.201, 0.201),  # BIGGER TO PREVENT FROM LEAKAGE
    #         scripted_position = tc.function13(position_function),
    #         scripted_rotation = tc.constant_function13(tc.Vector(90, 0, 0)),
    #         friction          = Friction,
    #         codimensional     = True,
    #         angular_damping   = 0,
    #         mesh_fn           = 'projects/mpm/data/confPlate.obj'
    # )

    #---------------------------------------------------------------------------
    mpm.simulate(
        clear_output_directory = True,
      # frame_update           = frame_update,
        print_profile_info     = True,
    )
