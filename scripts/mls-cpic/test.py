import taichi as tc
from math import *
if __name__ == '__main__':
    rx        = 301
    ry        = 301
    rz        = 301
    omega     = 0.025 * 180/pi # deg/sec
    finalTime = 22 # sec
    frameRate = 100 # Hz
    random    = 0

    mpm = tc.dynamics.MPM(
    res                    = (rx, ry, rz),
    # delta_x              = 0.003,
    base_delta_t           = 0.0001,
    frame_dt               = 1/frameRate, # inverse frame rate (sec or 1/Hz)
    num_frames             = finalTime*frameRate,
    gravity                = (0, -9.81, 0),
    rigid_body_collision   = False,
    print_rigid_body_state = False
    # verbose_bgeo         = True, # output particle attributes other than position
    ## SHOULD BE CHECKED:
    # penalty              = 1e4,
    # num_frames           = 200,
    # sand_crawler         = True,
    # pushing_force        = 0
    )

    ## add levelset
    levelset = mpm.create_levelset()
    #levelset.add_sphere(tc.Vector(0.3, 0.2, 0.3), 0.1, True)
    levelset.add_cylinder(tc.Vector(0.3, random, 0.3), 0.2, True)
    levelset.add_plane(tc.Vector(0, 1, 0), -0.1)
    levelset.set_friction(0.4)
    mpm.set_levelset(levelset, False)

    ## add soil
    tex = tc.Texture(
          'mesh',
          resolution = (rx, ry, rz),
          translate  = (0.3, 0.3, 0.3),
          scale      = (0.2, 0.1, 0.2),
          filename   = 'projects/mpm/data/tcSoilFine.obj') * 5
    print('Texture Density = {}'.format(tex.id))
    mpm.add_particles(
    type             = 'sand',
    # pd_periodic    = False,
    density_tex      = tex.id,
    initial_velocity = (0, 0, 0),
    density          = 2550,
    color            = (0.8, 0.7, 1.0),
    friction_angle   = 10
    )

    # ## inner cylinder
    # def rotation_function(t):
    #     return tc.Vector(0, omega*t, 0)
    # inner = mpm.add_particles(
    #         type              = 'rigid',
    #         density           = 40,
    #         scale             = (0.1999, 0.102, 0.1999),
    #         scripted_rotation = tc.function13(rotation_function),
    #         scripted_position = tc.constant_function13(tc.Vector(0.3, 0.15, 0.3)),
    #         friction          = 0.4,
    #         codimensional     = True,
    #         angular_damping   = 0,
    #         mesh_fn           = 'projects/mpm/data/cylinderInner.obj'
    #         )

    mpm.simulate(
    clear_output_directory = True,
    # frame_update         = frame_update,
    print_profile_info     = True
    )
