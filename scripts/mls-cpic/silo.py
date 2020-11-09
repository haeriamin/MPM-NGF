import taichi as tc
from math import *
from random import *


if __name__ == '__main__':
    r = 251 #301 #251
    dx = 1/r
    finalTime = 15  # sec
    frameRate = 120  # Hz
    Scale = 1  # problem scale
    oneG = 9.81
    micG = .02*oneG

    FrictionLS = 2

    offset = 0.2
    length_bin = 0.25 * Scale
    height_bin = 0.25 * Scale
    width_bin = length_bin/4
    orifice = 0.025 * min(1,Scale*2)  #0.01
    length_soil = length_bin*10
    height_soil = 0.125*10 * Scale
    width_soil = width_bin*10

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
    levelset.add_plane(tc.Vector(0, +1, 0), -offset)
    # # left wall:
    # levelset.add_plane(tc.Vector(+1, 0, 0), -offset)
    # level-set friction
    levelset.set_friction(FrictionLS)
    # dynamic level-set: False
    mpm.set_levelset(levelset, False)

    # --------------------------------------------------------------------------
    walls_pos = [
        # top floor
        tc.Vector(
            offset + length_bin/2 + orifice,
            offset + height_bin,
            offset + width_bin/2),
        # left wall
        tc.Vector(
            offset,
            offset + (height_bin+(height_soil/10))/2,
            offset + width_bin/2),
        # right wall
        tc.Vector(
            offset + length_bin,
            offset + height_bin+(height_soil/10)/2,
            offset + width_bin/2),
        # rear wall
        tc.Vector(
            offset + length_bin/2,
            offset + height_bin+(height_soil/10)/2,
            offset),
        # front wall
        tc.Vector(
            offset + length_bin/2,
            offset + height_bin+(height_soil/10)/2,
            offset + width_bin)
    ]

    walls_rot = [
        tc.Vector(0, 0, 0),   # top floor
        tc.Vector(0, 0, 90),  # left wall  
        tc.Vector(0, 0, 90),  # right wall
        tc.Vector(90, 0, 0),  # rear wall
        tc.Vector(90, 0, 0)   # front wall
    ]

    extra = 0.005
    walls_scale = [
        (length_bin, 1, (.125+extra)*Scale), 
        (.45*Scale, 1, .25*Scale),
        (.15*Scale, 1, .15*Scale),
        (length_bin+extra, 1, .15*Scale),
        (length_bin+extra, 1, .15*Scale),
    ]

    walls_friction = [
        -1,
        -2,  # left wall
        2,
        2,
        2,
    ]

    for i in range(5):
        pos = walls_pos[i]
        rot = walls_rot[i]
        scale = walls_scale[i]
        fric = walls_friction[i]
        mpm.add_particles(
            type='rigid',
            density=1e5,
            friction=fric,
            scale=scale,
            scripted_position=tc.constant_function13(pos),
            scripted_rotation=tc.constant_function13(rot),
            codimensional=True,
            mesh_fn='projects/mpm/data/flat_cutter_low_res.obj'
        )

    # --------------------------------------------------------------------------
    # add soil
    tex = tc.Texture(
        'mesh',
        scale=tc.Vector(length_soil, height_soil, width_soil),
        translate=(
            offset + (length_soil/10)/2,
            offset + height_bin + (height_soil/10)/2, #+dx
            offset + (width_soil/10)/2),
        resolution=(r, r, r),
        mesh_accuracy=3,
        filename='projects/mpm/data/cube_smooth.obj',
    ) * 8 #4 # more than 1 makes cell-crossing issue

    # Lower this looser sand
    # Up to some value (1e5Pa) due to large Me_tr!
    Y_mod = 15e4 #15e4 #1e6
    nu = 0.3
    grain_density = 2583 #2450
    packing_fraction = 0.67
    critical_density = packing_fraction * grain_density
    mpm.add_particles(
        type='nonlocal',
        pd=True,
        density_tex=tex.id,
        density=grain_density, 
        critical_density=critical_density,
        packing_fraction=packing_fraction,
        initial_velocity=(0, 0, 0),

        S_mod=Y_mod/2/(1+nu),
        B_mod=Y_mod/3/(1-2*nu),
        # S_mod=5.0e4,#5.0e4,  # in FEM was 5.0e6
        # B_mod=2.0e4,#2.0e4,  # in FEM was 2.0e7

        A_mat=.48,#0.48,#0.58,#0.51
        dia=0.0003*Scale,

        # These 3 parameters are corrolated [local, mu_2-mu_s=0.2616] (but not from [PhD nonlocal])
        # mu_s should be larger than sqrt(3)*(1-(2*nu))/(1+nu)
        mu_s=.3819,  # .7000 = 35deg
        mu_2=1.5, #.6435,  # .9616 = 44deg
        I_0=0.278,

        # Larger this larger oscillations & denser sand
        t_0=1e-4, #1e-1 #1e-3 #1e-4
    )

    # --------------------------------------------------------------------------
    mpm.simulate(
        clear_output_directory=True,
        print_profile_info=True)
