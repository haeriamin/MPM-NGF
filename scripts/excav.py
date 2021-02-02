## Excavation model

# MPM         EXP
# 1 u     ->  4 m
# 0.5000  <-  2.00  m (bin's length)  => (k = 4)
# 0.0375  <-  0.15  m (soil's height)
# 0.2450  <-  1.00  m (bin's width)
# 0.0763  <-  0.305 m (plate's length)
# 0.1142  <-  0.457 m (plate's width)
# 0.0015  <-  0.006 m (plate's thickness)
# 0.2000  <-  0.80  m (plate's forward dist)
# 0.0125  <-  0.05  m (plate's vertical dist)
# 0.0100  <-  0.04  m/s (plate's forward velocity)
# 0.0033  <-  4*(grid size)

import taichi as tc


if __name__ == '__main__':
    r = 301 #251 #301
    dx = 1/r

    dim_scale = 0.25
    time_scale = dim_scale*2
    velocity_scale = dim_scale*2

    verticalTime = 2*time_scale  # 2s for 5cm depth, 0.8s for 2cm depth, with 2.5cm/s mean vel in exp
    forwardTime = 21*time_scale  #19.38
    stopTime = 0*time_scale
    finalTime = forwardTime + 2*verticalTime + 2*stopTime
    frameRate = 60

    FrictionRB = 0.3 #0.3 or -1
    FrictionLS = -1

    G = 9.81

    offset = 0.05/2
    x_bin = 4  #3 #3.5 #4 #5
    y_bin = 0.375  #0.225 #0.25 #0.375 #0.414
    z_bin = 1.75  #1.70 #1.75 #2.0 #2.5
    bin_size = tc.Vector(x_bin, y_bin, z_bin)

    # --------------------------------------------------------------------------
    mpm = tc.dynamics.MPM(
        res=(r, r, r),
        base_delta_t=4e-5,  # can't be larger than 3.5e-4
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
        clean_boundary=False,  # def: true, clear boundary particles
        warn_particle_deletion=False,
        cdf_3d_modified=True,
        compute_particle_impulses=True,
        visualize_particle_impulses=False,
        affect_particle_impulses=False,
        particle_bc_at_levelset=False, 

        snapshots=False,
        verbose_bgeo=False,
        write_particle=False,
        write_rigid_body=False,
        write_partio=True,
        write_dataset=False,
    )

    # level-set ----------------------------------------------------------------
    levelset = mpm.create_levelset()
    levelset.add_plane(tc.Vector(1, 0, 0), -offset)
    levelset.add_plane(tc.Vector(0, 1, 0), -offset)
    levelset.add_plane(tc.Vector(0, 0, 1), -offset)
    levelset.add_plane(tc.Vector(-1, 0, 0), offset+x_bin/10)
    levelset.add_plane(tc.Vector(0, 0, -1), offset+z_bin/10)
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
    ) * 4 #1

    ## Check 1: dt < (1/res^2)*t0/2/A^2/d^2
    ## Check 2: dt < (1/res)*sqrt(rho/E)
    ## Check 3: mu_s > sqrt(3)*(1-(2*nu))/(1+nu)

    E_mod = 15e6*dim_scale
    nu = 0.3
    grain_density = 2583
    packing_fraction = 0.67
    critical_density = packing_fraction * grain_density
    mpm.add_particles(
        type='nonlocal',
        pd=True,
        density_tex=tex.id,
        density=grain_density,
        critical_density=critical_density,
        packing_fraction=packing_fraction,

        S_mod=E_mod/2/(1+nu),
        B_mod=E_mod/3/(1-2*nu),
        A_mat=0.48,
        dia=0.0003*dim_scale,  # verified

        # These 3 parameters are corrolated [local, mu_2-mu_s=0.2616] (but not from [PhD nonlocal])
        # mu_s should be larger than sqrt(3)*(1-(2*nu))/(1+nu)
        mu_s=.7000, # 35deg
        mu_2=.9616, # 44deg
        I_0=0.278,

        # Larger this larger oscillations & denser sand
        t_0=1e-4, #1e-4, #5e-5, 
    )

    # bucket -------------------------------------------------------------------
    def position_function(t):

        # Wrong run but correct results!!
        # if (t <= verticalTime):
        #     speed = 0.025 * velocity_scale  # m/sec
        #     return tc.Vector(offset+0.95*x_bin/10-speed*(t),
        #                      offset+(2.04*y_bin/10)-speed*(t),
        #                      offset+0.50*z_bin/10)
        # elif (t > forwardTime+verticalTime):
        #     speed = 0.025 * velocity_scale  # m/sec
        #     return tc.Vector(offset+0.95*x_bin/10-speed*(t),
        #                      offset+(2.04*y_bin/10)-speed*(verticalTime)+speed*(t-(forwardTime+verticalTime)),
        #                      offset+0.50*z_bin/10)
        # else:
        #     speed = 0.04 * velocity_scale  # m/sec
        #     return tc.Vector(offset+0.95*x_bin/10-speed*(t),
        #                      offset+(2.04*y_bin/10)-speed*(verticalTime),
        #                      offset+0.50*z_bin/10)

        offset_coeff = 0.9 #1.08 #1.07 #1.06 #0.9
        depth = 0#.0125+1*dx #4.7625*dx #=.0125/dx+1 #0
        depth_coeff = 2.04 #2.04 for 3.8 deg #1.925 for being in soil #was 2.5 or 2.7
        if (t <= verticalTime):
            speed = (0.0250*velocity_scale) # m/sec
            return tc.Vector(offset+offset_coeff*x_bin/10
                                -speed*(t),
                             offset+(depth_coeff*y_bin/10-depth)  
                                -speed*(t),
                             offset+0.50*z_bin/10)
        elif (t > verticalTime and t <= verticalTime+stopTime):
            return tc.Vector(offset+offset_coeff*x_bin/10
                                -((0.0250*velocity_scale))*(verticalTime),
                             offset+(depth_coeff*y_bin/10-depth)
                                -((0.0250*velocity_scale))*(verticalTime),
                             offset+0.50*z_bin/10)
        elif (t > verticalTime+stopTime and t <= verticalTime+stopTime+forwardTime):
            speed = 0.04*velocity_scale  # m/sec
            return tc.Vector(offset+offset_coeff*x_bin/10
                                -((0.0250*velocity_scale))*(verticalTime)
                                -speed*(t-(verticalTime+stopTime)),
                            offset+(depth_coeff*y_bin/10-depth)
                                -((0.0250*velocity_scale))*(verticalTime),
                            offset+0.50*z_bin/10)
        elif (t > verticalTime+stopTime+forwardTime and t <= verticalTime+stopTime+forwardTime+stopTime):
            return tc.Vector(offset+offset_coeff*x_bin/10
                                -((0.0250*velocity_scale))*(verticalTime)
                                -(0.04*velocity_scale)*(forwardTime),
                             offset+(depth_coeff*y_bin/10-depth)
                                -((0.0250*velocity_scale))*(verticalTime),
                             offset+0.50*z_bin/10)
        elif (t > verticalTime+stopTime+forwardTime+stopTime):
            speed = (0.0250*velocity_scale)  # m/sec
            return tc.Vector(offset+offset_coeff*x_bin/10
                                -((0.0250*velocity_scale))*(verticalTime)
                                -(0.04*velocity_scale)*(forwardTime)
                                -speed*(t-(verticalTime+stopTime+forwardTime+stopTime)),
                             offset+(depth_coeff*y_bin/10-depth)
                                -((0.0250*velocity_scale))*(verticalTime)
                                +speed*(t-(verticalTime+stopTime+forwardTime+stopTime)),
                             offset+0.50*z_bin/10)

    mpm.add_particles(
        type='rigid',
        density=1e5,
        packing_fraction=1,  
        rotation_axis=(0, 0, 1),
        friction=FrictionRB,
        scripted_position=tc.function13(position_function),
        scripted_rotation=tc.constant_function13(tc.Vector(0, 0, 90-3.8)), #3.8 #30.8
        scale=(0.0763, 0.0070, 0.1142),  # thickness should be at least equal to 2dx  // was 0.0050 and 1.5x
        # scale=(0.0763, dx, 0.1142),
        codimensional=False, #False
        mesh_fn='projects/mpm/data/plate_houdini.obj',
    )

    # --------------------------------------------------------------------------
    mpm.simulate(
        clear_output_directory=True,
        print_profile_info=True,
    )


    # Accelerated motion:
        # v = 0.04 * velocity_scale  # m/sec
        # a = 0.04  # m/sec2
        # v_if = 0.0283 * velocity_scale  # m/sec
        # a_if = 0.0283  # m/sec2
        # # Downward
        # if (t <= 1*time_scale):
        #     return tc.Vector(offset+0.95*x_bin/10
        #                         -(a_if*t**2)/2,
        #                      offset+(2.04*y_bin/10)
        #                         -(a_if*t**2)/2,
        #                      offset+0.50*z_bin/10)
        # elif (t > 1*time_scale and t <= 1.8*time_scale):
        #     return tc.Vector(offset+0.95*x_bin/10
        #                         -(a_if*(1*time_scale)**2)/2
        #                         -v_if*(t-1*time_scale),
        #                      offset+(2.04*y_bin/10)
        #                         -(a_if*(1*time_scale)**2)/2
        #                         -v_if*(t-1*time_scale),
        #                      offset+0.50*z_bin/10) 
        # elif (t > 1.8*time_scale and t <= verticalTime):
        #     return tc.Vector(offset+0.95*x_bin/10
        #                         -(a_if*(1*time_scale)**2)/2
        #                         -v_if*(.8*time_scale)
        #                         -(-a_if*(t**2-(1.8*time_scale)**2)/2 + a_if*(verticalTime)*(t-1.8*time_scale)),
        #                      offset+(2.04*y_bin/10)
        #                         -(a_if*(1*time_scale)**2)/2
        #                         -v_if*(.8*time_scale)
        #                         -(-a_if*(t**2-(1.8*time_scale)**2)/2 + a_if*(verticalTime)*(t-1.8*time_scale)),
        #                      offset+0.50*z_bin/10)
        # # Forward
        # elif (t > verticalTime and t <= verticalTime+1*time_scale):
        #     return tc.Vector(offset+0.95*x_bin/10
        #                         -(a_if*(1*time_scale)**2)/2
        #                         -v_if*(.8*time_scale)
        #                         -(-a_if*((verticalTime)**2-(1.8*time_scale)**2)/2 + a_if*(verticalTime)*((verticalTime)-1.8*time_scale))
        #                         -(a*(t-verticalTime)**2)/2,
        #                      offset+(2.04*y_bin/10)
        #                         -(a_if*(1*time_scale)**2)/2
        #                         -v_if*(.8*time_scale)
        #                         -(-a_if*((verticalTime)**2-(1.8*time_scale)**2)/2 + a_if*(verticalTime)*((verticalTime)-1.8*time_scale)),
        #                      offset+0.50*z_bin/10)
        # elif (t > verticalTime+1*time_scale and t <= verticalTime+forwardTime-1*time_scale):
        #     return tc.Vector(offset+0.95*x_bin/10
        #                         -(a_if*(1*time_scale)**2)/2
        #                         -v_if*(.8*time_scale)
        #                         -(-a_if*((verticalTime)**2-(1.8*time_scale)**2)/2 + a_if*(verticalTime)*((verticalTime)-1.8*time_scale))
        #                         -(a*(1*time_scale)**2)/2
        #                         -v*(t-(verticalTime+1*time_scale)),
        #                      offset+(2.04*y_bin/10)
        #                         -(a_if*(1*time_scale)**2)/2
        #                         -v_if*(.8*time_scale)
        #                         -(-a_if*((verticalTime)**2-(1.8*time_scale)**2)/2 + a_if*(verticalTime)*((verticalTime)-1.8*time_scale)),
        #                      offset+0.50*z_bin/10)
        # elif (t > verticalTime+forwardTime-1*time_scale and t <= verticalTime+forwardTime):
        #     return tc.Vector(offset+0.95*x_bin/10
        #                         -(a_if*(1*time_scale)**2)/2
        #                         -v_if*(.8*time_scale)
        #                         -(-a_if*((verticalTime)**2-(1.8*time_scale)**2)/2 + a_if*(verticalTime)*((verticalTime)-1.8*time_scale))
        #                         -(a*(1*time_scale)**2)/2
        #                         -v*(forwardTime-2*time_scale)
        #                         -(-a*(t**2-(verticalTime+forwardTime-1*time_scale)**2)/2 + a*(verticalTime+forwardTime)*(t-(verticalTime+forwardTime-1*time_scale))),
        #                      offset+(2.04*y_bin/10)
        #                         -(a_if*(1*time_scale)**2)/2
        #                         -v_if*(.8*time_scale)
        #                         -(-a_if*((verticalTime)**2-(1.8*time_scale)**2)/2 + a_if*(verticalTime)*((verticalTime)-1.8*time_scale)),
        #                      offset+0.50*z_bin/10)
        # # Upward
        # elif (t > verticalTime+forwardTime and t <= verticalTime+forwardTime+1*time_scale):
        #     return tc.Vector(offset+0.95*x_bin/10
        #                         -(a_if*(1*time_scale)**2)/2
        #                         -v_if*(.8*time_scale)
        #                         -(-a_if*((verticalTime)**2-(1.8*time_scale)**2)/2 + a_if*(verticalTime)*((verticalTime)-1.8*time_scale))
        #                         -(a*(1*time_scale)**2)/2
        #                         -v*(forwardTime-2*time_scale)
        #                         -(-a*((verticalTime+forwardTime)**2-(verticalTime+forwardTime-1*time_scale)**2)/2 + a*(verticalTime+forwardTime)*((verticalTime+forwardTime)-(verticalTime+forwardTime-1*time_scale)))
        #                         -(a_if*(t-(verticalTime+forwardTime))**2)/2,
        #                      offset+(2.04*y_bin/10)
        #                         -(a_if*(1*time_scale)**2)/2
        #                         -v_if*(.8*time_scale)
        #                         -(-a_if*((verticalTime)**2-(1.8*time_scale)**2)/2 + a_if*(verticalTime)*((verticalTime)-1.8*time_scale))
        #                         +(a_if*(t-(verticalTime+forwardTime))**2)/2,
        #                      offset+0.50*z_bin/10)
        # elif (t > verticalTime+forwardTime+1*time_scale and t <= verticalTime+forwardTime+1.8*time_scale):
        #     return tc.Vector(offset+0.95*x_bin/10
        #                         -(a_if*(1*time_scale)**2)/2
        #                         -v_if*(.8*time_scale)
        #                         -(-a_if*((verticalTime)**2-(1.8*time_scale)**2)/2 + a_if*(verticalTime)*((verticalTime)-1.8*time_scale))
        #                         -(a*(1*time_scale)**2)/2
        #                         -v*(forwardTime-2*time_scale)
        #                         -(-a*((verticalTime+forwardTime)**2-(verticalTime+forwardTime-1*time_scale)**2)/2 + a*(verticalTime+forwardTime)*((verticalTime+forwardTime)-(verticalTime+forwardTime-1*time_scale)))
        #                         -(a_if*(1*time_scale)**2)/2
        #                         -v_if*(t-(verticalTime+forwardTime+1*time_scale)),
        #                      offset+(2.04*y_bin/10)
        #                         -(a_if*(1*time_scale)**2)/2
        #                         -v_if*(.8*time_scale)
        #                         -(-a_if*((verticalTime)**2-(1.8*time_scale)**2)/2 + a_if*(verticalTime)*((verticalTime)-1.8*time_scale))
        #                         +(a_if*(1*time_scale)**2)/2
        #                         +v_if*(t-(verticalTime+forwardTime+1*time_scale)),
        #                      offset+0.50*z_bin/10)
        # elif (t > verticalTime+forwardTime+1.8*time_scale):
        #     return tc.Vector(offset+0.95*x_bin/10
        #                         -(a_if*(1*time_scale)**2)/2
        #                         -v_if*(.8*time_scale)
        #                         -(-a_if*((verticalTime)**2-(1.8*time_scale)**2)/2 + a_if*(verticalTime)*((verticalTime)-1.8*time_scale))
        #                         -(a*(1*time_scale)**2)/2
        #                         -v*(forwardTime-2*time_scale)
        #                         -(-a*((verticalTime+forwardTime)**2-(verticalTime+forwardTime-1*time_scale)**2)/2 + a*(verticalTime+forwardTime)*((verticalTime+forwardTime)-(verticalTime+forwardTime-1*time_scale)))
        #                         -(a_if*(1*time_scale)**2)/2
        #                         -v_if*(.8*time_scale)
        #                         -(-a_if*(t**2-(verticalTime+forwardTime+1.8*time_scale)**2)/2 + a_if*(verticalTime+forwardTime+verticalTime)*(t-(verticalTime+forwardTime+1.8*time_scale))),                             
        #                      offset+(2.04*y_bin/10)
        #                         -(a_if*(1*time_scale)**2)/2
        #                         -v_if*(.8*time_scale)
        #                         -(-a_if*((verticalTime)**2-(1.8*time_scale)**2)/2 + a_if*(verticalTime)*((verticalTime)-1.8*time_scale))
        #                         +(a_if*(1*time_scale)**2)/2
        #                         +v_if*(.8*time_scale)
        #                         +(-a_if*(t**2-(verticalTime+forwardTime+1.8*time_scale)**2)/2 + a_if*(verticalTime+forwardTime+verticalTime)*(t-(verticalTime+forwardTime+1.8*time_scale))),
        #                      offset+0.50*z_bin/10)