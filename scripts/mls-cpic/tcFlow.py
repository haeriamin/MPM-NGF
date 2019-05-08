# Taylor-Couette Flow, by Amin Haeri [email](ahaeri92@gmail.com)
# # from /ppo: models
# from ppo.mpm import *
# from ppo.trainer import Trainer # Trainer is a class in trainer.py located in /ppo folder
import taichi as tc
import math
if __name__ == '__main__':
    r            = 301
    dx           = 1/r
    dt           = .0001
    Omega        = 0.05     # rad/sec
    finalTime    = 22       # sec
    frameRate    = 50       # Hz
    Scale        = 1        # problem scale
    Friction     = -1       # rigidBody or levelset friction (not particle's), if not -1, causes segmentation fault:
    frictionAng  = 40;      # deg
    E            = 353700   # Young's modulus, def: 353700 Pa
    nu           = .3       # Poisson ratio, def: .3
    grainDensity = 2550     # kg/m^3
    oneG         = -9.81
    micG         = .02*oneG

    #---------------------------------------------------------------------------
    mpm = tc.dynamics.MPM(
    res                       = (r, r, r),
    delta_x                   = dx,
    base_delta_t              = dt,
    frame_dt                  = 1/frameRate, # inverse frame rate (sec or 1/Hz)
    num_frames                = finalTime*frameRate,
    gravity                   = (0, oneG, 0),
    particle_gravity          = True,
    rigidBody_gravity         = False,
    num_threads               = -1,
    rigid_body_collision      = False,
    rigid_body_levelset_collision = False, # apply impulse on rB if its particle's phi<0
    particle_collision        = True,  # update particle pos and vel if its phi<0
    pushing_force             = 0,     # G2P, add normal-to-boundary velocity to boundary particle, "transfer.cpp"
    penalty                   = 1e4,   # G2P, add normal-to-boundary velocity to boundary particle, typical: 1e3, "transfer.cpp", due to numerical advection error
  # rigid_penalty             = 0,     # def: 1e3, for rigidBody collision
    print_rigid_body_state    = False, # print pos, vel and angVel on terminal
    dirichlet_boundary_radius = 0,     # 1/0 means dirichletBoundary is ON/OFF (apply velocity!)
    clean_boundary            = True,  # def: true, clear boundary particles
    warn_particle_deletion    = True,
    verbose_bgeo              = False, # write other particles attributes, "visualize.cpp"
    Houdini                   = False,
    snapshots                 = False,
    Hardening                 = False,
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
    levelset.set_friction(Friction)
    # dynamic level-set: False
    mpm.set_levelset(levelset, False)

    #---------------------------------------------------------------------------
    ## add soil
    tex = tc.Texture(
        'mesh',
        resolution       = (2*r, 2*r, 2*r),
        mesh_accuracy    = 3,
        translate        = (.3, .15, .3),
        scale            = (.2*Scale, .1*Scale, .2*Scale),
        filename         = 'projects/mpm/data/tcSoilFine.obj'
    ) * 8 # number of particles per cell (ppc, max sample density)
    mpm.add_particles(
        type             = 'sand',
        minDistAH        = -1,    # .325*dx, -1 means use default
        pd_periodic      = True, # periodic domain
        density_tex      = tex.id,
        initial_velocity = (0, 0, 0),
        lambda_0         = E*nu/((1+nu)*(1-2*nu)), # first Lame parameter
        mu_0             = E/(2*(1+nu)), # shear modulus (G)
        friction_angle   = frictionAng,
        cohesion         = 0,
        density          = grainDensity,
        #color           = (.8, .7, 1),
    )
    #---------------------------------------------------------------------------
    ## inner cylinder
    inner = mpm.add_particles(
        type              = 'rigid',
        density           = 1e6,
        scale             = (.0999*2*Scale, .1*Scale, .0999*2*Scale), # ri=coef/2
        scripted_position = tc.constant_function13(tc.Vector(.3, .15, .3)),
        rotation_axis     = (0, 1, 0),
        initial_angular_velocity = Omega,
        friction          = Friction,
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
    mpm.simulate(
        clear_output_directory = True,
        print_profile_info     = True,
        frame_update           = frame_update,
    )

## Singularity 2.6.1:
# singularity pull shub://haeriamin/singularity-ubuntu:1604
# mv haeriamin-singularity-ubuntu-master-1604.simg ubuntu.simg
# singularity shell -H /private/k/kskoniec/a_hae ubuntu.simg
# singularity exec -H /private/k/kskoniec/a_hae ubuntu.simg /bin/bash
# singularity exec ubuntu.simg /bin/bash python install.py

## install taichi & mpm:
# python3 install.py # customized
# source ~/.bashrc
# modify main.py
# ti install mpm

## post processing:
# ffmpeg -framerate 50 -i untitled%d.jpg -c:v libx264 -profile:v high -crf 18 -pix_fmt yuv420p output.mp4

## ssh:
# rsync -v -e ssh ~/Desktop a_hae@speed-submit.encs.concordia.ca:/private/k/kskoniec/a_hae/install.py
# rsync -v -e ssh ~/Desktop/install.py a_hae@speed-submit.encs.concordia.ca:/private/k/kskoniec/a_hae
# rsync -avzhe ssh /Volumes/Z/taichi a_hae@speed-submit.encs.concordia.ca:/private/k/kskoniec/a_hae/taichi
# rsync -avzhe ssh a_hae@speed-submit.encs.concordia.ca:/private/k/kskoniec/a_hae/taichi/outputs/mpm/tcFlow /Volumes/X/
# to zip on mac: tar cvfz mpmCode.tar.gz mpmCode
# to unzip on server: tar -zxvf modules-4.2.3.tar.gz
# to delete on server: rm -f modules-4.2.3.tar.gz // -r for dirs
# to rename: mv pavel-demin-singularity-ubuntu-master-1804.simg ubuntu.simg
# /private/k/kskoniec/a_hae/opt/python-3.6.2/bin/python3 install.py
# vim ~/.bashrc
# see the used space: du -msh a_hae // du -m a_hae
# module avail
# module load python/3.5.1/default
# module load gcc/5.1/default
# check os version: cat /etc/*release
# copy: cp /home/n/nul-uge/tcsh.sh /private/k/kskoniec/a_hae

## vim
# delete all: dG
# paste: :set paste OR :set nopaste
# save and exit: :wq
# exit: :q

## run:
# ssh a_hae@speed-submit.encs.concordia.ca
# cd /private/k/kskoniec/a_hae
# singularity shell -H /private/k/kskoniec/a_hae ubuntu.simg
# source ~/.bashrc
# cd /private/k/kskoniec/a_hae/taichi/projects/mpm/scripts/mls-cpic
# python3 tcFlow.py

## modify:
# rsync -avzhe ssh /Volumes/Z/taichi/projects/mpm/data/ a_hae@speed-submit.encs.concordia.ca:/private/k/kskoniec/a_hae/taichi/projects/mpm/data
# rsync -avzhe ssh /Volumes/Z/taichi/projects/mpm/scripts/mls-cpic/ a_hae@speed-submit.encs.concordia.ca:/private/k/kskoniec/a_hae/taichi/projects/mpm/scripts/mls-cpic
# rsync -avzhe ssh /Volumes/Z/taichi/projects/mpm/src/ a_hae@speed-submit.encs.concordia.ca:/private/k/kskoniec/a_hae/taichi/projects/mpm/src
# rsync -avzhe ssh /Volumes/Z/taichi/projects/mpm/scripts/mls-cpic/tcFlow.py a_hae@speed-submit.encs.concordia.ca:/private/k/kskoniec/a_hae/taichi/projects/mpm/scripts/mls-cpic
# cd /private/k/kskoniec/a_hae/taichi/projects/mpm/src/ ## a file in this folder should be changed using vim not rsync
# cd /private/k/kskoniec/a_hae/taichi/python/taichi
# ti build

## results:
# cd /private/k/kskoniec/a_hae/taichi/outputs/mpm/
# rsync -avzhe ssh a_hae@speed-submit.encs.concordia.ca:/private/k/kskoniec/a_hae/taichi/outputs/mpm/tcFlow /Volumes/X
