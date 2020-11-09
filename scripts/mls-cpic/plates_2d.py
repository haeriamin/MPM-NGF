import taichi as tc
from math import *



if __name__ == '__main__':
  r = 301
  dx = 1/r
  finalTime = 10  # 10
  frameRate = 60 # 60

  FrictionRB = 0.3
  FrictionLS = 0.4

  oneG = 9.81

  mpm = tc.dynamics.MPM(
    res = (r, r),
    base_delta_t=1e-4,
    frame_dt=1/frameRate,
    num_frames=finalTime*frameRate,
    num_threads=-1,
    gravity=(0, -oneG,),
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
    # penalty=1e4,
    # sand_crawler=True,
    # rigid_penalty=0,  # def: 1e3, for rigidBody collision
    print_rigid_body_state=False,  # print pos, vel and angVel on terminal
    clean_boundary=True,  # def: true, clear boundary particles
    warn_particle_deletion=False,
    verbose_bgeo=False,  # write other particles attributes, "visualize.cpp"
    Houdini=True,
    snapshots=False,
    # visualize_cdf=True,
    # visualize_particle_cdf=True,
  )

  levelset = mpm.create_levelset()
  levelset.add_plane(tc.Vector(0, 1,), -0.1)
  levelset.set_friction(FrictionLS)
  mpm.set_levelset(levelset, False)

  def position_function(t):
    speed = 0.01  # m/sec
    return tc.Vector(0.9-speed*t, .15,)

  mesh = tc.SegmentMesh()
  segment = ((0.9, 0.1), (0.9+dx, 0.1), )
  mesh.add_segment(segment)
  segment = ((0.9+dx, 0.1),(0.9+dx, 0.2),)
  mesh.add_segment(segment)
  segment = ((0.9+dx, 0.2),(0.9, 0.2),)
  mesh.add_segment(segment)
  segment = ((0.9, 0.2),(0.9, 0.1),)
  mesh.add_segment(segment)
  mpm.add_particles(
    type='rigid',
    codimensional=True,
    density=400,
    rotation_axis=(0, 0, 1),
    scripted_position=tc.function12(position_function),
    scripted_rotation=tc.constant_function((0, 0)),
    scale=(1, 1),
    friction=FrictionRB,
    segment_mesh=mesh,
  )
  # mpm.add_articulation(
  #   type='rotation',
  #   obj0=rigids[0],
  #   obj1=rigids[1]
  # )

  tex = tc.Texture(
    'ring',
    inner=0.,
    outer=0.05,
    translate=(0.5, 0.15)
  ) * 4

  mpm.add_particles(
    type='nonlocal',
    pd=True,
    density_tex=tex.id,
    S_mod=1.3e5,
    B_mod=9.5e3,
    A_mat=0.48,
    dia=0.0003,
    rho_s=2583,
    mu_s=0.7, #0.3819,  # 20 deg
    mu_2=1.5, #0.6435,  # 32 deg, so the friction coefficient is varying
    I_0=0.278,
    t_0=1e-1,  # should be calibrated depending on soil properties, larger this looser sand
  )

  mpm.simulate(
    clear_output_directory=True,
    print_profile_info=True,
  )
