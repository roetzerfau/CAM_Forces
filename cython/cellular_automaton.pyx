#distutils : language = c++

cdef class PythonClassName: 
  cdef CythonClassName* thisptr #hold a C++ instance which we're wrapping 
  def __cinit__(self, jump_param_composites = "default" , subaggregate_threshold = "default"):
    if jump_param_composites   == "default" : jump_param_composites = 5
    if subaggregate_threshold ==  "default": subaggregate_threshold= 0.01
    self.thisptr = new CythonClassName(jump_param_composites, subaggregate_threshold) 
  def __dealloc__(self): 
    del self.thisptr 
  def print_array(self): 
    self.thisptr.print_array() 
  def fields(self): 
    return self.thisptr.fields() 
  def place_single_cell_bu_randomly(self, _jump_parameter, _porosity, random_seed): 
    if _jump_parameter == "default": _jump_parameter = 1.0 
    if _porosity == "default" : _porosity = 0.5 
    if random_seed == "default" : random_seed = 0 
    self.thisptr.place_single_cell_bu_randomly(_jump_parameter,_porosity,  random_seed) 
  def place_sphere(self, _jump_parameter, _radius, _position, _face_charge): 
    if _jump_parameter == "default" : _jump_parameter = 1.0 
    if _radius == "default" : _radius = 1.0 
    if _position == "default" : _position = -1 
    if _face_charge == "default": _face_charge = [ 0 ]* 2
    return self.thisptr.place_sphere(_jump_parameter, _radius, _position, _face_charge) 
  def place_plane(self, _jump_parameter, _extent, _position, _face_charge): 
    if _jump_parameter == "default" : _jump_parameter = 1.0 
    if _extent == "default" : _extent = [ 0 ]* 3 
    if _position == "default": _position = -1 
    if _face_charge == "default": _face_charge = [ 0 ]* 2 
    return self.thisptr.place_plane(_jump_parameter, _extent, _position, _face_charge) 
  def place_particles(self): 
    self.thisptr.place_particles() 
  def do_cam(self) : 
    self.thisptr.do_cam()

  def eval_measures(self):
    return self.thisptr.eval_measures()


  def average_particle_size_d(self):
    return self.thisptr.average_particle_size_d()
  def particle_size_distribution_d(self) : 
    return self.thisptr.particle_size_distribution_d()

  def average_particle_size(self,vec):
    return self.thisptr.average_particle_size(vec)
  def particle_size_distribution(self, vec) : 
    return self.thisptr.particle_size_distribution(vec)

  def bulk_distance(self,vec_a, vec_b):
    return self.thisptr.bulk_distance (vec_a, vec_b)
  def __getstate__(self):
    return None
