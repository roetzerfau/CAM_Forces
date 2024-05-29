from libcpp.vector cimport vector 
from libcpp.string cimport string 
from libcpp cimport bool

IncludeFiles


cdef extern from "<CAM/cam_interface.hxx>" : 
  cdef cppclass CythonClassName C++ClassName : 
    CythonClassName(const double jump_param_composite, const double subaggregate_threshold) except + 
    void print_array() 
    const vector[unsigned int] &fields() 
    void do_cam() 
    void place_single_cell_bu_randomly( double _jump_parameter, double _porosity, unsigned int random_seed)
    bool place_sphere(double _jump_parameter, double _radius, int _position, vector[double] _face_charge) 
    bool place_plane(double _jump_parameter, vector[unsigned int] _extent, int _position, vector[double] _face_charge) 
    void place_particles()
    double average_particle_size_d()
    vector[unsigned int] particle_size_distribution_d()
    vector[double] eval_measures()

    double average_particle_size(const vector[ unsigned int]&)
    vector[unsigned int] particle_size_distribution(const vector[ unsigned int]&)
    unsigned int bulk_distance (const vector[ unsigned int]&, const vector[ unsigned int ]&)
    
