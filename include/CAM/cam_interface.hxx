/*!*********************************************************************************************
 * \file cam_interface.hxx
 * \brief Interface for Python
 **********************************************************************************************/
#pragma once
#include <CAM/building_units.hxx>
#include <CAM/cellular_automaton.hxx>
#include <CAM/composite.hxx>
#include <CAM/domain.hxx>
#include <CAM/evaluation.hxx>
#include <array>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
namespace CAM
{
template <auto nx,
          unsigned int const_stencil_size,
          typename fields_array_t = std::array<unsigned int, CAM::n_fields<nx>()>>
class CAMInterface
{
 private:
  CAM::Domain<nx, fields_array_t> domain;
  unsigned int rand_seed;

 public:
  CAMInterface(const double _jump_parameter_composites = 5,
               const double _subaggregate_threshold = 0.01)
  {
    CAM::jump_parameter_composite = _jump_parameter_composites;
    CAM::subaggregate_threshold = _subaggregate_threshold;
    CAM::stencils_precomputed = get_n_stencils<nx>(jump_parameter_composite);
  }
  /*!*********************************************************************************************
   * \brief Creates domain with building units containing only one cell
   *
   * \param _porosity The percentage of void space, not occupied by solid/bu.
   * \param _jump_parameter How far individual particles are allowed to jump.
   * \param random_seed If given, sets random seed to given seed.
   ************************************************************************************************/
  void place_single_cell_bu_randomly(double _jump_parameter = 1,
                                     double _porosity = 0.5,
                                     unsigned int _random_seed = 0)
  {
    if (_random_seed == 0)
    {
      rand_seed = std::chrono::system_clock::now().time_since_epoch().count();
      std::srand(rand_seed);
    }

    else
      std::srand(_random_seed);

    unsigned int n_particles = (1. - _porosity) * n_fields<nx>();
    unsigned int position = std::rand() % (n_fields<nx>());
    for (unsigned int i = 0; i < n_particles; ++i)
    {
      while (domain.domain_fields[position] != 0)
        position = std::rand() % (n_fields<nx>());
      domain.place_bu(CAM::create_single_cell_bu<nx>(_jump_parameter, position,
                                                     domain.building_units.size() + 1));
    }
  }
  /*!*********************************************************************************************
   * \brief Places HyperSphere (2D: Circle, 3D: Sphere) into domain
   *
   * \param _position field index of center point; no position is given (-1) -> random position
   * \param _radius all cells are completely within the radius
   * \param _jump_parameter How far sphere is allowed to jump.
   * \return true: sphere is placed
   * \return false: sphere could not be placed. all its cells are removed
   ************************************************************************************************/
  bool place_sphere(double _jump_parameter = 1,
                    double _radius = 1,
                    int _position = -1,
                    std::vector<double> _face_charge_v = std::vector<double>(2 * nx.size(), 0))
  {
    std::array<double, 2 * nx.size()> face_charge_a;
    for (unsigned int i = 0; i < _face_charge_v.size(); i++)
      face_charge_a[i] = _face_charge_v[i];

    const CAM::BuildingUnit<nx> bu = CAM::create_hyper_sphere<nx>(
      _jump_parameter, _radius, _position, domain.building_units.size() + 1, face_charge_a);
    return domain.place_bu(bu);
  }
  /*!*********************************************************************************************
   * \brief Place a limited hyper plane (2D: Rectangle, 3D: cuboid) into domain
   *
   * \param _position field index of center point; no position is given (-1) -> random position
   * \param _extent size of plane in each dimension
   * \param _jump_parameter How far hyper plane is allowed to jump.
   * \return true: plane is placed
   * \return false: plane could not be placed. all its cells are removed
   ************************************************************************************************/
  bool place_plane(double _jump_parameter = 1,
                   std::vector<unsigned int> _extent_v = std::vector<unsigned int>(nx.size(), 0),
                   int _position = -1,
                   std::vector<double> _face_charge_v = std::vector<double>(2 * nx.size(), 0))
  {
    std::array<unsigned int, nx.size()> extent_a;
    for (unsigned int i = 0; i < _extent_v.size(); i++)
      extent_a[i] = _extent_v[i];

    std::array<double, 2 * nx.size()> face_charge_a;
    for (unsigned int i = 0; i < _face_charge_v.size(); i++)
      face_charge_a[i] = _face_charge_v[i];

    const BuildingUnit<nx> bu = CAM::create_hyper_plane<nx>(
      _jump_parameter, extent_a, _position, domain.building_units.size() + 1, face_charge_a);
    return domain.place_bu(bu);
  }
  // TODO parameterize this function
  void place_particles()
  {
    rand_seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::srand(rand_seed);
    unsigned int random_point;

    // std::vector<unsigned int> shape1 = CAMget_shape_by_feret<nx>(0, 30,
    // 20, 100); std::vector<unsigned int> shape2 =
    // CAM::get_shape_by_feret<nx>(0, 20, 10, 200);

    std::vector<unsigned int> shape1 = CAM::get_shape_by_feret<nx>(0, 10, 5, 10);
    std::vector<unsigned int> shape2 = CAM::get_shape_by_feret<nx>(0, 5, 1, 1);
    unsigned int index = 0;
    while (index < 10)
    {
      random_point = 0;
      for (unsigned int i = 0; i < nx.size(); ++i)
      {
        random_point += std::rand() % nx[i] * direct_neigh<nx>(2 * i + 1);
      }
      index = index + (unsigned int)domain.place_bu(CAM::create_custom_bu<nx>(
                        5, shape1, random_point, domain.building_units.size() + 1));
    }
    while (index < 20)
    {
      random_point = 0;
      for (unsigned int i = 0; i < nx.size(); ++i)
      {
        random_point += std::rand() % nx[i] * direct_neigh<nx>(2 * i + 1);
      }
      index = index + (unsigned int)domain.place_bu(CAM::create_custom_bu<nx>(
                        5, shape2, random_point, domain.building_units.size() + 1));
    }
  }

  void print_array() { domain.print_array(); }

  void do_cam() { CAM::CellularAutomaton<nx, const_stencil_size, fields_array_t>::apply(domain); }

  const fields_array_t& fields()
  {
    // show composites
    // domain.find_composites_via_bu_boundary();
    // std::shuffle(domain.composites.begin(), domain.composites.end(),
    //              std::default_random_engine(std::rand()));
    // for (unsigned int i = 0; i < domain.composites.size(); i++)
    // {
    //   std::cout << i << std::endl;
    //   for (unsigned int j : domain.composites[i].field_indices)
    //   {
    //     domain.domain_fields[j] = i * 10;
    //   }
    // }

    return domain.domain_fields;
  }

  std::vector<double> eval_measures()
  {
    const std::array<double, 12> meas_a =
      CAM::Evaluation<nx, fields_array_t>::eval_measures(domain);
    std::vector<double> meas_v;
    for (unsigned int k = 0; k < 12; ++k)
    {
      // std::cout << "Meas[" << k << "] = " << meas_a[k] << std::endl;
      meas_v.push_back(meas_a[k]);
    }

    return meas_v;
  }

  double average_particle_size_d()
  {
    return CAM::Evaluation<nx, fields_array_t>::average_particle_size(domain);
  }
  double average_particle_size(const fields_array_t& _domain)
  {
    return CAM::Evaluation<nx, fields_array_t>::average_particle_size(_domain);
  }
  std::vector<unsigned int> particle_size_distribution_d()
  {
    return CAM::Evaluation<nx, fields_array_t>::particle_size_distribution(domain);
  }
  std::vector<unsigned int> particle_size_distribution(const fields_array_t& _domain)
  {
    return CAM::Evaluation<nx, fields_array_t>::particle_size_distribution(_domain);
  }
  unsigned int bulk_distance(const fields_array_t& domain_a, const fields_array_t& domain_b)
  {
    return CAM::Evaluation<nx, fields_array_t>::bulk_distance(domain_a, domain_b);
  }
};
}  // namespace CAM
