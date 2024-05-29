/*!*********************************************************************************************
 * \file domain.hxx
 *
 * \brief This class implements a domain and holds all its properties
 * TODO Mayber outsource all evaluation functions to own class
 *
 * \tparam  nx              The size of a row for each dimension of the matrix
 * \tparam  n_fields        Size of the domain
 *
 ************************************************************************************************/
#pragma once

#include <CAM/building_units.hxx>
#include <CAM/composite.hxx>
#include <CAM/utils.hxx>
#include <array>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <unordered_map>
#ifndef SUB_AGGREGATES
#define SUB_AGGREGATES false
#endif
namespace CAM
{
/*!***********************************************************************************************
 * \brief two subaggregates stick together above this value
 ************************************************************************************************/
static double subaggregate_threshold;
/*!***********************************************************************************************
 * \brief precomputed stencils
 ************************************************************************************************/
std::vector<std::vector<unsigned int>> stencils_precomputed;
/*!***********************************************************************************************
 * \brief   Describes one connected set of bulk cells/fields.
 * Information gathered only by analysing the domain independ of defintion of bu and composites
 ************************************************************************************************/
struct Particle
{
  Particle(const std::vector<unsigned int>& _field_indices,
           const std::vector<unsigned int>& _numbers,
           unsigned int _n_surfaces_solid_fluid,
           unsigned int _n_surfaces_solid_solid)
  {
    field_indices = _field_indices;
    numbers = _numbers;
    n_surfaces_solid_fluid = _n_surfaces_solid_fluid;
    n_surfaces_solid_solid = _n_surfaces_solid_solid;
  }
  /*!*********************************************************************************************
   * \brief   Indices of building units contained in a particle
   * min_size = 1, max_size = number of building units in domain
   **********************************************************************************************/
  std::vector<unsigned int> numbers;
  /*!*********************************************************************************************
   * \brief   Location of the particle.
   **********************************************************************************************/
  std::vector<unsigned int> field_indices;

  unsigned int n_surfaces_solid_fluid;
  unsigned int n_surfaces_solid_solid;
};

/*!*********************************************************************************************
 * \brief   Domain class
 * TODO: implement subaggregate_threshold as template parameter
 * clang: non-type template argument of type 'double' is not yet supported
 **********************************************************************************************/
// template <auto nx, double subaggregate_threshold, typename fields_array_t>
template <auto nx, typename fields_array_t>
class Domain
{
 public:
  /*!*********************************************************************************************
   * \brief Construct a new Domain object
   *
   * \param _jump_parameter_composites How far composites particles are allowed to jump.
   ************************************************************************************************/
  Domain(const double _jump_parameter_composites = -1)
  {
    if (_jump_parameter_composites != -1.)
      CAM::jump_parameter_composite = _jump_parameter_composites;
    else
      CAM::jump_parameter_composite = 5;  // double _jump_parameter_composites = 1.0
    CAM::stencils_precomputed = get_n_stencils<nx>(jump_parameter_composite);
    if constexpr (std::is_same<fields_array_t,
                               std::vector<typename fields_array_t::value_type>>::value)
      domain_fields.resize(n_fields<nx>(), 0);
    else
    {
      static_assert(
        std::is_same<fields_array_t,
                     std::array<typename fields_array_t::value_type, n_fields<nx>()>>::value,
        "The fields array has incorrect size");
      domain_fields.fill(0);
    }
    max_field_number = 0;
  }
  /*!*********************************************************************************************
   * \brief Place any building unit (bu) into the domain
   * Basis function for placement for all special BUs
   * \param _unit building unit
   * \return true: cells of bu are marked with individual index, bu is stored in bu vector
   * \return false: bu could not be placed. all its cells are removed
   ************************************************************************************************/
  bool constexpr place_bu(const CAM::BuildingUnit<nx>& _unit)
  {
    // check if the number is already in the domain
    // auto iter = std::find_if(building_units.begin(), building_units.end(),
    //                          [&](const CAM::BuildingUnit<nx>& bu)
    //                          { return bu.get_number() == _unit.get_number(); });
    // if (_unit.get_number() == 0 || iter != building_units.end())
    //   return false;
    if (_unit.get_number() > max_field_number)
      max_field_number = _unit.get_number();

    unsigned int field;
    for (unsigned int i = 0; i < _unit.get_shape().size(); i++)
    {
      field = CAM::aim<nx>(_unit.get_reference_field(), _unit.get_shape()[i]);
      if (domain_fields[field] != 0)
        return false;
    }

    for (unsigned int i = 0; i < _unit.get_shape().size(); i++)
    {
      field = CAM::aim<nx>(_unit.get_reference_field(), _unit.get_shape()[i]);
      domain_fields[field] = _unit.get_number();
    }
    building_units.push_back(_unit);
    return true;
  }
  /*!*********************************************************************************************
   * \brief  Finds composites (particles containing more then one bu) in domain and stores
   * Combining particles/subaggregates as long as the connections are above a minimum threshold
   *defining the strength of the connection \param threshold connections/nofCells
   ************************************************************************************************/
  void find_sub_composites(const CAM::Composite<nx>& _composite,
                           std::vector<CAM::Composite<nx>>& sub_composites,
                           std::vector<std::vector<unsigned int>> _connections,
                           double threshold)
  {
    for (CAM::BuildingUnit<nx>* bu : _composite.building_units)
    {
      CAM::Composite<nx> sub_comp;
      sub_comp.building_units.push_back(bu);
      for (unsigned int shape_field : bu->get_shape())
        sub_comp.field_indices.push_back(CAM::aim<nx>(bu->get_reference_field(), shape_field));
      sub_comp.jump_parameter = CAM::get_jump_range_composite<nx>(sub_comp.field_indices.size());
      sub_composites.push_back(sub_comp);
    }
    bool combining = true;
    while (combining)
    {
      combining = false;
      unsigned int max_connection_a = 0, max_connection_b = 0;
      double connection_value_a, connection_value_b, max_connection_value = 0;
      for (unsigned int a = 0; a < sub_composites.size(); a++)
      {
        for (unsigned int b = a; b < sub_composites.size(); b++)
        {
          if (_connections[a][b] == 0 || a == b)
            continue;
          // connection value, the bigger the more connections are between two subaggregates
          connection_value_a =
            ((double)_connections[a][b] / 2) / (double)sub_composites[a].field_indices.size();
          connection_value_b =
            ((double)_connections[a][b] / 2) / (double)sub_composites[b].field_indices.size();

          if ((threshold < connection_value_a) && (threshold < connection_value_b) &&
              (max_connection_value < (connection_value_b + connection_value_b)))
          {
            max_connection_b = b;
            max_connection_a = a;
            max_connection_value = connection_value_b + connection_value_b;
            combining = true;
          }
        }
      }
      if (combining)
      {
        CAM::Composite<nx> combinded_composite =
          sub_composites[max_connection_a] + sub_composites[max_connection_b];
        sub_composites[max_connection_a] = combinded_composite;
        sub_composites.erase(sub_composites.begin() + max_connection_b);
        for (unsigned int i = 0; i < _connections.size(); i++)
        {
          _connections[max_connection_a][i] += _connections[max_connection_b][i];
          _connections[i][max_connection_a] += _connections[i][max_connection_b];
        }
        _connections.erase(_connections.begin() + max_connection_b);
        for (std::vector<unsigned int>& arr : _connections)
          arr.erase(arr.begin() + +max_connection_b);
      }
    }
  }
  /*!*********************************************************************************************
   * \brief  Finds composites (particles containing more then one bu) in domain and stores
   * information in std::vector<CAM::Composite<nx>*> composites; std::vector<Particle> particles;
   * Faster version using the connected particles and their boundaries
   ************************************************************************************************/
  void find_composites_via_bu_boundary()
  {
    constexpr unsigned int dim = nx.size();
    unsigned int neigh_field, boundaries_size;
    unsigned int field_number, field_number_neigh;
    unsigned int n_surfaces_solid_fluid, n_surfaces_solid_solid;
    unsigned int index;
    std::vector<unsigned int> boundaries, found_solids, helper;
    particles.clear();
    composites.clear();
    std::vector<bool> is_bu_visited(building_units.size(), false);

    for (unsigned int i = 0; i < building_units.size(); i++)
    {
      if (is_bu_visited[building_units[i].get_number()] == true)
        continue;

      CAM::Composite<nx> new_composite;
      new_composite.building_units.push_back(&building_units[i]);
#if SUB_COMPOSITES
      unsigned int local_a, local_b;
      std::vector<std::vector<unsigned int>> connections;
      connections.resize(new_composite.building_units.size());
      for (std::vector<unsigned int> arr : connections)
        arr.resize(new_composite.building_units.size(), 1);

      std::unordered_map<unsigned int, unsigned int> global_field_number_2_local_index;
      std::pair<unsigned int, unsigned int> pair(building_units[i].get_number(),
                                                 new_composite.building_units.size() - 1);
      global_field_number_2_local_index.insert(pair);
#endif
      n_surfaces_solid_fluid = 0;
      n_surfaces_solid_solid = 0;
      boundaries.clear();
      for (unsigned int boundary_field : building_units[i].get_boundary())
        boundaries.push_back(CAM::aim<nx>(building_units[i].get_reference_field(), boundary_field));

      found_solids.clear();
      for (unsigned int shape_field : building_units[i].get_shape())
        found_solids.push_back(CAM::aim<nx>(building_units[i].get_reference_field(), shape_field));

      field_number = building_units[i].get_number();
      is_bu_visited[field_number] = true;
      boundaries_size = boundaries.size();

      std::vector<unsigned int> composite_components;
      composite_components.push_back(field_number);

      for (unsigned int j = 0; j < boundaries_size; j++, boundaries_size = boundaries.size())
      {
        field_number = domain_fields[boundaries[j]];
        for (unsigned int k = 0; k < 2 * dim; ++k)
        {
          neigh_field = aim<nx>(boundaries[j], direct_neigh<nx>(k));
          field_number_neigh = domain_fields[neigh_field];

          if (field_number_neigh != 0 && is_bu_visited[field_number_neigh] != true)
          {
            is_bu_visited[field_number_neigh] = true;
            index = field_number_2_index[field_number_neigh];

            new_composite.building_units.push_back(&(building_units[index]));

            composite_components.push_back(field_number_neigh);

            for (unsigned int boundary_field : (building_units[index]).get_boundary())
              boundaries.push_back(
                CAM::aim<nx>((building_units[index]).get_reference_field(), boundary_field));

            for (unsigned int shape_field : (building_units[index]).get_shape())
              found_solids.push_back(
                CAM::aim<nx>((building_units[index]).get_reference_field(), shape_field));
#if SUB_COMPOSITES
            // second index
            std::pair<unsigned int, unsigned int> pair(field_number_neigh,
                                                       new_composite.building_units.size() - 1);
            global_field_number_2_local_index.insert(pair);

            connections.resize(new_composite.building_units.size());
            for (std::vector<unsigned int>& arr : connections)
              arr.resize(new_composite.building_units.size(), 0);
#endif
          }

          if (field_number_neigh == 0)
            n_surfaces_solid_fluid += 1;
          else if (field_number_neigh != field_number)
          {
            // this values have to be divided by two
            n_surfaces_solid_solid += 1;
#if SUB_COMPOSITES
            local_a = (global_field_number_2_local_index.find(field_number))->second;
            local_b = (global_field_number_2_local_index.find(field_number_neigh))->second;
            connections[local_a][local_b]++;
            connections[local_b][local_a]++;
#endif
          }
        }
      }
      // complete particle/composite is found
      if (composite_components.size() > 1)
      {
        new_composite.field_indices = found_solids;
        new_composite.jump_parameter =
          CAM::get_jump_range_composite<nx>(new_composite.field_indices.size());

#if SUB_COMPOSITES
        std::vector<CAM::Composite<nx>> sub_composites;
        find_sub_composites(new_composite, sub_composites, connections,
                            CAM::subaggregate_threshold);
        for (CAM::Composite<nx> comp : sub_composites)
          composites.push_back(comp);

#else
        composites.push_back(new_composite);
#endif
      }
      particles.push_back(Particle(found_solids, composite_components, n_surfaces_solid_fluid,
                                   n_surfaces_solid_solid));
    }
  }
  /*!***********************************************************************************************
   * \brief   Returns Domain
   *
   * \retval  domain_fields      Information about  status of each cell in domain
   ************************************************************************************************/
  const fields_array_t& fields() const
  {
    return domain_fields;
  }
  /*!***********************************************************************************************
   * \brief   Array of particle locations.
   ************************************************************************************************/
  fields_array_t domain_fields;
  /*!***********************************************************************************************
   * \brief   Vector of particles.
   ************************************************************************************************/
  std::vector<CAM::BuildingUnit<nx>> building_units;
  std::vector<CAM::Composite<nx>> composites;

  unsigned int max_field_number;
  std::vector<unsigned int> field_number_2_index;

  /*!*********************************************************************************************
   * \brief vector of particles (connected components)
   * contains
   ************************************************************************************************/
  std::vector<Particle> particles;

  // -------------------------------------------------------------------------------------------------
  // PRINTING SECTION STARTS HERE
  // -------------------------------------------------------------------------------------------------

  /*!*************************************************************************************************
   * \brief   Prints array of domain. Only used in print_n_dim().
   *
   * \param   fields      Particle index and size
   * \param   nx             Nx dimension and size
   * \param   init_index     Temporary index in current slice
   **************************************************************************************************/
  void print_1d(const fields_array_t& _fields, unsigned int init_index = 0)
  {
    const unsigned int len_numbers = std::log10(_fields.size() - 1) + 1;
    unsigned int index;

    for (unsigned int x = 0; x < nx[0]; ++x)
    {
      index = init_index + x;
      if (_fields[index] == 0)
        std::cout << std::setw(len_numbers) << std::setfill('0') << _fields[index] << "  ";
      else
        std::cout << "\033[0;31m" << std::setw(len_numbers) << std::setfill('0') << _fields[index]
                  << "\033[0m  ";
    }
    std::cout << std::endl;
  }
  /*!*************************************************************************************************
   * \brief   Prints array of fields in n dimensions. Only used in print_array().
   *
   * \tparam  n_dim          Temporary dimension of current slice
   * \param   fields      Particle index and size
   * \param   nx             Nx dimension and size
   * \param   init_index     Temporary index in current slice
   **************************************************************************************************/
  template <unsigned int n_dim>
  void print_n_dim(const fields_array_t& _fields, unsigned int init_index = 0)
  {
    if constexpr (n_dim == 1)
      print_1d(_fields, init_index);
    else
    {
      if (n_dim == 2 && nx.size() > 2)
      {
        unsigned int coord_i_helper = init_index;
        std::cout << std::endl;
        for (unsigned int i = 3; i < nx.size() + 1; ++i)
        {
          coord_i_helper = coord_i_helper / nx[i - 2];
          std::cout << i;
          if (i % 10 == 1 && i != 11)
            std::cout << "st";
          else if (i % 10 == 2 && i != 12)
            std::cout << "nd";
          else if (i % 10 == 3 && i != 13)
            std::cout << "rd";
          else
            std::cout << "th";
          std::cout << " coord: " << coord_i_helper % nx[i - 1] << "  ";
        }
        std::cout << std::endl;
      }
      for (unsigned int i = 0; i < nx[n_dim - 1]; ++i)
        print_n_dim<n_dim - 1>(_fields, (init_index + i) * nx[n_dim - 2]);
    }
  }
  /*!*************************************************************************************************
   * \brief   Runs print_n_dim() function
   *
   * \param   fields      Particle index and size
   * \param   nx             Nx dimension and size
   **************************************************************************************************/

  void print_array()
  {
    unsigned int sum = 0;
    std::for_each(domain_fields.begin(), domain_fields.end(), [&](unsigned int t) { sum += t; });
    print_n_dim<nx.size()>(domain_fields);
  }
};
}  // namespace CAM
