
#pragma once
#include <CAM/building_units.hxx>
#include <algorithm>
#include <cmath>
#include <iterator>
#include <vector>
namespace CAM
{
static double jump_parameter_composite;
template <auto nx>
static constexpr double get_jump_range_composite(const unsigned int _comp_size)
{
  return jump_parameter_composite /
         std::pow(
           _comp_size,
           1.0 / (double)nx.size());  //(jump_parameter_composite / std::sqrt(_comp_size) + 0.5);
}
/*!*********************************************************************************************
 * \brief Stores information about composites
 * \param building_units what bus is the composite made of
 * \param jump_parameter How far composite is allowed to jump.
 * \param field_indices what field indices is the composite made of (information also stored in bus)
 * \tparam nx
 **********************************************************************************************/
template <auto nx>
struct Composite
{
  std::vector<CAM::BuildingUnit<nx>*> building_units;
  double jump_parameter;
  std::vector<unsigned int> field_indices;

  Composite operator+(const Composite& comp) const
  {
    Composite composite;
    composite.building_units.resize(this->building_units.size() + comp.building_units.size());
    std::set_union(this->building_units.begin(), this->building_units.end(),
                   comp.building_units.begin(), comp.building_units.end(),
                   composite.building_units.begin());
    composite.field_indices.resize(this->field_indices.size() + comp.field_indices.size());
    std::set_union(this->field_indices.begin(), this->field_indices.end(),
                   comp.field_indices.begin(), comp.field_indices.end(),
                   composite.field_indices.begin());
    composite.jump_parameter = CAM::get_jump_range_composite<nx>(composite.field_indices.size());
    return composite;
  }
};

}  // namespace CAM
