import configparser
from .cpp import find_definition, extract_classname
from .names import cython_from_cpp
from .paths import main_dir

## \brief   Object that comprises all information for HyperHDG to create a problem.
class config:
  nx                  = []
  ca_settings        = []
  const_stencil_size = 0
  cython_replacements = []
  cpp_code            = ""
  include_files       = []
  debug_mode          = False
  print_message       = False

## \brief   Check that config is consistent.
def consistent(conf):
  assert isinstance(conf, config)
  if not (isinstance(conf.nx, list) and conf.nx != ""):
     return False
  for entry in conf.nx:
    if not isinstance(entry, int) or entry < 0:
      return False
  for entry in conf.ca_settings:
    if not isinstance(entry, bool):
      return False
  if not isinstance(conf.const_stencil_size, int):
    return False
  if not isinstance(conf.cpp_code, str):
    return False
  if not isinstance(conf.debug_mode, bool):
    return False
  if not isinstance(conf.print_message, bool):
    return False
  if not isinstance(conf.cython_replacements, list):
    return False
  for entry in conf.cython_replacements:
    if not isinstance(entry, str):
      return False
  if not isinstance(conf.include_files, list):
    return False
  for entry in conf.include_files:
    if not isinstance(entry, str):
      return False
  return True

## \brief   Add the files that need to be included for defining the problem to include list.
def extract_includes(conf, cpp_inc=""):
  assert isinstance(conf, config) and consistent(conf)
  helper = find_definition(".", "CAMInterface")
  if helper not in conf.include_files:
    conf.include_files.append(helper)
  include_string = ""
  for x in conf.include_files:
    include_string += "cdef extern from \"<" + x + ">\": pass\n"
  if cpp_inc != "":
    include_string += "cdef extern from \"<" + cpp_inc + ">\": pass\n"
  return include_string

## \brief   Evaluate config file that needs to be present for all .pyx/.pxd files.
def generate_cy_replace(conf):
  assert isinstance(conf, config)# and consistent(conf)
  config_file = configparser.ConfigParser()
  config_file.read(main_dir() + "/cython/" + "cellular_automaton" + ".cfg")
  n_replacements = int(config_file['default']['n_replacements'])
  cy_replace = []
  for i in range(n_replacements):
    if i < len(conf.cython_replacements) and conf.cython_replacements[i] != "":
      cy_replace.append(conf.cython_replacements[i])
    else:
      assert config_file.has_option('default','replacement' + '%02d' % (i+1))
      cy_replace.append(config_file['default']['replacement' + '%02d' % (i+1)])
  assert n_replacements == len(cy_replace)
  return cy_replace
