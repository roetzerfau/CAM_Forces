import hashlib, re

## \brief   Transform C++ class name to cython file name.
def cython_from_cpp(name):
  helper = re.sub('(.)([A-Z][a-z]+)', r'\1_\2', name)
  return re.sub('([a-z0-9])([A-Z])', r'\1_\2', helper).lower()

## \brief   Generate names of auxiliary classes that will be constructed (internal use only).
def files(conf):
  nx_string    = "std::array<unsigned int, " + str(len(conf.nx)) + ">({"
  for nx_dim in conf.nx:
    nx_string = nx_string + str(nx_dim) + ","
  nx_string    = nx_string[:-1] + "})"
  cpp_class    = "CAM::CAMInterface< " + nx_string + "," + str(conf.const_stencil_size) + ", std::vector<unsigned int> >"
  cython_file  = cpp_class + "_deb" if conf.debug_mode else cpp_class + "_rel"
  cython_file  = re.sub(' ', '', cython_file)
  cython_file  = "mod" + str(hashlib.sha256(cython_file.encode('utf-8')).hexdigest())
  cython_file  = re.sub('\+|\-|\*|\/|<|>|\,|\:', '_', cython_file)
  cython_class = cython_file + "CP"
  return nx_string, cpp_class, cython_file, cython_class
