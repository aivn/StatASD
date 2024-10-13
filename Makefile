name=LL3
headers=LL3.hpp
modules=LL3.cpp magnons.cpp

iheader='%include "std_vector.i"'
ifinish='%template(std_vector_float) std::vector<float>;  %template(std_vector_double) std::vector<double>; %pythoncode %{ std_vector_double._racs_pull_lock = std_vector_float._racs_pull_lock = True %}'

include /home/aiv/aiwlib/include/aiwlib/user.mk
