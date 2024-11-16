#.PHONY: LL3_SC LL3_BCC LL3_FCC

#all: LL3_SC LL3_BCC LL3_FCC;

#LL3_SC:; $(MAKE) name=LL3_SC CXXOPT="$(CXXOPT) -DSC"
#LL3_BCC:; $(MAKE) name=LL3_BCC CXXOPT="$(CXXOPT) -DBCC"
#LL3_FCC:; $(MAKE) name=LL3_FCC CXXOPT="$(CXXOPT) -DVCC"
# ввести понятие префикса при сборке для выполнения ODR?


name=LL3
headers=LL3.hpp
modules=LL3.cpp magnons.cpp

#iheader='%include "std_vector.i"'
#ifinish='%template(std_vector_float) std::vector<float>;  %template(std_vector_double) std::vector<double>; %pythoncode %{ std_vector_double._racs_pull_lock = std_vector_float._racs_pull_lock = True %}'

include /home/aiv/aiwlib/include/aiwlib/user.mk
