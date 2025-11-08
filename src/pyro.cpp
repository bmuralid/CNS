#include "pyro.H"

namespace pyro {

template <typename _DataTypeT, typename _ContainerT>
std::string pyro<_DataTypeT, _ContainerT>::get_species_name(int index) {
    return species_names[index];
}

template <typename _DataTypeT, typename _ContainerT>
std::string pyro<_DataTypeT, _ContainerT>::get_element_name(int index) {
    return element_names[index];
}

template <typename _DataTypeT, typename _ContainerT>
int pyro<_DataTypeT, _ContainerT>::get_species_index(const std::string& name) {
    if (name == "O2") return 0;
    if (name == "N2") return 1;
    return -1;
}

template <typename _DataTypeT, typename _ContainerT>
int pyro<_DataTypeT, _ContainerT>::get_element_index(const std::string& name) {
    if (name == "O") return 0;
    if (name == "N") return 1;
    return -1;
}


// Static member definitions
template <typename _DataTypeT, typename _ContainerT>
const char* pyro<_DataTypeT, _ContainerT>::species_names[] = {
    "O2", "N2"
};

template <typename _DataTypeT, typename _ContainerT>
const char* pyro<_DataTypeT, _ContainerT>::element_names[] = {
    "O", "N"
};

// Explicit instantiations for common types
template struct pyro<double, double>;
template struct pyro<float, float>;

} // namespace pyro
