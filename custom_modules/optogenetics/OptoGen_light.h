#ifndef OPTOGENETICS_LIGHT
#define OPTOGENETICS_LIGHT


namespace Opto::Light {

struct LightSource {
    // Parameters when light would pass through air directly onto surface of aparatus
    double free_intensity = 0;
    double free_frequency = 0;

    // TODO implement units for comparison later
    // std::string intensity_units = "stuff/area";
    // std::string frequency_units = "c/nm";

    // Could be used in the future to account for effects where light in 3D
    // media is reduced/altered
    // TODO implement later
    /* virtual double calculate_spatially_reduced_intensity(const double& distance, const double& density) {return 0.0;};
    virtual double calculate_spatially_altered_dispersion(const double& distance, const double& density) {return 0.0;};*/
};

}

#endif