#ifndef OPTOGENETICS_LIGHT
#define OPTOGENETICS_LIGHT


namespace Opto::Light {

class LightSource {
    private:

    public:
        // Parameters when light would pass through air directly onto surface of aparatus
        double free_intensity{};
        double free_frequency{};
        double free_dispersion{};

        // Could be used in the future to account for effects where light in 3D
        // media is reduced/altered
        // TODO implement later
        double calculate_spatially_reduced_intensity(const double& distance, const double& density) {return 0.0;};
        double calculate_spatially_altered_dispersion(const double& distance, const double& density) {return 0.0;};
};

}

#endif