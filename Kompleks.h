//Antoni Antoszek
#ifndef KOMPLEKS_H
#define KOMPLEKS_H
#include <iostream>
#include "WolnyModul.h"
#include "Sympleks.h"

namespace algebra {
    template<class S, unsigned d, unsigned p>
    class Kompleks : public WolnyModul<Sympleks<S, d>, p> {
    private:
        using BaseType = WolnyModul<Sympleks<S, d>, p>;

    public:
        // Constructors
        Kompleks() {
            WolnyModul<Sympleks<S,d>, p>();}

        explicit Kompleks(const Sympleks<S, d>& simplex) : BaseType(simplex) {}

        Kompleks(const Sympleks<S,d>& x) {
            this->generators.push_back(x);
            this->coefficients.push_back(ZMod<p>(1));}

        // Getters
        unsigned getDimension() const { return d; }
        unsigned getCharacteristic() const { return p; }

        // Assignment operators
        Kompleks& operator=(const Kompleks& other) {
            if (this != &other) {
                BaseType::operator=(other);
            }
            return *this;
        }

        Kompleks& operator=(const Sympleks<S, d>& simplex) {
            BaseType::operator=(simplex);
            return *this;
        }

        // Boundary computation
        Kompleks<S, d-1, p> brzeg() const {
            return computeBoundary();
        }

        Kompleks<S, d-1, p> boundary() const {
            return computeBoundary();
        }

    private:
        Kompleks<S, d-1, p> computeBoundary() const {
            Kompleks<S, d-1, p> result;

            // Special case: boundary of 1-simplex is always empty
            if (d == 1) {
                return result;
            }

            const auto& generators = this->getGenerators();
            const auto& coefficients = this->getCoefficients();

            for (size_t i = 0; i < generators.size(); ++i) {
                const Sympleks<S, d>& simplex = generators[i];
                const ZMod<p>& simplex_coeff = coefficients[i];

                WolnyModul<Sympleks<S, d-1>, d> boundary_module = simplex.boundary();

                for (auto it = boundary_module.begin(); it != boundary_module.end(); ++it) {
                    const Sympleks<S, d-1>& boundary_simplex = it->generator();
                    int coeff_val = static_cast<int>(it->wspolczynnik());

                    // Adjust for negatives if d > 1
                    if (d > 1 && coeff_val > static_cast<int>(d)/2) {
                        coeff_val -= d;
                    }

                    result.addGenerator(boundary_simplex, ZMod<p>(coeff_val) * simplex_coeff);
                }
            }

            return result;
        }

    public:
        // Composition methods
        Kompleks& operator+=(const Kompleks& other) {
            BaseType::operator+=(other);
            return *this;
        }

        friend Kompleks operator+(const Kompleks& lhs, const Kompleks& rhs) {
            Kompleks result = lhs;
            result += rhs;
            return result;
        }

        friend Kompleks operator-(const Kompleks& operand) {
            Kompleks result;
            const auto& generators = operand.getGenerators();
            const auto& coefficients = operand.getCoefficients();

            std::vector<Sympleks<S, d>> neg_generators = generators;
            std::vector<ZMod<p>> neg_coefficients;
            neg_coefficients.reserve(coefficients.size());

            for (const auto& coeff : coefficients) {
                neg_coefficients.push_back(-coeff);
            }

            result = Kompleks(neg_generators, neg_coefficients);
            return result;
        }

        friend Kompleks operator*(int scalar, const Kompleks& operand) {
            Kompleks result;
            const auto& generators = operand.getGenerators();
            const auto& coefficients = operand.getCoefficients();

            std::vector<Sympleks<S, d>> scalar_generators = generators;
            std::vector<ZMod<p>> scalar_coefficients;
            scalar_coefficients.reserve(coefficients.size());

            for (const auto& coeff : coefficients) {
                scalar_coefficients.push_back(ZMod<p>(scalar) * coeff);
            }

            result = Kompleks(scalar_generators, scalar_coefficients);
            return result;
        }
    };

    // Specialization for 0-dimensional complexes
    template<class S, unsigned p>
    class Kompleks<S, 0, p> : public WolnyModul<Sympleks<S, 0>, p> {
    private:
        using BaseType = WolnyModul<Sympleks<S, 0>, p>;

    public:
        // Constructors
        Kompleks() : BaseType() {}

        explicit Kompleks(const Sympleks<S, 0>& simplex) : BaseType(simplex) {}

        Kompleks(const Kompleks& other) : BaseType(other) {}

        Kompleks(const std::vector<Sympleks<S, 0>>& generators,
                const std::vector<ZMod<p>>& coefficients)
            : BaseType(generators, coefficients) {}

        // Destructor
        ~Kompleks() override = default;

        // Getters
        unsigned getDimension() const { return 0; }
        unsigned getCharacteristic() const { return p; }

        // Assignment operators
        Kompleks& operator=(const Kompleks& other) {
            if (this != &other) {
                BaseType::operator=(other);
            }
            return *this;
        }

        Kompleks& operator=(const Sympleks<S, 0>& simplex) {
            BaseType::operator=(simplex);
            return *this;
        }

        // Boundary computation - 0-simplices have empty boundary
        Kompleks<S, 0, p> brzeg() const {
            return Kompleks<S, 0, p>();
        }

        Kompleks<S, 0, p> boundary() const {
            return Kompleks<S, 0, p>();
        }

        // Composition methods
        Kompleks& operator+=(const Kompleks& other) {
            BaseType::operator+=(other);
            return *this;
        }

        friend Kompleks operator+(const Kompleks& lhs, const Kompleks& rhs) {
            Kompleks result = lhs;
            result += rhs;
            return result;
        }

        friend Kompleks operator-(const Kompleks& operand) {
            Kompleks result;
            const auto& generators = operand.getGenerators();
            const auto& coefficients = operand.getCoefficients();

            std::vector<Sympleks<S, 0>> neg_generators = generators;
            std::vector<ZMod<p>> neg_coefficients;
            neg_coefficients.reserve(coefficients.size());

            for (const auto& coeff : coefficients) {
                neg_coefficients.push_back(-coeff);
            }

            result = Kompleks(neg_generators, neg_coefficients);
            return result;
        }

        friend Kompleks operator*(int scalar, const Kompleks& operand) {
            Kompleks result;
            const auto& generators = operand.getGenerators();
            const auto& coefficients = operand.getCoefficients();

            std::vector<Sympleks<S, 0>> scalar_generators = generators;
            std::vector<ZMod<p>> scalar_coefficients;
            scalar_coefficients.reserve(coefficients.size());

            for (const auto& coeff : coefficients) {
                scalar_coefficients.push_back(ZMod<p>(scalar) * coeff);
            }

            result = Kompleks(scalar_generators, scalar_coefficients);
            return result;
        }

    };
}

#endif //KOMPLEKS_H
