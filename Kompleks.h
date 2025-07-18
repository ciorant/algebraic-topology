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
        Kompleks() : BaseType() {}

        explicit Kompleks(const Sympleks<S, d>& simplex) : BaseType(simplex) {}

        Kompleks(const Kompleks& other) : BaseType(other) {}

        Kompleks(const std::vector<Sympleks<S, d>>& generators,
                const std::vector<ZMod<p>>& coefficients)
            : BaseType(generators, coefficients) {}

        // Destructor
        ~Kompleks() override = default;

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
        // Utility methods for chain complex operations
        bool isZeroDimensional() const { return d == 0; }
        bool isOneDimensional() const { return d == 1; }

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

        // Chain complex specific operations
        bool isCycle() const {
            if (d == 0) return true;
            return boundary().getNonZeroCount() == 0;
        }

        bool isBoundary() const {
            // This would require computing d+1 dimensional chains
            // For now, just return false as it's a complex computation
            return false;
        }

        // Access methods that preserve encapsulation
        size_t getNumberOfSimplices() const {
            return this->getGenerators().size();
        }

        bool hasSimplices() const {
            return !this->getGenerators().empty();
        }

        void addSimplex(const Sympleks<S, d>& simplex, int coefficient = 1) {
            this->addGenerator(simplex, ZMod<p>(coefficient));
        }

        void removeSimplex(const Sympleks<S, d>& simplex) {
            this->setCoefficient(simplex, 0);
        }

        int getSimplexCoefficient(const Sympleks<S, d>& simplex) const {
            return this->getCoefficient(simplex);
        }

        void setSimplexCoefficient(const Sympleks<S, d>& simplex, int coefficient) {
            this->setCoefficient(simplex, coefficient);
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

        // Utility methods
        bool isZeroDimensional() const { return true; }
        bool isOneDimensional() const { return false; }

        // 0-dimensional complexes are always cycles
        bool isCycle() const { return true; }

        // 0-dimensional complexes are boundaries only if they're empty
        bool isBoundary() const {
            return this->getNonZeroCount() == 0;
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

        // Access methods that preserve encapsulation
        size_t getNumberOfPoints() const {
            return this->getGenerators().size();
        }

        bool hasPoints() const {
            return !this->getGenerators().empty();
        }

        void addPoint(const Sympleks<S, 0>& point, int coefficient = 1) {
            this->addGenerator(point, ZMod<p>(coefficient));
        }

        void removePoint(const Sympleks<S, 0>& point) {
            this->setCoefficient(point, 0);
        }

        int getPointCoefficient(const Sympleks<S, 0>& point) const {
            return this->getCoefficient(point);
        }

        void setPointCoefficient(const Sympleks<S, 0>& point, int coefficient) {
            this->setCoefficient(point, coefficient);
        }
    };

    // Type aliases for common usage patterns
    template<class S, unsigned d>
    using KompleksZ = Kompleks<S, d, 0>;  // Complex over integers

    template<class S, unsigned d>
    using KompleksZ2 = Kompleks<S, d, 2>;  // Complex over Z/2Z

    template<class S, unsigned d>
    using KompleksZ3 = Kompleks<S, d, 3>;  // Complex over Z/3Z

    // Utility functions for complex operations
    template<class S, unsigned d, unsigned p>
    bool areHomologous(const Kompleks<S, d, p>& chain1, const Kompleks<S, d, p>& chain2) {
        Kompleks<S, d, p> difference = chain1 + (-chain2);
        return difference.isCycle();
    }

    template<class S, unsigned d, unsigned p>
    Kompleks<S, d, p> createSimplexChain(const std::vector<Sympleks<S, d>>& simplices) {
        Kompleks<S, d, p> result;
        for (const auto& simplex : simplices) {
            result.addSimplex(simplex);
        }
        return result;
    }

    template<class S, unsigned d, unsigned p>
    std::vector<Kompleks<S, d-1, p>> getAllBoundaryComponents(const Kompleks<S, d, p>& complex) {
        std::vector<Kompleks<S, d-1, p>> boundaries;
        const auto& generators = complex.getGenerators();

        for (const auto& simplex : generators) {
            Kompleks<S, d, p> single_simplex(simplex);
            boundaries.push_back(single_simplex.boundary());
        }

        return boundaries;
    }
}

#endif //KOMPLEKS_H