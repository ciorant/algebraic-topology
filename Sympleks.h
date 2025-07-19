// Antoni Antoszek
#ifndef SYMPLEKS_H
#define SYMPLEKS_H
#include <iostream>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include "WolnyModul.h"

namespace algebra {
    template<class S, unsigned d>
    class Sympleks {
    private:
        std::vector<S> sequence_;

        void validateIndex(unsigned index) const {
            if (index < 1 || index > sequence_.size()) {
                throw std::out_of_range("Index out of bounds");
            }
        }

    public:
        // Constructors
        Sympleks() {
            sequence_.reserve(d);
            for (unsigned i = 0; i < d; ++i) {
                sequence_.push_back(S());
            }
        }

        Sympleks(const Sympleks& other) : sequence_(other.sequence_) {}

        explicit Sympleks(const S sequence[]) {
            sequence_.reserve(d);
            for (unsigned i = 0; i < d; ++i) {
                sequence_.push_back(sequence[i]);
            }
        }

        // Getters
        const std::vector<S>& getSequence() const { return sequence_; }

        const S& at(unsigned index) const {
            validateIndex(index);
            return sequence_[index - 1];
        }

        // Setters
        void setElement(unsigned index, const S& value) {
            validateIndex(index);
            sequence_[index - 1] = value;
        }

        void setSequence(const std::vector<S>& new_sequence) {
            if (new_sequence.size() != d) {
                throw std::invalid_argument("Sequence size must match dimension");
            }
            sequence_ = new_sequence;
        }

        // Assignment operator
        Sympleks& operator=(const Sympleks& other) {
            if (this != &other) {
                sequence_ = other.sequence_;
            }
            return *this;
        }

        // Comparison operators
        bool operator<(const Sympleks& other) const {
            return std::lexicographical_compare(sequence_.begin(), sequence_.end(),
                other.sequence_.begin(), other.sequence_.end());
        }

        bool operator<=(const Sympleks& other) const {
            return *this < other || *this == other;
        }

        bool operator!=(const Sympleks& other) const {
            return sequence_ != other.sequence_;
        }

        bool operator==(const Sympleks& other) const {
            return sequence_ == other.sequence_;
        }


        // Subscript operators (1-indexed for mathematical consistency)
        S& operator[](unsigned index) {
            validateIndex(index);
            return sequence_[index - 1];
        }

        const S& operator[](unsigned index) const {
            validateIndex(index);
            return sequence_[index - 1];
        }

        // Boundary computation
        WolnyModul<Sympleks<S, d-1>, d> boundary() const {
            if (d == 0) {
                return WolnyModul<Sympleks<S, d-1>, d>();
            }
            if (d == 1) {
                return WolnyModul<Sympleks<S, d-1>, d>();
            }

            std::vector<Sympleks<S, d-1>> generators;
            std::vector<ZMod<d>> coefficients;

            for (unsigned i = 0; i < d; ++i) {
                Sympleks<S, d-1> boundary_simplex;
                unsigned target_index = 0;

                for (unsigned j = 0; j < d; ++j) {
                    if (j != i) {
                        boundary_simplex.sequence_[target_index] = this->sequence_[j];
                        target_index++;
                    }
                }

                generators.push_back(boundary_simplex);
                ZMod<d> sign = (i % 2 == 0) ? ZMod<d>(1) : -ZMod<d>(1);
                coefficients.push_back(sign);
            }

            return WolnyModul<Sympleks<S, d-1>, d>(generators, coefficients);
        }

        // Stream operator
        friend std::ostream& operator<<(std::ostream& out, const Sympleks& simplex) {
            if (simplex.sequence_.empty()) {
                out << "()";
            } else {
                out << '(';
                for (size_t i = 0; i < simplex.sequence_.size() - 1; ++i) {
                    out << simplex.sequence_[i] << ',';
                }
                out << simplex.sequence_[simplex.sequence_.size() - 1] << ')';
            }
            return out;
        }
    };
}

#endif //SYMPLEKS_H
