//Antoni Antoszek
#ifndef WOLNYMODUL_H
#define WOLNYMODUL_H
#include "ZMod.h"
#include <algorithm>
#include <vector>

namespace algebra {
    template<class S, unsigned p>
    class WolnyModul {
    private:
        mutable std::vector<S> generators_;
        mutable std::vector<ZMod<p>> coefficients_;
        mutable bool is_normalized_;

        void normalize() const {
            if (is_normalized_) return;
            is_normalized_ = true;

            std::vector<std::pair<S, ZMod<p>>> pairs;
            pairs.reserve(generators_.size());
            for (size_t i = 0; i < generators_.size(); ++i) {
                pairs.emplace_back(generators_[i], coefficients_[i]);
            }

            std::sort(pairs.begin(), pairs.end(),
                [](const std::pair<S, ZMod<p>>& a, const std::pair<S, ZMod<p>>& b) {
                    return a.first < b.first;
                });

            generators_.clear();
            coefficients_.clear();

            if (!pairs.empty()) {
                S current_gen = pairs[0].first;
                ZMod<p> current_coeff = pairs[0].second;
                for (size_t i = 1; i < pairs.size(); ++i) {
                    if (pairs[i].first == current_gen) {
                        current_coeff = current_coeff + pairs[i].second;
                    } else {
                        if (current_coeff != ZMod<p>(0)) {
                            generators_.push_back(current_gen);
                            coefficients_.push_back(current_coeff);
                        }
                        current_gen = pairs[i].first;
                        current_coeff = pairs[i].second;
                    }
                }
                if (current_coeff != ZMod<p>(0)) {
                    generators_.push_back(current_gen);
                    coefficients_.push_back(current_coeff);
                }
            }
        }

    public:
        // Constructors
        WolnyModul() : is_normalized_(true) {}

        explicit WolnyModul(const S& generator) : is_normalized_(false) {
            generators_.push_back(generator);
            coefficients_.push_back(ZMod<p>(1));
        }

        WolnyModul(const WolnyModul& other)
            : generators_(other.generators_), coefficients_(other.coefficients_),
              is_normalized_(other.is_normalized_) {}

        WolnyModul(const std::vector<S>& generators, const std::vector<ZMod<p>>& coefficients)
            : generators_(generators), coefficients_(coefficients), is_normalized_(false) {}

        // Getters
        const std::vector<S>& getGenerators() const {
            normalize();
            return generators_;
        }

        const std::vector<ZMod<p>>& getCoefficients() const {
            normalize();
            return coefficients_;
        }

        bool isNormalized() const { return is_normalized_; }

        unsigned getNonZeroCount() const {
            WolnyModul temp = *this;
            temp.normalize();
            return temp.coefficients_.size();
        }

        int getCoefficient(const S& generator) const {
            WolnyModul temp = *this;
            temp.normalize();
            auto it = std::find(temp.generators_.begin(), temp.generators_.end(), generator);
            if (it != temp.generators_.end()) {
                size_t index = std::distance(temp.generators_.begin(), it);
                return static_cast<int>(temp.coefficients_[index]);
            }
            return 0;
        }

        // Setters
        void setCoefficient(const S& generator, int coefficient) {
            auto it = std::find(generators_.begin(), generators_.end(), generator);
            if (it != generators_.end()) {
                size_t index = std::distance(generators_.begin(), it);
                coefficients_[index] = ZMod<p>(coefficient);
            } else {
                generators_.push_back(generator);
                coefficients_.push_back(ZMod<p>(coefficient));
            }
            is_normalized_ = false;
        }

        void addGenerator(const S& generator, const ZMod<p>& coefficient = ZMod<p>(1)) {
            generators_.push_back(generator);
            coefficients_.push_back(coefficient);
            is_normalized_ = false;
        }

        void clear() {
            generators_.clear();
            coefficients_.clear();
            is_normalized_ = true;
        }

        // Assignment operators
        WolnyModul& operator=(const WolnyModul& other) {
            if (this != &other) {
                generators_ = other.generators_;
                coefficients_ = other.coefficients_;
                is_normalized_ = other.is_normalized_;
            }
            return *this;
        }

        WolnyModul& operator=(const S& generator) {
            clear();
            generators_.push_back(generator);
            coefficients_.push_back(ZMod<p>(1));
            is_normalized_ = false;
            return *this;
        }

        // Compound assignment operators
        WolnyModul& operator+=(const WolnyModul& other) {
            generators_.insert(generators_.end(), other.generators_.begin(), other.generators_.end());
            coefficients_.insert(coefficients_.end(), other.coefficients_.begin(), other.coefficients_.end());
            is_normalized_ = false;
            return *this;
        }

        WolnyModul& operator,(const S& generator) {
            *this += WolnyModul(generator);
            return *this;
        }

        // Arithmetic operators
        friend WolnyModul operator+(const WolnyModul& lhs, const WolnyModul& rhs) {
            WolnyModul result = lhs;
            result += rhs;
            return result;
        }

        friend WolnyModul operator-(const WolnyModul& operand) {
            WolnyModul result;
            result.generators_ = operand.generators_;
            result.coefficients_.reserve(operand.coefficients_.size());
            for (const auto& coeff : operand.coefficients_) {
                result.coefficients_.push_back(-coeff);
            }
            result.is_normalized_ = operand.is_normalized_;
            return result;
        }

        friend WolnyModul operator*(int scalar, const WolnyModul& operand) {
            WolnyModul result;
            result.generators_ = operand.generators_;
            result.coefficients_.reserve(operand.coefficients_.size());
            for (const auto& coeff : operand.coefficients_) {
                result.coefficients_.push_back(ZMod<p>(scalar) * coeff);
            }
            result.is_normalized_ = operand.is_normalized_;
            return result;
        }

        // Stream operator
        friend std::ostream& operator<<(std::ostream& out, const WolnyModul<S, p>& module) {
            WolnyModul result = module;
            result.normalize();
            if (result.coefficients_.empty()) {
                out << "[]";
            } else {
                out << "[";
                for (size_t i = 0; i < result.generators_.size() - 1; ++i) {
                    out << '(' << result.coefficients_[i] << ',' << result.generators_[i] << "),";
                }
                out << '(' << result.coefficients_[result.generators_.size() - 1]
                    << ',' << result.generators_[result.generators_.size() - 1] << ')';
                out << "]";
            }
            return out;
        }

        // Iterator support
        class IteratorHelper {
        private:
            const S* generator_ptr_;
            const ZMod<p>* coefficient_ptr_;

        public:
            IteratorHelper(const S* gen_ptr, const ZMod<p>* coeff_ptr)
                : generator_ptr_(gen_ptr), coefficient_ptr_(coeff_ptr) {}

            IteratorHelper() : generator_ptr_(nullptr), coefficient_ptr_(nullptr) {}

            const S& generator() const { return *generator_ptr_; }
            const S& operator*() const { return *generator_ptr_; }
            const ZMod<p>& wspolczynnik() const { return *coefficient_ptr_; }
        };

        class const_iterator {
        private:
            friend class WolnyModul;
            const std::vector<S>* generators_;
            const std::vector<ZMod<p>>* coefficients_;
            size_t index_;
            mutable IteratorHelper helper_;

        public:
            const_iterator(const std::vector<S>* generators, const std::vector<ZMod<p>>* coefficients, size_t index)
                : generators_(generators), coefficients_(coefficients), index_(index) {}

            IteratorHelper operator*() const {
                return IteratorHelper(&(*generators_)[index_], &(*coefficients_)[index_]);
            }

            IteratorHelper* operator->() const {
                helper_ = IteratorHelper(&(*generators_)[index_], &(*coefficients_)[index_]);
                return &helper_;
            }

            const_iterator& operator++() {
                ++index_;
                return *this;
            }

            const_iterator operator++(int) {
                const_iterator temp = *this;
                ++index_;
                return temp;
            }

            bool operator==(const const_iterator& other) const {
                return index_ == other.index_;
            }

            bool operator!=(const const_iterator& other) const {
                return index_ != other.index_;
            }
        };

        const_iterator begin() const {
            normalize();
            return const_iterator(&generators_, &coefficients_, 0);
        }

        const_iterator end() const {
            return const_iterator(&generators_, &coefficients_, generators_.size());
        }
    };
}

#endif //WOLNYMODUL_H
