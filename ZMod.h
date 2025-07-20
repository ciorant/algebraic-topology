// Antoni Antoszek
#ifndef ZMOD_H
#define ZMOD_H
#include <iostream>

namespace algebra {
    template<unsigned p>
    class ZMod {
    private:
        unsigned value_;

        unsigned normalize(int x) const {
            int result = x % static_cast<int>(p);
            if (result < 0) {
                result += static_cast<int>(p);
            }
            return static_cast<unsigned>(result);
        }

    public:
        // Constructors
        ZMod() : value_(0) {}
        explicit ZMod(int x) : value_(normalize(x)) {}
        ZMod(const ZMod& other) : value_(other.value_) {}

        // Getters
        unsigned getValue() const { return value_; }

        // Setters
        void setValue(int x) { value_ = normalize(x); }

        // Type conversion
        explicit operator int() const {
            return static_cast<int>(value_);
        }

        // Assignment operators
        ZMod& operator=(const ZMod& other) {
            if (this != &other) {
                value_ = other.value_;
            }
            return *this;
        }

        ZMod& operator=(int x) {
            setValue(x);
            return *this;
        }

        // Comparison operators
        bool operator==(const ZMod& other) const {
            return value_ == other.value_;
        }

        bool operator==(int x) const {
            return static_cast<int>(value_) == x;
        }

        bool operator!=(const ZMod& other) const {
            return !(*this == other);
        }

        bool operator!=(int x) const {
            return !(*this == x);
        }

        // Arithmetic operators
        friend ZMod operator+(const ZMod& lhs, const ZMod& rhs) {
            return ZMod(static_cast<int>((lhs.value_ + rhs.value_) % p));
        }

        friend ZMod operator*(const ZMod& lhs, const ZMod& rhs) {
            return ZMod(static_cast<int>((lhs.value_ * rhs.value_) % p));
        }

        friend ZMod operator-(const ZMod& x) {
            return ZMod(static_cast<int>(p - x.value_));
        }

        // Stream operator
        friend std::ostream& operator<<(std::ostream& out, const ZMod& x) {
            out << static_cast<int>(x.value_);
            return out;
        }
    };

    // Specialization for p=0 (integers)
    template<>
    class ZMod<0> {
    private:
        int value_;

    public:
        // Constructors
        ZMod() : value_(0) {}
        explicit ZMod(int x) : value_(x) {}
        ZMod(const ZMod& other) : value_(other.value_) {}

        // Getters
        int getValue() const { return value_; }

        // Setters
        void setValue(int x) { value_ = x; }

        // Type conversion
        explicit operator int() const {
            return value_;
        }

        // Assignment operators
        ZMod& operator=(const ZMod& other) {
            if (this != &other) {
                value_ = other.value_;
            }
            return *this;
        }

        ZMod& operator=(int x) {
            setValue(x);
            return *this;
        }

        // Comparison operators
        bool operator==(const ZMod& other) const {
            return value_ == other.value_;
        }

        bool operator==(int x) const {
            return value_ == x;
        }

        bool operator!=(const ZMod& other) const {
            return !(*this == other);
        }

        bool operator!=(int x) const {
            return !(*this == x);
        }

        // Arithmetic operators
        friend ZMod operator+(const ZMod& lhs, const ZMod& rhs) {
            return ZMod(lhs.value_ + rhs.value_);
        }

        friend ZMod operator*(const ZMod& lhs, const ZMod& rhs) {
            return ZMod(lhs.value_ * rhs.value_);
        }

        friend ZMod operator-(const ZMod& x) {
            return ZMod(-x.value_);
        }

        // Stream operator
        friend std::ostream& operator<<(std::ostream& out, const ZMod& x) {
            out << x.value_;
            return out;
        }
    };
}

#endif //ZMOD_H
