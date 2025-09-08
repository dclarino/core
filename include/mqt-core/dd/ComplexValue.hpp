/*
 * Copyright (c) 2023 - 2025 Chair for Design Automation, TUM
 * Copyright (c) 2025 Munich Quantum Software Company GmbH
 * All rights reserved.
 *
 * SPDX-License-Identifier: MIT
 */

#pragma once

#include <array>
#include <cstdint>
#include <iostream>
#include <string>
#include <utility>
#include <cmath>
#include <algorithm>
#include <complex>

namespace dd {

constexpr int OmegaDegree = 7;

// Integer-coefficient algebraic omega complex value
struct ComplexValue {
    using Coeffs = std::array<int64_t, OmegaDegree + 1>;
    Coeffs coeffs{};

    double r = 0.0;
    double i = 0.0;

    ComplexValue() = default;

    ComplexValue(int64_t real) {
        coeffs.fill(0);
        coeffs[0] = real;
        updateComponents();
    }
    ComplexValue(double real) {
        coeffs.fill(0);
        coeffs[0] = static_cast<int64_t>(std::round(real));
        updateComponents();
    }
    ComplexValue(int64_t real, int64_t imag) {
        coeffs.fill(0);
        coeffs[0] = real;
        coeffs[2] = imag;
        updateComponents();
    }
    ComplexValue(const Coeffs& cfs) : coeffs(cfs) {
        updateComponents();
    }
    ComplexValue(const std::complex<double>& c) {
        coeffs.fill(0);
        coeffs[0] = static_cast<int64_t>(std::round(c.real()));
        coeffs[2] = static_cast<int64_t>(std::round(c.imag()));
        updateComponents();
    }
    ComplexValue(double real, double imag) {
      coeffs.fill(0);
      coeffs[0] = static_cast<int64_t>(std::round(real));
      coeffs[2] = static_cast<int64_t>(std::round(imag));
      r = real;
      i = imag;
      // updateComponents(); // if needed for compatibility
    }
    // Factory for a complex value with floating-point components
    static ComplexValue fromRealImag(double real, double imag) {
        ComplexValue val;
        val.r = real;
        val.i = imag;
        // Optionally: fill coeffs for compatibility if needed
        return val;
    }

    static ComplexValue omega(int k);

    [[nodiscard]] bool operator==(const ComplexValue& other) const noexcept;
    [[nodiscard]] bool operator!=(const ComplexValue& other) const noexcept;
    [[nodiscard]] bool exactlyZero() const noexcept;
    [[nodiscard]] bool exactlyOne() const noexcept;
    [[nodiscard]] bool approximatelyEquals(const ComplexValue& other, int64_t tol = 0) const noexcept;
    [[nodiscard]] bool approximatelyZero() const noexcept;

    void writeBinary(std::ostream& os) const;
    void readBinary(std::istream& is);

    void fromString(const std::string& realStr, std::string imagStr);

    static std::string toString(const double& real, const double& imag, bool formatted = true, int precision = -1);
    std::string toString(bool formatted = true, int precision = -1) const;

    ComplexValue& operator+=(const ComplexValue& rhs) noexcept;
    ComplexValue& operator*=(int64_t s) noexcept;
    [[nodiscard]] ComplexValue operator*(const ComplexValue& other) const noexcept;
    [[nodiscard]] ComplexValue operator*(int64_t s) const noexcept;
    [[nodiscard]] ComplexValue operator-() const noexcept;

    [[nodiscard]] double real() const noexcept { return r; }
    [[nodiscard]] double imag() const noexcept { return i; }
    [[nodiscard]] double mag2() const noexcept;
    [[nodiscard]] double mag() const noexcept;

    [[nodiscard]] std::complex<double> asFloat() const noexcept;

    static std::pair<std::uint64_t, std::uint64_t> getLowestFraction(double x, std::uint64_t maxDenominator = 1U << 10);
    static void printFormatted(std::ostream& os, double num, bool imaginary = false);

    void updateComponents();
};

ComplexValue operator+(const ComplexValue& c1, const ComplexValue& c2);
ComplexValue operator*(const ComplexValue& c1, int64_t s);
ComplexValue operator*(int64_t s, const ComplexValue& c1);
ComplexValue operator*(const ComplexValue& c1, double s);
ComplexValue operator/(const ComplexValue& c1, int64_t s);
ComplexValue operator/(const ComplexValue& c1, const ComplexValue& c2);

std::ostream& operator<<(std::ostream& os, const ComplexValue& c);

} // namespace dd

namespace std {
template <>
struct hash<dd::ComplexValue> {
    size_t operator()(const dd::ComplexValue& c) const noexcept;
};
}
