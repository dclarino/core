#include "dd/ComplexValue.hpp"

#include <cmath>
#include <cstring>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <cassert>

namespace dd {

void ComplexValue::updateComponents() {
    r = 0.0;
    i = 0.0;
    for (int k = 0; k <= OmegaDegree; ++k) {
        r += coeffs[k] * std::cos(k * M_PI / 4.0);
        i += coeffs[k] * std::sin(k * M_PI / 4.0);
    }
}

ComplexValue ComplexValue::omega(int k) {
    Coeffs c{};
    c[(k % 8 + 8) % 8] = 1;
    ComplexValue val(c);
    val.updateComponents();
    return val;
}

bool ComplexValue::operator==(const ComplexValue& other) const noexcept {
    return coeffs == other.coeffs;
}

bool ComplexValue::operator!=(const ComplexValue& other) const noexcept {
    return !(*this == other);
}

bool ComplexValue::exactlyZero() const noexcept {
    for (auto c : coeffs) if (c) return false;
    return true;
}

bool ComplexValue::exactlyOne() const noexcept {
    return coeffs[0] == 1 &&
        std::all_of(coeffs.begin() + 1, coeffs.end(), [](int64_t v){ return v == 0; });
}

bool ComplexValue::approximatelyEquals(const ComplexValue& other, int64_t tol) const noexcept {
    for (size_t i = 0; i < coeffs.size(); ++i)
        if (std::abs(coeffs[i] - other.coeffs[i]) > tol)
            return false;
    return true;
}

bool ComplexValue::approximatelyZero() const noexcept {
    for (auto c : coeffs) if (c != 0) return false;
    return true;
}

void ComplexValue::writeBinary(std::ostream& os) const {
    for (auto c : coeffs)
        os.write(reinterpret_cast<const char*>(&c), sizeof(int64_t));
}

void ComplexValue::readBinary(std::istream& is) {
    for (auto& c : coeffs)
        is.read(reinterpret_cast<char*>(&c), sizeof(int64_t));
    updateComponents();
}

void ComplexValue::fromString(const std::string& realStr, std::string imagStr) {
    coeffs.fill(0);
    coeffs[0] = realStr.empty() ? 0 : std::stoll(realStr);
    coeffs[2] = imagStr.empty() ? 0 : std::stoll(imagStr);
    updateComponents();
}

std::string ComplexValue::toString(const double& real, const double& imag, bool formatted, int precision) {
    std::ostringstream oss;
    if (precision >= 0)
        oss << std::fixed << std::setprecision(precision);
    if (formatted) {
        oss << real;
        if (imag >= 0) oss << "+";
        oss << imag << "i";
    } else {
        oss << real << "," << imag;
    }
    return oss.str();
}

std::string ComplexValue::toString(bool formatted, int precision) const {
    return toString(r, i, formatted, precision);
}

ComplexValue& ComplexValue::operator+=(const ComplexValue& rhs) noexcept {
    for (size_t i = 0; i < coeffs.size(); ++i)
        coeffs[i] += rhs.coeffs[i];
    updateComponents();
    return *this;
}

ComplexValue& ComplexValue::operator*=(int64_t s) noexcept {
    for (size_t i = 0; i < coeffs.size(); ++i)
        coeffs[i] *= s;
    updateComponents();
    return *this;
}

ComplexValue ComplexValue::operator*(const ComplexValue& other) const noexcept {
    ComplexValue res;
    for (size_t i = 0; i < coeffs.size(); ++i)
        for (size_t j = 0; j < coeffs.size(); ++j)
            res.coeffs[(i + j) % coeffs.size()] += coeffs[i] * other.coeffs[j];
    res.updateComponents();
    return res;
}

ComplexValue ComplexValue::operator*(int64_t s) const noexcept {
    ComplexValue res;
    for (size_t i = 0; i < coeffs.size(); ++i)
        res.coeffs[i] = coeffs[i] * s;
    res.updateComponents();
    return res;
}

ComplexValue ComplexValue::operator-() const noexcept {
    ComplexValue res;
    for (size_t i = 0; i < coeffs.size(); ++i)
        res.coeffs[i] = -coeffs[i];
    res.updateComponents();
    return res;
}

double ComplexValue::mag2() const noexcept {
    return r*r + i*i;
}
double ComplexValue::mag() const noexcept {
    return std::sqrt(mag2());
}

std::complex<double> ComplexValue::asFloat() const noexcept {
    return std::complex<double>(r, i);
}

// Non-member operator+ (only this now, fixes ambiguity!)
ComplexValue operator+(const ComplexValue& c1, const ComplexValue& c2) {
    ComplexValue res;
    for (size_t i = 0; i < c1.coeffs.size(); ++i)
        res.coeffs[i] = c1.coeffs[i] + c2.coeffs[i];
    res.updateComponents();
    return res;
}
ComplexValue operator*(const ComplexValue& c1, int64_t s) { return c1.operator*(s); }
ComplexValue operator*(int64_t s, const ComplexValue& c1) { return c1.operator*(s); }
ComplexValue operator*(const ComplexValue& c1, double s) { return c1.operator*(static_cast<int64_t>(s)); }
ComplexValue operator/(const ComplexValue& c1, int64_t s) {
    ComplexValue res;
    for (size_t i = 0; i < c1.coeffs.size(); ++i) {
        if (c1.coeffs[i] % s != 0)
            throw std::runtime_error("Omega division: coefficient not divisible by scalar");
        res.coeffs[i] = c1.coeffs[i] / s;
    }
    res.updateComponents();
    return res;
}
ComplexValue operator/(const ComplexValue& /*c1*/, const ComplexValue& /*c2*/) {
    throw std::logic_error("Division for omega complex values is nontrivial and not implemented.");
}

std::ostream& operator<<(std::ostream& os, const ComplexValue& c) {
    bool first = true;
    for (int i = 0; i <= OmegaDegree; ++i) {
        if (c.coeffs[i] != 0) {
            if (!first) os << " + ";
            os << c.coeffs[i];
            if (i > 0) os << "*omega^" << i;
            first = false;
        }
    }
    if (first) os << "0";
    return os;
}

std::pair<std::uint64_t, std::uint64_t> ComplexValue::getLowestFraction(double x, std::uint64_t maxDenominator) {
    assert(x >= 0.);

    std::pair<std::uint64_t, std::uint64_t> lowerBound{0U, 1U};
    std::pair<std::uint64_t, std::uint64_t> upperBound{1U, 0U};

    while ((lowerBound.second <= maxDenominator) &&
           (upperBound.second <= maxDenominator)) {
        auto num = lowerBound.first + upperBound.first;
        auto den = lowerBound.second + upperBound.second;
        auto median = static_cast<double>(num) / static_cast<double>(den);
        if (std::abs(x - median) <= 1e-12) {
            if (den <= maxDenominator) {
                return std::pair{num, den};
            }
            if (upperBound.second > lowerBound.second) {
                return upperBound;
            }
            return lowerBound;
        }
        if (x > median) {
            lowerBound = {num, den};
        } else {
            upperBound = {num, den};
        }
    }
    if (lowerBound.second > maxDenominator) {
        return upperBound;
    }
    return lowerBound;
}

void ComplexValue::printFormatted(std::ostream& os, double num, bool imaginary) {
    if (std::signbit(num)) {
        os << "-";
        num = -num;
    } else if (imaginary) {
        os << "+";
    }

    if (std::abs(num) < 1e-12) {
        os << "0" << (imaginary ? "i" : "");
        return;
    }

    const auto absnum = std::abs(num);
    auto fraction = getLowestFraction(absnum);
    auto approx = static_cast<double>(fraction.first) / static_cast<double>(fraction.second);

    if (std::abs(absnum - approx) < 1e-12) {
        if (fraction.first == 1U && fraction.second == 1U) {
            os << (imaginary ? "i" : "1");
        } else if (fraction.second == 1U) {
            os << fraction.first << (imaginary ? "i" : "");
        } else if (fraction.first == 1U) {
            os << (imaginary ? "i" : "1") << "/" << fraction.second;
        } else {
            os << fraction.first << (imaginary ? "i" : "") << "/" << fraction.second;
        }
        return;
    }

    const double SQRT2_2 = std::sqrt(2.0) / 2.0;
    const double PI = M_PI;
    const auto abssqrt = absnum / SQRT2_2;
    fraction = getLowestFraction(abssqrt);
    approx = static_cast<double>(fraction.first) / static_cast<double>(fraction.second);
    if (std::abs(abssqrt - approx) < 1e-12) {
        if (fraction.first == 1U && fraction.second == 1U) {
            os << (imaginary ? "i" : "1") << "/√2";
        } else if (fraction.second == 1U) {
            os << fraction.first << (imaginary ? "i" : "") << "/√2";
        } else if (fraction.first == 1U) {
            os << (imaginary ? "i" : "1") << "/(" << fraction.second << "√2)";
        } else {
            os << fraction.first << (imaginary ? "i" : "") << "/(" << fraction.second << "√2)";
        }
        return;
    }

    const auto abspi = absnum / PI;
    fraction = getLowestFraction(abspi);
    approx = static_cast<double>(fraction.first) / static_cast<double>(fraction.second);
    if (std::abs(abspi - approx) < 1e-12) {
        const std::string imagUnit = imaginary ? "i" : "";

        if (fraction.first == 1U && fraction.second == 1U) {
            os << "π" << imagUnit;
        } else if (fraction.second == 1U) {
            os << fraction.first << "π" << imagUnit;
        } else if (fraction.first == 1U) {
            os << "π" << imagUnit << "/" << fraction.second;
        } else {
            os << fraction.first << "π" << imagUnit << "/" << fraction.second;
        }
        return;
    }

    if (imaginary) {
        os << num << "i";
    } else {
        os << num;
    }
}

} // namespace dd

namespace std {
size_t hash<dd::ComplexValue>::operator()(const dd::ComplexValue& c) const noexcept {
    size_t h = 0;
    for (auto coeff : c.coeffs)
        h ^= std::hash<int64_t>{}(coeff) + 0x9e3779b9 + (h << 6) + (h >> 2);
    return h;
}
}
