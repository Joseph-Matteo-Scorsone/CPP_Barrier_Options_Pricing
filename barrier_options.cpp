#include <cmath>
#include <algorithm>
#include <iostream>

// BarrierOption class to model and price barrier options
class BarrierOption {
public:
    // Enum for option type: Call or Put
    enum class OptionType { Call, Put };
    // Enum for barrier type
    enum class BarrierType { UpAndOut, DownAndOut, UpAndIn, DownAndIn };

    // Constructor with parameters for option specification
    BarrierOption(double S,        // Spot price
                  double K,        // Strike price
                  double H,        // Barrier level
                  double r,        // Risk-free interest rate
                  double q,        // Dividend yield
                  double sigma,    // Volatility
                  double T,        // Time to maturity
                  OptionType type,
                  BarrierType barrierType,
                  double rebate = 0.0) // Rebate in case of knock-out
        : S(S), K(K), H(H), r(r), q(q), sigma(sigma), T(T),
          type(type), barrierType(barrierType), rebate(rebate) {}
    
    // Main pricing function based on barrier type
    double price() const {
        switch (barrierType) {
            case BarrierType::UpAndOut:
                return _price_up_and_out();
            case BarrierType::DownAndOut:
                return _price_down_and_out();
            case BarrierType::UpAndIn:
                return priceVanilla() - _price_up_and_out();
            case BarrierType::DownAndIn:
                return priceVanilla() - _price_down_and_out();
            default:
                throw std::runtime_error("Unsupported barrier type");
        }
    }

    // Calculate price of equivalent vanilla option (no barrier)
    double priceVanilla() const {
        return (type == OptionType::Call)
            ? priceVanillaCall()
            : priceVanillaPut();
    }

    // Black-Scholes formula for vanilla call
    double priceVanillaCall() const {
        double d1 = (std::log(S / K) + (r - q + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));
        double d2 = d1 - sigma * std::sqrt(T);
        return S * std::exp(-q * T) * norm_cdf(d1) - K * std::exp(-r * T) * norm_cdf(d2);
    }

    // Black-Scholes formula for vanilla put
    double priceVanillaPut() const {
        double d1 = (std::log(S / K) + (r - q + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));
        double d2 = d1 - sigma * std::sqrt(T);
        return K * std::exp(-r * T) * norm_cdf(-d2) - S * std::exp(-q * T) * norm_cdf(-d1);
    }

private:
    double S, K, H, r, q, sigma, T;
    OptionType type;
    BarrierType barrierType;
    double rebate;

    // Standard normal cumulative distribution function
    static double norm_cdf(double x) {
        return 0.5 * std::erfc(-x / std::sqrt(2));
    }

    // Pricing for up-and-out options (both calls and puts)
    double _price_up_and_out() const {
        double lambda = (r - q + 0.5 * sigma * sigma) / (sigma * sigma);
        double z = std::log(H / S) / (sigma * std::sqrt(T)) + lambda * sigma * std::sqrt(T);
        double x1 = std::log(S / K) / (sigma * std::sqrt(T)) + lambda * sigma * std::sqrt(T);
        double y1 = std::log(H * H / (S * K)) / (sigma * std::sqrt(T)) + lambda * sigma * std::sqrt(T);

        if (type == OptionType::Call) {
            if (K >= H) {
                // Immediately knocked out
                return rebate * std::exp(-r * T);
            }

            // Core pricing components for up-and-out call
            double A = S * std::exp(-q * T) * norm_cdf(x1) - K * std::exp(-r * T) * norm_cdf(x1 - sigma * std::sqrt(T));
            double B = S * std::exp(-q * T) * std::pow(H / S, 2 * lambda) * (norm_cdf(y1) - norm_cdf(z));
            double C = K * std::exp(-r * T) * std::pow(H / S, 2 * lambda - 2) * (norm_cdf(y1 - sigma * std::sqrt(T)) - norm_cdf(z - sigma * std::sqrt(T)));

            double knockOutPrice = A - B + C;
            return std::max(knockOutPrice, 0.0) + rebate * std::exp(-r * T);
        } else { // Put option
            if (K <= H) {
                // Strike below barrier - option worthless
                return rebate * std::exp(-r * T);
            }

            // Core pricing components for up-and-out put
            double A = K * std::exp(-r * T) * norm_cdf(-x1 + sigma * std::sqrt(T)) - S * std::exp(-q * T) * norm_cdf(-x1);
            double B = K * std::exp(-r * T) * std::pow(H / S, 2 * lambda - 2) * norm_cdf(-y1 + sigma * std::sqrt(T));
            double C = S * std::exp(-q * T) * std::pow(H / S, 2 * lambda) * norm_cdf(-y1);

            double knockOutPrice = A - B + C;
            return std::max(knockOutPrice, 0.0) + rebate * std::exp(-r * T);
        }
    }

    // Pricing for down-and-out options (both calls and puts)
    double _price_down_and_out() const {
        double lambda = (r - q + 0.5 * sigma * sigma) / (sigma * sigma);
        double z = std::log(H / S) / (sigma * std::sqrt(T)) + lambda * sigma * std::sqrt(T);
        double x1 = std::log(S / K) / (sigma * std::sqrt(T)) + lambda * sigma * std::sqrt(T);
        double y1 = std::log(H * H / (S * K)) / (sigma * std::sqrt(T)) + lambda * sigma * std::sqrt(T);

        if (type == OptionType::Put) {
            if (K <= H) {
                // Immediately knocked out
                return rebate * std::exp(-r * T);
            }

            // Core pricing components for down-and-out put
            double A = K * std::exp(-r * T) * norm_cdf(-x1 + sigma * std::sqrt(T)) - S * std::exp(-q * T) * norm_cdf(-x1);
            double B = K * std::exp(-r * T) * std::pow(H / S, 2 * lambda - 2) * (norm_cdf(-y1 + sigma * std::sqrt(T)) - norm_cdf(-z + sigma * std::sqrt(T)));
            double C = S * std::exp(-q * T) * std::pow(H / S, 2 * lambda) * (norm_cdf(-y1) - norm_cdf(-z));

            double knockOutPrice = A - B + C;
            return std::max(knockOutPrice, 0.0) + rebate * std::exp(-r * T);
        } else { // Call option
            if (K >= H) {
                // Strike above barrier - option worthless
                return rebate * std::exp(-r * T);
            }

            // Core pricing components for down-and-out call
            double A = S * std::exp(-q * T) * norm_cdf(x1) - K * std::exp(-r * T) * norm_cdf(x1 - sigma * std::sqrt(T));
            double B = S * std::exp(-q * T) * std::pow(H / S, 2 * lambda) * norm_cdf(y1);
            double C = K * std::exp(-r * T) * std::pow(H / S, 2 * lambda - 2) * norm_cdf(y1 - sigma * std::sqrt(T));

            double knockOutPrice = A - B + C;
            return std::max(knockOutPrice, 0.0) + rebate * std::exp(-r * T);
        }
    }
};

int main() {
    try {
        
        // Up-and-Out Call
        BarrierOption upOutCall(100, 95, 120, 0.05, 0.02, 0.2, 1.0,
            BarrierOption::OptionType::Call,
            BarrierOption::BarrierType::UpAndOut,
            2.0);

        // Up-and-In Call
        BarrierOption upInCall(100, 95, 120, 0.05, 0.02, 0.2, 1.0,
            BarrierOption::OptionType::Call,
            BarrierOption::BarrierType::UpAndIn,
            2.0);

        // Up-and-Out Put
        BarrierOption upOutPut(100, 105, 120, 0.05, 0.02, 0.2, 1.0,
            BarrierOption::OptionType::Put,
            BarrierOption::BarrierType::UpAndOut,
            1.5);

        // Up-and-In Put
        BarrierOption upInPut(100, 105, 120, 0.05, 0.02, 0.2, 1.0,
            BarrierOption::OptionType::Put,
            BarrierOption::BarrierType::UpAndIn,
            1.5);

        // Down-and-Out Put
        BarrierOption downOutPut(100, 105, 80, 0.05, 0.02, 0.2, 1.0,
            BarrierOption::OptionType::Put,
            BarrierOption::BarrierType::DownAndOut,
            1.5);

        // Down-and-In Put
        BarrierOption downInPut(100, 105, 80, 0.05, 0.02, 0.2, 1.0,
            BarrierOption::OptionType::Put,
            BarrierOption::BarrierType::DownAndIn,
            1.5);

        // Down-and-Out Call
        BarrierOption downOutCall(100, 95, 80, 0.05, 0.02, 0.2, 1.0,
            BarrierOption::OptionType::Call,
            BarrierOption::BarrierType::DownAndOut,
            2.0);

        // Down-and-In Call
        BarrierOption downInCall(100, 95, 80, 0.05, 0.02, 0.2, 1.0,
            BarrierOption::OptionType::Call,
            BarrierOption::BarrierType::DownAndIn,
            2.0);

        // Output results
        std::cout << "=== CALL OPTIONS ===\n";
        std::cout << "Vanilla Call Price: " << upOutCall.priceVanillaCall() << "\n";
        std::cout << "Up-and-Out Call: " << upOutCall.price() << "\n";
        std::cout << "Up-and-In Call: " << upInCall.price() << "\n";
        std::cout << "Down-and-Out Call: " << downOutCall.price() << "\n";
        std::cout << "Down-and-In Call: " << downInCall.price() << "\n\n";

        std::cout << "=== PUT OPTIONS ===\n";
        std::cout << "Vanilla Put Price: " << downOutPut.priceVanillaPut() << "\n";
        std::cout << "Up-and-Out Put: " << upOutPut.price() << "\n";
        std::cout << "Up-and-In Put: " << upInPut.price() << "\n";
        std::cout << "Down-and-Out Put: " << downOutPut.price() << "\n";
        std::cout << "Down-and-In Put: " << downInPut.price() << "\n\n";

        // In + Out should equal Vanilla
        std::cout << "==================\n";
        std::cout << "Call: Up-In + Up-Out = " << upInCall.price() + upOutCall.price() 
                  << " vs Vanilla = " << upOutCall.priceVanillaCall() << "\n";
        std::cout << "Call: Down-In + Down-Out = " << downInCall.price() + downOutCall.price() 
                  << " vs Vanilla = " << downOutCall.priceVanillaCall() << "\n";
        std::cout << "Put: Up-In + Up-Out = " << upInPut.price() + upOutPut.price() 
                  << " vs Vanilla = " << upOutPut.priceVanillaPut() << "\n";
        std::cout << "Put: Down-In + Down-Out = " << downInPut.price() + downOutPut.price() 
                  << " vs Vanilla = " << downInPut.priceVanillaPut() << "\n";

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
    }

    return 0;
}
