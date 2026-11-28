#include "weather.h"
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <set>
#include <sstream>
#include <vector>

// ============================================================================
// Weather scenario files (see weather.h for the format)
// ============================================================================

namespace {
    // ------------------------------------------------------------------------
    // The keys: name, units / meaning, allowed range and where the value goes
    // ------------------------------------------------------------------------
    // A key is either a number (stored through setNumber, checked against [lo, hi]) or a boolean (setFlag).
    struct Key {
        const char* name;
        const char* help;                                      // Units and meaning, for weatherKeysHelp()
        bool isFlag;
        double lo, hi;                                         // Allowed range (numbers only)
        std::function<void(FireWeather&, double)> setNumber;   // Numbers
        std::function<void(FireWeather&, bool)> setFlag;       // Booleans
    };

    // Percent in the file, fraction in FuelMoisture (6 -> 0.06: the quotient of two exact values, so it is the
    // same double as the literal 0.06 in rothermel.h)
    double fromPercent(double pct) { return pct / 100.0; }

    const std::vector<Key>& keys() {
        static const std::vector<Key> k = {
            { "wind_20ft_kmh", "open wind 20 ft (6.1 m) above the vegetation, km/h", false, 0.0, 150.0,
              [](FireWeather& w, double v) { w.wind20ftKmh = static_cast<float>(v); }, nullptr },
            { "wind_from_deg", "direction the wind blows FROM, degrees clockwise from north (225 = SW)", false, 0.0, 360.0,
              [](FireWeather& w, double v) { w.windFromDeg = static_cast<float>(v); }, nullptr },
            { "wind_adjustment_factor", "midflame / 20-ft wind; 0 = worked out from fuel depth and canopy (Andrews 2012)",
              false, 0.0, 1.0, [](FireWeather& w, double v) { w.windAdjustmentFactor = static_cast<float>(v); }, nullptr },
            { "wind_limit", "cap the effective wind at Rothermel's maximum reliable wind (RMRS-GTR-371 3.2.7), true/false",
              true, 0, 0, nullptr, [](FireWeather& w, bool v) { w.windLimit = v; } },
            { "wind_stability_alpha", "terrain wind stability alpha_H / alpha_V: 1 neutral, < 1 stable, > 1 unstable",
              false, 0.05, 10.0, [](FireWeather& w, double v) { w.windStabilityAlpha = static_cast<float>(v); }, nullptr },
            { "crown_fire", "allow torching and crown fire where LANDFIRE maps canopy, true/false", true, 0, 0, nullptr,
              [](FireWeather& w, bool v) { w.crownFire = v; } },
            { "foliar_moisture_pct", "live canopy foliage moisture, % (normal range 75-150)", false, 30.0, 300.0,
              [](FireWeather& w, double v) { w.foliarMoisturePct = static_cast<float>(v); }, nullptr },
            { "moisture_1h_pct", "dead 1-h fuel moisture, %", false, 1.0, 60.0,
              [](FireWeather& w, double v) { w.moisture.dead1h = fromPercent(v); }, nullptr },
            { "moisture_10h_pct", "dead 10-h fuel moisture, %", false, 1.0, 60.0,
              [](FireWeather& w, double v) { w.moisture.dead10h = fromPercent(v); }, nullptr },
            { "moisture_100h_pct", "dead 100-h fuel moisture, %", false, 1.0, 60.0,
              [](FireWeather& w, double v) { w.moisture.dead100h = fromPercent(v); }, nullptr },
            { "moisture_live_herb_pct", "live herbaceous fuel moisture, % (30 = fully cured)", false, 30.0, 300.0,
              [](FireWeather& w, double v) { w.moisture.liveHerb = fromPercent(v); }, nullptr },
            { "moisture_live_woody_pct", "live woody fuel moisture, %", false, 30.0, 300.0,
              [](FireWeather& w, double v) { w.moisture.liveWoody = fromPercent(v); }, nullptr },
        };
        return k;
    }

    // Text without leading and trailing spaces and tabs
    std::string trim(const std::string& s) {
        const size_t a = s.find_first_not_of(" \t\r");
        if (a == std::string::npos) return "";
        const size_t b = s.find_last_not_of(" \t\r");
        return s.substr(a, b - a + 1);
    }

    // "<name>:<line>: <message>" as a WeatherFileError
    [[noreturn]] void fail(const std::string& name, int line, const std::string& message) {
        throw WeatherFileError(name + ":" + std::to_string(line) + ": " + message);
    }

    // A whole token as a finite number, else false
    bool parseNumber(const std::string& text, double& out) {
        if (text.empty()) return false;
        char* end = nullptr;
        out = std::strtod(text.c_str(), &end);
        return end == text.c_str() + text.size() && std::isfinite(out);
    }

    // A boolean word, else false
    bool parseFlag(const std::string& text, bool& out) {
        std::string t;
        for (char c : text) t += static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        if (t == "true" || t == "yes" || t == "on" || t == "1") { out = true; return true; }
        if (t == "false" || t == "no" || t == "off" || t == "0") { out = false; return true; }
        return false;
    }
}

// ----------------------------------------------------------------------------
// Defaults
// ----------------------------------------------------------------------------
FireWeather defaultWeather() { return FireWeather(); }   // FireWeather's member defaults are config.h's values

// ----------------------------------------------------------------------------
// Parse a scenario: one "key = value" per line on top of `base`
// ----------------------------------------------------------------------------
FireWeather parseWeather(std::istream& in, const std::string& name, const FireWeather& base) {
    FireWeather w = base;
    std::set<std::string> seen;   // Keys given so far, to refuse a key given twice
    std::string raw;
    int lineNo = 0;
    while (std::getline(in, raw)) {
        ++lineNo;
        // Comments run from '#' or ';' to the end of the line
        const size_t hash = raw.find_first_of("#;");
        const std::string line = trim(hash == std::string::npos ? raw : raw.substr(0, hash));
        if (line.empty()) continue;

        const size_t eq = line.find('=');
        if (eq == std::string::npos) fail(name, lineNo, "expected 'key = value', got '" + line + "'");
        const std::string key = trim(line.substr(0, eq)), value = trim(line.substr(eq + 1));
        if (key.empty()) fail(name, lineNo, "missing key before '='");

        const Key* k = nullptr;
        for (const Key& candidate : keys()) if (key == candidate.name) k = &candidate;
        if (!k) fail(name, lineNo, "unknown key '" + key + "' (known keys: see scenarios/default.ini)");
        if (!seen.insert(key).second) fail(name, lineNo, "'" + key + "' is given twice");

        if (k->isFlag) {
            bool v = false;
            if (!parseFlag(value, v)) fail(name, lineNo, "'" + key + "' needs true or false, got '" + value + "'");
            k->setFlag(w, v);
        }
        else {
            double v = 0.0;
            if (!parseNumber(value, v)) fail(name, lineNo, "'" + key + "' needs a number, got '" + value + "'");
            if (v < k->lo || v > k->hi) {
                std::ostringstream range;
                range << "'" << key << "' = " << value << " is outside " << k->lo << " .. " << k->hi;
                fail(name, lineNo, range.str());
            }
            k->setNumber(w, v);
        }
    }
    return w;
}

// ----------------------------------------------------------------------------
// Read a scenario file
// ----------------------------------------------------------------------------
FireWeather loadWeatherFile(const std::string& path, const FireWeather& base) {
    std::ifstream in(path);
    if (!in) throw WeatherFileError(path + ": cannot open the weather file");
    return parseWeather(in, path, base);
}

// ----------------------------------------------------------------------------
// The keys, for help text
// ----------------------------------------------------------------------------
std::string weatherKeysHelp() {
    std::ostringstream out;
    for (const Key& k : keys()) {
        out << "  " << k.name << ": " << k.help;
        if (!k.isFlag) out << " [" << k.lo << " .. " << k.hi << "]";
        out << "\n";
    }
    return out.str();
}
