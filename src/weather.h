#pragma once
#include <istream>
#include <stdexcept>
#include <string>
#include "fire_model.h"   // FireWeather

// ============================================================================
// Weather scenario: the fire-weather inputs, from config.h's defaults or from a file
//
// The weather is not in the LANDFIRE data (wind, fuel moisture, foliar moisture, atmospheric stability), so the
// app runs one fixed scenario. defaultWeather() is config.h's placeholder scenario; a scenario file (see
// scenarios/default.ini) replaces any of its values without a rebuild:
//
//     # comment                      ; also a comment
//     wind_20ft_kmh = 20             key = value, one per line; blank lines are ignored
//     crown_fire    = true           true/false, yes/no, on/off or 1/0
//
// Keys left out keep their default. Unknown keys, a key given twice, a line without '=', a value that is not
// a number (or a boolean where one is expected) and values out of range are all errors, reported with the
// file name and line number (WeatherFileError), so a typo never silently runs a different scenario.
// Moistures are written in percent in the file and stored as fractions in FireWeather (6 -> 0.06).
// ============================================================================

// Error in a scenario file: the message starts with "<file>:<line>: "
class WeatherFileError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

// config.h's placeholder scenario (the values FireWeather starts with)
FireWeather defaultWeather();

// Reads a scenario from `in` on top of `base`. `name` is used in error messages (normally the file name).
// Throws WeatherFileError.
FireWeather parseWeather(std::istream& in, const std::string& name, const FireWeather& base = defaultWeather());

// Reads the scenario file at `path` on top of `base`. Throws WeatherFileError (including when the file cannot
// be opened).
FireWeather loadWeatherFile(const std::string& path, const FireWeather& base = defaultWeather());

// The keys a scenario file may use, one per line with units, for the --help text and the README
std::string weatherKeysHelp();
