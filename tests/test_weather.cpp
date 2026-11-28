// ============================================================================
// Tests: weather scenarios (weather.h): defaults, the scenario file format and its error messages
// Part of wf_tests (see test_framework.h for how to run one file or one test).
// ============================================================================
#include <cstdio>
#include <cstring>
#include <sstream>
#include <string>

#include "fire_model.h"
#include "weather.h"
#include "test_framework.h"
#include "test_support.h"

namespace {
    // True if every field of two scenarios is identical to the last bit
    bool sameWeather(const FireWeather& a, const FireWeather& b) {
        return std::memcmp(&a.moisture, &b.moisture, sizeof(FuelMoisture)) == 0 && a.wind20ftKmh == b.wind20ftKmh
            && a.windFromDeg == b.windFromDeg && a.windAdjustmentFactor == b.windAdjustmentFactor
            && a.windLimit == b.windLimit && a.crownFire == b.crownFire && a.foliarMoisturePct == b.foliarMoisturePct
            && a.windStabilityAlpha == b.windStabilityAlpha;
    }

    // Parses `text` as a scenario named "test.ini"; returns the error message, or "" if it parsed
    std::string parseError(const std::string& text) {
        std::istringstream in(text);
        try { parseWeather(in, "test.ini"); }
        catch (const WeatherFileError& e) { return e.what(); }
        return "";
    }
}

// defaultWeather() is config.h's scenario
TEST(default_weather_is_the_config_scenario) {
    const FireWeather w = defaultWeather();
    CHECK(w.wind20ftKmh == WIND_SPEED_20FT_KMH && w.windFromDeg == WIND_FROM_DEG);
    CHECK(w.windAdjustmentFactor == WIND_ADJUSTMENT_FACTOR && w.windLimit == APPLY_WIND_LIMIT);
    CHECK(w.crownFire == ENABLE_CROWN_FIRE && w.foliarMoisturePct == FOLIAR_MOISTURE_PCT);
    CHECK(w.windStabilityAlpha == WIND_STABILITY_ALPHA);
    CHECK(sameWeather(w, FireWeather()));
}

// The shipped scenarios/default.ini holds exactly the built-in defaults, so running with it changes nothing
TEST(default_scenario_file_matches_the_defaults) {
    const FireWeather w = loadWeatherFile(std::string(WF_SOURCE_DIR) + "/scenarios/default.ini");
    CHECK(sameWeather(w, defaultWeather()));
    // ...and a fire set up from it has the same rates as one set up from the defaults
    std::vector<Cell> cells(N);
    const FireSimulation a(G, cells, w), b(G, cells, defaultWeather());
    bool same = true;
    for (int i = 0; i < N && same; i += 1009)
        for (int d = 0; d < NUM_DIRS; ++d) same = same && a.directionalRos(i, d) == b.directionalRos(i, d);
    CHECK(same);
}

// A sample file setting every key: each field gets its value (percent -> fraction for the moistures),
// comments, blank lines, spacing and the different boolean spellings are accepted
TEST(weather_file_sets_every_field) {
    std::istringstream in(
        "# A hot, dry, windy afternoon\n"
        "\n"
        "wind_20ft_kmh = 35.5        ; trailing comment\n"
        "   wind_from_deg=270\n"
        "wind_adjustment_factor = 0.4\n"
        "wind_limit = off\n"
        "wind_stability_alpha = 2.5\n"
        "crown_fire = NO\n"
        "foliar_moisture_pct = 85\n"
        "moisture_1h_pct = 3\n"
        "moisture_10h_pct = 4\n"
        "moisture_100h_pct = 5\n"
        "moisture_live_herb_pct = 30\n"
        "moisture_live_woody_pct = 70\n");
    const FireWeather w = parseWeather(in, "sample.ini");
    CHECK(w.wind20ftKmh == 35.5f);
    CHECK(w.windFromDeg == 270.0f);
    CHECK(w.windAdjustmentFactor == 0.4f);
    CHECK(!w.windLimit);
    CHECK(w.windStabilityAlpha == 2.5f);
    CHECK(!w.crownFire);
    CHECK(w.foliarMoisturePct == 85.0f);
    CHECK(w.moisture.dead1h == 0.03 && w.moisture.dead10h == 0.04 && w.moisture.dead100h == 0.05);
    CHECK(w.moisture.liveHerb == 0.30 && w.moisture.liveWoody == 0.70);
}

// Keys left out keep the base scenario's values
TEST(weather_file_keeps_unset_keys) {
    std::istringstream in("wind_20ft_kmh = 5\n");
    FireWeather base = defaultWeather();
    base.foliarMoisturePct = 120.0f;
    const FireWeather w = parseWeather(in, "partial.ini", base);
    FireWeather expected = base;
    expected.wind20ftKmh = 5.0f;
    CHECK(sameWeather(w, expected));
}

// Every kind of mistake is an error naming the file and the line
TEST(weather_file_errors_give_the_line) {
    const std::string unknown = parseError("wind_20ft_kmh = 10\n# fine\nwind_speed = 10\n");
    std::printf("  (%s)\n", unknown.c_str());
    CHECK(unknown.find("test.ini:3:") == 0 && unknown.find("unknown key 'wind_speed'") != std::string::npos);
    CHECK(parseError("\n\nwind_from_deg 225\n").find("test.ini:3: expected 'key = value'") == 0);
    CHECK(parseError("wind_20ft_kmh = fast\n").find("test.ini:1: 'wind_20ft_kmh' needs a number") == 0);
    CHECK(parseError("wind_20ft_kmh = 10 km/h\n").find("test.ini:1:") == 0);   // Units belong in comments
    CHECK(parseError("wind_from_deg = 400\n").find("test.ini:1: 'wind_from_deg' = 400 is outside 0 .. 360") == 0);
    CHECK(parseError("moisture_1h_pct = 0.06\n").find("test.ini:1:") == 0);    // A fraction, not a percent
    CHECK(parseError("crown_fire = maybe\n").find("test.ini:1: 'crown_fire' needs true or false") == 0);
    CHECK(parseError("wind_limit = on\nwind_limit = off\n").find("test.ini:2: 'wind_limit' is given twice") == 0);
    CHECK(parseError(" = 3\n").find("test.ini:1: missing key") == 0);
    CHECK(parseError("wind_20ft_kmh = nan\n").find("test.ini:1:") == 0);
    // A missing file
    bool threw = false;
    try { loadWeatherFile("no/such/file.ini"); }
    catch (const WeatherFileError& e) { threw = std::string(e.what()).find("cannot open") != std::string::npos; }
    CHECK(threw);
    // An empty file (or only comments) is the base scenario
    CHECK(parseError("# nothing here\n\n").empty());
}

// The help text lists every key the parser accepts
TEST(weather_help_lists_every_key) {
    const std::string help = weatherKeysHelp();
    for (const char* key : { "wind_20ft_kmh", "wind_from_deg", "wind_adjustment_factor", "wind_limit",
                             "wind_stability_alpha", "crown_fire", "foliar_moisture_pct", "moisture_1h_pct",
                             "moisture_10h_pct", "moisture_100h_pct", "moisture_live_herb_pct",
                             "moisture_live_woody_pct" }) {
        CHECK(help.find(key) != std::string::npos);
        std::istringstream in(std::string(key) + (std::strcmp(key, "wind_limit") == 0 || std::strcmp(key, "crown_fire") == 0
                                                   ? " = true\n" : " = 50\n"));
        bool accepted = true;
        try { parseWeather(in, "key.ini"); }
        catch (const WeatherFileError& e) {
            // Some keys don't allow 50; any error must then be about the range, not an unknown key
            accepted = std::string(e.what()).find("outside") != std::string::npos;
        }
        CHECK(accepted);
    }
}
