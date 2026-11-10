#pragma once
#include <cmath>
#include <cstdint>
#include "config.h"

// Plain 2D vector in grid units (cell centres are at +0.5). Deliberately not sf::Vector2f, so the
// simulation and planning code builds and runs without SFML.
struct Vec2 {
    float x = 0.0f, y = 0.0f;
    constexpr Vec2() = default;
    constexpr Vec2(float x_, float y_) : x(x_), y(y_) {}
    constexpr Vec2 operator+(Vec2 o) const { return { x + o.x, y + o.y }; }
    constexpr Vec2 operator-(Vec2 o) const { return { x - o.x, y - o.y }; }
    constexpr Vec2 operator*(float s) const { return { x * s, y * s }; }
    constexpr Vec2 operator/(float s) const { return { x / s, y / s }; }
    Vec2& operator+=(Vec2 o) { x += o.x; y += o.y; return *this; }
};

inline float sqDist(Vec2 a, Vec2 b) {
    const float dx = a.x - b.x, dy = a.y - b.y;
    return dx * dx + dy * dy;
}

inline Vec2 cellCentre(int idx) { return Vec2(idx % GRID_WIDTH + 0.5f, idx / GRID_WIDTH + 0.5f); }

enum class CellState : std::uint8_t { Unburned, Burning, Burned };

// One 30 m patch of landscape: the LANDFIRE layers for that patch plus its fire state.
struct Cell {
    float elevation = 0.0f;  // Metres (LANDFIRE Elev)
    float slopeDeg = 0.0f;   // Slope steepness in degrees (LANDFIRE SlpD, or derived from elevation)
    float upslopeX = 0.0f;   // Unit vector pointing uphill, grid frame (x = east, y = south);
    float upslopeY = 0.0f;   //   (0, 0) on flat ground. Derived from elevation (the aspect).
    float burnTime = 0.0f;   // Minutes spent burning
    std::int16_t fuelModel = DEFAULT_FUEL_MODEL; // FBFM40 fuel model number (Scott & Burgan); 91-99 = nonburnable
    std::uint8_t roadClass = 0;                  // LANDFIRE road class 20-23, 0 = no road
    CellState state = CellState::Unburned;
};
