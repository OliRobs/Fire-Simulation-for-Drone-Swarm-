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

struct Cell {
    float elevation = 0.0f; // Altitude in metres
    float slope_deg = 0.0f; // slop in degrees 
    int fuel_model = 0;      // Scott & Burgan code 
    bool is_road = flase; // Operational road barrier flag 

    //Existing Fire state properties 
    uint8_t state = 0; // 0: unburned, 1: burning, 2: burned 
    float burnTime = 0.0f;  // Ticks spent burning
};
