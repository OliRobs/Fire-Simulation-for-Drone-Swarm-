#pragma once
#include <cstddef>
#include <stdexcept>
#include <string>
#include "config.h"
#include "types.h"

// ============================================================================
// GridSpec: the size of the simulation grid, as a value passed at run time
//
// The grid is width x height cells, stored row by row (index = y * width + x), each cellSizeX x cellSizeY
// metres on the ground. Everything in wf_core that is sized by the grid takes a GridSpec (or reads it from an
// object built with one, such as a DistanceTransform), so tests can run on small or non-square grids and a
// different area needs no recompile. config.h's GRID_WIDTH x GRID_HEIGHT is only the app's default size
// (main.cpp, wf_bench).
//
// The cell size is here for the modules that work in metres (fire spread, terrain, wind). The planning and
// drone distances in config.h are converted to cells with cellsFromMetres from SIM_CELL_SIZE_M at compile
// time; a GridSpec with another cell size does not change them.
// ============================================================================
struct GridSpec {
    int width = 0, height = 0;                 // Cells
    float cellSizeX = SIM_CELL_SIZE_M;         // Ground size of one cell, metres (east-west)
    float cellSizeY = SIM_CELL_SIZE_M;         // (north-south)

    // Row-major index of cell (x, y)
    int index(int x, int y) const { return y * width + x; }
    // Column and row of a cell index
    int xOf(int idx) const { return idx % width; }
    int yOf(int idx) const { return idx / width; }
    // True if (x, y) is a cell of the grid
    bool inside(int x, int y) const { return x >= 0 && x < width && y >= 0 && y < height; }
    // Number of cells
    std::size_t count() const { return static_cast<std::size_t>(width) * static_cast<std::size_t>(height); }
    // Centre of cell idx in grid units (cell centres are at +0.5)
    Vec2 centre(int idx) const { return Vec2(idx % width + 0.5f, idx / width + 0.5f); }
    // The whole grid as a rectangle
    GridRect rect() const { return GridRect{ 0, 0, width, height }; }

    bool operator==(const GridSpec& o) const {
        return width == o.width && height == o.height && cellSizeX == o.cellSizeX && cellSizeY == o.cellSizeY;
    }
    bool operator!=(const GridSpec& o) const { return !(*this == o); }
};

// Throws std::invalid_argument unless `values` has exactly one entry per cell of `spec`. Public entry points
// that take a per-cell vector call this, so a vector built for another grid fails loudly instead of being read
// out of bounds. `what` names the argument in the message.
template <class Vector>
void requireGridSize(const Vector& values, const GridSpec& spec, const char* what) {
    if (values.size() != spec.count())
        throw std::invalid_argument(std::string(what) + ": " + std::to_string(values.size()) + " values for a "
                                    + std::to_string(spec.width) + " x " + std::to_string(spec.height) + " grid");
}
