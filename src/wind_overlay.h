#pragma once
#include <cstdint>
#include <vector>
#include "viewport.h"
#include "wind_solver.h"

// ============================================================================
// Wind arrow overlay: which arrows to draw for the terrain wind field
//
// Pure maths in map ("world") units, like viewport.h, so it is part of wf_core and testable without SFML; the
// renderer only turns the arrows into lines.
//
// - Where: on a square lattice fixed to the map, so arrows don't swim about while panning. Its spacing is a
//   power-of-two number of cells chosen so the arrows sit about `spacingPx` screen pixels apart at any zoom
//   (between spacingPx and 2 x spacingPx): zooming in shows more arrows, each sampling a smaller patch.
// - Which way: the direction the 20-ft wind blows TOWARDS in the cell under the lattice point.
// - How long: in proportion to that cell's 20-ft speed, the fastest wind on the map filling 85 % of the
//   lattice spacing, so lengths compare across the map. Very slow wind still gets a short stub (15 %).
// - What colour: the speed against the open wind from the weather, so you see at a glance where the terrain
//   speeds the wind up or slows it down. White = same as the open wind, shading to blue at 30 % slower or
//   less, and to magenta at 30 % faster or more. (Blue and magenta stay distinct from the red and orange of
//   the fire.)
// ============================================================================

constexpr float WIND_ARROW_SPACING_PX = 48.0f;   // Target gap between arrows on screen, pixels

struct WindArrow {
    float x0 = 0, y0 = 0;        // Tail (world units)
    float x1 = 0, y1 = 0;        // Tip
    float hx0 = 0, hy0 = 0;      // Ends of the two arrowhead strokes, each drawn to the tip
    float hx1 = 0, hy1 = 0;
    std::uint8_t r = 255, g = 255, b = 255;   // Colour: speed against the open wind (see above)
    float speedKmh = 0;          // The cell's 20-ft wind speed
};

// Colour for a 20-ft wind of speedKmh when the open wind is openKmh (white / blue / magenta, see above)
void windArrowColour(float speedKmh, float openKmh, std::uint8_t& r, std::uint8_t& g, std::uint8_t& b);

// Lattice spacing in world units for a view showing `perPixel` world units per screen pixel: CELL_SIZE x the
// smallest power of two cells that is at least spacingPx screen pixels (at most the grid's longer side)
float windArrowSpacing(const GridSpec& grid, float perPixel, float spacingPx = WIND_ARROW_SPACING_PX);

// The arrows inside `view` for `field` (on field.grid, one cell = CELL_SIZE world units).
// perPixel: world units per screen pixel. openKmh: the weather's open 20-ft wind. maxKmh: the fastest speed in
// the field (the full-length arrow). Empty when the field is empty.
std::vector<WindArrow> buildWindArrows(const WindField& field, const Viewport& view, float perPixel,
                                       float openKmh, float maxKmh, float spacingPx = WIND_ARROW_SPACING_PX);
