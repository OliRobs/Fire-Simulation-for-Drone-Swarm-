#include "terrain.h"

Terrain::Terrain(const std::vector<Cell>& grid, int width, int height)
    : grid_(grid), width_(width), height_(height) {}

sf::Color Terrain::GetCellColor(int x, int y) const {
    int idx = y * width_ + x;
    const Cell& cell = grid_[idx];

    // Priority 1: Burning state visualizer
    if (cell.state == 1) return sf::Color(255, 69, 0);   // Burning (Orange/Red)
    if (cell.state == 2) return sf::Color(50, 50, 50);   // Burned (Charcoal)

    // Priority 2: Road fuel break visualizer
    if (cell.is_road) return sf::Color(30, 30, 30);      // Road (Dark Gray)

    // Priority 3: Non-burnable fuel categories
    if (cell.fuel_model == 98) return sf::Color(0, 105, 148); // Water (Blue)
    if (cell.fuel_model == 99) return sf::Color(160, 160, 160); // Barren (Light Gray)

    // Priority 4: Vegetation landscape elevation shading
    uint8_t shade = static_cast<uint8_t>(std::clamp(cell.elevation / 10.0f, 40.0f, 200.0f));
    return sf::Color(34, shade, 34); // Forest/Grassland Green
}
