#pragma once
#include <SFML/Graphics.hpp>
#include <cstdint>
#include <string>
#include <vector>
#include "distance_transform.h"
#include "drone.h"
#include "mission_types.h"
#include "types.h"
#include "viewport.h"
#include "wind_solver.h"

// All SFML drawing. One RGBA texture for the grid (only changed pixels are rewritten) plus a second
// transparent texture for the mission overlay, with drones, routes and points of interest drawn on top.
class Renderer {
public:
    // A window showing grid `spec` (one pixel per cell at zoom 1 when CELL_SIZE is 1), shaded by terrainShade
    // (one value per cell)
    Renderer(const GridSpec& spec, const std::vector<std::uint8_t>& terrainShade);

    sf::RenderWindow& window() { return window_; }

    bool showPerimeter = true;  // P
    bool showMission = true;    // V: Voronoi overlay, points of interest and routes
    bool showRoads = false;     // R: LANDFIRE roads drawn over the elevation shading (call paintAll after toggling).
                                // Hidden by default; roads still act as fuel breaks either way.
    bool showWind = false;      // W: terrain wind arrows (wind_overlay.h). Hidden by default; needs setWindField.

    void paintAll(const std::vector<Cell>& grid);
    void paintCells(const std::vector<Cell>& grid, const std::vector<int>& cells);

    // Replaces the cyan perimeter ring (erases the old one first). Pass the new ring even when hidden.
    void setRing(const std::vector<int>& ring, const std::vector<Cell>& grid);

    // Rebuilds the sector / Voronoi / standoff overlay, the lanes and the route polylines after a re-plan, and
    // the fleet status line
    void setMissionVisuals(const MissionPlan& plan, const std::vector<Drone>& drones,
                           const DistanceTransform& liveEdt);

    // The terrain wind field to draw as arrows when showWind is on, and the weather's open 20-ft wind the
    // arrow colours compare against. The field is not copied: it must outlive the renderer's use of it.
    // nullptr = no field (USE_TERRAIN_WIND off): nothing is drawn.
    void setWindField(const WindField* field, float openKmh);
    bool hasWindField() const { return wind_ != nullptr; }
    float windMaxKmh() const { return windMaxKmh_; }   // Fastest 20-ft wind in the field (the longest arrow)
    float windMinKmh() const { return windMinKmh_; }   // Slowest

    void draw(const MissionPlan& plan, const std::vector<Drone>& drones);

    // Drone counts from the last setMissionVisuals (sweep with longest lap, ACO with revisit time, drones
    // returning to base) and a warning when the fleet is at its maximum but over a time limit. For the title bar.
    const std::string& fleetStatus() const { return fleetStatus_; }

    // --- Zoom and pan ---
    // Handles mouse wheel (zoom at cursor), right/middle drag (pan), +/- (zoom), arrow keys (pan) and
    // 0/Home (whole map). Returns true if the event was used, so the caller can skip it.
    bool handleViewEvent(const sf::Event& event);
    // Zooms to fit the non-zero cells of mask (e.g. the fire), with a margin of a few cells
    void frameMask(const std::vector<std::uint8_t>& mask);
    float zoom() const { return view_.zoom(); }

private:
    void applyView();   // Pushes view_ to the window as an sf::View
    void setColor(int idx, std::uint8_t r, std::uint8_t g, std::uint8_t b);
    void paintBase(int idx, const std::vector<Cell>& grid);

    const int W_, H_;           // Grid size, cells (plain ints: safe to set before the window)
    sf::RenderWindow window_;   // Must come before the textures: they need the GL context it creates
    std::vector<std::uint8_t> terrainShade_;
    std::vector<sf::Uint8> pixels_;
    sf::Texture texture_;
    sf::Sprite sprite_;

    std::vector<sf::Uint8> overlayPixels_;
    sf::Texture overlayTex_;
    sf::Sprite overlaySprite_;
    std::vector<sf::VertexArray> paths_;
    sf::VertexArray lanes_;     // Sweep lanes of the current plan (pairs of points)
    std::string fleetStatus_;

    // Rows of pixels_ changed since the last upload, [dirtyY0_, dirtyY1_) (empty when dirtyY1_ <= dirtyY0_).
    // draw() sends only these rows to the GPU instead of the whole image (4.8 MB on a 1200 x 1000 grid) every frame.
    int dirtyY0_ = 0, dirtyY1_ = 0;
    GridRect overlayDrawn_;   // Part of overlayPixels_ the last setMissionVisuals() may have written

    void drawWindArrows(float perPixel);   // The wind arrows inside the view (showWind)

    const WindField* wind_ = nullptr;   // Terrain wind field to draw (owned by main)
    float windOpenKmh_ = 0.0f;          // Open 20-ft wind: arrows this fast are white
    float windMinKmh_ = 0.0f, windMaxKmh_ = 0.0f;   // Speed range of the field

    std::vector<int> ringCells_;
    sf::CircleShape marker_, sensor_, body_;

    Viewport view_;
    bool dragging_ = false;          // Right or middle mouse button held: panning
    sf::Vector2i lastMouse_;         // Mouse position at the previous drag event, window pixels
};
