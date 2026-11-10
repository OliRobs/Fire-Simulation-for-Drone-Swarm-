#pragma once
#include <SFML/Graphics.hpp>
#include <cstdint>
#include <vector>
#include "distance_transform.h"
#include "drone.h"
#include "mission_types.h"
#include "types.h"
#include "viewport.h"

// All SFML drawing. One RGBA texture for the grid (only changed pixels are rewritten) plus a second
// transparent texture for the mission overlay, with drones, routes and points of interest drawn on top.
class Renderer {
public:
    explicit Renderer(const std::vector<std::uint8_t>& terrainShade);

    sf::RenderWindow& window() { return window_; }

    bool showPerimeter = true;  // P
    bool showMission = true;    // V: Voronoi overlay, points of interest and routes
    bool showRoads = false;     // R: LANDFIRE roads drawn over the elevation shading (call paintAll after toggling).
                                // Hidden by default; roads still act as fuel breaks either way.

    void paintAll(const std::vector<Cell>& grid);
    void paintCells(const std::vector<Cell>& grid, const std::vector<int>& cells);

    // Replaces the cyan perimeter ring (erases the old one first). Pass the new ring even when hidden.
    void setRing(const std::vector<int>& ring, const std::vector<Cell>& grid);

    // Rebuilds the Voronoi / standoff overlay and the route polylines after a re-plan
    void setMissionVisuals(const MissionPlan& plan, const std::vector<Drone>& drones,
                           const DistanceTransform& liveEdt);

    void draw(const MissionPlan& plan, const std::vector<Drone>& drones);

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

    sf::RenderWindow window_;   // Must be first: textures need the GL context it creates
    std::vector<std::uint8_t> terrainShade_;
    std::vector<sf::Uint8> pixels_;
    sf::Texture texture_;
    sf::Sprite sprite_;

    std::vector<sf::Uint8> overlayPixels_;
    sf::Texture overlayTex_;
    sf::Sprite overlaySprite_;
    std::vector<sf::VertexArray> paths_;

    std::vector<int> ringCells_;
    sf::CircleShape marker_, sensor_, body_;

    Viewport view_;
    bool dragging_ = false;          // Right or middle mouse button held: panning
    sf::Vector2i lastMouse_;         // Mouse position at the previous drag event, window pixels
};
