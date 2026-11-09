#include <SFML/Graphics.hpp>
#include <iostream>
#include "raster_io.h"
#include "fire_model.h"
#include "terrain.h"

int main() {
    int width = 0;
    int height = 0;
    std::vector<Cell> grid;

    std::cout << "Loading LANDFIRE GeoTIFF layers..." << std::endl;
    if (!LoadLandfireDataset("elevation.tif", "slope.tif", "fuel.tif", "roads.tif", grid, width, height)) {
        std::cerr << "Failed to load landscape layers." << std::endl;
        return -1;
    }

    std::cout << "Loaded dataset dimensions: " << width << "x" << height << std::endl;

    FireSimulation fire_sim(grid, width, height);
    Terrain terrain(grid, width, height);

    sf::RenderWindow window(sf::VideoMode(1280, 720), "Wildfire Simulation - LANDFIRE Integration");
    sf::Clock clock;

    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();

            if (event.type == sf::Event::MouseButtonPressed && event.mouseButton.button == sf::Mouse::Left) {
                // Ignite clicked cell position
                int mx = event.mouseButton.x;
                int my = event.mouseButton.y;
                fire_sim.Ignite(mx, my);
            }
        }

        float dt = clock.restart().asSeconds();
        fire_sim.Update(dt, 5.0f, 45.0f); // Update with 5 m/s wind at 45 deg

        window.clear(sf::Color::Black);
        // Drawing logic ...
        window.display();
    }

    return 0;
}
