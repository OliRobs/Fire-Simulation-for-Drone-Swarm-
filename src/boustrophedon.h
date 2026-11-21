#pragma once
#include <cstdint>
#include <vector>
#include "distance_transform.h"
#include "mission_types.h"

// ============================================================================
// Sweep lanes ("boustrophedon", as the ox ploughs) for the inner band
//
// Lanes are seeded every LANE_SPACING (at most) along the perimeter line and run straight in from it towards
// the fire, down the gradient of the committed-fire distance field, i.e. at right angles to the nearest fire
// edge. Each lane keeps only the parts that lie in the inner band (zone == 1), so it normally ends at the
// live-fire standoff. A sweep drone flies a contiguous run of lanes, alternately inward and outward.
// ============================================================================

// A cell counts as seen by a lane when its centre is within this distance of the lane: the sensor radius
// plus half a cell, i.e. the sensor footprint reaches the cell. Lanes 120 m (2 sensor radii) apart then leave
// no cell unseen between them on straight stretches; without the half cell, cells exactly midway would be
// missed by rounding and gap fill would add pointless lanes all along the perimeter line.
constexpr float LANE_COVER_RADIUS = SENSOR_RADIUS + 0.5f;

// Lays one lane per seed along each perimeter-line loop, in order round each loop (lane.loop = loop index).
// zone: per-cell zone from the planner (1 = inner band). Seeds whose straight line never meets the band (the
// fire has grown over that stretch) give no lane.
std::vector<SweepLane> layLanes(const std::vector<std::vector<Vec2>>& loops, const DistanceTransform& perimEdt,
                                const std::vector<std::uint8_t>& zone);

// Gap fill. Where the fire edge is concave (bays), lanes seeded 120 m apart on the perimeter line spread apart
// as they go in, leaving inner cells that no lane sees. For each patch of such cells this adds a short lane
// through the patch, pointing the same way (down the distance-field gradient), clipped to the band, and
// inserts it in the sequence next to the main lane it sits closest to. Repeats until every inner cell inside
// `area` is seen (at most a few passes). Returns the number of lanes added.
int fillCoverageGaps(std::vector<SweepLane>& lanes, const DistanceTransform& perimEdt,
                     const std::vector<std::uint8_t>& zone, const GridRect& area);

// Waypoints of the sector made of lanes[first..last] (inclusive, in sequence order): the first lane flown
// inward, the next outward, and so on, each step to the next lane flown from whichever end the drone is at.
// The closing leg back to the start is not included (safeLoop adds it).
std::vector<Vec2> sectorWaypoints(const std::vector<SweepLane>& lanes, int first, int last);

// Estimated length (cells) of one closed lap of that sector: the lanes, the straight steps between them and
// the straight return leg from the last lane's end to the first lane's start. Fire detours are not included.
float estimateSectorLength(const std::vector<SweepLane>& lanes, int first, int last);

// Builds the estimated lap of a sector one lane at a time, so code that grows a sector lane by lane (the greedy
// sector fill in mission_planner.cpp) pays O(1) per added lane instead of re-adding every lane each time.
// It does exactly the same floating-point additions, in the same order, as estimateSectorLength, so
// lap() == estimateSectorLength(lanes, first, last()) bit for bit.
class SectorLapBuilder {
public:
    // Starts a sector holding lanes[first] alone
    SectorLapBuilder(const std::vector<SweepLane>& lanes, int first);
    int last() const { return last_; }                    // Last lane in the sector so far
    float lap() const;                                    // Estimated closed lap of lanes[first..last()]
    float lapWithNext() const;                            // Same, if lanes[last() + 1] were added (it must exist)
    void addNext();                                       // Adds lanes[last() + 1] to the sector

private:
    // Open length (no return leg) and end point after adding lanes[i] to a sector whose open length is `len`
    // and which ends at `end`
    void extended(int i, float& len, Vec2& end) const;

    const std::vector<SweepLane>& lanes_;
    int first_, last_;
    float len_ = 0.0f;        // Lanes plus the steps between them, without the return leg (cells)
    Vec2 start_, end_;        // Where the lap starts (first lane's entry) and where it currently ends
};

// Sector of the nearest lane for every inner-band cell (inside `area`), -1 elsewhere. Found by a breadth-first
// search through the band from every lane's cells, so each cell goes to the lane it can be reached from first.
// The second form fills `owner`, which must already be GRID_WIDTH x GRID_HEIGHT and -1 everywhere (e.g. a buffer
// recycled from an earlier plan and cleared), instead of allocating a new whole-grid array.
void buildOwnerMap(const std::vector<SweepLane>& lanes, const std::vector<std::uint8_t>& zone,
                   const GridRect& area, std::vector<int>& owner);
std::vector<int> buildOwnerMap(const std::vector<SweepLane>& lanes, const std::vector<std::uint8_t>& zone,
                               const GridRect& area);
