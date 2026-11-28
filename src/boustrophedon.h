#pragma once
#include <cstdint>
#include <memory>
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
// Marks every cell of `area` (covered is a local bitmap, (y - area.y0) * width + (x - area.x0)) whose centre is
// within LANE_COVER_RADIUS of one of the lane's samples (every RAY_STEP along each run). Used by the gap fill;
// public so the tests can check it against the original disc-per-sample stamp.
void stampLaneCoverage(const SweepLane& lane, const GridRect& area, std::vector<std::uint8_t>& covered);

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

// Remembers every estimated sector lap worked out during one re-plan, so the sector split (which asks for the
// same laps again and again: a binary search of ~20 greedy fills per split, and many splits per re-plan while
// the fleet is sized) works each one out only once.
//
// For each first lane it keeps one SectorLapBuilder and the laps it has produced so far, lap(first, first),
// lap(first, first + 1), ..., extended only as far as a question needs. Every value is the builder's own
// result, so lap() == estimateSectorLength() bit for bit, and greedyLast() returns exactly what the lane-by-lane
// greedy fill (grow while lapWithNext() <= limit) returned. The lanes must not change while the table is in use
// (call reset() after they do, e.g. after loops are joined).
class SectorLapTable {
public:
    explicit SectorLapTable(const std::vector<SweepLane>& lanes);

    // Estimated lap of lanes[first..last] (cells), == estimateSectorLength(lanes, first, last) bit for bit
    float lap(int first, int last);

    // Last lane of a sector that starts at `first` and greedily takes lanes while its lap stays <= limit,
    // stopping at lane `bound` at the latest. The first lane is always taken, so the result is >= first.
    int greedyLast(int first, int bound, float limit);

    // Forgets every lap (the lanes have changed)
    void reset();

private:
    // Laps of the sectors that start at one lane, as far as they have been worked out
    struct Row {
        explicit Row(const std::vector<SweepLane>& lanes, int first) : builder(lanes, first) {}
        SectorLapBuilder builder;          // Holds lanes[first..first + laps.size() - 1]
        std::vector<float> laps;           // laps[k] = lap(first, first + k)
        std::vector<float> runMax;         // runMax[k] = largest of laps[1..k] (runMax[0] unused): never falls,
                                           // so the first lap over a limit can be found by binary search
    };
    Row& row(int first);                   // Creates the row on first use
    void grow(Row& r);                     // Adds the next lane to a row

    const std::vector<SweepLane>& lanes_;
    std::vector<std::unique_ptr<Row>> rows_;   // One per first lane, created lazily
};

// Sector of the nearest lane for every inner-band cell (inside `area`), -1 elsewhere. Found by a breadth-first
// search through the band from every lane's cells, so each cell goes to the lane it can be reached from first.
// `zone` and `owner` have one entry per cell of `grid`. The second form fills `owner`, which must already be -1
// everywhere (e.g. a buffer recycled from an earlier plan and cleared), instead of allocating a new array.
// bandCells: the number of zone-1 cells inside `area` if the caller already knows it (the planner counts them
// while classifying), else -1. When the search reaches that many cells every band cell has an owner, so the last
// stage (a scan of the whole area for band cells cut off from every lane) is skipped: it would change nothing.
void buildOwnerMap(const GridSpec& grid, const std::vector<SweepLane>& lanes, const std::vector<std::uint8_t>& zone,
                   const GridRect& area, std::vector<int>& owner, long long bandCells = -1);
std::vector<int> buildOwnerMap(const GridSpec& grid, const std::vector<SweepLane>& lanes,
                               const std::vector<std::uint8_t>& zone, const GridRect& area, long long bandCells = -1);
// The same map as the first form, without copying `zone`: the flood marks claimed band cells in `zone` as it
// goes (so it reads one byte per neighbour instead of two arrays) and puts every mark back before returning, so
// `zone` is unchanged on return. Used by the planner, which owns its zone map; the forms above call it on a copy.
void buildOwnerMapUsingZone(const GridSpec& grid, const std::vector<SweepLane>& lanes, std::vector<std::uint8_t>& zone,
                            const GridRect& area, std::vector<int>& owner, long long bandCells = -1);
