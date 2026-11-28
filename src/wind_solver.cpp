#include "wind_solver.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include "fire_model.h"    // windTowardGrid
#include "fuel_models.h"
#include "rothermel.h"     // FT_PER_M

// ============================================================================
// Mass-conserving terrain wind solver. The method, equations and sources are in wind_solver.h; this file
// follows its six steps. Units inside: metres, seconds, m/s. Axes: X east, Y south (the fire grid's frame,
// so a wind vector (u, v) is directly a grid-frame vector), Z up.
// ============================================================================

namespace {
    constexpr double WIND_REF_HEIGHT_M = 6.096;   // 20 ft: the reference wind height above the vegetation
    constexpr double VON_KARMAN = 0.4;            // As in WindNinja's mass solver (initialize.cpp)
    constexpr double MIN_ROUGHNESS_M = 0.01;      // z0 of short grass / bare ground (WindNinja's "grass" value)
    constexpr double KMH_PER_MS = 3.6;

    // 2 x 2 x 2 Gauss points sit at +-1/sqrt(3) in each natural coordinate, all with weight 1
    const double GAUSS = 1.0 / std::sqrt(3.0);

    // ------------------------------------------------------------------------
    // One mesh column: ground elevation and the vegetation's roughness
    // ------------------------------------------------------------------------
    struct Column {
        double ground = 0.0;    // Mean ground elevation of the fire cells in the column, m
        double vegH = 0.0;      // Mean vegetation height h (canopy, else fuel bed depth), m
        double z0 = MIN_ROUGHNESS_M;   // Roughness length 0.13 h, at least MIN_ROUGHNESS_M
        double d = 0.0;         // Zero-plane displacement 0.64 h
        double zLow = 0.0;      // Below this height above ground the profile is linear down to 0 (d + 7 z0)
        double zTop = 0.0;      // Above this height the speed stops increasing (d + boundary-layer height)
    };

    // ------------------------------------------------------------------------
    // Shape of the neutral log wind profile, S(agl) (step 2)
    // ------------------------------------------------------------------------
    // The wind at height agl above the ground is U(agl) = U20 * S(agl) / S(6.1 m + h), which is exactly U20
    // at the reference height. S is the log law ln((agl - d) / z0), with two patches WindNinja also makes
    // (windProfile.cpp): straight down to 0 at the ground below d + 7 z0 (the AERMOD convention, the log law
    // is meaningless inside the roughness), and constant above the top of the boundary layer.
    double profileShape(const Column& c, double agl) {
        if (agl <= 0.0) return 0.0;
        const double a = std::min(agl, c.zTop);
        if (a < c.zLow) return std::log((c.zLow - c.d) / c.z0) * a / c.zLow;
        return std::log((a - c.d) / c.z0);
    }

    // ------------------------------------------------------------------------
    // Trilinear hexahedron geometry at one point (steps 4 and 5)
    // ------------------------------------------------------------------------
    // Corner a (0..7) of an element sits at offsets (a & 1, (a >> 1) & 1, (a >> 2) & 1) from its lowest
    // corner, i.e. natural coordinates xi_a, eta_a, zeta_a = -1 or +1.
    struct HexPoint {
        double N[8];            // Shape functions N_a at the point
        double dN[8][3];        // dN_a/dX, dN_a/dY, dN_a/dZ (physical derivatives)
        double detJ = 0.0;      // Jacobian determinant: volume of the element per unit natural volume
        double pos[3] = { 0.0, 0.0, 0.0 };   // Physical position of the point
    };

    // Shape functions and their derivatives in natural coordinates at one point (xi, eta, zeta). They depend
    // only on the point, not on the element, so the solver works them out once for each of the 8 Gauss points.
    struct NaturalPoint {
        double N[8];            // N_a
        double dNdXi[8][3];     // dN_a/dxi, dN_a/deta, dN_a/dzeta
    };

    NaturalPoint naturalPoint(double xi, double eta, double zeta) {
        NaturalPoint np;
        for (int a = 0; a < 8; ++a) {
            const double sx = (a & 1) ? 1.0 : -1.0, sy = ((a >> 1) & 1) ? 1.0 : -1.0, sz = ((a >> 2) & 1) ? 1.0 : -1.0;
            const double fx = 1.0 + sx * xi, fy = 1.0 + sy * eta, fz = 1.0 + sz * zeta;
            np.N[a] = 0.125 * fx * fy * fz;
            np.dNdXi[a][0] = 0.125 * sx * fy * fz;
            np.dNdXi[a][1] = 0.125 * fx * sy * fz;
            np.dNdXi[a][2] = 0.125 * fx * fy * sz;
        }
        return np;
    }

    // The 2 x 2 x 2 Gauss points of an element, g = 0..7: xi, eta, zeta = +-GAUSS from bits 0, 1, 2 of g
    struct GaussTable {
        NaturalPoint p[8];
        GaussTable() {
            for (int g = 0; g < 8; ++g)
                p[g] = naturalPoint((g & 1) ? GAUSS : -GAUSS, ((g >> 1) & 1) ? GAUSS : -GAUSS,
                                    ((g >> 2) & 1) ? GAUSS : -GAUSS);
        }
    };

    // One physical coordinate (x, y or z) of the 8 corners mapped through a point's shape functions: the
    // coordinate of the point, sum N_a c_a, and its row of the Jacobian, d c / d xi_k = sum c_a dN_a/dxi_k.
    // Each sum runs over the corners in order a = 0..7, starting from 0.
    struct JacobianRow {
        double pos = 0.0;
        double J[3] = { 0.0, 0.0, 0.0 };
    };

    JacobianRow jacobianRow(const NaturalPoint& np, const double coord[8]) {
        JacobianRow row;
        for (int a = 0; a < 8; ++a) {
            row.pos += np.N[a] * coord[a];
            for (int c = 0; c < 3; ++c) row.J[c] += coord[a] * np.dNdXi[a][c];
        }
        return row;
    }

    // Completes hp at one point from the point's natural shape functions and the three Jacobian rows (x, y, z).
    // The Jacobian J[r][c] = d x_r / d xi_c; physical derivatives are dN/dx = J^-T dN/dxi.
    void finishHex(const NaturalPoint& np, const JacobianRow& rx, const JacobianRow& ry, const JacobianRow& rz,
                   HexPoint& hp) {
        const double(&dNdXi)[8][3] = np.dNdXi;
        for (int a = 0; a < 8; ++a) hp.N[a] = np.N[a];
        hp.pos[0] = rx.pos;
        hp.pos[1] = ry.pos;
        hp.pos[2] = rz.pos;
        const double J[3][3] = { { rx.J[0], rx.J[1], rx.J[2] }, { ry.J[0], ry.J[1], ry.J[2] }, { rz.J[0], rz.J[1], rz.J[2] } };
        // Inverse of J by cofactors
        const double det = J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1])
                         - J[0][1] * (J[1][0] * J[2][2] - J[1][2] * J[2][0])
                         + J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);
        hp.detJ = det;
        const double inv = 1.0 / det;
        double Ji[3][3];
        Ji[0][0] = (J[1][1] * J[2][2] - J[1][2] * J[2][1]) * inv;
        Ji[0][1] = (J[0][2] * J[2][1] - J[0][1] * J[2][2]) * inv;
        Ji[0][2] = (J[0][1] * J[1][2] - J[0][2] * J[1][1]) * inv;
        Ji[1][0] = (J[1][2] * J[2][0] - J[1][0] * J[2][2]) * inv;
        Ji[1][1] = (J[0][0] * J[2][2] - J[0][2] * J[2][0]) * inv;
        Ji[1][2] = (J[0][2] * J[1][0] - J[0][0] * J[1][2]) * inv;
        Ji[2][0] = (J[1][0] * J[2][1] - J[1][1] * J[2][0]) * inv;
        Ji[2][1] = (J[0][1] * J[2][0] - J[0][0] * J[2][1]) * inv;
        Ji[2][2] = (J[0][0] * J[1][1] - J[0][1] * J[1][0]) * inv;
        // dN/dx_r = sum over c of dN/dxi_c * (J^-1)[c][r]
        for (int a = 0; a < 8; ++a)
            for (int r = 0; r < 3; ++r)
                hp.dN[a][r] = dNdXi[a][0] * Ji[0][r] + dNdXi[a][1] * Ji[1][r] + dNdXi[a][2] * Ji[2][r];
    }

    // Fills hp for the point (xi, eta, zeta) of the element whose corner coordinates are xyz[a][0..2]
    void evalHex(const double xyz[8][3], double xi, double eta, double zeta, HexPoint& hp) {
        const NaturalPoint np = naturalPoint(xi, eta, zeta);
        double c[3][8];
        for (int a = 0; a < 8; ++a)
            for (int r = 0; r < 3; ++r) c[r][a] = xyz[a][r];
        finishHex(np, jacobianRow(np, c[0]), jacobianRow(np, c[1]), jacobianRow(np, c[2]), hp);
    }

    // ------------------------------------------------------------------------
    // Gauss-point geometry of every element, x and y parts shared down each column
    // ------------------------------------------------------------------------
    // The mesh columns stand on a regular x-y grid, so the x and y rows of an element's Jacobian (and the x
    // and y of its Gauss points) depend only on the element's i (for x) or j (for y), never on k or on the
    // terrain. They are worked out once per i and per j (8 Gauss points each) and reused for every element;
    // only the z row is computed per element. Each row is the same sum as evalHex computes, so every HexPoint
    // comes out bit-identical to evalHex at the same Gauss point. This halves the cost of the two passes that
    // visit every Gauss point (assembly and gradient smoothing).
    struct ElementGeometry {
        GaussTable gauss;
        std::vector<JacobianRow> xRows, yRows;   // [i * 8 + g], [j * 8 + g]
    };

    // ------------------------------------------------------------------------
    // The wind mesh: columns of terrain-following nodes (step 1)
    // ------------------------------------------------------------------------
    struct Mesh {
        int nx = 0, ny = 0, nz = 0;   // Nodes along X, Y and up each column
        double hx = 0.0, hy = 0.0;    // Column spacing, m
        std::vector<Column> cols;     // nx * ny columns, row-major (j * nx + i)
        std::vector<double> z;        // Node elevations, nx * ny * nz (index node(i, j, k))

        std::size_t node(int i, int j, int k) const {
            return (static_cast<std::size_t>(k) * ny + j) * nx + i;
        }
        std::size_t count() const { return static_cast<std::size_t>(nx) * ny * nz; }
        // Phi = 0 (Dirichlet) on the four sides and the top; the ground (k = 0) is left natural
        bool boundary(int i, int j, int k) const {
            return i == 0 || j == 0 || i == nx - 1 || j == ny - 1 || k == nz - 1;
        }
        // Corner coordinates of element (i, j, k), corner order as in evalHex
        void elementCorners(int i, int j, int k, double xyz[8][3], std::size_t ids[8]) const {
            for (int a = 0; a < 8; ++a) {
                const int ci = i + (a & 1), cj = j + ((a >> 1) & 1), ck = k + ((a >> 2) & 1);
                ids[a] = node(ci, cj, ck);
                xyz[a][0] = (ci + 0.5) * hx;
                xyz[a][1] = (cj + 0.5) * hy;
                xyz[a][2] = z[ids[a]];
            }
        }
    };

    // Builds the shared x and y Jacobian rows for mesh m (see ElementGeometry). Corner x coordinates are
    // (ci + 0.5) hx with ci = i + (a & 1), as Mesh::elementCorners sets them; y likewise with j.
    ElementGeometry elementGeometry(const Mesh& m) {
        ElementGeometry eg;
        eg.xRows.resize(static_cast<std::size_t>(m.nx) * 8);
        eg.yRows.resize(static_cast<std::size_t>(m.ny) * 8);
        for (int i = 0; i + 1 < m.nx; ++i) {
            double cx[8];
            for (int a = 0; a < 8; ++a) cx[a] = (i + (a & 1) + 0.5) * m.hx;
            for (int g = 0; g < 8; ++g) eg.xRows[static_cast<std::size_t>(i) * 8 + g] = jacobianRow(eg.gauss.p[g], cx);
        }
        for (int j = 0; j + 1 < m.ny; ++j) {
            double cy[8];
            for (int a = 0; a < 8; ++a) cy[a] = (j + ((a >> 1) & 1) + 0.5) * m.hy;
            for (int g = 0; g < 8; ++g) eg.yRows[static_cast<std::size_t>(j) * 8 + g] = jacobianRow(eg.gauss.p[g], cy);
        }
        return eg;
    }

    // HexPoint at Gauss point g of element (i, j, *) whose corners have z coordinates xyz[a][2]: the same
    // values evalHex(xyz, +-GAUSS, +-GAUSS, +-GAUSS, hp) gives, from the shared x and y rows
    void gaussHex(const ElementGeometry& eg, int i, int j, int g, const double xyz[8][3], HexPoint& hp) {
        double cz[8];
        for (int a = 0; a < 8; ++a) cz[a] = xyz[a][2];
        const NaturalPoint& np = eg.gauss.p[g];
        finishHex(np, eg.xRows[static_cast<std::size_t>(i) * 8 + g], eg.yRows[static_cast<std::size_t>(j) * 8 + g],
                  jacobianRow(np, cz), hp);
    }

    // Vegetation height of one fire cell: the tree canopy where LANDFIRE maps one, else the surface fuel bed
    double vegetationHeightM(const Cell& c) {
        if (c.canopyCoverPct > 0.0f && c.canopyHeightM > 0.0f) return c.canopyHeightM;
        if (const FuelModel* m = findFuelModel(c.fuelModel)) return m->depthFt / FT_PER_M;
        return 0.0;   // Nonburnable (rock, water, urban, ...) or no data
    }

    // Groups the fire grid into columns, sets their roughness and stretches the node layers between each
    // column's ground and a flat top above the highest summit (WindNinja's Mesh::get_z / compute_domain_height).
    // wind20Ms is used only for the boundary-layer height (0.2 u* / f, Blackadar & Tennekes 1968).
    Mesh buildMesh(const GridSpec& spec, const std::vector<Cell>& grid, double wind20Ms, const WindSolverSettings& s) {
        Mesh m;
        const float cellSizeX = spec.cellSizeX, cellSizeY = spec.cellSizeY;
        const int W = spec.width, H = spec.height;
        const int fx = std::max(1, static_cast<int>(std::lround(s.columnSizeM / cellSizeX)));   // Fire cells per column
        const int fy = std::max(1, static_cast<int>(std::lround(s.columnSizeM / cellSizeY)));
        m.nx = (W + fx - 1) / fx;
        m.ny = (H + fy - 1) / fy;
        m.nz = std::max(3, s.layers);
        m.hx = static_cast<double>(fx) * cellSizeX;
        m.hy = static_cast<double>(fy) * cellSizeY;
        m.cols.assign(static_cast<std::size_t>(m.nx) * m.ny, Column());

        // Column means of ground elevation and vegetation height (the last row/column of blocks may be partial)
        double maxGround = -1e30, maxVeg = 0.0;
        for (int j = 0; j < m.ny; ++j) {
            for (int i = 0; i < m.nx; ++i) {
                double sumZ = 0.0, sumH = 0.0;
                int n = 0;
                for (int y = j * fy; y < std::min((j + 1) * fy, H); ++y)
                    for (int x = i * fx; x < std::min((i + 1) * fx, W); ++x) {
                        const Cell& c = grid[static_cast<std::size_t>(y) * W + x];
                        sumZ += c.elevation;
                        sumH += vegetationHeightM(c);
                        ++n;
                    }
                Column& col = m.cols[static_cast<std::size_t>(j) * m.nx + i];
                col.ground = sumZ / n;
                col.vegH = sumH / n;
                col.z0 = std::max(0.13 * col.vegH, MIN_ROUGHNESS_M);
                col.d = 0.64 * col.vegH;
                col.zLow = col.d + 7.0 * col.z0;
                // Neutral friction velocity from the 20-ft wind, then the neutral boundary-layer height
                const double uStar = VON_KARMAN * wind20Ms / std::log((WIND_REF_HEIGHT_M + col.vegH - col.d) / col.z0);
                col.zTop = std::max(col.d + 0.2 * uStar / s.coriolisPerS, col.zLow);
                maxGround = std::max(maxGround, col.ground);
                maxVeg = std::max(maxVeg, col.vegH);
            }
        }

        // Domain top: topAboveSummitM above the highest column. Never lower than WindNinja's own top (first layer
        // column width / maxAspectRatio thick over the summit, layers growing geometrically) or 3 x (20-ft
        // height + tallest vegetation), so the sampling height always lies well inside the mesh.
        const double r = s.layerGrowth;
        const double firstLayer = 0.5 * (m.hx + m.hy) / s.maxAspectRatio;
        double above = firstLayer * (std::pow(r, m.nz) - 1.0) / (r - 1.0);
        above = std::max({ above, static_cast<double>(s.topAboveSummitM), 3.0 * (WIND_REF_HEIGHT_M + maxVeg) });
        const double top = maxGround + above;

        // Layer k of a column sits at ground + (top - ground) (r^k - 1) / (r^(nz-1) - 1): k = 0 on the ground,
        // k = nz - 1 at the top, layers thin near the ground and thicker upwards (terrain-following)
        m.z.assign(m.count(), 0.0);
        const double denom = std::pow(r, m.nz - 1) - 1.0;
        for (int k = 0; k < m.nz; ++k) {
            const double frac = (std::pow(r, k) - 1.0) / denom;
            for (int j = 0; j < m.ny; ++j)
                for (int i = 0; i < m.nx; ++i) {
                    const double g = m.cols[static_cast<std::size_t>(j) * m.nx + i].ground;
                    m.z[m.node(i, j, k)] = g + (top - g) * frac;
                }
        }
        return m;
    }

    // ------------------------------------------------------------------------
    // The 27-point stencil matrix (step 4)
    // ------------------------------------------------------------------------
    // Row n holds the coefficients coupling node n to node n + offset[s], s = (dk+1)*9 + (dj+1)*3 + (di+1)
    // for di, dj, dk in {-1, 0, 1}. s = 13 is the diagonal; s < 13 are lower-numbered nodes, s > 13 higher.
    // Rows of boundary nodes (Phi fixed at 0) are identity rows, flagged in `fixed` and never read through the
    // stencil: their neighbours can lie off the mesh. Every other node has i, j inside the mesh and k below
    // the top, so all its neighbours exist except the layer below a ground node (k = 0): those s < 9 entries
    // stay 0 and are skipped.
    struct StencilMatrix {
        std::vector<double> a;                // 27 coefficients per node
        std::vector<char> fixed;              // 1 = boundary node: identity row, Phi = 0
        std::ptrdiff_t offset[27] = {};       // Node-index offset of each stencil entry
        int layerSize = 0;                    // nx * ny: node-index step between layers
        int nx = 0, ny = 0, nz = 0;           // Mesh nodes along X, Y and up (used by the row-wise kernels)

        double* row(std::size_t n) { return &a[n * 27]; }
        const double* row(std::size_t n) const { return &a[n * 27]; }
        // First stencil entry to read in a free row: skip the layer below on the ground
        int firstEntry(std::size_t n) const { return n < static_cast<std::size_t>(layerSize) ? 9 : 0; }
        std::size_t neighbour(std::size_t n, int s) const {
            return static_cast<std::size_t>(static_cast<std::ptrdiff_t>(n) + offset[s]);
        }
    };

    // ------------------------------------------------------------------------
    // Reference kernels: node by node, in node order, exactly as first written
    // ------------------------------------------------------------------------
    // Kept so the tests can check the row-wise kernels below against them bit for bit
    // (WindSolverSettings::referenceKernels, which also makes the assembly and smoothing call evalHex at every
    // Gauss point instead of using the shared ElementGeometry). The app always uses the row-wise kernels.

    // y = A x
    void multiplyReference(const StencilMatrix& A, const std::vector<double>& x, std::vector<double>& y) {
        const std::size_t total = x.size();
        for (std::size_t n = 0; n < total; ++n) {
            if (A.fixed[n]) { y[n] = x[n]; continue; }   // Identity row
            const double* row = A.row(n);
            double sum = 0.0;
            for (int s = A.firstEntry(n); s < 27; ++s) sum += row[s] * x[A.neighbour(n, s)];
            y[n] = sum;
        }
    }

    // Symmetric Gauss-Seidel preconditioner (SSOR with omega = 1, as WindNinja's preconditioner.cpp):
    // z = M^-1 r with M = (D + L) D^-1 (D + U). The forward sweep solves (D + L) t = r, the backward sweep
    // solves (D + U) z = D t, both in place in z. Identity rows just copy r.
    void preconditionReference(const StencilMatrix& A, const std::vector<double>& r, std::vector<double>& z) {
        const std::size_t total = r.size();
        for (std::size_t n = 0; n < total; ++n) {   // Forward: lower-numbered neighbours are already final
            if (A.fixed[n]) { z[n] = r[n]; continue; }
            const double* row = A.row(n);
            double sum = r[n];
            for (int s = A.firstEntry(n); s < 13; ++s) sum -= row[s] * z[A.neighbour(n, s)];
            z[n] = sum / row[13];
        }
        for (std::size_t n = total; n-- > 0;) {     // Backward: higher-numbered neighbours are already final
            if (A.fixed[n]) continue;
            const double* row = A.row(n);
            double sum = 0.0;
            for (int s = 14; s < 27; ++s) sum += row[s] * z[A.neighbour(n, s)];
            z[n] -= sum / row[13];
        }
    }

    double dot(const std::vector<double>& a, const std::vector<double>& b) {
        double s = 0.0;
        for (std::size_t i = 0; i < a.size(); ++i) s += a[i] * b[i];
        return s;
    }

    // ------------------------------------------------------------------------
    // Row-wise kernels (what the solver runs): same arithmetic, faster order
    // ------------------------------------------------------------------------
    // Every node's value is computed with exactly the same operations, in the same order, as in the reference
    // kernels above, so the results are identical to the last bit. What changes is only the order in which
    // nodes are visited and some bookkeeping:
    //  - The mesh is walked row by row (fixed j and k, i running), so the boundary tests and the "skip the
    //    layer below the ground" test are made once per row instead of once per node and stencil entry.
    //  - The Gauss-Seidel sweeps work on bands of SWEEP_BAND rows at once in a skewed wavefront (below).

    // Stencil entries of one row of coefficients times the values round node n, in stencil order:
    // sum (starting from `start`) combined with row[s] * v[n + offset[s]] for s = First .. Last - 1.
    // Sub = true subtracts each product (forward sweep), false adds it (multiply, backward sweep).
    // First and Last are template parameters so the compiler unrolls the loop completely.
    template <int First, int Last, bool Sub>
    inline double stencilSum(double start, const double* row, const double* v, const std::ptrdiff_t* offset) {
        double sum = start;
        for (int s = First; s < Last; ++s) {
            if (Sub) sum -= row[s] * v[offset[s]];
            else sum += row[s] * v[offset[s]];
        }
        return sum;
    }

    // y = A x, and returns dot(x, y) = x . A x (summed in node order, as dot() would)
    //
    // Rows of fixed nodes copy x. Free nodes sum all 27 entries (18 on the ground layer, which has no layer
    // below), each from a sum that starts at 0, so y[n] is bit-identical to multiplyReference. The dot product
    // CG needs next (p . A p) is accumulated in the same loop, in increasing node order, which is the order
    // dot() adds in, so it too is bit-identical; this saves one pass over two vectors per iteration.
    double multiplyAndDot(const StencilMatrix& A, const std::vector<double>& x, std::vector<double>& y) {
        const int nx = A.nx, ny = A.ny, nz = A.nz;
        const double* xs = x.data();
        double* ys = y.data();
        double dotXY = 0.0;
        for (int k = 0; k < nz; ++k)
            for (int j = 0; j < ny; ++j) {
                const std::size_t n0 = (static_cast<std::size_t>(k) * ny + j) * nx;   // First node of the row
                if (k == nz - 1 || j == 0 || j == ny - 1) {   // Whole row on the boundary: identity rows
                    for (int i = 0; i < nx; ++i) { ys[n0 + i] = xs[n0 + i]; dotXY += xs[n0 + i] * ys[n0 + i]; }
                    continue;
                }
                // i = 0 and i = nx - 1 are side-boundary nodes (identity); i = 1 .. nx - 2 are free
                ys[n0] = xs[n0];
                dotXY += xs[n0] * ys[n0];
                for (int i = 1; i < nx - 1; ++i) {
                    const std::size_t n = n0 + i;
                    const double* row = A.row(n);
                    ys[n] = k == 0 ? stencilSum<9, 27, false>(0.0, row, xs + n, A.offset)
                                   : stencilSum<0, 27, false>(0.0, row, xs + n, A.offset);
                    dotXY += xs[n] * ys[n];
                }
                const std::size_t last = n0 + nx - 1;
                ys[last] = xs[last];
                dotXY += xs[last] * ys[last];
            }
        return dotXY;
    }

    // Rows swept together by the Gauss-Seidel kernels. Each node's update is a chain of 13 dependent
    // multiply-adds and a division, and in the backward sweep that chain starts with the node just finished,
    // so node after node the processor waits on the previous result. Rows of a band are interleaved so it
    // always has independent chains to work on. On a 2-core Xeon this cut the preconditioner's time by about a
    // quarter; bands of 2, 4 and 8 rows gave much the same gain (within the machine's timing noise).
    constexpr int SWEEP_BAND = 4;

    // Symmetric Gauss-Seidel preconditioner, z = M^-1 r, bit-identical to preconditionReference.
    //
    // Why the order can change without changing a single bit: in the forward sweep node (i, j, k) reads only
    // lower-numbered neighbours, which must already hold their final values: the layer below, row j - 1 at
    // i - 1 .. i + 1, and its own row at i - 1. In node order those are all done first, but they are also all
    // done if row j runs two nodes behind row j - 1. So a band of rows is swept as a skewed wavefront: at step
    // t, row q of the band does node i = 1 + t - 2q. When row q reaches i, row q - 1 finished i + 1 one step
    // earlier. Every node therefore reads exactly the values it reads in node order, and does the same
    // operations on them in the same order, so z comes out bit-identical. The backward sweep is the mirror
    // image (rows from the top of the layer down, i from nx - 2 down, reading only higher-numbered nodes).
    // Fixed nodes (z = r) are set first; nothing else writes them, so every read of them sees r as before.
    void precondition(const StencilMatrix& A, const std::vector<double>& r, std::vector<double>& z) {
        const int nx = A.nx, ny = A.ny, nz = A.nz;
        const double* rs = r.data();
        double* zs = z.data();
        auto node = [&](int i, int j, int k) { return (static_cast<std::size_t>(k) * ny + j) * nx + i; };

        // Fixed nodes: z = r (the top layer, rows j = 0 and ny - 1, columns i = 0 and nx - 1)
        for (int k = 0; k < nz; ++k)
            for (int j = 0; j < ny; ++j) {
                const std::size_t n0 = node(0, j, k);
                if (k == nz - 1 || j == 0 || j == ny - 1) {
                    for (int i = 0; i < nx; ++i) zs[n0 + i] = rs[n0 + i];
                }
                else {
                    zs[n0] = rs[n0];
                    zs[n0 + nx - 1] = rs[n0 + nx - 1];
                }
            }
        const int freeRows = ny - 2;   // Free rows j = 1 .. ny - 2 in each layer below the top
        const int freeCols = nx - 2;   // Free nodes i = 1 .. nx - 2 in each free row
        if (freeRows <= 0 || freeCols <= 0) return;

        // Forward sweep, (D + L) z = r: layers upwards, bands of rows northwards (increasing j)
        for (int k = 0; k < nz - 1; ++k)
            for (int jb = 1; jb <= ny - 2; jb += SWEEP_BAND) {
                const int rows = std::min(SWEEP_BAND, ny - 1 - jb);
                const int steps = freeCols + 2 * (rows - 1);   // Until the last row of the band reaches nx - 2
                for (int t = 0; t < steps; ++t)
                    for (int q = 0; q < rows; ++q) {
                        const int i = 1 + t - 2 * q;
                        if (i < 1 || i > nx - 2) continue;   // This row has not started yet, or has finished
                        const std::size_t n = node(i, jb + q, k);
                        const double* row = A.row(n);
                        const double sum = k == 0 ? stencilSum<9, 13, true>(rs[n], row, zs + n, A.offset)
                                                  : stencilSum<0, 13, true>(rs[n], row, zs + n, A.offset);
                        zs[n] = sum / row[13];
                    }
            }

        // Backward sweep, (D + U) z = D t: layers downwards, bands of rows southwards (decreasing j)
        for (int k = nz - 2; k >= 0; --k)
            for (int jb = ny - 2; jb >= 1; jb -= SWEEP_BAND) {
                const int rows = std::min(SWEEP_BAND, jb);
                const int steps = freeCols + 2 * (rows - 1);
                for (int t = 0; t < steps; ++t)
                    for (int q = 0; q < rows; ++q) {
                        const int i = nx - 2 - t + 2 * q;
                        if (i < 1 || i > nx - 2) continue;
                        const std::size_t n = node(i, jb - q, k);
                        const double* row = A.row(n);
                        const double sum = stencilSum<14, 27, false>(0.0, row, zs + n, A.offset);
                        zs[n] -= sum / row[13];
                    }
            }
    }

    // RMS of div(u, v, w) over the mesh, each element evaluated at its centre and weighted by its volume.
    // The two lowest layers of elements are left out. The smoothed wind at a ground node is a "slip" velocity
    // (the first guess is 0 there but the gradient of Phi is not), and the smoothing at the first nodes above
    // the ground mixes in gradients from the ground layer, so div(u) in those two thin layers measures the
    // step-5 smoothing, not the solve. The solve itself conserves mass in the weak (finite element) sense
    // everywhere; from the third layer up the smoothed wind shows it too (div(u) ~20 times smaller than the
    // first guess's on the test hill).
    double rmsDivergence(const Mesh& m, const std::vector<double>& u, const std::vector<double>& v,
                         const std::vector<double>& w) {
        double sum = 0.0, vol = 0.0;
        double xyz[8][3];
        std::size_t ids[8];
        HexPoint hp;
        for (int k = 2; k < m.nz - 1; ++k)
            for (int j = 0; j < m.ny - 1; ++j)
                for (int i = 0; i < m.nx - 1; ++i) {
                    m.elementCorners(i, j, k, xyz, ids);
                    evalHex(xyz, 0.0, 0.0, 0.0, hp);
                    double div = 0.0;
                    for (int a = 0; a < 8; ++a)
                        div += hp.dN[a][0] * u[ids[a]] + hp.dN[a][1] * v[ids[a]] + hp.dN[a][2] * w[ids[a]];
                    const double elemVol = 8.0 * hp.detJ;   // Natural volume of the cube is 2 x 2 x 2
                    sum += div * div * elemVol;
                    vol += elemVol;
                }
        return vol > 0.0 ? std::sqrt(sum / vol) : 0.0;
    }

    // The open wind everywhere: used when there is nothing to adjust (no wind, or a map too small for a mesh)
    WindField uniformField(const GridSpec& spec, Vec2 toward, double wind20Kmh) {
        WindField f;
        f.grid = spec;
        f.toward.assign(spec.count(), toward);
        f.speedKmh.assign(spec.count(), static_cast<float>(wind20Kmh));
        f.converged = true;
        return f;
    }
}

// ----------------------------------------------------------------------------
// Solve the terrain wind (steps 1-6 of wind_solver.h)
// ----------------------------------------------------------------------------
WindField solveTerrainWind(const GridSpec& spec, const std::vector<Cell>& grid, float wind20Kmh, float windFromDeg,
                           const WindSolverSettings& s) {
    requireGridSize(grid, spec, "solveTerrainWind cells");
    const float cellSizeX = spec.cellSizeX, cellSizeY = spec.cellSizeY;   // Ground size of a fire cell, m
    const Vec2 toward = windTowardGrid(windFromDeg);
    const double wind20Ms = wind20Kmh / KMH_PER_MS;
    if (wind20Ms <= 0.0) return uniformField(spec, toward, 0.0);

    // --- Step 1: mesh ---
    const Mesh m = buildMesh(spec, grid, wind20Ms, s);
    if (m.nx < 3 || m.ny < 3) return uniformField(spec, toward, wind20Kmh);   // No interior nodes to solve for
    const std::size_t total = m.count();

    // --- Step 2: first guess, the open wind spread up each column with the log profile ---
    std::vector<double> u0(total, 0.0), v0(total, 0.0), w0(total, 0.0);
    for (int j = 0; j < m.ny; ++j)
        for (int i = 0; i < m.nx; ++i) {
            const Column& col = m.cols[static_cast<std::size_t>(j) * m.nx + i];
            const double sRef = profileShape(col, WIND_REF_HEIGHT_M + col.vegH);
            for (int k = 1; k < m.nz; ++k) {   // k = 0 is the ground: no wind there
                const std::size_t n = m.node(i, j, k);
                const double speed = wind20Ms * profileShape(col, m.z[n] - col.ground) / sRef;
                u0[n] = speed * toward.x;
                v0[n] = speed * toward.y;
            }
        }

    // --- Step 4: assemble K Phi = Q (step 3's equation in weak form) ---
    // K_ab = integral of Rx dNa/dX dNb/dX + Ry dNa/dY dNb/dY + Rz dNa/dZ dNb/dZ;  Q_a = integral of N_a div(u0)
    const double Rx = 0.5, Ry = 0.5, Rz = 0.5 * s.alpha * s.alpha;
    StencilMatrix A;
    A.a.assign(total * 27, 0.0);
    A.layerSize = m.nx * m.ny;
    A.nx = m.nx;
    A.ny = m.ny;
    A.nz = m.nz;
    for (int dk = -1; dk <= 1; ++dk)
        for (int dj = -1; dj <= 1; ++dj)
            for (int di = -1; di <= 1; ++di)
                A.offset[(dk + 1) * 9 + (dj + 1) * 3 + (di + 1)] =
                    static_cast<std::ptrdiff_t>(dk) * A.layerSize + static_cast<std::ptrdiff_t>(dj) * m.nx + di;
    std::vector<double> rhs(total, 0.0);
    const ElementGeometry geom = elementGeometry(m);
    {
        double xyz[8][3];
        std::size_t ids[8];
        HexPoint hp;
        for (int k = 0; k < m.nz - 1; ++k)
            for (int j = 0; j < m.ny - 1; ++j)
                for (int i = 0; i < m.nx - 1; ++i) {
                    m.elementCorners(i, j, k, xyz, ids);
                    double S[8][8] = {};
                    double Q[8] = {};
                    for (int g = 0; g < 8; ++g) {
                        if (s.referenceKernels)
                            evalHex(xyz, (g & 1) ? GAUSS : -GAUSS, ((g >> 1) & 1) ? GAUSS : -GAUSS,
                                    ((g >> 2) & 1) ? GAUSS : -GAUSS, hp);
                        else
                            gaussHex(geom, i, j, g, xyz, hp);   // Bit-identical to the evalHex call above
                        double div0 = 0.0;   // div(u0) at this Gauss point
                        for (int b = 0; b < 8; ++b)
                            div0 += hp.dN[b][0] * u0[ids[b]] + hp.dN[b][1] * v0[ids[b]] + hp.dN[b][2] * w0[ids[b]];
                        for (int a = 0; a < 8; ++a) {
                            Q[a] += hp.N[a] * div0 * hp.detJ;
                            for (int b = 0; b < 8; ++b)
                                S[a][b] += (Rx * hp.dN[a][0] * hp.dN[b][0] + Ry * hp.dN[a][1] * hp.dN[b][1]
                                          + Rz * hp.dN[a][2] * hp.dN[b][2]) * hp.detJ;
                        }
                    }
                    // Scatter into the stencil rows: corner b is at offset (bx - ax, by - ay, bz - az) from a
                    for (int a = 0; a < 8; ++a) {
                        rhs[ids[a]] += Q[a];
                        double* row = A.row(ids[a]);
                        for (int b = 0; b < 8; ++b) {
                            const int di = (b & 1) - (a & 1), dj = ((b >> 1) & 1) - ((a >> 1) & 1),
                                      dk = ((b >> 2) & 1) - ((a >> 2) & 1);
                            row[(dk + 1) * 9 + (dj + 1) * 3 + (di + 1)] += S[a][b];
                        }
                    }
                }
    }

    // Boundary conditions: Phi = 0 on the sides and top. Their rows become identity rows with a 0 right-hand
    // side, and every other row drops its coupling to them (the known Phi = 0 moves to the right-hand side as
    // 0). This keeps the matrix symmetric. The ground needs nothing: the natural condition is no flux.
    A.fixed.assign(total, 0);
    for (int k = 0; k < m.nz; ++k)
        for (int j = 0; j < m.ny; ++j)
            for (int i = 0; i < m.nx; ++i)
                A.fixed[m.node(i, j, k)] = m.boundary(i, j, k) ? 1 : 0;
    for (std::size_t n = 0; n < total; ++n) {
        double* row = A.row(n);
        if (A.fixed[n]) {
            std::fill(row, row + 27, 0.0);
            row[13] = 1.0;
            rhs[n] = 0.0;
            continue;
        }
        for (int s2 = A.firstEntry(n); s2 < 27; ++s2)
            if (s2 != 13 && A.fixed[A.neighbour(n, s2)]) row[s2] = 0.0;
    }

    // --- Step 4: preconditioned conjugate gradients ---
    WindField field;
    field.grid = spec;
    field.columnsX = m.nx;
    field.columnsY = m.ny;
    field.layers = m.nz;
    std::vector<double> phi(total, 0.0);
    // Nothing to correct when the first guess is already mass-conserving (flat ground with even roughness):
    // the right-hand side is then only rounding noise, far below the size of a real term (wind x column area)
    const double rhsNorm = std::sqrt(dot(rhs, rhs));
    const double rhsScale = wind20Ms * m.hx * m.hy * std::sqrt(static_cast<double>(total));
    if (rhsNorm > 1e-10 * rhsScale) {
        std::vector<double> r = rhs, z(total, 0.0), p(total, 0.0), Ap(total, 0.0);
        // The row-wise kernels give bit-identical results to the reference ones (tests check this)
        auto applyM = [&](const std::vector<double>& in, std::vector<double>& out) {
            if (s.referenceKernels) preconditionReference(A, in, out);
            else precondition(A, in, out);
        };
        applyM(r, z);
        p = z;
        double rz = dot(r, z);
        int it = 0;
        double rel = 1.0;
        while (it < s.maxIterations) {
            double pAp;
            if (s.referenceKernels) { multiplyReference(A, p, Ap); pAp = dot(p, Ap); }
            else pAp = multiplyAndDot(A, p, Ap);
            const double alpha = rz / pAp;
            // Step along p, and ||r||^2 in the same pass (same additions, in node order, as dot(r, r))
            double rr = 0.0;
            if (s.referenceKernels) {
                for (std::size_t n = 0; n < total; ++n) {
                    phi[n] += alpha * p[n];
                    r[n] -= alpha * Ap[n];
                }
                rr = dot(r, r);
            }
            else {
                for (std::size_t n = 0; n < total; ++n) {
                    phi[n] += alpha * p[n];
                    r[n] -= alpha * Ap[n];
                    rr += r[n] * r[n];
                }
            }
            ++it;
            rel = std::sqrt(rr) / rhsNorm;
            if (rel < s.tolerance) break;
            applyM(r, z);
            const double rzNew = dot(r, z);
            const double beta = rzNew / rz;
            rz = rzNew;
            for (std::size_t n = 0; n < total; ++n) p[n] = z[n] + beta * p[n];
        }
        field.iterations = it;
        field.relResidual = rel;
        field.converged = rel < s.tolerance;
    }
    else {
        field.converged = true;   // The first guess is already mass-conserving (e.g. flat ground)
    }
    A.a.clear();
    A.a.shrink_to_fit();   // Free the matrix before the next step allocates

    // --- Step 5: wind at the nodes, u = u0 + Rx dPhi/dx etc., gradients smoothed onto the nodes ---
    // Each Gauss point's gradient goes to the element's 8 corners, weighted by 1 / distance
    std::vector<double> gx(total, 0.0), gy(total, 0.0), gz(total, 0.0), wsum(total, 0.0);
    {
        double xyz[8][3];
        std::size_t ids[8];
        HexPoint hp;
        for (int k = 0; k < m.nz - 1; ++k)
            for (int j = 0; j < m.ny - 1; ++j)
                for (int i = 0; i < m.nx - 1; ++i) {
                    m.elementCorners(i, j, k, xyz, ids);
                    for (int g = 0; g < 8; ++g) {
                        if (s.referenceKernels)
                            evalHex(xyz, (g & 1) ? GAUSS : -GAUSS, ((g >> 1) & 1) ? GAUSS : -GAUSS,
                                    ((g >> 2) & 1) ? GAUSS : -GAUSS, hp);
                        else
                            gaussHex(geom, i, j, g, xyz, hp);   // Bit-identical to the evalHex call above
                        double grad[3] = { 0.0, 0.0, 0.0 };
                        for (int a = 0; a < 8; ++a)
                            for (int c = 0; c < 3; ++c) grad[c] += hp.dN[a][c] * phi[ids[a]];
                        for (int a = 0; a < 8; ++a) {
                            const double dx = xyz[a][0] - hp.pos[0], dy = xyz[a][1] - hp.pos[1], dz = xyz[a][2] - hp.pos[2];
                            const double wgt = 1.0 / std::sqrt(dx * dx + dy * dy + dz * dz);
                            gx[ids[a]] += wgt * grad[0];
                            gy[ids[a]] += wgt * grad[1];
                            gz[ids[a]] += wgt * grad[2];
                            wsum[ids[a]] += wgt;
                        }
                    }
                }
    }
    std::vector<double> u(total), v(total), w(total);
    for (std::size_t n = 0; n < total; ++n) {
        u[n] = u0[n] + Rx * gx[n] / wsum[n];
        v[n] = v0[n] + Ry * gy[n] / wsum[n];
        w[n] = w0[n] + Rz * gz[n] / wsum[n];
    }
    field.divergenceBefore = rmsDivergence(m, u0, v0, w0);
    field.divergenceAfter = rmsDivergence(m, u, v, w);

    // --- Step 6a: horizontal wind 6.1 m above the vegetation in each column (WindNinja's interp_uvw) ---
    std::vector<double> colU(static_cast<std::size_t>(m.nx) * m.ny), colV(colU.size());
    for (int j = 0; j < m.ny; ++j)
        for (int i = 0; i < m.nx; ++i) {
            const std::size_t c = static_cast<std::size_t>(j) * m.nx + i;
            const Column& col = m.cols[c];
            const double hRef = WIND_REF_HEIGHT_M + col.vegH;   // Sampling height above the ground
            int k = 1;
            while (k < m.nz - 1 && m.z[m.node(i, j, k)] - col.ground < hRef) ++k;
            const std::size_t nk = m.node(i, j, k);
            if (k == 1) {
                // Inside the lowest mesh cell: scale the first node's wind down the log profile
                const double ratio = profileShape(col, hRef) / profileShape(col, m.z[nk] - col.ground);
                colU[c] = u[nk] * ratio;
                colV[c] = v[nk] * ratio;
            }
            else {
                // Between nodes k - 1 and k: interpolate wind / S(height) linearly, then multiply by S at the
                // sampling height. Plain linear interpolation of a log profile reads ~1 % low between nodes;
                // this way an undisturbed log profile (flat ground) comes back exactly.
                const std::size_t nb = m.node(i, j, k - 1);
                const double h1 = m.z[nb] - col.ground, h2 = m.z[nk] - col.ground;
                const double t = std::clamp((hRef - h1) / (h2 - h1), 0.0, 1.0);
                const double s1 = profileShape(col, h1), s2 = profileShape(col, h2), sRef = profileShape(col, hRef);
                colU[c] = sRef * ((1.0 - t) * u[nb] / s1 + t * u[nk] / s2);
                colV[c] = sRef * ((1.0 - t) * v[nb] / s1 + t * v[nk] / s2);
            }
        }

    // --- Step 6b: bilinear interpolation from the column centres to every fire cell ---
    field.toward.resize(spec.count());
    field.speedKmh.resize(field.toward.size());
    for (int y = 0; y < spec.height; ++y) {
        const double ty = std::clamp((y + 0.5) * cellSizeY / m.hy - 0.5, 0.0, m.ny - 1.0);
        const int j0 = std::min(static_cast<int>(ty), m.ny - 2);
        const double fy = ty - j0;
        for (int x = 0; x < spec.width; ++x) {
            const double tx = std::clamp((x + 0.5) * cellSizeX / m.hx - 0.5, 0.0, m.nx - 1.0);
            const int i0 = std::min(static_cast<int>(tx), m.nx - 2);
            const double fx = tx - i0;
            const std::size_t c00 = static_cast<std::size_t>(j0) * m.nx + i0, c10 = c00 + 1;
            const std::size_t c01 = c00 + m.nx, c11 = c01 + 1;
            const double cu = (1 - fy) * ((1 - fx) * colU[c00] + fx * colU[c10]) + fy * ((1 - fx) * colU[c01] + fx * colU[c11]);
            const double cv = (1 - fy) * ((1 - fx) * colV[c00] + fx * colV[c10]) + fy * ((1 - fx) * colV[c01] + fx * colV[c11]);
            const double speed = std::sqrt(cu * cu + cv * cv);
            const std::size_t idx = static_cast<std::size_t>(y) * spec.width + x;
            field.speedKmh[idx] = static_cast<float>(speed * KMH_PER_MS);
            // Direction of the solved wind; keep the open wind's direction where the wind has died away
            field.toward[idx] = speed > 1e-6 ? Vec2(static_cast<float>(cu / speed), static_cast<float>(cv / speed)) : toward;
        }
    }
    return field;
}
