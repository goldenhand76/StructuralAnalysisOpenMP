#include "DelaunayTriangulationOpenMP.h"
#include <omp.h>  // Include OpenMP header

// Point class methods (same as before)
Point::Point(float x, float y) : x(x), y(y) {}

std::vector<float> Point::pos() const {
    return { x, y };
}

bool Point::isEqual(const Point& other) const {
    return (x == other.x && y == other.y);
}

std::string Point::pointToStr() const {
    return "(" + std::to_string(x) + ", " + std::to_string(y) + ")";
}

// Edge class methods (same as before)
Edge::Edge(const Point& a, const Point& b) : a(a), b(b) {}

bool Edge::isEqual(const Edge& other) const {
    return (a.isEqual(other.a) && b.isEqual(other.b)) || (a.isEqual(other.b) && b.isEqual(other.a));
}

float Edge::length() const {
    return std::sqrt(std::pow(b.x - a.x, 2) + std::pow(b.y - a.y, 2));
}

bool Edge::edgeIntersection(const Edge& other) const {
    float x1 = a.x, y1 = a.y, x2 = b.x, y2 = b.y;
    float x3 = other.a.x, y3 = other.a.y, x4 = other.b.x, y4 = other.b.y;

    float t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) /
        ((x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4));
    float u = ((x2 - x1) * (y1 - y3) - (y2 - y1) * (x1 - x3)) /
        ((x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4));

    return (t >= 0 && t <= 1 && u >= 0 && u <= 1);
}

// Triangle class methods (same as before)
Triangle::Triangle(const Point& a, const Point& b, const Point& c) : a(a), b(b), c(c) {}

bool Triangle::isEqual(const Triangle& other) const {
    return (a.isEqual(other.a) || a.isEqual(other.b) || a.isEqual(other.c)) &&
        (b.isEqual(other.a) || b.isEqual(other.b) || b.isEqual(other.c)) &&
        (c.isEqual(other.a) || c.isEqual(other.b) || c.isEqual(other.c));
}

// Graph class methods with OpenMP added to `addPoint`
std::vector<float> Graph::circumcircle(const Triangle& tri) const {
    float D = (tri.a.x - tri.c.x) * (tri.b.y - tri.c.y) - (tri.b.x - tri.c.x) * (tri.a.y - tri.c.y);

    float centerX = (((tri.a.x - tri.c.x) * (tri.a.x + tri.c.x) + (tri.a.y - tri.c.y) * (tri.a.y + tri.c.y)) / 2 * (tri.b.y - tri.c.y) -
        ((tri.b.x - tri.c.x) * (tri.b.x + tri.c.x) + (tri.b.y - tri.c.y) * (tri.b.y + tri.c.y)) / 2 * (tri.a.y - tri.c.y)) / D;

    float centerY = (((tri.b.x - tri.c.x) * (tri.b.x + tri.c.x) + (tri.b.y - tri.c.y) * (tri.b.y + tri.c.y)) / 2 * (tri.a.x - tri.c.x) -
        ((tri.a.x - tri.c.x) * (tri.a.x + tri.c.x) + (tri.a.y - tri.c.y) * (tri.a.y + tri.c.y)) / 2 * (tri.b.x - tri.c.x)) / D;

    float radius = std::sqrt(std::pow(tri.c.x - centerX, 2) + std::pow(tri.c.y - centerY, 2));

    return { centerX, centerY, radius };
}

bool Graph::pointInCircle(const Point& point, const std::vector<float>& circle) const {
    float d = std::sqrt(std::pow(point.x - circle[0], 2) + std::pow(point.y - circle[1], 2));
    return (d < circle[2]);
}

void Graph::addPoint(const Point& point) {
    // Check if point is unique before adding
    for (const auto& p : points) {
        if (p.isEqual(point)) {
            return;  // Equivalent point already exists
        }
    }
    points.push_back(point);

    // Clear previous triangles and edges
    triangles.clear();
    edges.clear();

    size_t numPoints = points.size();
    std::vector<Triangle> tempTriangles;

    // Parallelize triangle creation using OpenMP
#pragma omp parallel for collapse(3) shared(tempTriangles)
    for (int i = 0; i < numPoints; ++i) {
        for (int j = i + 1; j < numPoints; ++j) {
            for (size_t k = j + 1; k < numPoints; ++k) {
                if (!areCollinear(points[i], points[j], points[k])) {
                    Triangle newTri(points[i], points[j], points[k]);
                    if (triangleIsDelaunay(newTri)) {
#pragma omp critical
                        tempTriangles.push_back(newTri);  // Avoid race condition
                    }
                }
            }
        }
    }

    triangles = tempTriangles;  // Assign the calculated triangles

    // Update edges for valid triangles
    for (const auto& tri : triangles) {
        edges.push_back(Edge(tri.a, tri.b));
        edges.push_back(Edge(tri.b, tri.c));
        edges.push_back(Edge(tri.c, tri.a));
    }
}

// Remaining methods are the same...


bool Graph::triangleIsDelaunay(const Triangle& triangle) const {
    std::vector<float> cc = circumcircle(triangle);
    for (const auto& p : points) {
        if (!p.isEqual(triangle.a) && !p.isEqual(triangle.b) && !p.isEqual(triangle.c)) {
            if (pointInCircle(p, cc)) {
                return false;
            }
        }
    }
    return true;
}

void Graph::printTriangles() const {
    for (const auto& tri : triangles) {
        std::cout << "[" << tri.a.pointToStr() << ", " << tri.b.pointToStr() << ", " << tri.c.pointToStr() << "]," << std::endl;
    }
}

bool Graph::areCollinear(const Point& a, const Point& b, const Point& c) {
    return (b.y - a.y) * (c.x - b.x) == (c.y - b.y) * (b.x - a.x);
}

int Graph::getNodeIndex(const Point& point) const {
    for (size_t i = 0; i < points.size(); ++i) {
        if (points[i].isEqual(point)) {
            return i;  // Return the index if the point matches
        }
    }
    return -1;  // Return -1 if the point is not found (should not happen)
}

// Divide points into domains and process each in parallel
void Graph::addPointsWithDomainDecomposition(const std::vector<Point>& points, int numDomains) {
    size_t numPoints = points.size();
    this->points = points;

    // Partition points into domains
    std::vector<std::vector<Point>> domains(numDomains);
    partitionPoints(points, domains, numDomains);

    // Store temporary triangles from each domain
    std::vector<std::vector<Triangle>> domainTriangles(numDomains);

    // Perform triangulation in each domain in parallel
#pragma omp parallel for shared(domainTriangles)
    for (int d = 0; d < numDomains; ++d) {
        std::vector<Triangle> tempTriangles;

        size_t numDomainPoints = domains[d].size();
        for (int i = 0; i < numDomainPoints; ++i) {
            for (int j = i + 1; j < numDomainPoints; ++j) {
                for (int k = j + 1; k < numDomainPoints; ++k) {
                    if (!areCollinear(domains[d][i], domains[d][j], domains[d][k])) {
                        Triangle newTri(domains[d][i], domains[d][j], domains[d][k]);
                        if (triangleIsDelaunay(newTri)) {
                            tempTriangles.push_back(newTri);
                        }
                    }
                }
            }
        }
        domainTriangles[d] = tempTriangles;
    }

    // Merge all triangles from each domain
    for (const auto& triSet : domainTriangles) {
        triangles.insert(triangles.end(), triSet.begin(), triSet.end());
    }

    // Update edges based on the merged triangles
    updateEdges();
}

// Partition points into subdomains
void Graph::partitionPoints(const std::vector<Point>& points, std::vector<std::vector<Point>>& domains, int numDomains) {
    int domainSize = points.size() / numDomains;
    for (size_t i = 0; i < points.size(); ++i) {
        int domainIdx = i / domainSize;
        domainIdx = std::min(domainIdx, numDomains - 1);  // Ensure last domain gets any extra points
        domains[domainIdx].push_back(points[i]);
    }
}

void Graph::updateEdges() {
    for (const auto& tri : triangles) {
        edges.push_back(Edge(tri.a, tri.b));
        edges.push_back(Edge(tri.b, tri.c));
        edges.push_back(Edge(tri.c, tri.a));
    }
}