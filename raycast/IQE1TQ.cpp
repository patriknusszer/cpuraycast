#include "framework.h"
#include <vector>

class Intersectable;

class Hit {
private:
    vec3 intersection;
    vec3 radiance;
    vec3 from;
    bool success;
    vec3 normal;
    Intersectable* i;
public:
    Hit(Intersectable* i, const vec3& from, const vec3&at, const vec3& normal) {
        intersection = at;
        this->from = from;
        this->normal = normal;
        success = true;
        this->i = i;
    }
    
    Hit() {
        success = false;
    }
    
    Intersectable* getObject() {
        return i;
    }
    
    vec3 getIntersection() const {
        return intersection;
    }
    
    vec3 getNormal() const {
        return normal;
    }
    
    vec3 getRadiance() const {
        return radiance;
    }
    
    double getDistance() const {
        return length(from - intersection);
    }
    
    bool getSuccess() const {
        return success;
    }
};

class Ray {
private:
    vec3 start, dir;
public:
    Ray(vec3 _start, vec3 _dir) { start = _start; dir = normalize(_dir); }
    
    vec3 getStart() const { return start; }
    vec3 getDirection() const { return  dir; }
};

class Intersectable {
protected:
    vec3 colour;
public:
    virtual bool transparentFromOutside() { return false; }
    virtual Hit intersect(const Ray& ray) = 0;
    virtual vec3 insideColour() {
        return vec3(0.2, 0.2, 0.2);
    }
    virtual vec3 outsideColour() {
        return vec3(0.2, 0.2, 0.2);
    }
};

inline bool operator<(const Hit& one, const Hit& other) {
    if (!one.getSuccess() && !other.getSuccess()) return false;
    else if (!one.getSuccess() && other.getSuccess())
        return false;
    else if (one.getSuccess() && other.getSuccess() && one.getDistance() > other.getDistance())
        return false;
    else return true;
}

namespace tools {
    const double delta = 0.0001;

    bool close(double a, double b) {
        return fabs(a - b) < delta;
    }
}

class Triangle : Intersectable {
private:
    std::vector<vec3> vertices;
    double A, B, C, D;
    double AREA;
    
    void setSurfaceNormal(const vec3& origin) {
        double discriminator = A * A + B * B + C * C;
        double lambda = (-1 * D - A * origin.x - B * origin.y - C * origin.z) / discriminator;
        vec3 intersection = origin + lambda * vec3(A, B, C);
        vec3 snorm = intersection - origin;
        A = snorm.x;
        B = snorm.y;
        C = snorm.z;
        vec3 normalised = normalize(vec3(A, B, C));
        A = normalised.x;
        B = normalised.y;
        C = normalised.z;
        D = -1 * A * vertices[0].x - B * vertices[0].y - C * vertices[0].z;
    }
    
    void init(const std::vector<vec3>& vertices, const vec3& origin) {
        this->vertices = vertices;
        vec3 side1 = vertices[0] - vertices[1];
        vec3 side2 = vertices[0] - vertices[2];
        vec3 normal = cross(side1, side2);
        A = normal.x;
        B = normal.y;
        C = normal.z;
        D = -1 * A * vertices[0].x - B * vertices[0].y - C * vertices[0].z;
        setSurfaceNormal(origin);
        double sideA = length(vertices[0] - vertices[2]);
        double sideB = length(vertices[0] - vertices[1]);
        double sideC = length(vertices[2] - vertices[1]);
        double S = (sideA + sideB + sideC) / 2;
        AREA = sqrt(S * (S - sideA) * (S - sideB) * (S - sideC));
    }
public:
    Triangle(const std::vector<vec3>& vertices, const vec3& origin) {
        init(vertices, origin);
    }
    
    Triangle(vec3 first, vec3 second, vec3 third, const vec3& origin) {
        std::vector<vec3> verts;
        verts.push_back(first);
        verts.push_back(second);
        verts.push_back(third);
        init(verts, origin);
    }
    
    Hit intersect(const Ray& ray) {
        double scalar = dot(ray.getDirection(), vec3(A,B,C));
        if (scalar == 0) return Hit();
        double t = dot(vertices[0] - ray.getStart(), vec3(A,B,C)) / scalar;
        if (t <= 0) return Hit();
        vec3 p = ray.getStart() + t * ray.getDirection();
        if (!tools::close(dot(p - vertices[0], vec3(A,B,C)), 0)) return Hit();
        double s1 = dot(cross(vertices[1] - vertices[0], p - vertices[0]), vec3(A,B,C));
        double s2 = dot(cross(vertices[2] - vertices[1], p - vertices[1]), vec3(A,B,C));
        double s3 = dot(cross(vertices[0] - vertices[2], p - vertices[2]), vec3(A,B,C));
        if (s1 > 0 == s2 > 0 && s2 > 0 == s3 > 0) return Hit(this, ray.getStart(), p, vec3(A, B, C));
        else return Hit();
        
        /*double discriminator = A * ray.getDirection().x + B * ray.getDirection().y + C * ray.getDirection().z;
        
        // the ray is parallel to the triangle: either no solution, or infinitely many solutions (when line is inside plane)
        if (tools::close(discriminator, 0))
            return Hit();
        
        double lambda = (-1 * D - A * ray.getStart().x - B * ray.getStart().y - C * ray.getStart().z) / discriminator;
        
        if (lambda <= 0) return Hit();
        
        vec3 intersection = ray.getStart() + ray.getDirection() * lambda;
        
        // the ray intersects the plane of the triangle, when direction is elongated by lambda
        // lets check if intersection is within bounds
        double sideA = length(vertices[0] - vertices[1]);
        double sideB = length(intersection - vertices[0]);
        double sideC = length(intersection - vertices[1]);
        double S = (sideA + sideB + sideC) / 2;
        double area1 = sqrt(S * (S - sideA) * (S - sideB) * (S - sideC));
        sideA = length(vertices[1] - vertices[2]);
        sideB = length(intersection - vertices[2]);
        // double sideC = length(intersection - vertices[1]);
        S = (sideA + sideB + sideC) / 2;
        double area2 = sqrt(S * (S - sideA) * (S - sideB) * (S - sideC));
        sideA = length(vertices[0] - vertices[2]);
        // sideB = length(intersection - vertices[2]);
        sideC = length(intersection - vertices[0]);
        S = (sideA + sideB + sideC) / 2;
        double area3 = sqrt(S * (S - sideA) * (S - sideB) * (S - sideC));
        double area = area1 + area2 + area3;
        if (tools::close(AREA, area))
            return Hit(this, ray.getStart(), intersection, vec3(A, B, C));
        else
            return Hit();*/
    }
};

class Wall;

Wall* flr = 0;

class ConvexQuadrilateral : public Intersectable {
private:
    std::vector<vec3> vertices;
    double A, B, C, D; // parameters of plane
    double AREA; // area of face
    vec3 origin;
    
    void setSurfaceNormal(const vec3& origin) {
        this->origin = origin;
        double discriminator = A * A + B * B + C * C;
        double lambda = (-1 * D - A * origin.x - B * origin.y - C * origin.z) / discriminator;
        vec3 intersection = origin + lambda * vec3(A, B, C);
        vec3 snorm = intersection - origin;
        A = snorm.x;
        B = snorm.y;
        C = snorm.z;
        vec3 normalised = normalize(vec3(A, B, C));
        A = normalised.x;
        B = normalised.y;
        C = normalised.z;
        D = -1 * A * vertices[0].x - B * vertices[0].y - C * vertices[0].z;
    }
protected:
    void init(const std::vector<vec3>& vertices, const vec3& origin) {
        this->vertices = vertices;
        vec3 side1 = vertices[0] - vertices[1];
        vec3 side2 = vertices[0] - vertices[3];
        vec3 normal = cross(side1, side2);
        A = normal.x;
        B = normal.y;
        C = normal.z;
        D = -1 * A * vertices[0].x - B * vertices[0].y - C * vertices[0].z;
        setSurfaceNormal(origin);
        double sideA = length(vertices[0] - vertices[2]);
        double sideB = length(vertices[0] - vertices[1]);
        double sideC = length(vertices[2] - vertices[1]);
        double S = (sideA + sideB + sideC) / 2;
        double area1 = sqrt(S * (S - sideA) * (S - sideB) * (S - sideC));
        //  double sideA = length(vertices[0] - vertices[2]);
        sideB = length(vertices[0] - vertices[3]);
        sideC = length(vertices[2] - vertices[3]);
        S = (sideA + sideB + sideC) / 2;
        double area2 = sqrt(S * (S - sideA) * (S - sideB) * (S - sideC));
        AREA = area1 + area2;
    }
    
public:
    ConvexQuadrilateral(const std::vector<vec3>& vertices, const vec3& origin) {
        init(vertices, origin);
    }
    
    ConvexQuadrilateral(vec3 first, vec3 second, vec3 third, vec3 fourth, const vec3& origin) {
        std::vector<vec3> verts;
        verts.push_back(first);
        verts.push_back(second);
        verts.push_back(third);
        verts.push_back(fourth);
        init(verts, origin);
    }
    
    Hit intersect(const Ray& ray) {
        double scalar = dot(ray.getDirection(), vec3(A,B,C));
        if (scalar == 0) return Hit();
        double t = dot(vertices[0] - ray.getStart(), vec3(A,B,C)) / scalar;
        if (t <= 0) return Hit();
        vec3 p = ray.getStart() + t * ray.getDirection();
        if (!tools::close(dot(p - vertices[0], vec3(A,B,C)), 0)) return Hit();
        double s1 = dot(cross(vertices[1] - vertices[0], p - vertices[0]), vec3(A,B,C));
        double s2 = dot(cross(vertices[2] - vertices[1], p - vertices[1]), vec3(A,B,C));
        double s3 = dot(cross(vertices[3] - vertices[2], p - vertices[2]), vec3(A,B,C));
        double s4 = dot(cross(vertices[0] - vertices[3], p - vertices[3]), vec3(A,B,C));
        if (s1 < 0 == s2 < 0 && s2 < 0 == s3 < 0 && s3 < 0 == s4 < 0)
            return Hit(this, ray.getStart(), p, vec3(A, B, C));
        else return Hit();
        
        //printf("%lf:%lf:%lf %lf:%lf:%lf\n", ray.getStart().x, ray.getStart().y, ray.getStart().z,
               //ray.getDirection().x, ray.getDirection().y, ray.getDirection().z);
        /*double discriminator = A * ray.getDirection().x + B * ray.getDirection().y + C * ray.getDirection().z;
        
        // the ray is parallel to the face: either no solution, or infinitely many solutions (when line is inside plane)
        if (tools::close(discriminator, 0))
            return Hit();
        
        double lambda = (-1 * D - A * ray.getStart().x - B * ray.getStart().y - C * ray.getStart().z) / discriminator;
        if (lambda <= 0) return Hit();
        vec3 intersection = ray.getStart() + ray.getDirection() * lambda;
        // the ray intersects the plane of the face, when direction is elongated by lambda
        // lets check if intersection is within bounds
        double sideA = length(vertices[0] - vertices[1]);
        double sideB = length(intersection - vertices[0]);
        double sideC = length(intersection - vertices[1]);
        double S = (sideA + sideB + sideC) / 2;
        double area1 = sqrt(S * (S - sideA) * (S - sideB) * (S - sideC));
        sideA = length(vertices[0] - vertices[3]);
        // sideB = length(intersection - vertices[0]);
        sideC = length(intersection - vertices[3]);
        S = (sideA + sideB + sideC) / 2;
        double area2 = sqrt(S * (S - sideA) * (S - sideB) * (S - sideC));
        sideA = length(vertices[2] - vertices[3]);
        // sideC = length(intersection - vertices[3]);
        sideB = length(intersection - vertices[2]);
        S = (sideA + sideB + sideC) / 2;
        double area3 = sqrt(S * (S - sideA) * (S - sideB) * (S - sideC));
        sideA = length(vertices[2] - vertices[1]);
        sideC = length(intersection - vertices[1]);
        // sideB = length(intersection - vertices[2]);
        S = (sideA + sideB + sideC) / 2;
        double area4 = sqrt(S * (S - sideA) * (S - sideB) * (S - sideC));
        double area = area1 + area2 + area3 + area4;
        
        if (tools::close(AREA, area)) {
            //printf("i: (%lf, %lf, %lf)\n", intersection.x, intersection.y, intersection.z);
            return Hit(this, ray.getStart(), intersection, vec3(A, B, C));
        }
        else
            return Hit();*/
    }
};

class Tetrahedron : public Intersectable {
private:
    std::vector<Triangle> faces;
    double a, H;
    vec3 center;
public:
    // center of base triangle, dir to obtain a vertex of base triangle, the normal to the dir
    Tetrahedron(vec3 baseCenter, vec3 dir, vec3 normal) {
        vec3 bottom1 = baseCenter + dir;
        vec3 perpendicular = cross(dir, normal);
        vec3 dir2 = cos(M_PI * 0.6666) * normalize(dir) + sin(M_PI * 0.6666) * normalize(perpendicular);
        dir2 = dir2 * length(dir);
        vec3 bottom2 = baseCenter + dir2;
        perpendicular = cross(dir2, normal);
        vec3 dir3 = cos(M_PI * 0.6666) * normalize(dir2) + sin(M_PI * 0.6666) * normalize(perpendicular);
        dir3 = dir3 * length(dir);
        vec3 bottom3 = baseCenter + dir3;
        a = 2 * length(dir) * cos(M_PI * 0.16666);
        H = sqrt(6) / 3 * a;
        vec3 apex = baseCenter + normalize(normal) * H;
        center = (bottom1 + bottom2 + bottom3 + apex) / 4;
        faces.push_back(Triangle(bottom1, bottom2, bottom3, center));
        faces.push_back(Triangle(bottom1, bottom2, apex, center));
        faces.push_back(Triangle(bottom2, bottom3, apex, center));
        faces.push_back(Triangle(bottom3, bottom1, apex, center));
    }
    
    Hit intersect(const Ray& ray) {
        Hit best = Hit();
        
        for (int i = 0; i < (int)faces.size(); i++) {
            Hit other = faces[i].intersect(ray);
            
            if (other < best)
                best = other;
        }
            
        return best;
    }
};

class Octahedron : public Intersectable {
private:
    std::vector<Triangle> faces;
    vec3 center;
    vec3 normal;
    vec3 dir;
    double a, H;
public:
    Octahedron(vec3 center, vec3 dir, vec3 normal) {
        this->center = center;
        this->dir = dir;
        this->normal = normal;
        vec3 bottom1 = center + dir;
        vec3 perpendicular = cross(dir, normal);
        vec3 dir2 = cos(M_PI * 0.5) * normalize(dir) + sin(M_PI * 0.5) * normalize(perpendicular);
        dir2 = dir2 * length(dir);
        vec3 bottom2 = center + dir2;
        
        perpendicular = cross(dir2, normal);
        vec3 dir3 = cos(M_PI * 0.5) * normalize(dir2) + sin(M_PI * 0.5) * normalize(perpendicular);
        dir3 = dir3 * length(dir);
        vec3 bottom3 = center + dir3;
        
        perpendicular = cross(dir3, normal);
        vec3 dir4 = cos(M_PI * 0.5) * normalize(dir3) + sin(M_PI * 0.5) * normalize(perpendicular);
        dir4 = dir4 * length(dir);
        vec3 bottom4 = center + dir4;
        
        a = sin(M_PI / 4) * length(dir) * 2;
        H = sqrt(a * a / 2);
        vec3 apex1 = center + H * normalize(normal);
        vec3 apex2 = center - H * normalize(normal);
        
        faces.push_back(Triangle(bottom1, bottom2, apex1, center));
        faces.push_back(Triangle(bottom2, bottom3, apex1, center));
        faces.push_back(Triangle(bottom3, bottom4, apex1, center));
        faces.push_back(Triangle(bottom4, bottom1, apex1, center));
        faces.push_back(Triangle(bottom1, bottom2, apex2, center));
        faces.push_back(Triangle(bottom2, bottom3, apex2, center));
        faces.push_back(Triangle(bottom3, bottom4, apex2, center));
        faces.push_back(Triangle(bottom4, bottom1, apex2, center));
    }
    
    Hit intersect(const Ray& ray) {
        Hit best = Hit();
        
        for (int i = 0; i < (int)faces.size(); i++) {
            Hit other = faces[i].intersect(ray);
            
            if (other < best)
                best = other;
        }
            
        return best;
    }
};

class Cube : public Intersectable {
private:
    vec3 center;
    vec3 nose;
    vec3 ear;
    std::vector<ConvexQuadrilateral> faces;
    
    void genCube() {
        // variable naames assume left ear provided
        vec3 noseDir = nose - center;
        vec3 leftEarDir = ear - center;
        // either forehead, or chin, depending on whether right or left ear was provided on construction
        vec3 foreHeadDir = cross(leftEarDir, noseDir) * 1 / length(noseDir);
        vec3 foreHead = nose + foreHeadDir;
        vec3 rightEarDir = -1 * leftEarDir;
        vec3 topRightFrontCorner = foreHead + rightEarDir;
        vec3 topLeftFrontCorner = foreHead + leftEarDir;
        vec3 tailDir = -1 * noseDir;
        vec3 tail = center + tailDir;
        vec3 topRightBackCorner = tail + foreHeadDir + rightEarDir;
        vec3 topLeftBackCorner = tail + foreHeadDir + leftEarDir;
        vec3 chinDir = -1 * foreHeadDir;
        vec3 bottomRightFrontCorner = nose + chinDir + rightEarDir;
        vec3 bottomLeftFrontCorner = nose + chinDir + leftEarDir;
        vec3 bottomRightBackCorner = tail + chinDir + rightEarDir;
        vec3 bottomLeftBackCorner = tail + chinDir + leftEarDir;

        // face
        faces.push_back(ConvexQuadrilateral(topRightFrontCorner, bottomRightFrontCorner, bottomLeftFrontCorner, topLeftFrontCorner, center));
        // top
        faces.push_back(ConvexQuadrilateral(topRightFrontCorner, topRightBackCorner, topLeftBackCorner, topLeftFrontCorner, center));
        // rear
        faces.push_back(ConvexQuadrilateral(topRightBackCorner, bottomRightBackCorner, bottomLeftBackCorner, topLeftBackCorner, center));
        // bottom
        faces.push_back(ConvexQuadrilateral(bottomRightFrontCorner, bottomRightBackCorner, bottomLeftBackCorner, bottomLeftFrontCorner, center));
        // right
        faces.push_back(ConvexQuadrilateral(topRightFrontCorner, bottomRightFrontCorner, bottomRightBackCorner, topRightBackCorner, center));
        // left
        faces.push_back(ConvexQuadrilateral(topLeftFrontCorner, bottomLeftFrontCorner, bottomLeftBackCorner, topLeftBackCorner, center));
    }
public:
    //either left or right ear can be specified
    Cube(vec3 center, vec3 nose, vec3 ear) {
        this->center = center;
        this->nose = nose;
        this->ear = ear;
        genCube();
    }

    Hit intersect(const Ray& ray) {
        Hit best = Hit();
        
        for (int i = 0; i < (int)faces.size(); i++) {
            Hit other = faces[i].intersect(ray);
            
            if (other < best)
                best = other;
        }
            
        return best;
    }
};


class Camera {
    vec3 eye, lookat, right, up;
public:
    void set(vec3 _eye, vec3 _lookat, vec3 vup, float fov) {
        eye = _eye;
        lookat = _lookat;
        vec3 w = eye - lookat;
        float focus = length(w);
        right = normalize(cross(vup, w)) * focus * tanf(fov / 2);
        up = normalize(cross(w, right)) * focus * tanf(fov / 2);
    }
    
    Ray getRay(int X, int Y) {
        vec3 dir = lookat + right * (2.0f * (X + 0.5f) / windowWidth - 1) + up * (2.0f * (Y + 0.5f) / windowHeight - 1) - eye;
        return Ray(eye, dir);
    }
};

class Wall : public ConvexQuadrilateral {
public:
    Wall(vec3 first, vec3 second, vec3 third, vec3 fourth, const vec3& origin) : ConvexQuadrilateral(first, second,  third, fourth,origin){
    }
    
    bool transparentFromOutside() override {
        return true;
    }
};

class Lamp : public Intersectable {
private:
    vec3 colour;
    double a,b,c;
    vec3 n;
    vec3 position;
    double angle;
    double h;
public:
    vec3 getColour() {
        return colour;
    }
    
    vec3 getPosition() {
        return position;
    }
    
    vec3 getN() {
        
        return n;
    }
    
    double getAngle() {
        return angle;
    }
    
    vec3 insideColour() override {
        return colour;
    }
    
    virtual vec3 outsideColour() {
        return colour;
    }
    
    Lamp(vec3 colour, vec3 position, vec3 n, double angle, double h) {
        this->colour = colour;
        this->position = position;
        this->angle = angle;
        this->n = normalize(n);
        this->h = h;
    }
    
    Hit intersect(const Ray& ray) override {
        double a = pow(dot(ray.getDirection(), n), 2) - dot(ray.getDirection(), ray.getDirection()) * pow(cos(angle), 2);
        double b = 2 * dot(ray.getDirection(), n) * dot(n, ray.getStart() - position) -2 * dot(ray.getDirection(), ray.getStart() - position) * pow(cos(angle), 2);
        double c = pow(dot(n, ray.getStart() - position), 2) - dot(ray.getStart() - position, ray.getStart() - position) * pow(cos(angle), 2);
        
        double D = b * b - 4 * a * c;
        
        if (D < 0)
            return Hit();
        
        double t1 = (-1 * b + sqrt(D)) / (2 * a);
        double t2 = (-1 * b - sqrt(D)) / (2 * a);
        vec3 int1 = t1 * ray.getDirection() + ray.getStart();
        vec3 int2 = t2 * ray.getDirection() + ray.getStart();
        vec3 n1 = 2 * dot(int1 - position, n) * n - 2 * (int1 - position) * pow(cos(angle), 2);
        n1 = n1 *-1;
        vec3 n2 = 2 * dot(int2 - position, n) * n - 2 * (int2 - position) * pow(cos(angle), 2);
        n2 = n2 *-1;
        double s1 = dot(int1 - position, n);
        double s2 = dot(int2 - position, n);
        bool g1 = s1 >= 0 && s1 <= this->h && t1 > 0;
        bool g2 = s2 >= 0 && s2 <= this->h && t2 > 0;
        double dst1 = length(ray.getStart() - int1);
        double dst2 = length(ray.getStart() - int2);
        
        if (!g1 && !g2)
            return Hit();
        else if (g1 && g2)
            return dst1 < dst2 ? Hit(this, ray.getStart(), int1, n1) : Hit(this, ray.getStart(), int2, n2);
        else if (g1)
            return Hit(this, ray.getStart(), int1, n1);
        else if (g2)
            return Hit(this, ray.getStart(), int2, n2);
        else return Hit(this, ray.getStart(), int1, n1);
    }
};

float rnd() { return (float)rand() / RAND_MAX; }

class Scene {
    std::vector<Intersectable *> objects;
    std::vector<Lamp *> lights;
    Camera camera;
public:
    void onClick(int x, int y) {
        Ray r = camera.getRay(x, y);
        Hit h =  firstIntersect(r);
        if (!h.getSuccess()) return;
        
        for (int i = 0; i < (int)lights.size(); i++)
            if (h.getObject() == lights[i])
                return;
        
        Lamp* l = 0;
        
        for (int i = 0; i < (int)lights.size(); i++) {
            if (l == 0)
                l = lights[i];
            else if (length(lights[i]->getPosition() - h.getIntersection()) < length(l->getPosition() - h.getIntersection()))
                l = lights[i];
        }
        
        if (l == 0) return;
        
        for (auto it = lights.begin(); it != lights.end(); it++)
            if (*it == l) {
                lights.erase(it);
                break;
            }
        
        for (auto it = objects.begin(); it != objects.end(); it++)
            if (*it == l) {
                objects.erase(it);
                break;
            }
        Lamp* l2 = new Lamp(l->getColour(), h.getIntersection(), h.getNormal(), M_PI / 9, 0.4);
        delete l;
        objects.push_back(l2);
        lights.push_back(l2);
    }
    
    void build() {
        vec3 eye = vec3(3.7, 1, 5.4), vup = vec3(0, 1, 0), lookat = vec3(-0.8, 0.2, 0);
        float fov = 35 * M_PI / 180;
        camera.set(eye, lookat, vup, fov);
        
        //back
        objects.push_back(new Wall(
                                                  vec3(-1.2, -0.5, -0.5),
                                                  vec3(-1.2, 1.5, -0.5),
                                                  vec3(1.2, 1.5, -0.5),
                                                  vec3(1.2, -0.5, -0.5),
                                                  vec3(0, 0, -1)));
        //front
        objects.push_back(new Wall(
                                                  vec3(-1.2, -0.5, 2),
                                                  vec3(-1.2, 1.5, 2),
                                                  vec3(1.2, 1.5, 2),
                                                  vec3(1.2, -0.5, 2),
                                                  vec3(0, 0, 3)));
        //floor
        objects.push_back(new Wall(
                                                  vec3(-1.2, -0.5, -0.5),
                                                  vec3(1.2, -0.5, -0.5),
                                                  vec3(1.2, -0.5, 2),
                                                  vec3(-1.2, -0.5, 2),
                                                  vec3(0,-0.6, 0)));
        
        flr = (Wall*)objects[objects.size() - 1];
        // top
        objects.push_back(new Wall(
                                                  vec3(-1.2, 1.5, -0.5),
                                                  vec3(1.2, 1.5, -0.5),
                                                  vec3(1.2, 1.5, 2),
                                                  vec3(-1.2, 1.5, 2),
                                                  vec3(0,1.6, 0)));
        // right
        objects.push_back(new Wall(
                                                  vec3(1.2, -0.5, -0.5),
                                                  vec3(1.2, 1.5, -0.5),
                                                  vec3(1.2, 1.5, 2),
                                                  vec3(1.2, -0.5, 2),
                                                  vec3(1.3,0, 0)));
        // left
        objects.push_back(new Wall(
                                                  vec3(-1.2, -0.5, -0.5),
                                                  vec3(-1.2, 1.5, -0.5),
                                                  vec3(-1.2, 1.5, 2),
                                                  vec3(-1.2, -0.5, 2),
                                                  vec3(-1.3,0, 0)));
        objects.push_back(new Cube(vec3(-0.3, -0.3, -0.3), vec3(-0.3, -0.3, -0.1), vec3(-0.5, -0.3, -0.3)));
        objects.push_back(new Tetrahedron(vec3(0.7, -0.5, 0.4), vec3(0.2, 0, 0.4), vec3(0, 1, 0)));
        objects.push_back(new Octahedron(vec3(-0.3, 0, 1.5), vec3(0.2,0,0.3), vec3(0,1,0)));
        objects.push_back(new Lamp(vec3(0,1,0), vec3 (0, 1.5, 1), vec3(0, -1, 0), M_PI / 9, 0.4));
        lights.push_back((Lamp*)objects[objects.size() - 1]);
        objects.push_back(new Lamp(vec3(1,0,0), vec3 (0.5, 1.5, 1), vec3(0, -1, 0), M_PI / 9, 0.4));
        lights.push_back((Lamp*)objects[objects.size() - 1]);
        objects.push_back(new Lamp(vec3(0,0,1), vec3 (1, 1.5, 1), vec3(0, -1, 0), M_PI / 9, 0.4));
        lights.push_back((Lamp*)objects[objects.size() - 1]);
        
        objects.push_back(new Lamp(vec3(0.2,0.2,0.2), vec3 (0.4, 0.4, 1.3), vec3(0, -1, 0), M_PI / 10, 0.9));
       
        //for (int i = 0; i < 500; i++)
        //    objects.push_back(new Sphere(vec3(rnd() - 0.5f, rnd() - 0.5f, rnd() - //0.5f), rnd() * 0.1f, material));
    }

    void render(std::vector<vec4>& image) {
        for (int Y = 0; Y < (int)windowHeight; Y++) {
            for (int X = 0; X < (int)windowWidth; X++) {
                vec3 color = trace(camera.getRay(X, Y));
                image[Y * windowWidth + X] = vec4(color.x, color.y, color.z, 1);
            }
        }
    }

    Hit firstIntersect(Ray ray) {
        Hit bestHit;
        
        for (Intersectable * object : objects) {
            Hit hit = object->intersect(ray);
            
            if (hit < bestHit) {
                double s = dot(hit.getNormal(), ray.getStart() - hit.getIntersection());
                
                if (s >= 0)
                    bestHit = hit;
                else if (!hit.getObject()->transparentFromOutside())
                    bestHit = hit;
                
            }
        }
        
        return bestHit;
    }

    vec3 trace(Ray ray, int depth = 0) {
        Hit h = firstIntersect(ray);
        vec3 c = vec3(0.2, 0.2, 0.2);
        
        if (h.getSuccess())
        {
            
            for (int i = 0; i < (int)lights.size(); i++) {
                vec3 v = normalize(h.getIntersection() - lights[i]->getPosition());
                
                double s = acos(dot(v, lights[i]->getN()));
                
                if (s <= lights[i]->getAngle() - 0.01) {
                    Ray r2 = Ray(h.getIntersection() + h.getNormal() * 0.005, lights[i]->getPosition() - h.getIntersection());
                    Hit h2 = firstIntersect(r2);
       
                    if (length(h2.getIntersection() - lights[i]->getPosition()) < 0.005 || h2.getObject() == lights[i])
                        c = c + lights[i]->getColour() / pow(length(lights[i]->getPosition() - h.getIntersection()), 1.5);
                }
            }
            
            double d = dot(normalize(h.getNormal()), normalize(ray.getStart() - h.getIntersection()));
                        
            if (d < 0 && !h.getObject()->transparentFromOutside()) {
                d = fabs(d);
                c = c + d * h.getObject()->insideColour();
                return c;
            }
            else return c + d * h.getObject()->outsideColour();
        }
        else return vec3(45.0 / 255, 27.0 / 255, 88.0 / 255);
    }
};

GPUProgram gpuProgram; // vertex and fragment shaders
Scene scene;

// vertex shader in GLSL
const char *vertexSource = R"(
    #version 330
    precision highp float;

    layout(location = 0) in vec2 cVertexPosition;    // Attrib Array 0
    out vec2 texcoord;

    void main() {
        texcoord = (cVertexPosition + vec2(1, 1))/2;                            // -1,1 to 0,1
        gl_Position = vec4(cVertexPosition.x, cVertexPosition.y, 0, 1);         // transform to clipping space
    }
)";

// fragment shader in GLSL
const char *fragmentSource = R"(
    #version 330
    precision highp float;

    uniform sampler2D textureUnit;
    in  vec2 texcoord;            // interpolated texture coordinates
    out vec4 fragmentColor;        // output that goes to the raster memory as told by glBindFragDataLocation

    void main() {
        fragmentColor = texture(textureUnit, texcoord);
    }
)";

class FullScreenTexturedQuad {
    unsigned int vao;    // vertex array object id and texture id
    unsigned int vbo;
    Texture texture;
    int w, h;
public:
    FullScreenTexturedQuad(int windowWidth, int windowHeight, std::vector<vec4>& image)
        : texture(windowWidth, windowHeight, image)
    {
        w = windowWidth;
        h = windowHeight;
        glGenVertexArrays(1, &vao);    // create 1 vertex array object
        glBindVertexArray(vao);        // make it active

        //unsigned int vbo;        // vertex buffer objects
        glGenBuffers(1, &vbo);    // Generate 1 vertex buffer objects

        // vertex coordinates: vbo0 -> Attrib Array 0 -> vertexPosition of the vertex shader
        glBindBuffer(GL_ARRAY_BUFFER, vbo); // make it active, it is an array
        float vertexCoords[] = { -1, -1,  1, -1,  1, 1,  -1, 1 };    // two triangles forming a quad
        glBufferData(GL_ARRAY_BUFFER, sizeof(vertexCoords), vertexCoords, GL_STATIC_DRAW);       // copy to that part of the memory which is not modified
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, NULL);     // stride and offset: it is tightly packed
    }
    
    void update(std::vector<vec4>& image, int windowWidth, int windowHeight) {
        texture = Texture(windowWidth, windowHeight,image);
    }

    void Draw() {
        glBindVertexArray(vao);    // make the vao and its vbos active playing the role of the data source
        gpuProgram.setUniform(texture, "textureUnit");
        glDrawArrays(GL_TRIANGLE_FAN, 0, 4);    // draw two triangles forming a quad
    }
};

FullScreenTexturedQuad * fullScreenTexturedQuad;

// Initialization, create an OpenGL context
void onInitialization() {
    glViewport(0, 0, windowWidth, windowHeight);
    scene.build();

    //std::vector<vec4> image(windowWidth * windowHeight);
    //scene.render(image);

    // copy image to GPU as a texture
    //fullScreenTexturedQuad = new FullScreenTexturedQuad(windowWidth, windowHeight, image);

    // create program for the GPU
    gpuProgram.create(vertexSource, fragmentSource, "fragmentColor");
}

// Window has become invalid: Redraw
void onDisplay() {
    std::vector<vec4> image(windowWidth * windowHeight);
    scene.render(image);
    fullScreenTexturedQuad = new FullScreenTexturedQuad(windowWidth, windowHeight, image);
    fullScreenTexturedQuad->Draw();
    glutSwapBuffers();
    // exchange the two buffers
}

// Key of ASCII code pressed
void onKeyboard(unsigned char key, int pX, int pY) {
}

// Key of ASCII code released
void onKeyboardUp(unsigned char key, int pX, int pY) {

}

// Mouse click event
void onMouse(int button, int state, int pX, int pY) {
    if (state == GLUT_DOWN) {
        scene.onClick(pX, windowHeight - pY);
        glutPostRedisplay();
    }
}

// Move mouse with key pressed
void onMouseMotion(int pX, int pY) {
}

// Idle event indicating that some time elapsed: do animation here
void onIdle() {
}
