#ifndef CONTACT_H_
#define CONTACT_H_

#include "interaction.hpp"
#include "physics_state.hpp"
#include "simulable.hpp"

enum GEOMETRY_TYPE {
SPHERE = 0,
INFPLANE = 1,
FINPLANE = 2,
};

struct ContactGeometry {
    virtual double distance_point(const Vec3& point, bool& valid) const = 0;
    virtual Vec3 outward_direction(const Vec3& point) const = 0;
};

struct Sphere : public ContactGeometry {
        Sphere(Vec3 center, double r) : r(r), center(center) {}
        double r;
        Vec3 center;
        double distance_point(const Vec3& point, bool& valid) const override;
        Vec3 outward_direction(const Vec3& point) const override;
};

struct InfPlane : public ContactGeometry {
        InfPlane(Vec3 position, Vec3 normal) : center(position), normal(normal) {}
        Vec3 normal;
        Vec3 center;
        double distance_point(const Vec3& point, bool& valid) const override;
        Vec3 outward_direction(const Vec3& point) const override;
};

struct FinPlane : public ContactGeometry {
        Vec3 normal;
        Vec3 center;
        Vec3 up;
        double radius;
        double distance_point(const Vec3& point, bool& valid) const override;
        Vec3 outward_direction(const Vec3& point) const override;
};

class Contact : public Interaction {
    public:
        /* Contact class does NOT own the ContactGeometry pointer */
        Contact(const Simulable& simulable, const Sphere* sphere);
        Contact(const Simulable& simulable, const InfPlane* plane);
        Contact(const Simulable& simulable, const FinPlane* plane);

        ~Contact();

        void apply(const PhysicsState* state, EnergyDerivatives* f) override;

        inline void set_stiffness(double stiffness) { contact_stiffness = stiffness; }

    private:
        Vec3 force(const Vec3& direction, const double dist);

        Mat3 force_derivative(const PhysicsState* state, const Vec3& direction, const double dist);

        GEOMETRY_TYPE geometry_type;
        const ContactGeometry* geometry;
        double contact_stiffness = 100000.0f;
};
#endif // CONTACT_H_
