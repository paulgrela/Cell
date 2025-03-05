
#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <unordered_set>
#include <algorithm>

#include "CellEngineMolecularDynamicsSimulationForceField2.h"

struct Vec3 {
    double x, y, z;
    Vec3() : x(0), y(0), z(0) {}
    Vec3(double x, double y, double z) : x(x), y(y), z(z) {}
    Vec3 operator+(const Vec3& other) const { return {x+other.x, y+other.y, z+other.z}; }
    Vec3 operator-(const Vec3& other) const { return {x-other.x, y-other.y, z-other.z}; }
    Vec3 operator*(double s) const { return {x*s, y*s, z*s}; }
    Vec3 operator/(double s) const { return {x/s, y/s, z/s}; }
    Vec3& operator+=(const Vec3& other) { x += other.x; y += other.y; z += other.z; return *this; }
    Vec3& operator-=(const Vec3& other) { x -= other.x; y -= other.y; z -= other.z; return *this; }
    Vec3& operator*=(const Vec3& other) { x *= other.x; y *= other.y; z *= other.z; return *this; }
    Vec3& operator*=(const double s) { x *= s; y *= s; z *= s; return *this; }
    double dot(const Vec3& other) const { return x*other.x + y*other.y + z*other.z; }
    Vec3 cross(const Vec3& other) const
    {
        return {y*other.z - z*other.y, z*other.x - x*other.z, x*other.y - y*other.x};
    }
    double length() const { return std::sqrt(x*x + y*y + z*z); }
    Vec3 normalized() const { double l = length(); return {x/l, y/l, z/l}; }
};

class Atom {
public:
    int type; // 0=O, 1=N, 2=P, 3=C
    Vec3 pos, vel, force;
    double charge, mass, sigma, epsilon;

    Atom(int type, const Vec3& pos, double charge) : type(type), pos(pos), charge(charge) {
        switch(type) {
            case 0: mass = 16.0; sigma = 3.0; epsilon = 0.2; break;
            case 1: mass = 14.0; sigma = 3.2; epsilon = 0.1; break;
            case 2: mass = 31.0; sigma = 3.5; epsilon = 0.15; break;
            case 3: mass = 12.0; sigma = 3.4; epsilon = 0.05; break;
            default: mass = 12.0; sigma = 3.0; epsilon = 0.1;
        }
    }
};

//struct Bond { int i, j; double r0, k; };
//struct Angle { int i, j, k; double theta0, k_theta; };
//struct Dihedral { int i, j, k, l; double k_phi, n, phi0; };

struct Bond
{
int i, j; double r0; // equilibrium distance
double k; // spring constant

Bond(int i, int j, double r0, double k) : i(i), j(j), r0(r0), k(k) {}
};


struct Angle
{
int i, j, k; double theta0; // equilibrium angle in radians
double k_theta; // force constant

Angle(int i, int j, int k, double theta0, double k_theta) : i(i), j(j), k(k), theta0(theta0), k_theta(k_theta) {}
};


struct Dihedral
{
int i, j, k, l; double k_phi; // force constant
int n; // multiplicity
double phi0; // phase angle in radians

Dihedral(int i, int j, int k, int l, double k_phi, int n, double phi0) : i(i), j(j), k(k), l(l), k_phi(k_phi), n(n), phi0(phi0) {}
};


class Simulation
{
    std::vector<Atom> atoms;
    std::vector<Bond> bonds;
    std::vector<Angle> angles;
    std::vector<Dihedral> dihedrals;
    std::unordered_set<int> bonded_pairs;
    double box_size, temperature, time_step;
    const double k_coulomb = 332.0636, k_boltzmann = 0.0019872041;

    void computeBondForces()
    {
        for (const auto& b : bonds)
        {
            Atom &a = atoms[b.i], &b_atom = atoms[b.j];
            Vec3 delta = b_atom.pos - a.pos;
            double r = delta.length();
            double force_mag = -b.k * (r - b.r0) / r;
            Vec3 force = delta * force_mag;
            a.force += force;
            b_atom.force -= force;
        }
    }

    void computeAngleForces()
    {
        for (const auto& a : angles) {
            Atom &i = atoms[a.i], &j = atoms[a.j], &k = atoms[a.k];
            Vec3 r_ij = i.pos - j.pos, r_kj = k.pos - j.pos;
            double len_ij = r_ij.length(), len_kj = r_kj.length();
            double cos_theta = r_ij.dot(r_kj) / (len_ij * len_kj);
            double theta = std::acos(cos_theta);
            double dtheta = theta - a.theta0;
            double dE_dtheta = a.k_theta * dtheta;

            Vec3 dcos_di = (r_kj * len_ij - r_ij * (r_ij.dot(r_kj)/len_ij)) / (len_ij*len_ij*len_kj);
            Vec3 dcos_dk = (r_ij * len_kj - r_kj * (r_ij.dot(r_kj)/len_kj)) / (len_kj*len_kj*len_ij);
            Vec3 dcos_dj = dcos_di + dcos_dk;
            dcos_dj *= -1.0f;

            Vec3 force_i = dcos_di * (-dE_dtheta / std::sqrt(1 - cos_theta*cos_theta));
            Vec3 force_j = dcos_dj * (-dE_dtheta / std::sqrt(1 - cos_theta*cos_theta));
            Vec3 force_k = dcos_dk * (-dE_dtheta / std::sqrt(1 - cos_theta*cos_theta));

            i.force += force_i;
            j.force += force_j;
            k.force += force_k;
        }
    }

    void computeDihedralForces()
    {
        for (const auto& d : dihedrals)
        {
            Atom &i = atoms[d.i], &j = atoms[d.j], &k = atoms[d.k], &l = atoms[d.l];
            Vec3 r_ij = i.pos - j.pos, r_kj = k.pos - j.pos, r_kl = l.pos - k.pos;
            Vec3 n1 = r_ij.cross(r_kj), n2 = r_kj.cross(r_kl);
            double n1_len = n1.length(), n2_len = n2.length();
            if (n1_len == 0 || n2_len == 0) continue;

            double cos_phi = n1.dot(n2) / (n1_len * n2_len);
            double phi = std::acos(cos_phi);
            double sign = n1.cross(n2).dot(r_kj) < 0 ? -1 : 1;
            phi *= sign;

            double dV_dphi = d.k_phi * d.n * std::sin(d.n * phi - d.phi0);

            Vec3 f1 = (n1 * (d.n * d.k_phi / n1_len)).cross(r_kj);
            Vec3 f4 = (n2 * (d.n * d.k_phi / n2_len)).cross(r_kj);
            Vec3 f2 = (f1.cross(r_ij) - f4.cross(r_kl)) / r_kj.dot(r_kj);
//            Vec3 f3 = -f2 - f4;
            Vec3 f3 = (f2 + f4) * -1.0;

            i.force += f1 * dV_dphi;
//            j.force += (-f1 + f2) * dV_dphi;
            j.force += (f1 * -1.0 + f2) * dV_dphi;
//            k.force += (-f3 - f4) * dV_dphi;
            k.force += (f3 + f4) * -1.0 * dV_dphi;
            l.force += f4 * dV_dphi;
        }
    }

    void computeNonBondedForces()
    {
        size_t n = atoms.size();
        for (size_t i = 0; i < n; ++i)
        {
            for (size_t j = i+1; j < n; ++j)
            {
                if (bonded_pairs.count(i*n + j))
                   continue;

                Atom &a = atoms[i], &b = atoms[j];
                Vec3 delta = a.pos - b.pos;
                double r = delta.length();
                if (r == 0)
                   continue;

                // Coulomb interaction
                double coulomb = k_coulomb * a.charge * b.charge / (r*r*r);
                Vec3 f_coulomb = delta * coulomb;
                a.force += f_coulomb;
                b.force -= f_coulomb;

                // Lennard-Jones interaction
                double sigma_avg = (a.sigma + b.sigma)/2;
                double epsilon_avg = std::sqrt(a.epsilon * b.epsilon);
                double sr6 = std::pow(sigma_avg/r, 6);
                double lj = 24 * epsilon_avg * (2*sr6*sr6 - sr6) / (r*r);
                Vec3 f_lj = delta * lj;
                a.force += f_lj;
                b.force -= f_lj;
            }
        }
    }

public:
    Simulation(double box, double temp, double dt) : box_size(box), temperature(temp), time_step(dt) {}

    void initializeAtoms(int n)
    {
        std::mt19937 gen(std::random_device{}());
        std::uniform_real_distribution<> pos(0, box_size), charge(-1, 1);
        std::uniform_int_distribution<> type(0, 3);

        for (int i = 0; i < n; ++i)
        {
            atoms.emplace_back(type(gen), Vec3(pos(gen), pos(gen), pos(gen)), charge(gen));
            double scale = std::sqrt(3 * k_boltzmann * temperature / atoms.back().mass);
            std::normal_distribution<> vel(0, scale);
            atoms.back().vel = Vec3(vel(gen), vel(gen), vel(gen));
        }

        // Create molecular topology (example: linear chain)
        for (int i = 0; i < n/2; ++i)
        {
            bonds.push_back({i, i+1, 1.5, 100});
            bonded_pairs.insert(i*n + i+1);
            if (i < n/2 - 2)
            {
                angles.push_back({i, i+1, i+2, 1.5708, 100});
                dihedrals.push_back({i, i+1, i+2, i+3, 5, 2, 0});
            }
        }
    }

void initializeAtoms1(int num_atoms)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> pos_dist(0.0, box_size);
    std::uniform_int_distribution<> type_dist(0, 3); // 0=O, 1=N, 2=P, 3=C
    std::uniform_real_distribution<> charge_dist(-1.0, 1.0);

    for (int i = 0; i < num_atoms; ++i)
    {
        Vec3 pos(pos_dist(gen), pos_dist(gen), pos_dist(gen));
        int type = type_dist(gen);
        double charge = charge_dist(gen);
        atoms.emplace_back(type, pos, charge);
    }

    // Assign random velocities based on Maxwell-Boltzmann distribution
    std::normal_distribution<> vel_dist(0.0, std::sqrt(k_boltzmann * temperature));
    for (auto& atom : atoms)
    {
        atom.vel.x = vel_dist(gen);
        atom.vel.y = vel_dist(gen);
        atom.vel.z = vel_dist(gen);
        // Subtract center-of-mass velocity to avoid drift
        // (omitted for brevity)
    }

    // Generate random bonds (example: 500 bonds)
    std::uniform_int_distribution<> atom_dist(0, num_atoms - 1);

    for (int b = 0; b < 500; ++b)
    {
        int i = atom_dist(gen);
        int j = atom_dist(gen);
        if (i != j)
        {
            // Ensure i < j to avoid duplicates
            if (i > j)
               std::swap(i, j);
            // Check if bond already exists
            int key = i * num_atoms + j;
            if (bonded_pairs.find(key) == bonded_pairs.end())
            {
                bonded_pairs.insert(key);
                // Assign random r0 and k (example values)
                double r0 = 1.5; // Angstrom
                double k = 100.0; // kcal/(mol·Å²)
                //bonds.emplace_back(i, j, r0, k);
            }
        }
    }

    // Generate angles and dihedrals (example)
    // For angles, find triplets i-j-k where i-j and j-k are bonded
    // This is a simplified approach; realistically, need to find connected triplets
    // For demonstration, randomly select triplets

    for (int a = 0; a < 300; ++a)
    {
        int i = atom_dist(gen);
        int j = atom_dist(gen);
        int k = atom_dist(gen);
        if (i != j && j != k && i != k) {
            angles.emplace_back(i, j, k,
            1.5708 // 90 degree
            , 100.0 // k_theta
            );
        }
    }

    // Similarly for dihedrals
    for (int d = 0; d < 200; ++d)
    {
        int i = atom_dist(gen);
        int j = atom_dist(gen);
        int k = atom_dist(gen);
        int l = atom_dist(gen);
        if (i != j && j != k && k != l && i != k && i != l && j != l) {
            dihedrals.emplace_back(i, j, k, l,
            5.0 // k_phi
            ,
            2 // n
            ,
            0.0 // phi0
            );
        }
    }

}


    void computeForces()
    {
        for (auto& a : atoms)
            a.force = Vec3();

        computeBondForces();
        computeAngleForces();
        computeDihedralForces();
        computeNonBondedForces();
    }

    void integrate()
    {
        // Velocity half-step
        for (auto& a : atoms)
        {
            a.vel += a.force * (time_step / (2 * a.mass));
        }

        // Position update with periodic boundary conditions
        for (auto& a : atoms)
        {
            a.pos += a.vel * time_step;
            a.pos.x = fmod(a.pos.x + box_size, box_size);
            a.pos.y = fmod(a.pos.y + box_size, box_size);
            a.pos.z = fmod(a.pos.z + box_size, box_size);

            if (a.pos.x < 0) a.pos.x += box_size;
            if (a.pos.y < 0) a.pos.y += box_size;
            if (a.pos.z < 0) a.pos.z += box_size;
        }

        // Compute new forces
        computeForces();

        // Velocity second half-step
        for (auto& a : atoms)
        {
            a.vel += a.force * (time_step / (2 * a.mass));
        }
    }

    void applyBerendsenThermostat(double target_temp, double tau)
    {
        double current_temp = 0;
        for (const auto& a : atoms)
        {
            current_temp += a.mass * a.vel.dot(a.vel);
        }
        current_temp /= (3 * atoms.size() * k_boltzmann);

        double lambda = std::sqrt(1 + time_step/tau * (target_temp/current_temp - 1));
        for (auto& a : atoms)
        {
            a.vel = a.vel * lambda;
        }
    }
};

int main2()
{
    const double box_size = 100.0;    // Angstroms
    const double temperature = 303.15; // 30°C in Kelvin
    const double time_step = 0.001;    // picoseconds
    const int num_steps = 1000;
    const int output_interval = 100;

    Simulation sim(box_size, temperature, time_step);
    sim.initializeAtoms(1000);
    sim.computeForces();

    for (int step = 0; step < num_steps; ++step)
    {
        sim.integrate();
        sim.applyBerendsenThermostat(temperature, 0.1);

        if (step % output_interval == 0)
        {
            std::cout << "Step " << step << " completed" << std::endl;
        }
    }

    return 0;
}

