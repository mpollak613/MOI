/* Copyright 2020 Michael Pollak.
 *
 * Use of this source code is governed by an MIT-style
 * licence that can be found in the LICENSE file.
 */

#ifndef BODY_H
#define BODY_H

#include <cmath>
#include <string>
#include <valarray>
#include <vector>


namespace hyx
{
    template<typename T>
    inline T euclidNorm(std::valarray<T> A)
    {
        return std::sqrt((A * A).sum());
    }

    class Body
    {
    private:
        std::string _name;
        std::valarray<double> _pos;
        std::valarray<double> _vel;
        double _mass;
        double _rad;
        std::valarray<float> _color;

    public:

        // constuctors

        Body(std::string name = "temp",
            std::valarray<double> pos = { 0,0,0 },
            std::valarray<double> vel = { 0,0,0 },
            double mass = 0,
            double rad = 0,
            std::valarray<float> color = { 1,1,1 })
            : _name(name),
            _pos(pos),
            _vel(vel),
            _mass(mass),
            _rad(rad),
            _color(color)
        {
            // empty
        }

        Body(Body const& copy)
        {
            *this = copy;
        }

        // NO custom destructor


        // getters

        std::string getName() const
        {
            return _name;
        }

        std::valarray<double> getPos() const
        {
            return _pos;
        }

        std::valarray<double> getVel() const
        {
            return _vel;
        }

        double getMass() const
        {
            return _mass;
        }

        double getRad() const
        {
            return _rad;
        }

        std::valarray<float> getColor() const
        {
            return _color;
        }

        // (operator; =, ==, !=)

        Body& operator= (Body const& rhs)
        {
            this->_name = rhs.getName();
            this->_pos = rhs.getPos();
            this->_vel = rhs.getVel();
            this->_mass = rhs.getMass();
            this->_rad = rhs.getRad();
            this->_color = rhs.getColor();

            return *this;
        }

        friend bool operator== (Body const& lhs, Body const& rhs);
        friend bool operator!= (Body const& lhs, Body const& rhs);

        // helper functions

        friend void updateStats(std::vector<Body>& sys, std::vector<std::valarray<double>> const& dr, std::vector<std::valarray<double>> const& dv);
        friend std::vector<Body> copy(std::vector<Body> sys, double nval);
        friend std::vector<Body> copy2(std::vector<Body> sys, std::vector<std::valarray<double>> nval);
    };

    // (operator; ==, !=)

    inline bool operator== (Body const& lhs, Body const& rhs)
    {
        return (lhs.getName() == rhs.getName() &&
            (lhs.getPos() == rhs.getPos()).min() &&
            (lhs.getVel() == rhs.getVel()).min() &&
            lhs.getMass() == rhs.getMass() &&
            lhs.getRad() == rhs.getRad());
    }

    inline bool operator!= (Body const& lhs, Body const& rhs)
    {
        return !(lhs == rhs);
    }

    inline std::valarray<double> getForce(std::vector<Body> const& sys, size_t idx)
    {
        double G = 6.6742e-11; //m^3 * kg^-1 * s^-2
        std::valarray<double> acc = { 0,0,0 };
        std::valarray<double> delt;

        for (size_t i = 0; i < sys.size(); ++i)
        {
            if (i != idx)
            {
                delt = sys[i].getPos() - sys[idx].getPos();
                acc += G * sys[i].getMass() / std::pow(euclidNorm(delt), 3) * delt;
            }
        }

        return acc;
    }

    inline double getEnergy(std::vector<Body> const& sys)
    {
        double G = 6.6742e-11; //m^3 * kg^-1 * s^-2
        double K = 0;
        double U = 0;

        for (size_t i = 0; i < sys.size(); ++i)
        {
            for (size_t j = 0; j < i; ++j)
            {
                K -= G * sys[j].getMass() * sys[i].getMass() / euclidNorm<double>(sys[i].getPos() - sys[j].getPos());
            }

            U += std::pow(euclidNorm<double>(sys[i].getMass() * sys[i].getVel()), 2) / (2 * sys[i].getMass());
        }

        return K + U;
    }

    inline void updateStats(std::vector<Body>& sys, std::vector<std::valarray<double>> const& dr, std::vector<std::valarray<double>> const& dv)
    {
        for (size_t i = 0; i < sys.size(); ++i)
        {
            sys[i]._pos += dr[i];
            sys[i]._vel += dv[i];
        }
    }

    inline std::vector<Body> copy(std::vector<Body> sys, double nval)
    {

        for (size_t i = 0; i < sys.size(); ++i)
        {
            sys[i]._pos += nval * sys[i].getVel();
        }

        return sys;
    }

    inline std::vector<Body> copy2(std::vector<Body> sys, std::vector<std::valarray<double>> nval)
    {
        for (size_t i = 0; i < sys.size(); ++i)
        {
            sys[i]._pos += nval[i];
        }

        return sys;
    }

    inline void stepFWEuler(std::vector<Body>& sys, double dt)
    {
        std::vector<std::valarray<double>> dr(sys.size(), std::valarray<double>(3));
        std::vector<std::valarray<double>> dv(sys.size(), std::valarray<double>(3));

        for (size_t i = 0; i < sys.size(); ++i)
        {
            dr[i] = sys[i].getVel() * dt;
            dv[i] = getForce(sys, i) * dt;
        }

        updateStats(sys, dr, dv);
    }

    inline void stepBWEuler(std::vector<Body>& sys, double dt)
    {
        std::vector<std::valarray<double>> dr(sys.size(), std::valarray<double>(3));
        std::vector<std::valarray<double>> dv(sys.size(), std::valarray<double>(3));


        std::vector<Body> tmpsys = copy(sys, dt);
        for (size_t i = 0; i < sys.size(); ++i)
        {
            dv[i] = getForce(tmpsys, i) * dt;
            dr[i] = (dv[i] + sys[i].getVel()) * dt;
        }

        updateStats(sys, dr, dv);
    }

    inline void stepEulerCromer(std::vector<Body>& sys, double dt)
    {
        std::vector<std::valarray<double>> dr(sys.size(), std::valarray<double>(3));
        std::vector<std::valarray<double>> dv(sys.size(), std::valarray<double>(3));

        for (size_t i = 0; i < sys.size(); ++i)
        {
            dv[i] = getForce(sys, i) * dt;
            dr[i] = (dv[i] + sys[i].getVel()) * dt;
        }

        updateStats(sys, dr, dv);
    }

    inline void stepLeapFrog(std::vector<Body>& sys, double dt)
    {
        std::vector<std::valarray<double>> dr(sys.size(), std::valarray<double>(3));
        std::vector<std::valarray<double>> dv(sys.size(), std::valarray<double>(3));

        std::vector<Body> tmpsys = copy(sys, dt);
        for (size_t i = 0; i < sys.size(); ++i)
        {
            std::valarray<double> a = getForce(sys, i);
            dr[i] = sys[i].getVel() * dt + 0.5 * a * dt * dt;
            dv[i] = 0.5 * (a + getForce(tmpsys, i)) * dt;
        }

        updateStats(sys, dr, dv);
    }

    inline void stepRKN89(std::vector<Body>& sys, double dt)
    {
        // function obtained from https://ntrs.nasa.gov/api/citations/19730015887/downloads/19730015887.pdf
        double sqrt15 = std::sqrt(15);

        std::vector<double> alpha = {
            1. / 3,
            2. / 3,
            1. / 2,
            1. / 3,
            1.0,
            1. / 9,
            1. / 2,
            1. / 3,
            2. / 3,
            0.1 * (5 - sqrt15),
            0.1 * (5 + sqrt15),
            1.0,
            1.0,
        };

        // generate sparce matrix
        std::vector<std::valarray<double>> gamma(13);
        for (int i = 0; i < 13; ++i)
        {
            gamma[i] = std::valarray<double>(i + 1);
        }
        gamma[0] = {
            1. / 18
        };
        gamma[1] = {
            2. / 27,
            4. / 27
        };
        gamma[2] = {
            7. / 128,
            5. / 64,
            -1. / 128
        };
        gamma[3] = {
            89. / 3240,
            31. / 540,
            11. / 1080,
            -16. / 405
        };
        gamma[4] = {
            11. / 120,
            0.0,
            9. / 40,
            -4. / 15,
            9. / 20
        };
        gamma[5] = {
            33259. / 7085880,
            0.0,
            343. / 157464,
            -4708. / 885735,
            1879. / 393660,
            -139. / 885735
        };
        gamma[6] = {
            29. / 1920,
            0.0,
            0.0,
            0.0,
            99. / 2560,
            1. / 30720,
            729. / 10240
        };
        gamma[7] = {
            13. / 1215,
            0.0,
            0.0,
            0.0,
            1. / 144,
            1. / 77760,
            87. / 2240,
            -8. / 8505
        };
        gamma[8] = {
            22. / 1215,
            0.0,
            0.0,
            0.0,
            0.0,
            1. / 4860,
            3. / 28,
            256. / 8505,
            1. / 15
        };
        gamma[9] = {
            1. / 420000 * (7561 - 1454 * sqrt15),
            0.0,
            0.0,
            0.0,
            -9. / 800000 * (1373 + 45 * sqrt15),
            1. / 1344000 * (397 - 145 * sqrt15),
            729. / 78400000 * (6997 - 1791 * sqrt15),
            1. / 183750 * (999 - 473 * sqrt15),
            27. / 5600000 * (19407 - 3865 * sqrt15),
            297. / 700000 * (78 - 19 * sqrt15)
        };
        gamma[10] = {
            1. / 840000 * (12647 + 2413 * sqrt15),
            0.0,
            0.0,
            0.0,
            -9. / 800000 * (1373 - 45 * sqrt15),
            -1. / 336000 * (29 - 61 * sqrt15),
            729. / 19600000 * (14743 + 3789 * sqrt15),
            1. / 183750 * (999 + 143 * sqrt15),
            27. / 5600000 * (20157 + 4315 * sqrt15),
            27. / 1400000 * (1641 + 463 * sqrt15),
            -1. / 56 * (27 + 7 * sqrt15)
        };
        gamma[11] = { 9. / 280 - 35. / 6561 * 5. / 2,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            5. / 2,
            16. / 315 - 160. / 19683 * 5. / 2,
            243. / 1540 + 35. / 2673 * 5. / 2,
            243. / 3080 + 7. / 2673 * 5. / 2,
            25. / 1386 * (5 + sqrt15) - 3500. / 216513 * (31 + 8 * sqrt15) * 5. / 2,
            25. / 1386 * (5 - sqrt15) - 3500. / 216513 * (31 - 8 * sqrt15) * 5. / 2
        };
        gamma[12] = {
            9. / 280,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            16. / 315,
            243. / 1540,
            243. / 3080,
            25. / 1386 * (5 + sqrt15),
            25. / 1386 * (5 - sqrt15),
            0.0
        };

        std::valarray<double> cdot = {
            9. / 280,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            32. / 315,
            729. / 3080,
            729. / 3080,
            125. / 693,
            125. / 693,
            9. / 280
        };

        std::vector<std::vector<std::valarray<double>>> f(sys.size(), std::vector<std::valarray<double>>(13, std::valarray<double>(3)));
        std::vector<std::valarray<double>> dr(sys.size(), std::valarray<double>(3));
        std::vector<std::valarray<double>> dv(sys.size(), std::valarray<double>(3));
        for (size_t k = 0; k < 13; ++k)
        {
            std::vector<Body> tmpsys = copy2(sys, dr);
            for (size_t i = 0; i < sys.size(); ++i)
            {
                f[i][k] = getForce(tmpsys, i);
                // splitting dr like this gives +10% performance
                dr[i] = gamma[k][0] * f[i][0];
                for (size_t j = 1; j < k + 1; ++j)
                {
                    dr[i] += gamma[k][j] * f[i][j];
                }
                dr[i] *= dt * dt;
                dr[i] += dt * alpha[k] * sys[i].getVel();
            }
        }

        for (size_t i = 0; i < sys.size(); ++i)
        {
            for (size_t j = 0; j < 13; ++j)
            {
                dv[i] += cdot[j] * f[i][j];
            }
            dv[i] *= dt;
        }

        updateStats(sys, dr, dv);
    }

} // namespace hyx

#endif // BODY_H
