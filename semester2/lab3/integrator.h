#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <memory>
#include "def.h"
#include "distribution.h"

class integral {
protected:
    func Function;
    value_distribution DistrY;
    double Tollerance;

private:
    bool isSetGridX;

public:
    integral( void ) : Function(nullptr), isSetGridX(false) {}

    integral & setFunction( func F ) {
        Function = F;
        DistrY.setFunc(F);
        return *this;
    }

    integral & setBorders( double a, double b ) {
        DistrY.setBorders(a, b);
        return *this;
    }

protected:
    integral & setGrid( shared_ptr<distribution> NewDistrX ) {
        DistrY.setGridX(NewDistrX);
        isSetGridX = true;
        return *this;
    }

    integral & setGridStep( size_t Step ) {
        DistrY.setBorders(Step).eval();

        return *this;
    }

public:
    integral & setTollerance( double T ) {
        Tollerance = T;
        return *this;
    }

    /* eval integral */
    virtual double operator()( void ) { return 0; }
};

class trapezium_integral : public integral {
private:
    uint Fragmentation;
    func f;

public:
    trapezium_integral( void ) {
      setGrid(shared_ptr<uniform>(new uniform()));
    }

    /* eval integral */
    double operator()( void ) {
        double
                integralWithStep = 0,
                integralWithStepx2;
        Fragmentation = 2;

        auto eval = [this] ( void ) -> double {
            double res = 0;

            for (uint i = 1; i < DistrY.getNodeCount() - 1; i++)
                res += DistrY[i];

            res += (DistrY[0] + DistrY[DistrY.getNodeCount() - 1]) / 2.0;

            res *= (DistrY.getB() - DistrY.getA()) / (DistrY.getNodeCount() - 1);

            return res;
        };


        setGridStep(Fragmentation);
        integralWithStepx2 = eval();

        for(;fabs(integralWithStepx2 - integralWithStep) > 3 * Tollerance; Fragmentation <<= 1) {

            integralWithStep = integralWithStepx2;
            setGridStep(Fragmentation << 1);
            integralWithStepx2 = eval();
        }

        return integralWithStepx2;
    }

    uint getFragmentation( void ) const {
        return Fragmentation;
    }
};

class rado_integral_n2 : public integral {
private:
    uint Fragmentation;

public:
    rado_integral_n2( void ) { setGrid(shared_ptr<uniform>(new uniform())); }

    /* eval integral */
    double operator()( void ) {
        double
                B1 = 2 / 9.0,
                x1 = -0.289897948556636,
                x2 = 0.689897948556636,
                A1 = 1.024971652376843,
                A2 = 0.752806125400934;
        double
                integralWithStep = 0,
                integralWithStepx2, prevB = Function(DistrY.X()[0]);
        Fragmentation = 2;

        auto eval = [this, B1, x1, x2, A1, A2, &prevB] ( void ) -> double {
            double res = 0, resA1 = 0, resA2 = 0, resB = prevB;

            for (uint i = 1; i < DistrY.getNodeCount() - 1; i += 2)
                resB += Function(DistrY.X()[i]);

            prevB = resB;
            resB *= B1;

            double h = (DistrY.getB() - DistrY.getA()) / (DistrY.getNodeCount() - 1);

            for (uint i = 1; i < DistrY.getNodeCount(); i++) {
                resA1 += Function(h / 2.0 * x1 + (DistrY.X()[i - 1] + DistrY.X()[i]) / 2.0);
                resA2 += Function(h / 2.0 * x2 + (DistrY.X()[i - 1] + DistrY.X()[i]) / 2.0);
            }

            resA1 *= A1;
            resA2 *= A2;

            res = (resA1 + resA2 + resB) * h / 2.0;

            return res;
        };

        DistrY.setBorders(Fragmentation);
        integralWithStepx2 = eval();

        for(;fabs(integralWithStepx2 - integralWithStep) > 15 * Tollerance;) {

            integralWithStep = integralWithStepx2;
            DistrY.setBorders(Fragmentation  = ((Fragmentation - 1) << 1) + 1);
            integralWithStepx2 = eval();
        }

        return integralWithStepx2;
    }

    uint getFragmentation( void ) const {
        return Fragmentation;
    }
};

#endif // INTEGRATOR_H
