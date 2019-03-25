#ifndef INTEGRATOR_H
#define INTEGRATOR_H

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
        DistrY.setGridX().setBorders(a, b);
        return *this;
    }

    integral & setGrid( distribution const &NewDistrX ) {
        DistrY.setGridX(NewDistrX);
        isSetGridX = true;
        return *this;
    }

    integral & setTollerance( double T ) {
        Tollerance = T;
        return *this;
    }

    integral & setGridStep( size_t Step ) {
        DistrY.setGridX().setBorders(Step).eval();

        return *this;
    }

    /* eval integral */
    virtual double operator()( void ) = 0;
};

class trapezium_integral : public integral {
private:
    uint Fragmentation;

public:
    trapezium_integral( void ) {}

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

            res *= (DistrY.getB() - DistrY.getA()) / DistrY.getNodeCount();

            return res;
        };


        setGridStep(Fragmentation);
        integralWithStepx2 = eval();

        for(;fabs(integralWithStepx2 - integralWithStep) * 1 / 3.0 > Tollerance; Fragmentation <<= 1) {

            integralWithStep = integralWithStepx2;
            setGridStep(Fragmentation << 1);
            integralWithStepx2 = eval();
        }

        return integralWithStep;
    }

    uint getFragmentation( void ) const {
        return Fragmentation;
    }
};

#endif // INTEGRATOR_H
