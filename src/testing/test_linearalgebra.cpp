#include <math.h>
#include <stdlib.h>

#include "linearalgebra.h"

#include "gtest/gtest.h"

class LinearalgebraTest: public ::testing::Test {
protected:
    LinearalgebraTest() {
    };

    virtual ~LinearalgebraTest() {
    };
    
    /* Tolerance for comparing values */
    static const double tol_ = 1e-12;
};


TEST_F(LinearalgebraTest, dot) {
    double _A[9] = { 1.0, 0.0, 0.0,  0.0, 1.0, 0.0,  0.0, 0.0, 1.0 };
    double _x[3] = { 1.0, 2.0, 3.0 };
    double _y[3] = { 4.0, 5.0, 6.0 };

    mat<double> A(3, _A);
    vec<double> x(3, _x);
    vec<double> y(3, _y);

    y = A.dot(x);

    ASSERT_NEAR(y[0], x[0], tol_);
    ASSERT_NEAR(y[1], x[1], tol_);
    ASSERT_NEAR(y[2], x[2], tol_);
}
