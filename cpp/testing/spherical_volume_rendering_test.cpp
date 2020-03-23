#include "gtest/gtest.h"
#include "../spherical_volume_rendering_util.h"

// Utilizes the Google Test suite.
// To run:
//       1. Clone the google test suite to spherical-volume-rendering/cpp/testing
//          - It can be found at: https://github.com/google/googletest.git
//          - Currently uses the default folder name (googletest) in the CMake file, so there is no need to change this.
//       2. Load (or re-load) CMakeLists.txt
//
// For information on Google Test, see: https://github.com/google/googletest/blob/master/googletest/README.md
// For examples of Google Test, see: https://github.com/google/googletest/tree/master/googletest/samples

TEST(SampleTest, SimpleMacros) {
    EXPECT_EQ(1, 1);
    EXPECT_TRUE(1 == 1);
    EXPECT_FALSE(1 == 0);
    EXPECT_GT(1, 0);

    const double tolerance = 0.1;
    EXPECT_NEAR(1.0, 1.01, tolerance);
}
