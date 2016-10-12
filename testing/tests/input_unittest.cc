// Testing use of grvy for reading input from a file 
//
// Author: dsondak@ices.utexas.edu (David Sondak)


// Don't forget gtest.h, which declares the testing framework.
//#include <limits.h>
#include "grvy.h"
#include "read_input.h"
#include "gtest/gtest.h"


// See gtest.h for a complete list of test macros.
//
// You should not use underscore (_) in the test names.

// Tests using grvy to input

// Tests reading in integer from input file.
TEST(ReadInputTest, Integer) {
  EXPECT_EQ(3, ReadInt());

}

// Tests reading in double from input file.
TEST(ReadInputTest, Doubles) {
  EXPECT_EQ(17.0, ReadDouble());

}

// Tests reading in double from input file.
TEST(ReadInputTest, Strings) {
  EXPECT_EQ("hello", ReadString());

}

// Step 3. Call RUN_ALL_TESTS() in main().
//
// We do this by linking in src/gtest_main.cc file, which consists of
// a main() function which calls RUN_ALL_TESTS() for us.
//
// This runs all the tests you've defined, prints the result, and
// returns 0 if successful, or 1 otherwise.
//
// Did you notice that we didn't register the tests?  The
// RUN_ALL_TESTS() macro magically knows about all the tests we
// defined.  Isn't this convenient?
