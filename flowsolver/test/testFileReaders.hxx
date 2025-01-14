#ifndef TESTFILEREADERS_HPP_
#define TESTFILEREADERS_HPP_

#include "fileReaders.cxx"
#include "gtest/gtest.h"
#include "CrimsonGlobalArrayTransfer.h"



	// The fixture for testing class Foo.
	class testFileReaders : public ::testing::Test {
	 protected:
	  // You can remove any or all of the following functions if its body
	  // is empty.
	 RcrtReader* rcrtReader_instance;

	  testFileReaders() {
		rcrtReader_instance = RcrtReader::Instance();
		rcrtReader_instance->setFileName("rcrt_test.dat");
		rcrtReader_instance->readAndSplitMultiSurfaceInputFile();
	  }

	  virtual ~testFileReaders() {
	    // You can do clean-up work that doesn't throw exceptions here.
	  }

	  // If the constructor and destructor are not enough for setting up
	  // and cleaning up each test, you can define the following methods:

	  virtual void SetUp() {
	    // Code here will be called immediately after the constructor (right
	    // before each test).
	  }

	  virtual void TearDown() {
	    // Code here will be called immediately after each test (right
	    // before the destructor).
	    rcrtReader_instance->tearDown();
	    CrimsonGlobalArrayTransfer::Get()->tearDown();
	  }

	  // Objects declared here can be used by all tests in the test case for Foo.
	};

// Hack to force the compiler to link this test to the relevant main() for testing
	int PullInMyLibrary();

#endif
	