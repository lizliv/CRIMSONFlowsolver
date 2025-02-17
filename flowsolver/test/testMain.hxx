#ifndef TESTMAIN_HXX_
#define TESTMAIN_HXX_

#include "gtest/gtest.h"

#include "mpi.h"

#include <iostream>
#include <stdio.h>
#include "common_c.h"
#include "input.h"
#include "proces.h"
#include "itrdrv.h"
#include "itrPC.h"
#include "partition.h"
#include "input_fform.h"
#include "multidom.hxx"
#include "FortranBoundaryDataPointerManager.hxx"
#include "fileReaders.hxx"
#include "CrimsonGlobalArrayTransfer.h"

#include "debuggingToolsForCpp.hxx"
#include <boost/filesystem/path.hpp>

#ifdef intel
#include <direct.h>
#else
#include <unistd.h>
#endif



// The fixture for testing class Foo.
	class testMain : public ::testing::Test {
	 public:
		int rank;
        int numProcsTotal;
        MPI_Comm iNewComm_C;
	 protected:
	  // You can remove any or all of the following functions if its body
	  // is empty.
	 	boost::filesystem::path dirBinaryCalledFrom;

	 	// This constructor should just do exactly what main.cxx does in estimation/src/main.cxx
		testMain() {
			MPI_Barrier(MPI_COMM_WORLD);
			// dirBinaryCalledFrom = get_current_dir_name();
			// dirBinaryCalledFrom = boost::filesystem::current_path()
			dirBinaryCalledFrom = boost::filesystem::current_path();
			getRank();
		}

		void setSimDirectory(std::string dir)
		{
			int success = chdir(dir.c_str());
			if (success==0)
			{
				std::cout << "II: Changed to dir: " << dir << " rank was: " << rank << std::endl;
			}
			else
			{
				std::cout << "EE: Failed to change to dir: " << dir << " rank was: " << rank << " success was: " << success << std::endl;
				perror("EEEE: Failed to change dir in test");
			}
		}

		void clearOutOldFiles()
		{
			MPI_Barrier(iNewComm_C);
			if (rank == 0)
			{
				// Warning - this has the potential to delete multiple folders due to the *
				system("rm -rf *-procs-case");
			}
			MPI_Barrier(iNewComm_C);
		}

		void runSimulation()
		{
		   char pathToProcsCaseDir[100];

		   // Moved this to the gtest_main.cc
		   // MPI_Init(&fake_argc,(char***)&fake_argv);

		   // Moved this to the gtest_main.cc
		   // if(fake_argc > 2 ){
			  //  static volatile int debuggerPresent =0;
			  //  while (!debuggerPresent ); // assign debuggerPresent=1
		   // }

		   // Dont need this during tests:
		   // if ( argc < 2 ) {
		   //    if ( 0 == rank ) {
		   //       std::cout << "Usage: " << fake_argv[0] << "<.inp file> <optional debug flag> \n";
		   //    }
		   //    // return 0;
		   // }

		   // read configuration file
		   int errFlag = input_fform();
		   if (errFlag != 0)
		   {
		   	   throw std::runtime_error("EE: Failed during parsing of input files.");
		   }

		   // Preprocess data and run the problem
		   // Partition the problem to the correct number of processors

		   if( rank == 0 )
		   {
		      //cout << "number of procs " << numprocs_perparticle << endl;
		      Partition_Problem( numProcsTotal );
		   }

		   MPI_Barrier(iNewComm_C);
			
		   sprintf(pathToProcsCaseDir,"%d-procs-case",numProcsTotal);
		   boost::filesystem::path thisDir = boost::filesystem::current_path();
		   int errStat = chdir(pathToProcsCaseDir);

                   if (errStat != 0)
		   {
			std::cerr << "Failed to change to directory " << pathToProcsCaseDir << ". Rank is: " << rank << std::endl;
			perror("EEE");
			throw std::runtime_error("EEEEE");
	    	   }
		   else
		   {
			   std::cout << "changing directory to " << pathToProcsCaseDir << std::endl;			
		   }

		   input(&numProcsTotal, &rank);
		   proces();

		   {
		   	   // just initialise the time values that the AbstractBoundaryCondition needs (when it's called in multidom_initialise).
		   	   // This will be called again during itrdrv_init.
		   	   // This is really ugly, but a proper fix will take days - it's a BIG refactor.
			   int dummyInitialItseqValue=1;
			   callFortranSetupTimeParameters(dummyInitialItseqValue);
		   }

		   // initialise reduced order boundary conditions
		   multidom_initialise();
		   multidomSetupControlSystems();

		   itrdrv_init(); // initialize solver

		   // FortranBoundaryDataPointerManager* pointerManager;
		   // pointerManager = FortranBoundaryDataPointerManager::Get();
		   for (int kk = 1; kk <= inpdat.nstep[0]; kk++) {
		   	   multidom_iter_initialise();
			   itrdrv_iter_init();

			   itrdrv_iter_step();

			   itrdrv_iter_finalize();
			   multidom_iter_finalise();
		   }

		   itrdrv_finalize();
		   multidom_finalise();
	       MPI_Barrier(iNewComm_C);
		   // Moved this to the gtest_main.cc
		   // MPI_Finalize();

		   // return ierr;
		}

	  virtual ~testMain() {
	    // You can do clean-up work that doesn't throw exceptions here.
	    boost::filesystem::current_path(dirBinaryCalledFrom);
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
	    FortranBoundaryDataPointerManager::Get()->tearDown();
	    CrimsonGlobalArrayTransfer::Get()->tearDown();
	    // BoundaryConditionManager::Instance()->Term();
	  }
	private:
		void getRank()
		{
		   // save the communicator
		   iNewComm_C = MPI_COMM_WORLD;
		   newcom.iNewComm = MPI_Comm_c2f(iNewComm_C); // modifies newcom in fortran common block
	           MPI_Barrier(iNewComm_C);

		   MPI_Comm_size(iNewComm_C, &numProcsTotal);
		   MPI_Comm_rank(iNewComm_C, &rank);
		}

	  // Objects declared here can be used by all tests in the test case for Foo.
	};

// Hack to force the compiler to link this test to the relevant main() for testing
	int PullInMyLibraryTestMain();

#endif
