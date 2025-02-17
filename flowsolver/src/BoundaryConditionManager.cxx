#include "BoundaryConditionManager.hxx"
#include "RCR.hxx"
#include "ControlledCoronary.hxx"
#include "NetlistBoundaryCondition.hxx"
#include "ImpedanceBoundaryCondition.hxx"
#include "FortranBoundaryDataPointerManager.hxx"
#include "fileWriters.hxx"
#include "fileIOHelpers.hxx"
#include "../../estimation/src/CrimsonGlobalArrayTransfer.h"

// This file contains (and should continue to contain) all the tools needed to control the boundary conditions.
//
// This includes functions which can be called from Fortran, and should be the sole point of interface between Fortran and C++
// for the boundary conditions, as far as is possible.
//
// One thing that might be though of as an exception to this rule is the FortranBoundaryDataPointerManager class, but this
// is really a way of setting up the link between Fortran and C++, not so much a way of allowing Fortran to control the BCs.
//
// Another exception is that some of the global data is accessed directly by the C++ classes, such as in the constructor
// for the AbstractBoundaryCondition. This is not ideal, and should be phased out slowly so that we have fewer points
// of interface between the two languages.

// Static class static member variables:
BoundaryConditionManager* BoundaryConditionManager::instance = 0;
HistFileReader* BoundaryConditionManager::PHistReader = NULL;
bool BoundaryConditionManager::m_thisIsARestartedSimulation = false;

// Functions which affect features of the abstract class:
void BoundaryConditionManager::setNumberOfRCRSurfaces(const int numGRCRSrfs)
{
  assert(m_NumberOfRCRSurfaces == 0);
  m_NumberOfRCRSurfaces = numGRCRSrfs;
  m_numberOfBoundaryConditionsManaged += m_NumberOfRCRSurfaces;
  CrimsonGlobalArrayTransfer::Get()->initialiseForRCRFiltering(numGRCRSrfs);
}

void BoundaryConditionManager::setNumberOfControlledCoronarySurfaces(const int numControlledCoronarySrfs)
{
  assert(m_NumberOfControlledCoronarySurfaces == 0);
  m_NumberOfControlledCoronarySurfaces = numControlledCoronarySrfs;
  m_numberOfBoundaryConditionsManaged += numControlledCoronarySrfs;
}

void BoundaryConditionManager::setNumberOfNetlistSurfaces(const int numNetlistLPNSrfs)
{
  assert(m_NumberOfNetlistSurfaces == 0);
  m_NumberOfNetlistSurfaces = numNetlistLPNSrfs;
  m_numberOfBoundaryConditionsManaged += m_NumberOfNetlistSurfaces;
}

void BoundaryConditionManager::setNumberOfImpedanceSurfaces(const int numImpedanceSurfaces)
{
  assert(m_NumberOfImpedanceSurfaces == 0);
  m_NumberOfImpedanceSurfaces = numImpedanceSurfaces;
  m_numberOfBoundaryConditionsManaged += m_NumberOfImpedanceSurfaces;
}

void BoundaryConditionManager::setMasterControlScriptPresent(const int masterControlScriptPresent)
{
  if (masterControlScriptPresent == 1)
  {
    m_masterControlScriptPresent = true;
  }
    else
  {
    m_masterControlScriptPresent = false;
  }
}

void BoundaryConditionManager::setDelt(const double delt)
{
  m_delt = delt;
  m_deltHasBeenSet = true;
}

void BoundaryConditionManager::setHstep(const int hstep)
{
  m_hstep = hstep;
  m_hstepHasBeenSet = true;
}

void BoundaryConditionManager::setAlfi(const double alfi)
{
  m_alfi = alfi;
  m_alfiHasBeenSet = true;
}

void BoundaryConditionManager::setSimulationModePurelyZeroD(const int simulationIsPurelyZeroD)
{
  if (simulationIsPurelyZeroD == 1)
  {
    m_simulationIsPurelyZeroD = true;
  }
  else
  {
    m_simulationIsPurelyZeroD = false;
  }
  
}

// void BoundaryConditionManager::setLstep(const int currentTimestepIndex)
// {
//   m_currentTimestepIndex = currentTimestepIndex;
//   m_currentTimestepIndexHasBeenSet = true;
// }

void BoundaryConditionManager::setStartingTimestepIndex(const int startingTimestepIndex)
{
  assert(!m_startingTimestepIndexHasBeenSet);
  m_startingTimestepIndex = startingTimestepIndex;
  m_currentTimestepIndex = startingTimestepIndex;
  m_startingTimestepIndexHasBeenSet = true;
} 

void BoundaryConditionManager::incrementTimestepIndex()
{
  // assert(m_currentTimestepIndexHasBeenSet);

  // increment the internal timestep of the manager
  m_currentTimestepIndex++;

  // Tell all the bounddary conditions to increment too
  for (auto boundaryCondition = m_boundaryConditions.begin(); boundaryCondition != m_boundaryConditions.end(); boundaryCondition++)
  {
    (*boundaryCondition)->incrementTimestepIndex();
  }
}

void BoundaryConditionManager::setNtout(const int ntout)
{
  m_ntout = ntout;
  m_ntoutHasBeenSet = true;
}

void BoundaryConditionManager::setMaxsurf(const int maxsurf)
{
  m_maxsurf = maxsurf;
  m_maxsurfHasBeenSet = true;
}

void BoundaryConditionManager::setNstep(const int nstep)
{
  m_nstep = nstep;
  m_nstepHasBeenSet = true;
}

void BoundaryConditionManager::setNumLoopClosingnetlistCircuits(const int numLoopClosingCircuits)
{
  m_numLoopClosingNetlistCircuits = numLoopClosingCircuits;
  m_numLoopClosingNetlistCircuitsHasBeenSet = true;
}

void BoundaryConditionManager::checkIfThisIsARestartedSimulation()
{
  SimpleFileReader numstartReader("numstart.dat");

  bool successfullyReadNumstartDotDat = false;
  std::string numstartString = numstartReader.getNextDataSplitBySpacesOrEndOfLine(successfullyReadNumstartDotDat);
  assert(successfullyReadNumstartDotDat);

  int valueFromNumstartDotDat = boost::lexical_cast<int>(numstartString);

  setStartingTimestepIndex(valueFromNumstartDotDat);

  if (valueFromNumstartDotDat > 0)
  {
    m_thisIsARestartedSimulation = true;
    m_nextTimestepWrite_netlistBoundaries_start = valueFromNumstartDotDat + 1;
  }
  else
  {
    m_thisIsARestartedSimulation = false;
    m_nextTimestepWrite_netlistBoundaries_start = 0;
  }
}


void BoundaryConditionManager::giveBoundaryConditionsListsOfTheirAssociatedMeshNodes(const int* ndsurf_nodeToBoundaryAssociationArray, const int& lengthOfNodeToBoundaryAssociationArray)
{
  for (auto boundaryCondition=m_boundaryConditions.begin(); boundaryCondition!=m_boundaryConditions.end(); boundaryCondition++)
  {
    (*boundaryCondition)->setListOfMeshNodesAtThisBoundary(ndsurf_nodeToBoundaryAssociationArray, lengthOfNodeToBoundaryAssociationArray);
  }
}
// ---WRAPPED BY--->
extern "C" void callCPPGiveBoundaryConditionsListsOfTheirAssociatedMeshNodes(const int*& ndsurf_nodeToBoundaryAssociationArray, const int& lengthOfNodeToBoundaryAssociationArray)
{
  BoundaryConditionManager* boundaryConditionManager_instance = BoundaryConditionManager::Instance();
  boundaryConditionManager_instance->giveBoundaryConditionsListsOfTheirAssociatedMeshNodes(ndsurf_nodeToBoundaryAssociationArray, lengthOfNodeToBoundaryAssociationArray);
}


// RCR Boundary condition specific functions
void BoundaryConditionManager::setPressureFromFortran()
{
  // see the called funciton setPressureFromFortran comments for details of what this does.
  for (auto boundaryCondition = m_boundaryConditions.begin(); boundaryCondition!=m_boundaryConditions.end(); boundaryCondition++)
  {
    if (typeid(**boundaryCondition)==typeid(RCR))
    {
      // std::cout << "setting pressure for C++ RCRs" << std::endl;
      RCR* downcastRCR = dynamic_cast<RCR*>(boundaryCondition->get());
      downcastRCR->setPressureFromFortran();
    }
  }
}
// ---WRAPPED BY--->
extern "C" void callCPPSetPressureFromFortran()
{
  BoundaryConditionManager* boundaryConditionManager_instance = BoundaryConditionManager::Instance();
  boundaryConditionManager_instance->setPressureFromFortran();
}

void BoundaryConditionManager::getImplicitCoeff_rcr(double* const implicitCoeffs_toBeFilled)
{
  // This code is a bit tricky, becase FORTRAN/C++ interfacing doesn't yet support passing arrays which are sized
  // at run-time to C++ from FORTRAN. Therefore, I've had to just pass a pointer to the first entry, and then manage
  // dereferencing of that pointer manually to fill the whole array, but with the FORTRAN column-major array structure,
  // as opposed to the C++ row-major standard.
  int writeLocation = 0;
  
  for(auto iterator=m_boundaryConditions.begin(); iterator!=m_boundaryConditions.end(); iterator++)
  {
    if (typeid(**iterator)==typeid(RCR))
    {
      // std::cout << "RCR implicoeff: " << (*iterator)->getSurfaceIndex() << " " << (*iterator)->getdp_dq() << " " << (*iterator)->getHop() << std::endl;     
      implicitCoeffs_toBeFilled[writeLocation] = (*iterator)->getdp_dq();
      // +m_maxsurf+1 here to move to the next column of the array (the +1 is annoying, and is because of weird design decisions in old FORTRAN code)
      
      implicitCoeffs_toBeFilled[writeLocation + m_maxsurf + 1] = (*iterator)->getHop();
      
      writeLocation++;
    }
  }
  
}
// ---WRAPPED BY--->
extern "C" void callCppGetImplicitCoeff_rcr(double*& implicitCoeffs_toBeFilled)
{
  BoundaryConditionManager* boundaryConditionManager_instance = BoundaryConditionManager::Instance();
  boundaryConditionManager_instance->getImplicitCoeff_rcr(implicitCoeffs_toBeFilled);
}

void BoundaryConditionManager::updateAllRCRS_Pressure_n1_withflow()
{
  for(auto iterator=m_boundaryConditions.begin(); iterator!=m_boundaryConditions.end(); iterator++)
  {
    if (typeid(**iterator)==typeid(RCR) || typeid(**iterator)==typeid(NetlistBoundaryCondition))
    {
      (*iterator)->updpressure_n1_withflow();
    }
  }
}
// ---WRAPPED BY--->
extern "C" void callCPPUpdateAllRCRS_Pressure_n1_withflow()
{
  BoundaryConditionManager* boundaryConditionManager_instance = BoundaryConditionManager::Instance();
  boundaryConditionManager_instance->updateAllRCRS_Pressure_n1_withflow();
}

// void BoundaryConditionManager::storeAllBoundaryConditionFlowsAndPressuresAtStartOfTimestep()
// {
//   for (auto boundaryCondition = m_boundaryConditions.begin(); boundaryCondition != m_boundaryConditions.end(); boundaryCondition++)
//   {
//     (*boundaryCondition)->storeFlowAndPressureAtStartOfTimestep();
//   }
// }

void BoundaryConditionManager::setSurfaceList(const std::vector<std::pair<int,boundary_condition_t>> surfaceList)
{
  // Defensive:
  assert(m_deltHasBeenSet);
  assert(m_hstepHasBeenSet);
  assert(m_alfiHasBeenSet);
  // assert(m_currentTimestepIndexHasBeenSet);
  assert(m_ntoutHasBeenSet);
  assert(m_maxsurfHasBeenSet);
  assert(m_nstepHasBeenSet);
  assert(m_numLoopClosingNetlistCircuitsHasBeenSet);
  assert(m_startingTimestepIndexHasBeenSet);

  assert(!m_hasSurfaceList);
  m_hasSurfaceList = true;

  // Build a factory
  BoundaryConditionFactory factory(m_hstep, m_delt, m_alfi, m_maxsurf, m_nstep, m_numLoopClosingNetlistCircuits, m_simulationIsPurelyZeroD, m_startingTimestepIndex);

  factory.createNetlistLoopClosingCircuits(m_netlistDownstreamLoopClosingSubsections);

  for (auto iterator = surfaceList.begin(); iterator != surfaceList.end(); iterator++)
  {
    m_boundaryConditions.push_back(factory.createBoundaryCondition(iterator->first,iterator->second));
  }
}

void BoundaryConditionManager::markClosedLoopLinearSystemsForRebuilding()
{
  // Only do this if this simulation is using a closed loop:
  for (auto downstreamLoopClosingSubsection = m_netlistDownstreamLoopClosingSubsections.begin(); downstreamLoopClosingSubsection != m_netlistDownstreamLoopClosingSubsections.end(); downstreamLoopClosingSubsection++)
  {
    (*downstreamLoopClosingSubsection)->markLinearSystemAsNeedingBuildingAgain();
    (*downstreamLoopClosingSubsection)->markLinearSystemAsNeedingUpdatingAgain();
  }
}

void BoundaryConditionManager::setZeroDDomainReplacementPressuresAndFlows(double* zeroDDomainPressures, double* zeroDDomainFlows)
{
  for (auto boundaryCondition = m_boundaryConditions.begin(); boundaryCondition != m_boundaryConditions.end(); boundaryCondition++)
  {
    boost::shared_ptr<NetlistBoundaryCondition> downcastNetlist = boost::dynamic_pointer_cast<NetlistBoundaryCondition> (*boundaryCondition);
    if (downcastNetlist == NULL) {
      throw std::runtime_error("Can only use Netlist boundary conditions in pure zero-D simulations.");
    }

    // Get the (zero-indexed) Netlist index; this gives us the appropriate location of the pointers in the input variables
    int netlistIndex = downcastNetlist->getIndexAmongstNetlists();
    // Give the appropriate memory addresses of the pressures and flows to this NetlistBC:
    downcastNetlist->setPressureAndFlowPointers(&zeroDDomainPressures[netlistIndex], &zeroDDomainFlows[netlistIndex]);

    downcastNetlist->initialiseModel();
  }
}

void BoundaryConditionManager::ifRestartingLoadNecessaryData()
{
  if (m_thisIsARestartedSimulation)
  {
    // Load PHistRCR.dat, necessary for setting the pressure data in the 
    // LPN at the boundary when restarting
    if (m_NumberOfRCRSurfaces > 0)
    {
      PHistReader = new HistFileReader();
      PHistReader->setFileName("PHistRCR.dat");
      PHistReader->setNumColumns(m_NumberOfRCRSurfaces+1);
      PHistReader->readAndSplitMultiSurfaceRestartFile();
    }
  }
}

std::vector<boost::shared_ptr<AbstractBoundaryCondition>>* BoundaryConditionManager::getBoundaryConditions()
{
    return &m_boundaryConditions;
}

// FULLY DEFINED IN HEADER SO OTHER TRANSLATION UNITS CAN USE IT
// template <typename TemplateBoundaryConditionType>
// void BoundaryConditionManager::computeImplicitCoeff_solve(const int timestepNumber)
// {
//   for (auto&& boundaryCondition : m_boundaryConditions)
//   {
//     if (boost::dynamic_pointer_cast<TemplateBoundaryConditionType> (boundaryCondition))
//     {
//       boundaryCondition->computeImplicitCoeff_solve(timestepNumber);
//     }
//   }
// }
// ---WRAPPED BY--->
extern "C" void callCppComputeAllImplicitCoeff_solve(int& timestepNumber)
{
  BoundaryConditionManager* boundaryConditionManager_instance = BoundaryConditionManager::Instance();
  boundaryConditionManager_instance->computeImplicitCoeff_solve<AbstractBoundaryCondition>(timestepNumber);
}
// ---AND--->
extern "C" void callCppComputeAllNetlistImplicitCoeff_solve(int& timestepNumber)
{
  BoundaryConditionManager* boundaryConditionManager_instance = BoundaryConditionManager::Instance();
  boundaryConditionManager_instance->computeImplicitCoeff_solve<NetlistBoundaryCondition>(timestepNumber);
}
// ---AND--->
extern "C" void callCppComputeAllCoronaryImplicitCoeff_solve(int& timestepNumber)
{
  BoundaryConditionManager* boundaryConditionManager_instance = BoundaryConditionManager::Instance();
  boundaryConditionManager_instance->computeImplicitCoeff_solve<ControlledCoronary>(timestepNumber);
}
// ---AND--->
extern "C" void callCppComputeAllNumericalRCRImplicitCoeff_solve(int& timestepNumber)
{
  BoundaryConditionManager* boundaryConditionManager_instance = BoundaryConditionManager::Instance();
  boundaryConditionManager_instance->computeImplicitCoeff_solve<RCR>(timestepNumber);
}
// ---AND--->
extern "C" void callCppComputeAllImpedanceImplicitCoeff_solve(int& timestepNumber)
{
  BoundaryConditionManager* boundaryConditionManager_instance = BoundaryConditionManager::Instance();
  boundaryConditionManager_instance->computeImplicitCoeff_solve<ImpedanceBoundaryCondition>(timestepNumber);
}

// FULLY DEFINED IN HEADER SO OTHER TRANSLATION UNITS CAN USE IT
// template <typename TemplateBoundaryConditionType>
// void BoundaryConditionManager::computeImplicitCoeff_update(const int timestepNumber)
// {
//   for (auto&& boundaryCondition : m_boundaryConditions)
//   {
//     if (boost::dynamic_pointer_cast<TemplateBoundaryConditionType> (boundaryCondition))
//     {
//       boundaryCondition->computeImplicitCoeff_update(timestepNumber);
//     }
//   }
// }
// ---WRAPPED BY--->
extern "C" void callCppComputeAllImplicitCoeff_update(int& timestepNumber)
{
  BoundaryConditionManager* boundaryConditionManager_instance = BoundaryConditionManager::Instance();
  boundaryConditionManager_instance->computeImplicitCoeff_update<AbstractBoundaryCondition>(timestepNumber);
}
// ---AND--->
extern "C" void callCppComputeAllNetlistImplicitCoeff_update(int& timestepNumber)
{
  BoundaryConditionManager* boundaryConditionManager_instance = BoundaryConditionManager::Instance();
  boundaryConditionManager_instance->computeImplicitCoeff_update<NetlistBoundaryCondition>(timestepNumber);
}
// ---AND--->
extern "C" void callCppComputeAllCoronaryImplicitCoeff_update(int& timestepNumber)
{
  BoundaryConditionManager* boundaryConditionManager_instance = BoundaryConditionManager::Instance();
  boundaryConditionManager_instance->computeImplicitCoeff_update<ControlledCoronary>(timestepNumber);
}
// ---AND--->
extern "C" void callCppComputeAllNumericalRCRImplicitCoeff_update(int& timestepNumber)
{
  BoundaryConditionManager* boundaryConditionManager_instance = BoundaryConditionManager::Instance();
  boundaryConditionManager_instance->computeImplicitCoeff_update<RCR>(timestepNumber);
}
// ---AND--->
extern "C" void callCppComputeAllImpedanceImplicitCoeff_update(int& timestepNumber)
{
  BoundaryConditionManager* boundaryConditionManager_instance = BoundaryConditionManager::Instance();
  boundaryConditionManager_instance->computeImplicitCoeff_update<ImpedanceBoundaryCondition>(timestepNumber);
}

void BoundaryConditionManager::updateAllRCRS_setflow_n(const double* const flows)
{
  int readLocation = 0;
  for(auto iterator=m_boundaryConditions.begin(); iterator!=m_boundaryConditions.end(); iterator++)
  {
    if (typeid(**iterator)==typeid(RCR))
    {
      (*iterator)->setFlowN(flows[readLocation]);
      readLocation++;
    }
  }
}
// ---WRAPPED BY--->
extern "C" void callCPPUpdateAllRCRS_setflow_n(double*& flows)
{
  BoundaryConditionManager* boundaryConditionManager_instance = BoundaryConditionManager::Instance();
  boundaryConditionManager_instance->updateAllRCRS_setflow_n(flows); 
}


void BoundaryConditionManager::updateAllRCRS_setflow_n1(const double* const flows)
{
  int readLocation = 0;
  for(auto iterator=m_boundaryConditions.begin(); iterator!=m_boundaryConditions.end(); iterator++)
  {
    if (typeid(**iterator)==typeid(RCR))
    {
      (*iterator)->setFlowN1(flows[readLocation]);
      readLocation++;
    }
  }
}
// ---WRAPPED BY--->
extern "C" void callCPPUpdateAllRCRS_setflow_n1(double*& flows)
{
  BoundaryConditionManager* boundaryConditionManager_instance = BoundaryConditionManager::Instance();
  boundaryConditionManager_instance->updateAllRCRS_setflow_n1(flows); 
}

void BoundaryConditionManager::recordPressuresAndFlowsInHistoryArrays()
{
  for(auto iterator=m_boundaryConditions.begin(); iterator!=m_boundaryConditions.end(); iterator++)
  {
    (*iterator)->updatePressureAndFlowHistory();
  }
}
// ---WRAPPED BY--->
extern "C" void callCPPRecordPressuresAndFlowsInHistoryArrays()
{
  BoundaryConditionManager* boundaryConditionManager_instance = BoundaryConditionManager::Instance();
  boundaryConditionManager_instance->recordPressuresAndFlowsInHistoryArrays();
}

void BoundaryConditionManager::writePHistAndQHistRCR()
{
  // Open a file writer to append to Phist
  BasicFileWriter phistrcr_writer;
  phistrcr_writer.setFileName("PHistRCR.dat");

  // Open a file writer to append to Qhist
  BasicFileWriter qhistrcr_writer;
  qhistrcr_writer.setFileName("QHistRCR.dat");

  // Loop over all the updates since the last restart was written:
  for (int timestepToWrite = m_currentTimestepIndex - m_ntout + 1; timestepToWrite < m_currentTimestepIndex + 1; timestepToWrite++)
  {
    phistrcr_writer.writeStepIndex(timestepToWrite);
    qhistrcr_writer.writeStepIndex(timestepToWrite);

    // Loop the boundary conditions looking for the RCRs
    for(auto iterator=m_boundaryConditions.begin(); iterator!=m_boundaryConditions.end(); iterator++)
    {
      if (typeid(**iterator)==typeid(RCR))
      {
        // Write the pressure and flow for this timestep (indexed timestepToWrite)
        phistrcr_writer.writeToFile((*iterator)->getPressureHistoryValueByTimestepIndex(timestepToWrite));
        qhistrcr_writer.writeToFile((*iterator)->getFlowHistoryValueByTimestepIndex(timestepToWrite));
      }
    }
    phistrcr_writer.writeEndLine();
    qhistrcr_writer.writeEndLine();
  }
}
// ---WRAPPED BY--->
extern "C" void callCPPWritePHistAndQHistRCR()
{
  BoundaryConditionManager* boundaryConditionManager_instance = BoundaryConditionManager::Instance();
  boundaryConditionManager_instance->writePHistAndQHistRCR();
}


// =========== Controlled Coronary Block ===========

// void BoundaryConditionManager::setSurfacePressure_controlledCoronary(double* coronarySurfacePressures)
// {
//   int readLocation = int(0);
//   for(auto iterator=m_boundaryConditions.begin(); iterator!=m_boundaryConditions.end(); iterator++)
//   {
//     if (typeid(**iterator)==typeid(ControlledCoronary))
//     {
//      (*iterator)->setLPNInflowPressure(coronarySurfacePressures[readLocation]);
//      readLocation++;
//     }
//   }
// }
// // ---WRAPPED BY--->
// extern "C" void callCppSetSurfacePressure_controlledCoronary(double*& coronarySurfacePressures)
// {
//   BoundaryConditionManager* boundaryConditionManager_instance = BoundaryConditionManager::Instance();
//   boundaryConditionManager_instance->setSurfacePressure_controlledCoronary(coronarySurfacePressures);
// }

void BoundaryConditionManager::getImplicitCoeff_controlledCoronary(double* const implicitCoeffs_toBeFilled)
{
  // This code is a bit tricky, becase FORTRAN/C++ interfacing doesn't yet support passing arrays which are sized
  // at run-time to C++ from FORTRAN. Therefore, I've had to just pass a pointer to the first entry, and then manage
  // dereferencing of that pointer manually to fill the whole array, but with the FORTRAN column-major array structure,
  // as opposed to the C++ row-major standard.
  int writeLocation = 0;
  
  for(auto iterator=m_boundaryConditions.begin(); iterator!=m_boundaryConditions.end(); iterator++)
  {
    if (typeid(**iterator)==typeid(ControlledCoronary))
    {
      
      implicitCoeffs_toBeFilled[writeLocation] = (*iterator)->getdp_dq();
      // +MAXSURF+1 here to move to the next column of the array (the +1 is annoying, and is because of weird design decisions in old FORTRAN code)
      
      implicitCoeffs_toBeFilled[writeLocation+m_maxsurf+1] = (*iterator)->getHop();
      
      writeLocation++;
    }
  }
}
// ---WRAPPED BY--->
extern "C" void callCppGetImplicitCoeff_controlledCoronary(double*& implicitCoeffs_toBeFilled) 
{
  BoundaryConditionManager* boundaryConditionManager_instance = BoundaryConditionManager::Instance();
  boundaryConditionManager_instance->getImplicitCoeff_controlledCoronary(implicitCoeffs_toBeFilled);
}

void BoundaryConditionManager::updateAllControlledCoronaryLPNs()
{
  for(auto boundaryCondition=m_boundaryConditions.begin(); boundaryCondition!=m_boundaryConditions.end(); boundaryCondition++)
  {
    boost::shared_ptr<ControlledCoronary> downcastCoronary = boost::dynamic_pointer_cast<ControlledCoronary> (*boundaryCondition);
    if (downcastCoronary != NULL)
    {
      downcastCoronary->updateLPN();
    }
  }
}
// ---WRAPPED BY--->
extern "C" void callCppUpdateAllControlledCoronaryLPNs()
{
  BoundaryConditionManager* boundaryConditionManager_instance = BoundaryConditionManager::Instance();
  boundaryConditionManager_instance->updateAllControlledCoronaryLPNs();
}


void BoundaryConditionManager::finalizeLPNAtEndOfTimestep_controlledCoronary()
{
  for(auto boundaryCondition=m_boundaryConditions.begin(); boundaryCondition!=m_boundaryConditions.end(); boundaryCondition++)
  {
    boost::shared_ptr<ControlledCoronary> downcastCoronary = boost::dynamic_pointer_cast<ControlledCoronary> (*boundaryCondition);
    if (downcastCoronary != NULL)
    {
      downcastCoronary->finaliseAtEndOfTimestep();
    }
  }
}
// ---WRAPPED BY--->
extern "C" void callCppfinalizeLPNAtEndOfTimestep_controlledCoronary()
{
  BoundaryConditionManager* boundaryConditionManager_instance = BoundaryConditionManager::Instance();
  boundaryConditionManager_instance->finalizeLPNAtEndOfTimestep_controlledCoronary();
}

void BoundaryConditionManager::getImplicitCoeff_impedanceBoundaryConditions(double* const implicitCoeffs_toBeFilled)
{
  // This code is a bit tricky, becase FORTRAN/C++ interfacing doesn't yet support passing arrays which are sized
  // at run-time to C++ from FORTRAN. Therefore, I've had to just pass a pointer to the first entry, and then manage
  // dereferencing of that pointer manually to fill the whole array, but with the FORTRAN column-major array structure,
  // as opposed to the C++ row-major standard.
  int writeLocation = 0;

  for (auto& boundaryCondition : m_boundaryConditions)
  {
    boost::shared_ptr<ImpedanceBoundaryCondition> downcastImpedanceBC = boost::dynamic_pointer_cast<ImpedanceBoundaryCondition> (boundaryCondition);
    if (downcastImpedanceBC)
    {
      implicitCoeffs_toBeFilled[writeLocation] = downcastImpedanceBC->getdp_dq();
      implicitCoeffs_toBeFilled[writeLocation+m_maxsurf+1] = downcastImpedanceBC->getHop();
      writeLocation++;
    }
  }
}
// ---WRAPPED BY--->
extern "C" void callCPPGetImplicitCoeff_impedanceBoundaryConditions(double*& implicitCoeffs_toBeFilled)
{
  BoundaryConditionManager* boundaryConditionManager_instance = BoundaryConditionManager::Instance();
  boundaryConditionManager_instance->getImplicitCoeff_impedanceBoundaryConditions(implicitCoeffs_toBeFilled);
}

void BoundaryConditionManager::finalizeLPNAtEndOfTimestep_netlists()
{
  for(auto boundaryCondition=m_boundaryConditions.begin(); boundaryCondition!=m_boundaryConditions.end(); boundaryCondition++)
  {
    auto downcastNetlist = boost::dynamic_pointer_cast<NetlistBoundaryCondition> (*boundaryCondition);
    if (downcastNetlist != NULL)
    {
      downcastNetlist->finaliseAtEndOfTimestep();
    }
  }

  // Do the closed loop downstream subsections:
  for (auto downstreamLoopClosingSubsection = m_netlistDownstreamLoopClosingSubsections.begin(); downstreamLoopClosingSubsection != m_netlistDownstreamLoopClosingSubsections.end(); downstreamLoopClosingSubsection++)
  {
    (*downstreamLoopClosingSubsection)->finalizeLPNAtEndOfTimestep();
  }
}
// ---WRAPPED BY--->
extern "C" void callCppfinalizeLPNAtEndOfTimestep_netlists()
{
  BoundaryConditionManager* boundaryConditionManager_instance = BoundaryConditionManager::Instance();
  boundaryConditionManager_instance->finalizeLPNAtEndOfTimestep_netlists();
}

std::vector<double*> BoundaryConditionManager::getPointersToAllNetlistCapacitorNodalHistoryPressures() const
{
  std::vector<double*> capacitorNodalHistoryPressuresPointers;
  for (auto boundaryCondition : m_boundaryConditions)
  {
    auto downcastNetlist = boost::dynamic_pointer_cast<NetlistBoundaryCondition> (boundaryCondition);
    if (downcastNetlist != NULL)
    {
      std::vector<double*> nodalHistoryPressuresForThisNetlist = downcastNetlist->getCapacitorNodalHistoryPressurePointers();

      // prepend the values for this Netlist to the vector of values from all Netlists that will be returned to the caller:
      capacitorNodalHistoryPressuresPointers.insert(capacitorNodalHistoryPressuresPointers.begin(), nodalHistoryPressuresForThisNetlist.begin(), nodalHistoryPressuresForThisNetlist.end());
    }
  }
  return capacitorNodalHistoryPressuresPointers;
}


// void BoundaryConditionManager::updateAllControlledCoronaryLPNs_Pressure_n1_withflow()
// {
//   for(auto iterator=m_boundaryConditions.begin(); iterator!=m_boundaryConditions.end(); iterator++)
//   {
//     if (typeid(**iterator)==typeid(ControlledCoronary))
//     {
//       (*iterator)->updpressure_n1_withflow();
//     }
//   }
// }
// // ---WRAPPED BY--->
// extern "C" void callCPPUpdateAllControlledCoronaryLPNs_Pressure_n1_withflow()
// {
//   BoundaryConditionManager* boundaryConditionManager_instance = BoundaryConditionManager::Instance();
//   boundaryConditionManager_instance->updateAllRCRS_Pressure_n1_withflow();
// }

// ========== Controlled Coronary Block End =========

// ========== Netlist LPN Block Start =========
void BoundaryConditionManager::initialiseLPNAtStartOfTimestep_netlist()
{
  for (auto boundaryCondition=m_boundaryConditions.begin(); boundaryCondition!=m_boundaryConditions.end(); boundaryCondition++)
  {
    auto downcastNetlist = boost::dynamic_pointer_cast<NetlistBoundaryCondition> (*boundaryCondition);
    if (downcastNetlist != NULL)
    {
      downcastNetlist->initialiseAtStartOfTimestep();
    }
  }

  // Now initialise any closed loop downstream subsections for this timestep:
  for (auto downstreamCircuit = m_netlistDownstreamLoopClosingSubsections.begin(); downstreamCircuit != m_netlistDownstreamLoopClosingSubsections.end(); downstreamCircuit++)
  {
    (*downstreamCircuit)->initialiseAtStartOfTimestep();
  }
}
// ---WRAPPED BY--->
extern "C" void callCPPInitialiseLPNAtStartOfTimestep_netlist()
{
  BoundaryConditionManager* boundaryConditionManager_instance = BoundaryConditionManager::Instance();
  boundaryConditionManager_instance->initialiseLPNAtStartOfTimestep_netlist();
}


void BoundaryConditionManager::updateAllNetlistLPNs(const int timestepNumber)
{
  for(auto boundaryCondition=m_boundaryConditions.begin(); boundaryCondition!=m_boundaryConditions.end(); boundaryCondition++)
  {
    auto downcastNetlist = boost::dynamic_pointer_cast<NetlistBoundaryCondition> (*boundaryCondition);
    if (downcastNetlist != NULL)
    {
      downcastNetlist->updateLPN(timestepNumber);
    }
  }
}
// ---WRAPPED BY--->
extern "C" void callCPPUpdateAllNetlistLPNs(int& timestepNumber)
{
  BoundaryConditionManager* boundaryConditionManager_instance = BoundaryConditionManager::Instance();
  boundaryConditionManager_instance->updateAllNetlistLPNs(timestepNumber);
}

std::map<int,std::pair<double,double>> BoundaryConditionManager::getImplicitCoeff_netlistLPNs_toPassTo3DDomainReplacement()
{
  std::map<int,std::pair<double,double>> allNetlistImplicitCoefficients;
  
  for(auto iterator=m_boundaryConditions.begin(); iterator!=m_boundaryConditions.end(); iterator++)
  {
    if (typeid(**iterator)==typeid(NetlistBoundaryCondition))
    {
      std::pair<double,double> thisNetlistsImplicitCoefficients;
      
      thisNetlistsImplicitCoefficients.first = (*iterator)->getdp_dq();
      thisNetlistsImplicitCoefficients.second = (*iterator)->getHop();

      boost::shared_ptr<NetlistBoundaryCondition> downcastNetlist = boost::dynamic_pointer_cast<NetlistBoundaryCondition> (*iterator);

      allNetlistImplicitCoefficients.insert(std::make_pair(downcastNetlist->getIndexAmongstNetlists(), thisNetlistsImplicitCoefficients));
    }
  }

  return allNetlistImplicitCoefficients;
}

void BoundaryConditionManager::getImplicitCoeff_netlistLPNs(double* const implicitCoeffs_toBeFilled)
{
  // This code is a bit tricky, becase FORTRAN/C++ interfacing doesn't yet support passing arrays which are sized
  // at run-time to C++ from FORTRAN. Therefore, I've had to just pass a pointer to the first entry, and then manage
  // dereferencing of that pointer manually to fill the whole array, but with the FORTRAN column-major array structure,
  // as opposed to the C++ row-major standard.
  int writeLocation = 0;
  
  for(auto iterator=m_boundaryConditions.begin(); iterator!=m_boundaryConditions.end(); iterator++)
  {
    if (typeid(**iterator)==typeid(NetlistBoundaryCondition))
    {
      
      implicitCoeffs_toBeFilled[writeLocation] = (*iterator)->getdp_dq();
      // std::cout << "just got implicit dp_dq: " << implicitCoeffs_toBeFilled[writeLocation];

      // +m_maxsurf+1 here to move to the next column of the array (the +1 is annoying, and is because of weird design decisions in old FORTRAN code)
      implicitCoeffs_toBeFilled[writeLocation+m_maxsurf+1] = (*iterator)->getHop();
      // std::cout << " and H operator: " << implicitCoeffs_toBeFilled[writeLocation+m_maxsurf+1] << std::endl;
      // std::cout << "Netlist implicoeff: " << (*iterator)->getSurfaceIndex() << " " << (*iterator)->getdp_dq() << " " << (*iterator)->getHop() << std::endl;
      writeLocation++;
    }
  }
}
// ---WRAPPED BY--->
extern "C" void callCPPGetImplicitCoeff_netlistLPNs(double*& implicitCoeffs_toBeFilled) 
{
  BoundaryConditionManager* boundaryConditionManager_instance = BoundaryConditionManager::Instance();
  boundaryConditionManager_instance->getImplicitCoeff_netlistLPNs(implicitCoeffs_toBeFilled);
}

// The purpose of this function is to detect when flow across the 3D interface is blocked due
// to diode closure (or due to flows being prescribed -- still \todo).
// It provides the Fortran code with an array of zeros and ones, one for each boundary node in the mesh;
// a 1 indicates that the boundary codes should be left as-is (meaning flow is disallowed),
// whereas a 0 indicates that it should be Neumann.
void BoundaryConditionManager::getBinaryMaskToAdjustNodalBoundaryConditions(int* const binaryMask, const int binaryMaskLength)
{
  // Begin by setting the binary mask to all ones (i.e. flagging 1 to set Dirichlet boundary conditions at all nodes)
  // ... we will set zeros where we want Neumann conditions in a moment...
  for (int maskLocation=0; maskLocation<binaryMaskLength; maskLocation++)
  {
    binaryMask[maskLocation] = 1;
  }
  // Ask the boundary conditions to set zeros where they want Neumann conditions
  for (auto iterator=m_boundaryConditions.begin(); iterator!=m_boundaryConditions.end(); iterator++)
  {
      (*iterator)->setDirichletConditionsIfNecessary(binaryMask);
  }
}
// ---WRAPPED BY--->
extern "C" void callCPPGetBinaryMaskToAdjustNodalBoundaryConditions(int*& binaryMask, const int& binaryMaskLength)
{
  BoundaryConditionManager* boundaryConditionManager_instance = BoundaryConditionManager::Instance();
  boundaryConditionManager_instance->getBinaryMaskToAdjustNodalBoundaryConditions(binaryMask, binaryMaskLength);
}

void BoundaryConditionManager::getNumberOfBoundaryConditionsWhichCurrentlyDisallowFlow(int& numBCsWhichDisallowFlow)
{
  // Ensure we start from zero, before we count the surfaces which disallow flow due to closed valves (so we have to switch to Dirichlet)
  numBCsWhichDisallowFlow = 0;
  // ...do the counting (currently only netlists have valves...)
  for (auto boundaryCondition=m_boundaryConditions.begin(); boundaryCondition!=m_boundaryConditions.end(); boundaryCondition++)
  {
    if (typeid(**boundaryCondition)==typeid(NetlistBoundaryCondition))
    {
      NetlistBoundaryCondition* downcastNetlist = dynamic_cast<NetlistBoundaryCondition*>(boundaryCondition->get());
      if (!downcastNetlist->flowPermittedAcross3DInterface())
      {
        numBCsWhichDisallowFlow++;
      }
    }
  }
}
// ---WRAPPED BY--->
extern "C" void callCPPGetNumberOfBoundaryConditionsWhichCurrentlyDisallowFlow(int& numBCsWhichDisallowFlow)
{
  BoundaryConditionManager* boundaryConditionManager_instance = BoundaryConditionManager::Instance();
  boundaryConditionManager_instance->getNumberOfBoundaryConditionsWhichCurrentlyDisallowFlow(numBCsWhichDisallowFlow);
}

void BoundaryConditionManager::getNumberOfBoundaryConditionManagerBoundaryConditions_reference(int& totalNumberOfManagedBoundaryConditions) const
{
  totalNumberOfManagedBoundaryConditions = m_numberOfBoundaryConditionsManaged;
}
// ---WRAPPED BY--->
extern "C" void callCPPGetNumberOfCppManagedBoundaryConditions(int& totalNumberOfManagedBoundaryConditions)
{
  BoundaryConditionManager* boundaryConditionManager_instance = BoundaryConditionManager::Instance();
  boundaryConditionManager_instance->getNumberOfBoundaryConditionManagerBoundaryConditions_reference(totalNumberOfManagedBoundaryConditions);
}

void BoundaryConditionManager::getNumberOfNetlistBoundaryConditionsWhichCurrentlyAllowFlow(int& numBCsWhichAllowFlow)
{
  // Ensure we start from zero, before we count the surfaces which allow flow due to closed valves (so we have to switch to Dirichlet)
  numBCsWhichAllowFlow = 0;
  // ...do the counting (currently only netlists have valves...)
  for (auto boundaryCondition=m_boundaryConditions.begin(); boundaryCondition!=m_boundaryConditions.end(); boundaryCondition++)
  {
    if (typeid(**boundaryCondition)==typeid(NetlistBoundaryCondition))
    {
      NetlistBoundaryCondition* downcastNetlist = dynamic_cast<NetlistBoundaryCondition*>(boundaryCondition->get());
      if (downcastNetlist->flowPermittedAcross3DInterface())
      {
        numBCsWhichAllowFlow++;
      }
    }
  }
}
// ---WRAPPED BY--->
extern "C" void callCPPGetNumberOfNetlistsWhichCurrentlyAllowFlow(int& numBCsWhichAllowFlow)
{
  BoundaryConditionManager* boundaryConditionManager_instance = BoundaryConditionManager::Instance();
  boundaryConditionManager_instance->getNumberOfNetlistBoundaryConditionsWhichCurrentlyAllowFlow(numBCsWhichAllowFlow);
}

// Takes a surface index and a reference to an int - if flow is allowed across this surface,
// returns 1 in the referenced int, otherwise returns zero in that int.
void BoundaryConditionManager::discoverWhetherFlowPermittedAcrossSurface(const int& queriedSurfaceIndex, int& flowIsPermitted)
{
  // Begin by assuming flow is permitted; this will be changed below if flow is not permitted.
  flowIsPermitted = 1;
  // find the queried surface:
  for (auto boundaryCondition=m_boundaryConditions.begin(); boundaryCondition!=m_boundaryConditions.end(); boundaryCondition++)
  {
    bool thisIsANetlist = typeid(**boundaryCondition)==typeid(NetlistBoundaryCondition);
    if ((*boundaryCondition)->getSurfaceIndex() == queriedSurfaceIndex && thisIsANetlist)
    {
      // Discover whether we should report that flow is permitted or not:
      NetlistBoundaryCondition* downcastNetlist = static_cast<NetlistBoundaryCondition*>(boundaryCondition->get());
      if (!(downcastNetlist->flowPermittedAcross3DInterface()))
      {
        flowIsPermitted = 0;
      }
    }
  }
}
//---WRAPPED BY--->
extern "C" void callCPPDiscoverWhetherFlowPermittedAcrossSurface(const int& queriedSurfaceIndex, int& flowIsPermitted)
{
  BoundaryConditionManager* boundaryConditionManager_instance = BoundaryConditionManager::Instance();
  boundaryConditionManager_instance->discoverWhetherFlowPermittedAcrossSurface(queriedSurfaceIndex, flowIsPermitted);
}

void BoundaryConditionManager::haveBoundaryConditionTypesChanged(int& boundaryConditionTypesHaveChanged)
{
  // Warning: This thing is broken: the booleans it checks for are not reset correctly currently. Do not use without fixing first!
  assert(false);
  // Begin by assuming boundary conditions are as they were on the previous time-step; this will be changed below if the assumption is false.
  boundaryConditionTypesHaveChanged = 0;
  // find netlists
  for (auto boundaryCondition=m_boundaryConditions.begin(); boundaryCondition!=m_boundaryConditions.end(); boundaryCondition++)
  {
    if (typeid(**boundaryCondition)==typeid(NetlistBoundaryCondition))
    {
      // Discover whether we should report a change in boundary condition type (Neumann/Dirichlet):
      NetlistBoundaryCondition* downcastNetlist = dynamic_cast<NetlistBoundaryCondition*>(boundaryCondition->get());
      if (downcastNetlist->boundaryConditionTypeHasJustChanged())
      {
        boundaryConditionTypesHaveChanged = 1;
      }
    }
  }
}
//---WRAPPED BY--->
extern "C" void callCPPHaveBoundaryConditionTypesChanged(int& boundaryConditionTypesHaveChanged)
{
  BoundaryConditionManager* boundaryConditionManager_instance = BoundaryConditionManager::Instance();
  boundaryConditionManager_instance->haveBoundaryConditionTypesChanged(boundaryConditionTypesHaveChanged);
}

// void BoundaryConditionManager::setSurfacePressure_netlistLPNs(double* netlistSurfacePressures)
// {
//   int readLocation = int(0);
//   for(auto iterator=m_boundaryConditions.begin(); iterator!=m_boundaryConditions.end(); iterator++)
//   {
//     if (typeid(**iterator)==typeid(NetlistBoundaryCondition))
//     {
//      (*iterator)->setLPNInflowPressure(netlistSurfacePressures[readLocation]);
//      readLocation++;
//     }
//   }
// }
// // ---WRAPPED BY--->
// extern "C" void callCppSetSurfacePressure_netlistLPNs(double*& netlistSurfacePressures)
// {
//   BoundaryConditionManager* boundaryConditionManager_instance = BoundaryConditionManager::Instance();
//   boundaryConditionManager_instance->setSurfacePressure_netlistLPNs(netlistSurfacePressures);
// }

void BoundaryConditionManager::writeAllNetlistComponentFlowsAndNodalPressures()
{
  writeNetlistFlowsPressuresAndVolumes(m_boundaryConditions, m_netlistDownstreamLoopClosingSubsections, m_nextTimestepWrite_netlistBoundaries_start);
}
// ---WRAPPED BY--->
extern "C" void callCPPWriteAllNetlistComponentFlowsAndNodalPressures()
{
  BoundaryConditionManager* boundaryConditionManager_instance = BoundaryConditionManager::Instance();
  boundaryConditionManager_instance->writeAllNetlistComponentFlowsAndNodalPressures();
}

// void BoundaryConditionManager::loadAllNetlistComponentFlowsAndNodalPressures()
// {
//   assert(m_startingTimestepIndexHasBeenSet);
//   loadNetlistPressuresFlowsAndVolumesOnRestart(m_boundaryConditions, m_netlistDownstreamLoopClosingSubsections, m_startingTimestepIndex);
// }
// // ---WRAPPED BY--->
// extern "C" void callCPPLoadAllNetlistComponentFlowsAndNodalPressures()
// {
//   BoundaryConditionManager* boundaryConditionManager_instance = BoundaryConditionManager::Instance();
//   boundaryConditionManager_instance->loadAllNetlistComponentFlowsAndNodalPressures();
// }

// Control systems specific functions
void BoundaryConditionManager::updateBoundaryConditionControlSystems()
{
  if (m_controlSystemsPresent)
  {
    mp_controlSystemsManager->updateBoundaryConditionControlSystems();
  }
}
// ---WRAPPED BY--->
extern "C" void callCPPUpdateBoundaryConditionControlSystems()
{
  BoundaryConditionManager* boundaryConditionManager_instance = BoundaryConditionManager::Instance();
  boundaryConditionManager_instance->updateBoundaryConditionControlSystems();
}

void BoundaryConditionManager::resetStateUsingKalmanFilteredEstimate(const double flow, const double pressure, const int surfaceIndex, const int timestepNumber)
{
  for (auto boundaryCondition = m_boundaryConditions.begin(); boundaryCondition!=m_boundaryConditions.end(); boundaryCondition++)
  {
    if ((*boundaryCondition)->getSurfaceIndex() == surfaceIndex)
    {
      (*boundaryCondition)->resetStateUsingKalmanFilteredEstimate(flow, pressure, timestepNumber);
    }
  }
}
// ---WRAPPED BY--->
extern "C" void callCPPResetStateUsingKalmanFilteredEstimate(double& flow, double& pressure, int& surfaceIndex, int& timestepNumber)
{
  BoundaryConditionManager* boundaryConditionManager_instance = BoundaryConditionManager::Instance();
  boundaryConditionManager_instance->resetStateUsingKalmanFilteredEstimate(flow, pressure, surfaceIndex, timestepNumber);
}

void BoundaryConditionManager::debugPrintFlowPointerTarget_BCM()
{
  for (auto const &boundaryCondition : m_boundaryConditions) {
    boundaryCondition->debugPrintFlowPointerTarget();
  }  
}
// ---WRAPPED BY--->
extern "C" void callCPPDebugPrintFlowPointerTarget_BCM()
{
  BoundaryConditionManager* boundaryConditionManager_instance = BoundaryConditionManager::Instance();
  boundaryConditionManager_instance->debugPrintFlowPointerTarget_BCM();
}

void BoundaryConditionManager::createControlSystems()
{
  assert(m_startingTimestepIndexHasBeenSet);
  assert(m_ntoutHasBeenSet);
  m_controlSystemsPresent = true;
  // Instantiate the manager
  mp_controlSystemsManager = boost::shared_ptr<ControlSystemsManager>(new ControlSystemsManager(m_delt, m_masterControlScriptPresent, m_startingTimestepIndex, m_ntout));
  
  // Get the reader class for the netlist data file, and ask it for the control description data:
  NetlistXmlReader* netlistXmlReader_instance = NetlistXmlReader::Instance();
  
  // Get info for the components that need control (number of these, the component indices in the netlist, and the control types for each)
  // std::vector<int> numberOfComponentsWithControl = getNumberOfComponentsWithControl();
  std::map<int, std::map<int, ComponentControlSpecificationContainer>> mapsOfComponentControlTypes = netlistXmlReader_instance->getMapsOfComponentControlTypesForEachSurface();

  // Get info for the nodes that need control (number of these, the nodes indices in the netlist, and the control types for each)
  // std::vector<int> numberOfNodesWithControl = getNumberOfNodesWithControl();
  std::map<int, std::map<int,parameter_controller_t>> mapsOfNodeControlTypes = netlistXmlReader_instance->getMapsOfNodalControlTypesForEachSurface();


  // Check for the existence of netlists with input data setting up control of any of 
  // the components. If any are found, initialise the control appropriately.
  for (auto boundaryCondition = m_boundaryConditions.begin(); boundaryCondition!=m_boundaryConditions.end(); boundaryCondition++)
  {
    if (typeid(**boundaryCondition) == typeid(NetlistBoundaryCondition))
    {
      // Downcast to a shared_ptr to a Netlist:
      boost::shared_ptr<NetlistBoundaryCondition> currentNetlist = boost::dynamic_pointer_cast<NetlistBoundaryCondition> (*boundaryCondition);
      boost::shared_ptr<NetlistCircuit> currentNetlistCircuit = currentNetlist->getNetlistCircuit();
      // We now initialise all the controls which affect this netlist...
      int netlistIndex = currentNetlist->getIndexAmongstNetlists();
      // Create the controls for components by looping over the pairs which give the component index in the netlist, together with its prescribed control type from netlist_surfaces.dat:
      try {
        for (auto componentIndexAndControlTypes = mapsOfComponentControlTypes.at(netlistIndex).begin(); componentIndexAndControlTypes != mapsOfComponentControlTypes.at(netlistIndex).end(); componentIndexAndControlTypes++)
          {
            // The component may have multiple controllers attached (e.g. both unstressed volume and compliance), so we loop their names:
            for (int controlSpecificationIndex = 0; controlSpecificationIndex < componentIndexAndControlTypes->second.getNumberOfControlScripts(); controlSpecificationIndex++)
            {
              parameter_controller_t controlType = componentIndexAndControlTypes->second.getControlTypeByIndexLocalToComponent(controlSpecificationIndex);
              mp_controlSystemsManager->createParameterController(controlType, currentNetlistCircuit, componentIndexAndControlTypes->first);
            }
          }
          // Create the controls for nodes by looping over the pairs which give the component index in the netlist, together with its prescribed control type from netlist_surfaces.dat:
          for (auto nodeIndexAndControlType = mapsOfNodeControlTypes.at(netlistIndex).begin(); nodeIndexAndControlType != mapsOfNodeControlTypes.at(netlistIndex).end(); nodeIndexAndControlType++)
          {
            mp_controlSystemsManager->createParameterController(nodeIndexAndControlType->second, currentNetlistCircuit, nodeIndexAndControlType->first);
          }
      } catch (const std::exception& e) {
          std::cout << e.what() << " observed at line " << __LINE__ << " of " << __FILE__ << std::endl;
          throw;
      }
    }
  }

  // If there's a closed loop present:
  if (m_numLoopClosingNetlistCircuits > 0)
  {
    // Get the reader for the closed loop system
    // NetlistDownstreamCircuitReader* downstreamNetlistReader_instance = NetlistDownstreamCircuitReader::Instance();
    NetlistDownstreamXmlReader* downstreamNetlistReader_instance = NetlistDownstreamXmlReader::Instance();

    // Get info for the components that need control (number of these, the component indices in the netlist, and the control types for each)
    // std::vector<int> numberOfComponentsWithControl = getNumberOfComponentsWithControl();
    // std::vector<std::map<int,parameter_controller_t>> mapsOfComponentControlTypes_closedLoop = downstreamNetlistReader_instance->getMapsOfComponentControlTypesForEachSurface();
    std::map<int, std::map<int, ComponentControlSpecificationContainer>> mapsOfComponentControlTypes_closedLoop = downstreamNetlistReader_instance->getMapsOfComponentControlTypesForEachSurface();

    // Get info for the nodes that need control (number of these, the nodes indices in the netlist, and the control types for each)
    // std::vector<int> numberOfNodesWithControl = getNumberOfNodesWithControl();
    // std::vector<std::map<int,parameter_controller_t>> mapsOfNodeControlTypes_closedLoop = downstreamNetlistReader_instance->getMapsOfNodalControlTypesForEachSurface();
    std::map<int, std::map<int,parameter_controller_t>> mapsOfNodeControlTypes_closedLoop = downstreamNetlistReader_instance->getMapsOfNodalControlTypesForEachSurface();

    for (auto loopClosingCircuit = m_netlistDownstreamLoopClosingSubsections.begin(); loopClosingCircuit != m_netlistDownstreamLoopClosingSubsections.end(); loopClosingCircuit++)
    {
      const int closedLoopIndex = (*loopClosingCircuit)->getIndexOfClosedLoop_zeroIndexed();
      boost::shared_ptr<NetlistCircuit> currentNetlistCircuit = (*loopClosingCircuit)->getNetlistCircuit();

      try {
        // Create the controls for components by looping over the pairs which give the component index in the netlist, together with its prescribed control type from netlist_surfaces.dat:
          for (auto componentIndexAndControlTypes = mapsOfComponentControlTypes_closedLoop.at(closedLoopIndex).begin(); componentIndexAndControlTypes != mapsOfComponentControlTypes_closedLoop.at(closedLoopIndex).end(); componentIndexAndControlTypes++)
          {
            // The component may have multiple controllers attached (e.g. both unstressed volume and compliance), so we loop their names:
            for (int controlSpecificationIndex = 0; controlSpecificationIndex < componentIndexAndControlTypes->second.getNumberOfControlScripts(); controlSpecificationIndex++)
            {
              parameter_controller_t controlType = componentIndexAndControlTypes->second.getControlTypeByIndexLocalToComponent(controlSpecificationIndex);
              mp_controlSystemsManager->createParameterController(controlType, currentNetlistCircuit, componentIndexAndControlTypes->first);
            }
          }
          // Create the controls for nodes by looping over the pairs which give the component index in the netlist, together with its prescribed control type from netlist_surfaces.dat:
          for (auto nodeIndexAndControlType = mapsOfNodeControlTypes_closedLoop.at(closedLoopIndex).begin(); nodeIndexAndControlType != mapsOfNodeControlTypes_closedLoop.at(closedLoopIndex).end(); nodeIndexAndControlType++)
          {
            mp_controlSystemsManager->createParameterController(nodeIndexAndControlType->second, currentNetlistCircuit, nodeIndexAndControlType->first);
          }
      } catch (const std::exception& e) {
          std::cout << e.what() << " observed at line " << __LINE__ << " of " << __FILE__ << std::endl;
          throw;
      }
    }

  }

}

std::vector<std::pair<boundary_data_t,double>> BoundaryConditionManager::getBoundaryPressuresOrFlows_zeroDDomainReplacement()
{
  std::vector<std::pair<boundary_data_t,double>> pressuresOrFlowsAsAppropriate;
  for (auto boundaryCondition = m_boundaryConditions.begin(); boundaryCondition != m_boundaryConditions.end(); boundaryCondition++)
  {
    if (typeid(**boundaryCondition) == typeid(NetlistBoundaryCondition))
    {
      boost::shared_ptr<NetlistBoundaryCondition> downcastNetlist = boost::dynamic_pointer_cast<NetlistBoundaryCondition> (*boundaryCondition);
      pressuresOrFlowsAsAppropriate.push_back(downcastNetlist->computeAndGetFlowOrPressureToGiveToZeroDDomainReplacement());
    }
    else
    {
      std::stringstream errorMessage;
      errorMessage << "EE: You can only use a zero-D replacement for the 3D domain if all the boundary conditions are Netlists." << std::endl;
      throw std::runtime_error(errorMessage.str());
    }
  }
  assert(pressuresOrFlowsAsAppropriate.size() == m_NumberOfNetlistSurfaces);
  return pressuresOrFlowsAsAppropriate;
}

int BoundaryConditionManager::getNumberOfControlSystems() const
{
  return mp_controlSystemsManager->getNumberOfControlSystems();
}

// Put all the end-of-timestep cleanup work here
// (some of it is still located in Fortran at the time of writing 2016-11-22, but much of it can probably be moved here to reduce linking with Fortran)
void BoundaryConditionManager::finaliseOnTimeStep()
{
  for (auto& boundaryCondition : m_boundaryConditions){
    boost::shared_ptr<ImpedanceBoundaryCondition> downcastImpedanceBC = boost::dynamic_pointer_cast<ImpedanceBoundaryCondition> (boundaryCondition);
    if (downcastImpedanceBC)
    {
      downcastImpedanceBC->finaliseAtEndOfTimestep();
    }
  }
}