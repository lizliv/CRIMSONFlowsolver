#ifndef NETLISTCIRCUIT_HXX_
#define NETLISTCIRCUIT_HXX_

#include "gtest/gtest_prod.h"
#include "CircuitData.hxx"
#include <sstream>
#include "petscsys.h"
#include "petscmat.h"
#include "petscvec.h"
#include <set>
#include <boost/weak_ptr.hpp>
#include "fileReaders.hxx"
#include "NetlistXmlReader.hxx"

class NetlistCircuit
{
	friend class testMultidom;
	FRIEND_TEST(testMultidom,checkNetlistComponentNeighbourPointers);
	FRIEND_TEST(testMultidom, checkClosedDiodeWithRemainingOpenPathDetected);
	FRIEND_TEST(testMultidom, checkClosedDiodeWithoutRemainingOpenPathDetected);
public:
	NetlistCircuit(const int hstep, const int surfaceIndex, const int indexOfThisNetlistLPN, const bool thisIsARestartedSimulation, const double alfi, const double delt, const int startingTimestepIndex)
	: m_surfaceIndex(surfaceIndex),
	m_IndexOfThisNetlistLPNInInputFile(indexOfThisNetlistLPN),
	m_hstep(hstep),
	m_thisIsARestartedSimulation(thisIsARestartedSimulation),
	m_delt(delt),
	m_alfi(alfi),
	m_startingTimestepIndex(startingTimestepIndex)
	{
		initialisePetscArrayNames();

		// Fortran will report the wrong flow (zero) on the first timestep
		// of a restarted simulation. In this case, we trust the value
		// of flow that the boundary condition model already knows (it was
		// saved at the point-of-restart during the last simulation).
		//
		// Once this flag is used, we set it to false again (i.e. its one-shot)
		//
		//\todo fix this in the Fortran code instead.
		if (m_startingTimestepIndex > 0)
		{
			m_oneshotIgnoreIncorrectFortranFlow = true;
		}
		else
		{
			m_oneshotIgnoreIncorrectFortranFlow = false;
		}

		safetyCounterLimit = 1000;
		mp_circuitData = boost::shared_ptr<CircuitData> (new CircuitData(m_hstep));
		// mp_circuitDataWithoutDiodes = boost::shared_ptr<CircuitData> (new CircuitData(m_hstep));

		std::stringstream pressureFileNameBuilder;
		pressureFileNameBuilder << "netlistPressures_surface_" << m_surfaceIndex << ".dat";
		m_PressureHistoryFileName = pressureFileNameBuilder.str();

		std::stringstream flowFileNameBuilder;
		flowFileNameBuilder << "netlistFlows_surface_" << m_surfaceIndex << ".dat";
		m_FlowHistoryFileName = flowFileNameBuilder.str();

		std::stringstream volumeFileNameBuilder;
		volumeFileNameBuilder << "netlistVolumes_surface_" << m_surfaceIndex << ".dat";
		m_VolumeHistoryFileName = volumeFileNameBuilder.str();
	}

	virtual void initialiseCircuit();

	bool surfaceIndexMatches(const int surfaceIndexToTest) const;
	int getSurfaceIndex() const;

	bool flowPermittedAcross3DInterface() const;
	bool boundaryConditionTypeHasJustChanged();
	void closeAllDiodes();
	virtual void detectWhetherClosedDiodesStopAllFlowAt3DInterface();
	void switchDiodeStatesIfNecessary();
	void rebuildCircuitMetadata();

	void setPressureAndFlowPointers(double* pressurePointer, double* flowPointer);
	void recordPressuresFlowsAndVolumesInHistoryArrays();

	// void identifyAtomicSubcircuits();
	virtual void initialiseAtStartOfTimestep();
	void finalizeLPNAtEndOfTimestep();
	boost::shared_ptr<CircuitData> getCircuitDescription();
	int getNumberOfDegreesOfFreedom() const;
	std::vector<int> getNodesWithDeferredKirchoffEquations() const;
	int getIndexAmongstNetlists(){return m_IndexOfThisNetlistLPNInInputFile;}
	std::vector<std::pair<int,double*>> getComponentInputDataIndicesAndFlows() const;
	std::vector<std::pair<int,double*>> getNodeInputDataIndicesAndPressures() const;
	std::vector<std::pair<int,double*>> getVolumeTrackingComponentInputDataIndicesAndVolumes() const;
	
	bool hasPrescribedPressureAcross3DInterface() const;
	bool hasPrescribedFlowAcross3DInterface() const;

	void computeHistoryVariablesToMatchCurrentKalmanFilterParticle(const double alfi_delt);

	std::vector<double*> getCapacitorNodalHistoryPressurePointers() const;

	virtual void createCircuitDescription();
	virtual ~NetlistCircuit()
	{
		terminatePetscArrays();
	}

	// This can be used to give more than one pressure and one flow pointer to the netlist. Useful if this Netlist
	// has multiple interfaces with other domains (e.g. if this is a Netlist replacement for the 3D domain.)
	void setPointersToBoundaryPressuresAndFlows(double* const interfacePressures, double* const interfaceFlows, const int& numberOfPointers);

	void writePressuresFlowsAndVolumes(int& nextTimestepWrite_start);
	void loadPressuresFlowsAndVolumesOnRestart();

	virtual std::pair<double,double> computeImplicitCoefficients(const int timestepNumber, const double timeAtStepNplus1, const double alfi_delt);
	virtual void updateLPN(const int timestepNumber);

	virtual std::pair<boundary_data_t,double> computeAndGetFlowOrPressureToGiveToZeroDDomainReplacement();
	boost::shared_ptr<CircuitComponent> getComponentByInputDataIndex(const int componentIndex);
	boost::shared_ptr<CircuitPressureNode> getNodeByInputDataIndex(const int componentIndex);
	int getNumberOfHistoryPressures() const;
	double getInterfaceFlowSign() const;
protected:
	// Overload constructor for subclasses to call:
	NetlistCircuit(const int hstep, const int indexOfThisNetlistLPN, const bool thisIsARestartedSimulation, const double alfi, const double delt, const int startingTimestepIndex)
	: m_surfaceIndex(-1),
	m_IndexOfThisNetlistLPNInInputFile(indexOfThisNetlistLPN),
	m_hstep(hstep),
	m_thisIsARestartedSimulation(thisIsARestartedSimulation),
	m_delt(delt),
	m_alfi(alfi),
	m_startingTimestepIndex(startingTimestepIndex)
	{
		initialisePetscArrayNames();
	}
	std::string m_PressureHistoryFileName;
	std::string m_FlowHistoryFileName;
	std::string m_VolumeHistoryFileName;
	boost::shared_ptr<CircuitData> mp_circuitData;
	const int m_surfaceIndex;
	const int m_IndexOfThisNetlistLPNInInputFile;
	const int m_hstep;
	const bool m_thisIsARestartedSimulation;
	const double m_delt;
	const double m_alfi;
	const int m_startingTimestepIndex;
	std::vector<double*> pressure_n_ptrs;
	std::vector<double*> flow_n_ptrs;
	int m_NumberOfAtomicSubcircuits;

	NetlistReader* mp_netlistFileReader;
	NetlistXmlReader* mp_netlistXmlReader;

	void createBasicCircuitDescription();
	void createVectorsAndMatricesForCircuitLinearSystem();
	void createListOfNodesWithMultipleIncidentCurrents();
	void getMapOfPressHistoriesToCorrectPressNodes();
	void getMapOfFlowHistoriesToCorrectComponents();
	void getMapOfVolumeHistoriesToCorrectComponents();
	void getMapOfTrackedVolumesToCorrectComponents();
	void generateLinearSystemFromPrescribedCircuit(const double alfi_delt);
	void assembleRHS(const bool useHistoryHistoryPressure); // useHistoryHistoryPressure should usually be "false"
	void giveNodesTheirPressuresFromSolutionVector();
	void giveComponentsTheirFlowsFromSolutionVector();
	void giveComponentsTheirVolumesFromSolutionVector();
	void giveComponentsTheirProposedVolumesFromSolutionVector();
	std::vector<double> getVolumesFromSolutionVector();
	bool areThereNegativeVolumes(const double alfi_delt);
	void initialiseCircuit_common();
	void setInternalHistoryPressureFlowsAndVolumes();

	Mat m_systemMatrix;
	Mat m_inverseOfSystemMatrix;
	Mat m_identityMatrixForPetscInversionHack;
	Vec m_RHS;
	Vec m_solutionVector;

	std::vector<double> pressuresInSubcircuit;
	std::vector<double> historyPressuresInSubcircuit; // As pressuresInLPN, but for any nodes with histories. /Most/ of the entries in this array will never be used.
	std::vector<double> flowsInSubcircuit;            // Flow through each component in the LPN, in the order they appear in the netlist
	std::vector<double> historyFlowsInSubcircuit;	  // As flowsInLPN, but for any nodes with histories. /Most/ of the entries in this array will never be used.
	std::vector<double> volumesInSubcircuit;
	std::vector<double> historyVolumesInSubcircuit;
	// circuitData subcircuitInputData;
	std::map<int,int> nodeIndexToPressureHistoryNodeOrderingMap;
	std::map<int,int> componentIndexToFlowHistoryComponentOrderingMap;
	std::map<int,int> componentIndexToVolumeHistoryComponentOrderingMap;
	std::map<int,int> componentIndexToTrackedVolumeComponentOrderingMap;
	PetscInt m_numberOfSystemRows;
	PetscInt m_numberOfSystemColumns;
	std::vector<int> listOfNodesWithMultipleIncidentCurrents;
	int m_numberOfMultipleIncidentCurrentNodes;
	std::set<int> listOfHistoryPressures;            // generated from input data, listing pressure node indices and component flow indices where a history is needed (i.e. last time-step values for capacitors/inductors)
	std::set<int> listOfHistoryFlows;
	std::set<int> listOfHistoryVolumes;
	std::set<int> listOfTrackedVolumes;
	int numberOfPrescribedPressuresAndFlows;           // Just the sum of the previous two declared integers
	int m_numberOfHistoryPressures;
	int numberOfHistoryFlows;
	int numberOfHistoryVolumes;
	int m_numberOfTrackedVolumes;
	int m_numberOfVolumeTrackingPressureChambers;
	std::vector<int> columnMap;
	// int columnMapSize;//\todo check this is used
	std::vector<int> m_columnOf3DInterfacePrescribedPressureInLinearSystem;
	std::vector<int> m_locationOf3DInterfaceComputedPressureInSolutionVector;
	std::vector<int> m_columnOf3DInterfacePrescribedFlowInLinearSystem;
	std::vector<int> m_locationOf3DInterfaceComputedFlowInSolutionVector;

	std::vector<int> m_nodesWithKirchoffEquationsDeferredToClosedLoop;

	int safetyCounterLimit;

	PetscScalar m_interfaceFlow;
  	PetscScalar m_interfacePressure;
	
	void buildAndSolveLinearSystem(const double alfi_delt);
	void generateLinearSystemWithoutFactorisation(const double alfi_delt);

private:
	void initialisePetscArrayNames();
	void terminatePetscArrays();
	virtual void setupPressureNode(const int indexOfEndNodeInInputData, boost::shared_ptr<CircuitPressureNode>& node, boost::shared_ptr<CircuitComponent> component);
	virtual bool kirchoffEquationAtNodeDeferredToInterfacingCircuit(const int nodeIndex) const;
	void findLinearSystemIndicesOf3DInterfacePressureAndFlow();
	void setupCustomPythonControlSystems();
	void countVolumeTrackingPressureChambers();
	void solveLinearSystem();
	void buildAndSolveLinearSystemForUpdatingHistoryVariablesToMatchCurrentKalmanParticle(const double alfi_delt);
	void recordPressureHistory();
	void recordPressureHistoryHistory();
	// void createInitialCircuitDescriptionWithoutDiodes();
	// void assignComponentsToAtomicSubcircuits();

	// boost::shared_ptr<CircuitData> mp_circuitDataWithoutDiodes;
	std::vector<boost::shared_ptr<CircuitData>> m_activeSubcircuitCircuitData;
	std::vector<int> m_AtomicSubcircuitsComponentsBelongsTo; // This is indexed by component, as they appear in mp_circuitDataWithoutDiodes
	bool m_oneshotIgnoreIncorrectFortranFlow;
	std::vector<std::pair<int,double>> m_locationsInRHSForUnstressedVolumesAndTheirValues;

	// std::vector<double> m_PressuresInLPN;                       // Pressure at each LPN node, using the same node indexing as in the netlist
	// std::vector<double> m_HistoryPressuresInLPN;                // As m_PressuresInLPN, but for any nodes with histories. /Most/ of the entries in this array will never be used.
	// std::vector<double> m_FlowsInLPN;                           // Flow through each component in the LPN, in the order they appear in the netlist
	// std::vector<double> m_HistoryFlowsInLPN;					  // As m_FlowsInLPN, but for any nodes with histories. /Most/ of the entries in this array will never be used.

};

#endif