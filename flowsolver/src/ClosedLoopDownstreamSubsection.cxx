#include "ClosedLoopDownstreamSubsection.hxx"
#include <utility>
#include <stdexcept>
#include "debuggingToolsForCpp.hxx"
#include "indexShifters.hxx"

bool ClosedLoopDownstreamSubsection::boundaryConditionCircuitConnectsToThisDownstreamSubsection(const int boundaryConditionIndex) const
{
    bool returnValue = mp_NetlistCircuit->boundaryConditionCircuitConnectsToThisDownstreamSubsection(boundaryConditionIndex);
	return returnValue;
}

void ClosedLoopDownstreamSubsection::setPointerToNeighbouringBoundaryConditionCircuit(boost::shared_ptr<NetlistCircuit> upstreamBCCircuit)
{
    m_upstreamBoundaryConditionCircuits.push_back(upstreamBCCircuit);
    // m_numberOfUpstreamCircuits gets updated every time we add an upstreamBCCircuit
    // (although its value won't be used - at the time of writing this comment - until all upstreamBCCircuits have been set.)
    m_numberOfUpstreamCircuits = m_upstreamBoundaryConditionCircuits.size();
}

void ClosedLoopDownstreamSubsection::initialiseModel()
{
    // Get the input data
    mp_NetlistCircuit->createCircuitDescription();

    // Determine how many subcircuits are needed, and note which components belong to each subcircuit
    // mp_NetlistCircuit->identifyAtomicSubcircuits();

    // Initialise all diodes to their closed state, for stability
    //\todo change this if you're restarting and the diodes need to be open at restart!
    mp_NetlistCircuit->closeAllDiodes();

    mp_NetlistCircuit->initialiseCircuit();

    // count the diodes, and set up the AtomicSubcircuitConnectionManager, which is used it working out
    // what connections should be made when a diode/valve opens.
    // AtomicSubcircuitConnectionManager* toPassToSharedPtr = new AtomicSubcircuitConnectionManager(mp_CircuitDescription,m_CircuitDataForAtomicSubcircuits);
    //

    generateCircuitInterfaceNodeData();
}

void ClosedLoopDownstreamSubsection::initialiseAtStartOfTimestep()
{
    // Idetify and construct the appropriate subcircuits for this timestep
    mp_NetlistCircuit->initialiseAtStartOfTimestep();
}

void ClosedLoopDownstreamSubsection::giveNodesAndComponentsTheirUpdatedValues()
{
    // This method passes the entries in m_solutionVector which correspond to just the
    // solution of the linear system for the downstream circuit itself to that circuit,
    // and tells it to extract the solution values and give them to the components of
    // the downstream circuit (pressures, flows, volumes)

    // Get the solution entries to give to the downstream circuit
    std::vector<PetscScalar> solutionEntriesForDownstreamCircuit;
    solutionEntriesForDownstreamCircuit = extractContiguousRangeFromPetscVector(m_solutionVector, m_boundsOfDownstreamSolutionDataInSolutionVector.first, m_boundsOfDownstreamSolutionDataInSolutionVector.second);

    mp_NetlistCircuit->giveNodesAndComponentsTheirUpdatedValuesFromSolutionVector(solutionEntriesForDownstreamCircuit);
}

std::vector<PetscScalar> ClosedLoopDownstreamSubsection::extractContiguousRangeFromPetscVector(Vec vector, const int firstEntry, const int lastEntry) const
{
    // Get the raw data to return:
    PetscErrorCode errFlag;
    PetscScalar* rawDataInVector;
    errFlag = VecGetArray(vector, &rawDataInVector); CHKERRABORT(PETSC_COMM_SELF, errFlag);
    
    // Extract the specific range from the solution vector which belongs
    // to the surface with the requested surfaceIndex.
    //
    // Put the raw data in a vector for easy handling for return:
    std::vector<PetscScalar> returnData;
    for (int location = firstEntry; location <= lastEntry; location++)
    {
        returnData.push_back(rawDataInVector[location]);
    }

    errFlag = VecRestoreArray(vector, &rawDataInVector); CHKERRABORT(PETSC_COMM_SELF, errFlag);

    return returnData;
}

std::vector<PetscScalar> ClosedLoopDownstreamSubsection::getSolutionVectorEntriesCorrespondingToSurface(const int surfaceIndex) const
{
    std::pair<int,int> dataLocationRange;
    try {
        dataLocationRange = m_mapOfSurfaceIndicesToRangeOfEntriesInSolutionVector.at(surfaceIndex);
    } catch (const std::exception& e) {
        std::cout << e.what() << " observed at line " << __LINE__ << " of " << __FILE__ << std::endl;
        throw;
    }

    // Get the data to return:
    std::vector<PetscScalar> returnData = extractContiguousRangeFromPetscVector(m_solutionVector, dataLocationRange.first, dataLocationRange.second);

    return returnData;
}

void ClosedLoopDownstreamSubsection::buildAndSolveLinearSystem_internal(const double alfi_delt)
{
    PetscErrorCode errFlag;
    terminatePetscArrays(); // This does nothing if the arrays don't yet exist. \todo consider refactoring the matrix creation/termination.
    clearQueue(m_matrixContributionsFromUpstreamBoundaryConditions);
    clearQueue(m_rhsContributionsFromUpstreamBoundaryConditions);
    m_indicesOf3DInterfaceComputedFlowsInUpstreamSolutionVectors.clear();
    m_columnIndicesOf3DInterfacePrescribedFlowsInUpstreamLinearSystems.clear();
    // std::cout << "just cleared m_columnIndicesOf3DInterfacePrescribedFlowsInUpstreamLinearSystems." << std::endl;
    m_indicesOf3DInterfaceComputedPressuresInUpstreamSolutionVectors.clear();
    m_columnIndicesOf3DInterfacePrescribedPressuresInUpstreamLinearSystems.clear();

    // Call the upstream boundary conditions to ask for their contributions to the (closed loop)-type
    // linear system:
    int firstEntryBelongingToCurrentSurfaceInSolutionVector = 0;
    m_mapOfSurfaceIndicesToRangeOfEntriesInSolutionVector.clear();

    store3DInterfaceFlowSigns();

    m_systemSize = 0;
    for (auto upstreamBCCircuit = m_upstreamBoundaryConditionCircuits.begin(); upstreamBCCircuit != m_upstreamBoundaryConditionCircuits.end(); upstreamBCCircuit++)
    {
        Mat matrixContribution;
        Vec rhsContribution;

        boost::shared_ptr<NetlistBoundaryCircuitWhenDownstreamCircuitsExist> downcastCircuit = boost::dynamic_pointer_cast<NetlistBoundaryCircuitWhenDownstreamCircuitsExist> (*upstreamBCCircuit);
        
        downcastCircuit->getMatrixContribution(alfi_delt, matrixContribution);
        // Get some metadata (the index amongst netlists) with which to tag this circuit's matrixContribution, for convenience later.
        int indexAmongstNetlists = downcastCircuit->getIndexAmongstNetlists();
        // package the just-obtained data and push it onto the queue:
        std::pair<int, Mat> matrixContributionTaggedWithIndexAmongstNetlistBCs = std::make_pair(indexAmongstNetlists, matrixContribution);
        m_matrixContributionsFromUpstreamBoundaryConditions.push(matrixContributionTaggedWithIndexAmongstNetlistBCs);
        
        downcastCircuit->getRHSContribution(rhsContribution);
        m_rhsContributionsFromUpstreamBoundaryConditions.push(rhsContribution);

        // Record which entries in the eventual solution vector will belong to 
        // each upstream circuit:
        int surfaceIndex = downcastCircuit->getSurfaceIndex();
        PetscInt numberOfRows;
        PetscInt numberOfColumns;
        errFlag = MatGetSize(matrixContribution, &numberOfRows, &numberOfColumns); CHKERRABORT(PETSC_COMM_SELF,errFlag);
        std::pair<int,int> solutionVectorRange(firstEntryBelongingToCurrentSurfaceInSolutionVector, firstEntryBelongingToCurrentSurfaceInSolutionVector+toZeroIndexing(numberOfColumns));
        m_mapOfSurfaceIndicesToRangeOfEntriesInSolutionVector.insert(std::make_pair(surfaceIndex, solutionVectorRange));
        firstEntryBelongingToCurrentSurfaceInSolutionVector += numberOfColumns;

        m_systemSize += downcastCircuit->getNumberOfDegreesOfFreedom();
    }

    assert(m_systemSize > 0); // defensive
    
    // Get the final system size by adding in the number of degrees of freedom in the downstream closed loop subsection circuit.
    m_systemSize += mp_NetlistCircuit->getNumberOfDegreesOfFreedom();

    // Get the offsets which tell us where we will find the implicit coefficient conributions in 
    // the solved linear system.
    //
    // The first entry of the pair will be the boundary condition index (NOT the surface index from solver.inp). This
    // indexing is internal to the code, starts at zero, and numbers the netlist surfaces consecuitvely.
    for (auto upstreamBCCircuit = m_upstreamBoundaryConditionCircuits.begin(); upstreamBCCircuit != m_upstreamBoundaryConditionCircuits.end(); upstreamBCCircuit++)
    {
        boost::shared_ptr<NetlistBoundaryCircuitWhenDownstreamCircuitsExist> downcastCircuit = boost::dynamic_pointer_cast<NetlistBoundaryCircuitWhenDownstreamCircuitsExist> (*upstreamBCCircuit);
        int upstreamCircuitIndex = downcastCircuit->getCircuitIndex();

        int threeDInterfaceComputedFlowLocation = downcastCircuit->getLocationOf3DInterfaceComputedFlowInSolutionVector();
        m_indicesOf3DInterfaceComputedFlowsInUpstreamSolutionVectors.insert(std::make_pair(upstreamCircuitIndex,threeDInterfaceComputedFlowLocation));
        if (downcastCircuit->hasPrescribedFlowAcross3DInterface())
        {
            // std::cout << "inserting upstreamCircuitIndex " << upstreamCircuitIndex << " in m_columnIndicesOf3DInterfacePrescribedFlowsInUpstreamLinearSystems" << std::endl;
            int threeDInterfaceFlowPrescriptionLocationColumn = downcastCircuit->getColumnOf3DInterfacePrescribedFlowInLinearSystem();
            m_columnIndicesOf3DInterfacePrescribedFlowsInUpstreamLinearSystems.insert(std::make_pair(upstreamCircuitIndex,threeDInterfaceFlowPrescriptionLocationColumn));
        }
        
        int threeDInterfacePressureLocation = downcastCircuit->getLocationOf3DInterfaceComputedPressureInSolutionVector();
        m_indicesOf3DInterfaceComputedPressuresInUpstreamSolutionVectors.insert(std::make_pair(upstreamCircuitIndex,threeDInterfacePressureLocation));
        if (downcastCircuit->hasPrescribedPressureAcross3DInterface())
        {
            int threeDInterfacePressurePrescriptionLocationColumn = downcastCircuit->getColumnOf3DInterfacePrescribedPressureInLinearSystem();
            m_columnIndicesOf3DInterfacePrescribedPressuresInUpstreamLinearSystems.insert(std::make_pair(upstreamCircuitIndex,threeDInterfacePressurePrescriptionLocationColumn));
        }
    }

    createVectorsAndMatricesForCircuitLinearSystem();

    m_indicesOfFirstRowOfEachSubcircuitContributionInClosedLoopMatrix.clear();
    m_indicesOfFirstColumnOfEachSubcircuitContributionInClosedLoopMatrix.clear();
    // Tile the matrices to make the full closed loop system matrix
    {
        // these variables will be used during tiling to mark where the top-left corner of
        // the next matrix that we tile into the system goes:
        m_nextBlankSystemMatrixRow = 0;
        m_nextBlankSystemMatrixColumn = 0;

        // assert(m_indicesOfFirstRowOfEachSubcircuitContributionInClosedLoopMatrix.size() == 0);
        // assert(m_indicesOfFirstColumnOfEachSubcircuitContributionInClosedLoopMatrix.size() == 0);

        for (size_t upstreamCircuit = 0; upstreamCircuit < m_numberOfUpstreamCircuits; upstreamCircuit++)
        {
            std::pair<int, Mat> indexAmongstNetlistTagAndMatrix = m_matrixContributionsFromUpstreamBoundaryConditions.front();
            Mat nextMatrixToAddToSystem = indexAmongstNetlistTagAndMatrix.second;

            // std::cout << "System matrix for upstream circuit " << upstreamCircuit << ":" << std::endl;
            // errFlag = MatView(nextMatrixToAddToSystem,PETSC_VIEWER_STDOUT_WORLD); CHKERRABORT(PETSC_COMM_SELF,errFlag);

            m_matrixContributionsFromUpstreamBoundaryConditions.pop();
            PetscInt numberOfRows;
            PetscInt numberOfColumns;
            errFlag = MatGetSize(nextMatrixToAddToSystem, &numberOfRows, &numberOfColumns); CHKERRABORT(PETSC_COMM_SELF,errFlag);

            // Extract the actual data array in the petsc matrix nextMatrixToAddToSystem, so we can pass it to MatSetValues, for inclusion in m_closedLoopSystemMatrix:
            PetscScalar* rawDataInNextMatrixToAddToSystem;
            errFlag = MatDenseGetArray(nextMatrixToAddToSystem, &rawDataInNextMatrixToAddToSystem); CHKERRABORT(PETSC_COMM_SELF,errFlag);

            // Irritatingly, MatDenseGetArray returns a column-major array (in the case of matseqdense - the Petsc type of matrix that we have here),
            // but MatSetValues (below) expects a row-major array. So we need to transpose it. Le sigh.
            PetscScalar* transposedRawDataInNextMatrixToAddToSystem = new PetscScalar[numberOfColumns * numberOfRows];
            int writeLocation = 0;
            for (int row=0; row < numberOfRows; row++)
            {
                for (int column=0; column < numberOfColumns; column++)
                {
                    transposedRawDataInNextMatrixToAddToSystem[writeLocation] = rawDataInNextMatrixToAddToSystem[row + column*numberOfRows];
                    writeLocation++;
                }
            }
            assert(writeLocation = numberOfRows*numberOfColumns);

            // I'm not convinced this call is necessary given that we're about to destroy nextMatrixToAddToSystem anyway,
            // but the Petsc documentation says I MUST call it after MatDenseGetArray once access to an array is no longer needed,
            // and who am I to argue with a block-capital imperative?
            errFlag = MatDenseRestoreArray(nextMatrixToAddToSystem, &rawDataInNextMatrixToAddToSystem); CHKERRABORT(PETSC_COMM_SELF,errFlag);

            // std::cout << "first numRows entries in the underlying array for upstream circuit " << upstreamCircuit << ":" << std::endl;
            // for (int jj=0; jj<numberOfRows;jj++)
            // {
            //     for (int ii=0; ii<numberOfColumns; ii++)
            //     {
            //         std::cout << rawDataInNextMatrixToAddToSystem[jj + ii*numberOfRows] << " ";
            //     }
            //     std::cout << std::endl;
            // }

            // Create location data for the rows and columns where we will place nextMatrixToAddToSystem in m_closedLoopSystemMatrix:
            PetscInt globalRowIndices[numberOfRows];
            createContiguousIntegerRange(m_nextBlankSystemMatrixRow, numberOfRows, globalRowIndices);
            PetscInt globalColumnIndices[numberOfColumns];
            createContiguousIntegerRange(m_nextBlankSystemMatrixColumn, numberOfColumns, globalColumnIndices);


            // std::cout << "first numRows entries in the underlying array for upstream circuit " << upstreamCircuit << ":" << std::endl;
            // for (int jj=0; jj<numberOfRows;jj++)
            // {
            //     for (int ii=0; ii<numberOfColumns; ii++)
            //     {
            //         std::cout << transposedRawDataInNextMatrixToAddToSystem[ii + jj*numberOfColumns] << " ";
            //     }
            //     std::cout << std::endl;
            // }
            // std::cout << "numberOfRows: " << numberOfRows << std::endl;
            // for (int ii = 0; ii < numberOfRows; ii++)
            // {
            //     std::cout << globalRowIndices[ii] << " ";
            // }
            // std::cout << std::endl;
            // std::cout << "numberOfColumns: " << numberOfColumns << std::endl;
            // for (int ii = 0; ii < numberOfColumns; ii++)
            // {
            //     std::cout << globalColumnIndices[ii] << " ";
            // }
            // std::cout << std::endl;
            // PetscInt numberOfRows_tmp;
            // PetscInt numberOfColumns_tmp;
            // errFlag = MatGetSize(m_closedLoopSystemMatrix, &numberOfRows_tmp, &numberOfColumns_tmp); CHKERRABORT(PETSC_COMM_SELF,errFlag);
            // std::cout << "rows: " << numberOfRows_tmp << " columns: " << numberOfColumns_tmp << std::endl;

            errFlag = MatSetValues(m_closedLoopSystemMatrix, numberOfRows, globalRowIndices, numberOfColumns, globalColumnIndices, transposedRawDataInNextMatrixToAddToSystem, INSERT_VALUES);
            delete[] transposedRawDataInNextMatrixToAddToSystem;

            int indexAmongstNetlists = indexAmongstNetlistTagAndMatrix.first;            
            m_indicesOfFirstRowOfEachSubcircuitContributionInClosedLoopMatrix.insert(std::make_pair(indexAmongstNetlists, m_nextBlankSystemMatrixRow));
            m_nextBlankSystemMatrixRow += numberOfRows;

            m_indicesOfFirstColumnOfEachSubcircuitContributionInClosedLoopMatrix.insert(std::make_pair(indexAmongstNetlists, m_nextBlankSystemMatrixColumn));
            m_nextBlankSystemMatrixColumn += numberOfColumns;
        }

        // std::cout << "PARTIAL System matrix for closed loop " << m_index << ":" << std::endl;
        // errFlag = MatAssemblyBegin(m_closedLoopSystemMatrix,MAT_FINAL_ASSEMBLY); CHKERRABORT(PETSC_COMM_SELF,errFlag);
        // errFlag = MatAssemblyEnd(m_closedLoopSystemMatrix,MAT_FINAL_ASSEMBLY); CHKERRABORT(PETSC_COMM_SELF,errFlag);
        // errFlag = MatView(m_closedLoopSystemMatrix,PETSC_VIEWER_STDOUT_WORLD); CHKERRABORT(PETSC_COMM_SELF,errFlag);

        // This adds the location of the first row of the downstream closed loop circuit in m_closedLoopSystemMatrix.
        int minusOneToIndicateDownstreamCircuitNotBC = -1;
        m_indicesOfFirstRowOfEachSubcircuitContributionInClosedLoopMatrix.insert(std::make_pair(minusOneToIndicateDownstreamCircuitNotBC, m_nextBlankSystemMatrixRow));
        // This adds the location of the first column of the downstream closed loop circuit in m_closedLoopSystemMatrix.
        m_indicesOfFirstColumnOfEachSubcircuitContributionInClosedLoopMatrix.insert(std::make_pair(minusOneToIndicateDownstreamCircuitNotBC, m_nextBlankSystemMatrixColumn));

        // Scoping unit containing just the addition of the closed loop downstream circuit matrix to the full system matrix, m_closedLoopSystemMatrix:
        {
            // Add the closed loop downstream circuit's matrix:
            Mat matrixContribution;
            mp_NetlistCircuit->getMatrixContribution(alfi_delt, matrixContribution);

            // std::cout << "Location 2B: Downstream circuit for closed loop " << m_index << ":" << std::endl;
            // errFlag = MatView(matrixContribution,PETSC_VIEWER_STDOUT_WORLD); CHKERRABORT(PETSC_COMM_SELF,errFlag);

            PetscInt numberOfRows;
            PetscInt numberOfColumns;
            errFlag = MatGetSize(matrixContribution, &numberOfRows, &numberOfColumns); CHKERRABORT(PETSC_COMM_SELF,errFlag);

            // Extract the actual data array in the petsc matrix matrixContribution, so we can pass it to MatSetValues, for inclusion in m_closedLoopSystemMatrix:
            PetscScalar* rawDataInMatrix;
            errFlag = MatDenseGetArray(matrixContribution, &rawDataInMatrix); CHKERRABORT(PETSC_COMM_SELF,errFlag);

            // Irritatingly, MatDenseGetArray returns a column-major array (in the case of matseqdense - the Petsc type of matrix that we have here),
            // but MatSetValues (below) expects a row-major array. So we need to transpose it. Le sigh.
            PetscScalar* transposedRawDataInMatrix = new PetscScalar[numberOfColumns * numberOfRows];
            int writeLocation = 0;
            for (int row=0; row < numberOfRows; row++)
            {
                for (int column=0; column < numberOfColumns; column++)
                {
                    transposedRawDataInMatrix[writeLocation] = rawDataInMatrix[row + column*numberOfRows];
                    writeLocation++;
                }
            }
            assert(writeLocation = numberOfRows*numberOfColumns);

            // as noted above, I'm not sure this is needed... but just in case...
            errFlag = MatDenseRestoreArray(matrixContribution, &rawDataInMatrix); CHKERRABORT(PETSC_COMM_SELF,errFlag);

            // Create location data for the rows and columns where we will place matrixContribution in m_closedLoopSystemMatrix:
            PetscInt globalRowIndices[numberOfRows];
            createContiguousIntegerRange(m_nextBlankSystemMatrixRow, numberOfRows, globalRowIndices);
            PetscInt globalColumnIndices[numberOfColumns];
            createContiguousIntegerRange(m_nextBlankSystemMatrixColumn, numberOfColumns, globalColumnIndices);

            // Record where the closed loop downstream circuit's entries are in terms of rows in the system matrix
            // (and thus, also in terms of the m_solutionVector):
            m_boundsOfDownstreamSolutionDataInSolutionVector.first = m_nextBlankSystemMatrixColumn;
            m_boundsOfDownstreamSolutionDataInSolutionVector.second = m_nextBlankSystemMatrixColumn + toZeroIndexing(numberOfColumns);

            // std::cout << "location 2: first numRows entries in the underlying array for closed loop " << ":" << std::endl;
            // for (int jj=0; jj<numberOfRows;jj++)
            // {
            //     for (int ii=0; ii<numberOfColumns; ii++)
            //     {
            //         std::cout << transposedRawDataInMatrix[ii + jj*numberOfRows] << " ";
            //     }
            //     std::cout << std::endl;
            // }
            // std::cout << "numberOfRows: " << numberOfRows << std::endl;
            // for (int ii = 0; ii < numberOfRows; ii++)
            // {
            //     std::cout << globalRowIndices[ii] << " ";
            // }
            // std::cout << std::endl;
            // std::cout << "numberOfColumns: " << numberOfColumns << std::endl;
            // for (int ii = 0; ii < numberOfColumns; ii++)
            // {
            //     std::cout << globalColumnIndices[ii] << " ";
            // }
            // std::cout << std::endl;
            // PetscInt numberOfRows_tmp;
            // PetscInt numberOfColumns_tmp;
            // errFlag = MatGetSize(m_closedLoopSystemMatrix, &numberOfRows_tmp, &numberOfColumns_tmp); CHKERRABORT(PETSC_COMM_SELF,errFlag);
            // std::cout << "rows: " << numberOfRows_tmp << " columns: " << numberOfColumns_tmp << std::endl;

            errFlag = MatSetValues(m_closedLoopSystemMatrix, numberOfRows, globalRowIndices, numberOfColumns, globalColumnIndices, transposedRawDataInMatrix, INSERT_VALUES);
            delete[] transposedRawDataInMatrix;

            // std::cout << "Location 3: System matrix for closed loop " << m_index << ":" << std::endl;
            // errFlag = MatAssemblyBegin(m_closedLoopSystemMatrix,MAT_FINAL_ASSEMBLY); CHKERRABORT(PETSC_COMM_SELF,errFlag);
            // errFlag = MatAssemblyEnd(m_closedLoopSystemMatrix,MAT_FINAL_ASSEMBLY); CHKERRABORT(PETSC_COMM_SELF,errFlag);
            // errFlag = MatView(m_closedLoopSystemMatrix,PETSC_VIEWER_STDOUT_WORLD); CHKERRABORT(PETSC_COMM_SELF,errFlag);

            m_nextBlankSystemMatrixRow += numberOfRows;
            m_nextBlankSystemMatrixColumn += numberOfColumns;
        }
    }

    // By this stage, the matrices for all the subcircuits of the closed loop
    // are tiled into the big matrix, m_closedLoopSystemMatrix. If this assert
    // fails, it means there are no rows left to actually add equations which
    // couple the circuit together.
    assert(m_nextBlankSystemMatrixRow < m_systemSize);

    // We should have used all the columns by now, because we've now tiled
    // all the matrices. All that remains is to add more /rows/ for the
    // interfaces between the circuits (Kirchoff and pressure equality equations)
    assert(m_nextBlankSystemMatrixColumn == m_systemSize);

    // Add the Kirchoff laws for the connecting nodes
    appendKirchoffLawsAtInterfacesBetweenCircuits();

    // std::cout << "Location 4: System matrix for closed loop " << m_index << ":" << std::endl;
    // errFlag = MatAssemblyBegin(m_closedLoopSystemMatrix,MAT_FINAL_ASSEMBLY); CHKERRABORT(PETSC_COMM_SELF,errFlag);
    // errFlag = MatAssemblyEnd(m_closedLoopSystemMatrix,MAT_FINAL_ASSEMBLY); CHKERRABORT(PETSC_COMM_SELF,errFlag);
    // errFlag = MatView(m_closedLoopSystemMatrix,PETSC_VIEWER_STDOUT_WORLD); CHKERRABORT(PETSC_COMM_SELF,errFlag);

    // At the interface between an "upstream" boundary condition and the "downstream" 
    // closed loop circuit, the interfacing nodes are duplicated (as they belong to both
    // the upstream and the downstream circuit). We must enforce equality of pressure
    // at such copies of the same node. Do that now.
    enforcePressureEqualityBetweenDuplicatedNodes();

    // Finalise the matrix construction
    errFlag = MatAssemblyBegin(m_closedLoopSystemMatrix,MAT_FINAL_ASSEMBLY); CHKERRABORT(PETSC_COMM_SELF,errFlag);
    errFlag = MatAssemblyEnd(m_closedLoopSystemMatrix,MAT_FINAL_ASSEMBLY); CHKERRABORT(PETSC_COMM_SELF,errFlag);

    // std::cout << "System matrix for closed loop " << m_index << ":" << std::endl;
    // errFlag = MatView(m_closedLoopSystemMatrix,PETSC_VIEWER_STDOUT_WORLD); CHKERRABORT(PETSC_COMM_SELF,errFlag);

    // LU factor m_closedLoopSystemMatrix
    errFlag = MatLUFactor(m_closedLoopSystemMatrix,NULL,NULL,NULL);CHKERRABORT(PETSC_COMM_SELF,errFlag);

    // Tile the RHS contributions into our closed loop RHS:
    {
        // To track where the start of the next rhs contribution goes in m_closedLoopRHS:
        m_nextBlankRhsRow = 0;

        // Zero out the RHS ready for the construction on this timestep:
        errFlag = VecZeroEntries(m_closedLoopRHS);CHKERRABORT(PETSC_COMM_SELF,errFlag);

        for (size_t upstreamCircuit = 0; upstreamCircuit < m_numberOfUpstreamCircuits; upstreamCircuit++)
        {
            Vec nextVectorToAddToClosedLoopRHS = m_rhsContributionsFromUpstreamBoundaryConditions.front();
            m_rhsContributionsFromUpstreamBoundaryConditions.pop();
            PetscInt numberOfRows;
            errFlag = VecGetSize(nextVectorToAddToClosedLoopRHS, &numberOfRows); CHKERRABORT(PETSC_COMM_SELF, errFlag);

            // Get the raw data to add to the RHS:
            PetscScalar* rawDataInNextVectorToAddToClosedLoopRHS;
            errFlag = VecGetArray(nextVectorToAddToClosedLoopRHS, &rawDataInNextVectorToAddToClosedLoopRHS); CHKERRABORT(PETSC_COMM_SELF, errFlag);
            
            PetscInt globalRowIndices[numberOfRows];
            createContiguousIntegerRange(m_nextBlankRhsRow, numberOfRows, globalRowIndices);
            errFlag = VecSetValues(m_closedLoopRHS, numberOfRows, globalRowIndices, rawDataInNextVectorToAddToClosedLoopRHS, INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF, errFlag);

            errFlag = VecRestoreArray(nextVectorToAddToClosedLoopRHS, &rawDataInNextVectorToAddToClosedLoopRHS); CHKERRABORT(PETSC_COMM_SELF, errFlag);

            m_nextBlankRhsRow += numberOfRows;
        }

        // Scoping unit just to add the downstream closed loop subsection to m_closedLoopRHS:
        {
            // Add the closed loop downstream circuit's RHS:
            Vec rhsContribution;
            mp_NetlistCircuit->getRHSContribution(rhsContribution);

            PetscInt numberOfRows;
            errFlag = VecGetSize(rhsContribution, &numberOfRows); CHKERRABORT(PETSC_COMM_SELF, errFlag);

            // Get the raw data to add to the RHS:
            PetscScalar* rawDataToAddToClosedLoopRHS;
            errFlag = VecGetArray(rhsContribution, &rawDataToAddToClosedLoopRHS); CHKERRABORT(PETSC_COMM_SELF, errFlag);

            PetscInt globalRowIndices[numberOfRows];
            createContiguousIntegerRange(m_nextBlankRhsRow, numberOfRows, globalRowIndices);
            errFlag = VecSetValues(m_closedLoopRHS, numberOfRows, globalRowIndices, rawDataToAddToClosedLoopRHS, INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF, errFlag);
            errFlag = VecAssemblyBegin(m_closedLoopRHS); CHKERRABORT(PETSC_COMM_SELF,errFlag);
            errFlag = VecAssemblyEnd(m_closedLoopRHS); CHKERRABORT(PETSC_COMM_SELF,errFlag);

            errFlag = VecRestoreArray(rhsContribution, &rawDataToAddToClosedLoopRHS); CHKERRABORT(PETSC_COMM_SELF, errFlag);

            m_nextBlankRhsRow += numberOfRows;
        }
    }

    // By this stage, the matrices for all the subcircuits of the closed loop
    // are tiled into the big matrix, m_closedLoopSystemMatrix. If this assert
    // fails, it means there are no rows left to actually add equations which
    // couple the circuit together.
    assert(m_nextBlankRhsRow < m_systemSize);

    // std::cout << "m_identityMatrixForPetscInversionHack" << std::endl;
    // errFlag = MatView(m_identityMatrixForPetscInversionHack,PETSC_VIEWER_STDOUT_WORLD); CHKERRABORT(PETSC_COMM_SELF,errFlag);
    // std::cout << "m_closedLoopSystemMatrix" << std::endl;
    // errFlag = MatView(m_closedLoopSystemMatrix,PETSC_VIEWER_STDOUT_WORLD); CHKERRABORT(PETSC_COMM_SELF,errFlag);

    // get the inverse of the system matrix:
    errFlag = MatMatSolve(m_closedLoopSystemMatrix,m_identityMatrixForPetscInversionHack,m_inverseOfClosedLoopMatrix); CHKERRABORT(PETSC_COMM_SELF,errFlag);
    // Release the m_systemMatrix so we can edit it again on the next iteration (we only need the just-computed m_inverseOfClosedLoopMatrix for computations on this step now.)
    errFlag = MatSetUnfactored(m_closedLoopSystemMatrix); CHKERRABORT(PETSC_COMM_SELF,errFlag);

    // std::cout << "m_inverseOfClosedLoopMatrix" << std::endl;
    // errFlag = MatView(m_inverseOfClosedLoopMatrix,PETSC_VIEWER_STDOUT_WORLD); CHKERRABORT(PETSC_COMM_SELF,errFlag);

    // std::cout << "m_inverseOfClosedLoopMatrix" << std::endl;
    // errFlag = MatView(m_inverseOfClosedLoopMatrix,PETSC_VIEWER_STDOUT_WORLD); CHKERRABORT(PETSC_COMM_SELF,errFlag);
    // std::cout << "m_closedLoopRHS: " << std::endl;
    // errFlag = VecView(m_closedLoopRHS,PETSC_VIEWER_STDOUT_WORLD); CHKERRABORT(PETSC_COMM_SELF,errFlag);
    // errFlag = VecView(m_solutionVector,PETSC_VIEWER_STDOUT_WORLD); CHKERRABORT(PETSC_COMM_SELF,errFlag);
    // Solve the system
    assert(m_inverseOfClosedLoopMatrix != PETSC_NULL);
    assert(m_closedLoopRHS != PETSC_NULL);
    assert(m_solutionVector != PETSC_NULL);
    errFlag = MatMult(m_inverseOfClosedLoopMatrix,m_closedLoopRHS,m_solutionVector); CHKERRABORT(PETSC_COMM_SELF,errFlag);
    // std::cout << "m_solutionVector (just computed): " << std::endl;
    // errFlag = VecView(m_solutionVector,PETSC_VIEWER_STDOUT_WORLD); CHKERRABORT(PETSC_COMM_SELF,errFlag);
}

void ClosedLoopDownstreamSubsection::store3DInterfaceFlowSigns()
{
    for (auto upstreamCircuitSharedPtr : m_upstreamBoundaryConditionCircuits)
    {
        m_signForPrescribed3DInterfaceFlow.insert(std::make_pair(upstreamCircuitSharedPtr->getIndexAmongstNetlists(), upstreamCircuitSharedPtr->getInterfaceFlowSign()));
    }
}

void ClosedLoopDownstreamSubsection::buildAndSolveLinearSystemIfNotYetDone(const double alfi_delt)
{
    // Check whether the linear system still needs to be built and solved; if not, do nothing.
    // std::cout << "called buildAndSolveLinearSystemIfNotYetDone 1 " << m_linearSystemAlreadyBuiltAndSolvedOnThisTimestep<<std::endl;
    if (!m_linearSystemAlreadyBuiltAndSolvedOnThisTimestep)
    {
        buildAndSolveLinearSystem_internal(alfi_delt);
        // Set the "done for this timestep" flag. Reset this with ClosedLoopDownstreamSubsection::markLinearSystemAsNeedingBuildingAgain() (from BC manager for now)
        m_linearSystemAlreadyBuiltAndSolvedOnThisTimestep = true;
    }
}

void ClosedLoopDownstreamSubsection::buildAndSolveLinearSystemForUpdateIfNotYetDone(const double delt)
{
    // Check whether the linear system still needs to be built and solved; if not, do nothing.
    // std::cout << "called buildAndSolveLinearSystemIfNotYetDone 2 " << m_linearSystemAlreadyUpdatedOnThisTimestep << std::endl;
    if (!m_linearSystemAlreadyUpdatedOnThisTimestep)
    {
        buildAndSolveLinearSystem_internal(delt);
        // Set the "done for this timestep" flag. Reset this with ClosedLoopDownstreamSubsection::markLinearSystemAsNeedingUpdatingAgain() (from BC manager for now)
        m_linearSystemAlreadyUpdatedOnThisTimestep = true;
    }
}

void ClosedLoopDownstreamSubsection::createVectorsAndMatricesForCircuitLinearSystem()
{
    PetscErrorCode errFlag;
    // Create a vector to hold the solution
    errFlag = VecCreate(PETSC_COMM_SELF,&m_solutionVector);CHKERRABORT(PETSC_COMM_SELF,errFlag);
    errFlag = VecSetType(m_solutionVector,VECSEQ);CHKERRABORT(PETSC_COMM_SELF,errFlag); // Make m_solutionVector a VECSEQ sequential vector
    errFlag = VecSetSizes(m_solutionVector,m_systemSize,m_systemSize); CHKERRABORT(PETSC_COMM_SELF,errFlag);
    errFlag = VecZeroEntries(m_solutionVector);CHKERRABORT(PETSC_COMM_SELF,errFlag);
    errFlag = VecAssemblyBegin(m_solutionVector); CHKERRABORT(PETSC_COMM_SELF,errFlag);
    errFlag = VecAssemblyEnd(m_solutionVector); CHKERRABORT(PETSC_COMM_SELF,errFlag);

    // Create a vector to hold the RHS (m_closedLoopRHS)
    errFlag = VecCreate(PETSC_COMM_SELF,&m_closedLoopRHS);CHKERRABORT(PETSC_COMM_SELF,errFlag);
    errFlag = VecSetType(m_closedLoopRHS,VECSEQ);CHKERRABORT(PETSC_COMM_SELF,errFlag); // Make m_solutionVector a VECSEQ sequential vector
    errFlag = VecSetSizes(m_closedLoopRHS,m_systemSize,m_systemSize); CHKERRABORT(PETSC_COMM_SELF,errFlag);
    errFlag = VecZeroEntries(m_closedLoopRHS);CHKERRABORT(PETSC_COMM_SELF,errFlag);
    errFlag = VecAssemblyBegin(m_closedLoopRHS); CHKERRABORT(PETSC_COMM_SELF,errFlag);
    errFlag = VecAssemblyEnd(m_closedLoopRHS); CHKERRABORT(PETSC_COMM_SELF,errFlag);

    // Create the vector to hold the system matrix
    errFlag = MatCreateSeqDense(PETSC_COMM_SELF,m_systemSize,m_systemSize,NULL,&m_closedLoopSystemMatrix);CHKERRABORT(PETSC_COMM_SELF,errFlag);
    errFlag = MatZeroEntries(m_closedLoopSystemMatrix);CHKERRABORT(PETSC_COMM_SELF,errFlag);

    // Create a matrix to store the inverse, m_inverseOfClosedLoopMatrix:
    errFlag = MatCreateSeqDense(PETSC_COMM_SELF,m_systemSize,m_systemSize,NULL,&m_inverseOfClosedLoopMatrix);CHKERRABORT(PETSC_COMM_SELF,errFlag);
    errFlag = MatZeroEntries(m_inverseOfClosedLoopMatrix);CHKERRABORT(PETSC_COMM_SELF,errFlag);

    // Create an identity matrix for use when inverting the system matrix:
    errFlag = MatCreateSeqDense(PETSC_COMM_SELF,m_systemSize,m_systemSize,NULL,&m_identityMatrixForPetscInversionHack);CHKERRABORT(PETSC_COMM_SELF,errFlag);
    errFlag = MatZeroEntries(m_identityMatrixForPetscInversionHack);CHKERRABORT(PETSC_COMM_SELF,errFlag);
    // Fill the diagonal with ones:
    for (int ii=0; ii<m_systemSize; ii++)
    {
        errFlag = MatSetValue(m_identityMatrixForPetscInversionHack,ii,ii,1.0,INSERT_VALUES);CHKERRABORT(PETSC_COMM_SELF,errFlag);
    }
    errFlag = MatAssemblyBegin(m_identityMatrixForPetscInversionHack,MAT_FINAL_ASSEMBLY); CHKERRABORT(PETSC_COMM_SELF,errFlag);
    errFlag = MatAssemblyEnd(m_identityMatrixForPetscInversionHack,MAT_FINAL_ASSEMBLY); CHKERRABORT(PETSC_COMM_SELF,errFlag);
}

void ClosedLoopDownstreamSubsection::initialisePetscArrayNames()
{
    m_closedLoopRHS = PETSC_NULL;
    m_solutionVector = PETSC_NULL;
    m_closedLoopSystemMatrix = PETSC_NULL;
    m_inverseOfClosedLoopMatrix = PETSC_NULL;
    m_identityMatrixForPetscInversionHack = PETSC_NULL;
}

void ClosedLoopDownstreamSubsection::terminatePetscArrays()
{
    PetscErrorCode errFlag;
    if (m_closedLoopRHS)
    {
        errFlag = VecDestroy(&m_closedLoopRHS); CHKERRABORT(PETSC_COMM_SELF,errFlag);
    }
    if (m_solutionVector)
    {
        errFlag = VecDestroy(&m_solutionVector); CHKERRABORT(PETSC_COMM_SELF,errFlag);
    }
    if (m_closedLoopSystemMatrix)
    {
        errFlag = MatDestroy(&m_closedLoopSystemMatrix); CHKERRABORT(PETSC_COMM_SELF,errFlag);
    }
    if (m_inverseOfClosedLoopMatrix)
    {
        errFlag = MatDestroy(&m_inverseOfClosedLoopMatrix); CHKERRABORT(PETSC_COMM_SELF,errFlag);
    }
    if (m_identityMatrixForPetscInversionHack)
    {
        errFlag = MatDestroy(&m_identityMatrixForPetscInversionHack); CHKERRABORT(PETSC_COMM_SELF,errFlag);
    }
}

void ClosedLoopDownstreamSubsection::generateCircuitInterfaceNodeData()
{
    std::vector<int> downstreamNodeIndices;
    std::vector<int> upstreamNodeIndices;
    std::vector<int> upstreamSurfaceIndices;
    mp_NetlistCircuit->getSharedNodeDownstreamAndUpstreamAndCircuitUpstreamIndices(downstreamNodeIndices, upstreamNodeIndices, upstreamSurfaceIndices);

    for (size_t interfaceConnectionIndex = 0; interfaceConnectionIndex < upstreamNodeIndices.size(); interfaceConnectionIndex++)
    {
        int downstreamNode;
        int upstreamSurface;
        try {
            downstreamNode = downstreamNodeIndices.at(interfaceConnectionIndex);
            upstreamSurface = upstreamSurfaceIndices.at(interfaceConnectionIndex);
        } catch (const std::exception& e) {
            std::cout << e.what() << " observed at line " << __LINE__ << " of " << __FILE__ << std::endl;
            throw;
        }
        // The map will throw an exception if the key does not exist; but we don't want this --
        // we want to create the key if it does not already exist, so we handle out_of_range
        // exceptions appropriately:
        try
        {
            m_mapOfSurfaceIndicesConnectedToEachDownstreamInterfaceNode.at(downstreamNode).insert(upstreamSurface);
        }
        catch (const std::out_of_range& outOfRange)
        {
            std::set<int> toInsert;
            toInsert.insert(upstreamSurface);
            m_mapOfSurfaceIndicesConnectedToEachDownstreamInterfaceNode.insert(std::make_pair(downstreamNode, toInsert));
        }
    }
}


void ClosedLoopDownstreamSubsection::appendKirchoffLawsAtInterfacesBetweenCircuits()
{
    std::vector<int> downstreamNodeIndices;
    std::vector<int> upstreamNodeIndices;
    std::vector<int> upstreamSurfaceIndices;
    mp_NetlistCircuit->getSharedNodeDownstreamAndUpstreamAndCircuitUpstreamIndices(downstreamNodeIndices, upstreamNodeIndices, upstreamSurfaceIndices);

    std::set<int> nodesInterfacingWithUpstreamCircuits = mp_NetlistCircuit->getPressureNodesConnectingToUpstreamCircuits();

    // Loop the downstream circuit nodes which connect to upstream circuits
    for (auto downstreamInterfaceNode = nodesInterfacingWithUpstreamCircuits.begin(); downstreamInterfaceNode != nodesInterfacingWithUpstreamCircuits.end(); downstreamInterfaceNode++)
    {
        // Loop the upstream BC circuits, checking whether they connect to this downstream node
        for (auto upstreamCircuit = m_upstreamBoundaryConditionCircuits.begin(); upstreamCircuit != m_upstreamBoundaryConditionCircuits.end(); upstreamCircuit++)
        {
            // If it does, add the contribution to the Kirchoff equation that we're currently writing
            int upstreamSurfaceIndex = (*upstreamCircuit)->getSurfaceIndex();
            bool thisUpstreamCircuitConnectsToThisDownstreamNode;
            try {
                thisUpstreamCircuitConnectsToThisDownstreamNode = (m_mapOfSurfaceIndicesConnectedToEachDownstreamInterfaceNode.at(*downstreamInterfaceNode).count(upstreamSurfaceIndex) == 1);
            } catch (const std::exception& e) {
                std::cout << e.what() << " observed at line " << __LINE__ << " of " << __FILE__ << std::endl;
                throw;
            }
            if (thisUpstreamCircuitConnectsToThisDownstreamNode)
            {
                boost::shared_ptr<NetlistBoundaryCircuitWhenDownstreamCircuitsExist> downcastUpstreamCircuit = boost::dynamic_pointer_cast<NetlistBoundaryCircuitWhenDownstreamCircuitsExist> (*upstreamCircuit);
                
                // Get the circuit data for this boundary (so we can check the sign for
                // the Kirchoff equation by seeing if a node is the start node or end
                // node of a component):
                boost::shared_ptr<CircuitData> upstreamCircuitData = downcastUpstreamCircuit->getCircuitDescription();
                int numberOfHistoryPressures_upstreamCircuit = downcastUpstreamCircuit->getNumberOfHistoryPressures();
                // Get the column where the data block for the current upstream boundarycondition circuit
                // begins in the big matrix, m_closedLoopSystemMatrix
                int upstreamCircuitIndex = (*upstreamCircuit)->getIndexAmongstNetlists();
                int columnOffsetOfCurrentUpstreamCircuit;
                try {
                    columnOffsetOfCurrentUpstreamCircuit = m_indicesOfFirstColumnOfEachSubcircuitContributionInClosedLoopMatrix.at(upstreamCircuitIndex);
                } catch (const std::exception& e) {
                    std::cout << e.what() << " observed at line " << __LINE__ << " of " << __FILE__ << std::endl;
                    throw;
                }

                // Write the part of the Kirchoff equation for this node (multipleIncidentCurrentNode)
                // which corresponds to the components incident at multipleIndidentCurrentNode in
                // the upstream boundary condition:
                int multipleIncidentCurrentNode = mp_NetlistCircuit->convertInterfaceNodeIndexFromDownstreamToUpstreamCircuit(upstreamSurfaceIndex,*downstreamInterfaceNode);
                writePartOfKirchoffEquationIntoClosedLoopSysteMatrix(upstreamCircuitData, multipleIncidentCurrentNode, m_nextBlankSystemMatrixRow, numberOfHistoryPressures_upstreamCircuit, columnOffsetOfCurrentUpstreamCircuit);
            }
        }

        // Write the part of the Kirchoff equation for this node (multipleIncidentCurrentNode)
        // which corresponds to the components incident at multipleIndidentCurrentNode in
        // the downstream closed loop subsection circuit:
        int numberOfHistoryPressures_downstreamCircuit = mp_NetlistCircuit->getNumberOfHistoryPressures();
        int minusOneToIndicateDownstreamCircuitNotBC = -1;
        int columnOffsetOfDownstreamClosedLoopCircuit;
        try {
            columnOffsetOfDownstreamClosedLoopCircuit = m_indicesOfFirstColumnOfEachSubcircuitContributionInClosedLoopMatrix.at(minusOneToIndicateDownstreamCircuitNotBC);
        } catch (const std::exception& e) {
            std::cout << e.what() << " observed at line " << __LINE__ << " of " << __FILE__ << std::endl;
            throw;
        }
        boost::shared_ptr<CircuitData> downstreamCircuitData = mp_NetlistCircuit->getCircuitDescription();
        writePartOfKirchoffEquationIntoClosedLoopSysteMatrix(downstreamCircuitData, *downstreamInterfaceNode, m_nextBlankSystemMatrixRow, numberOfHistoryPressures_downstreamCircuit, columnOffsetOfDownstreamClosedLoopCircuit);
        
        // Move to the next available row in the matrix that isn't used for an equation yet:
        m_nextBlankSystemMatrixRow++;

    }
}

void ClosedLoopDownstreamSubsection::writePressuresFlowsAndVolumes(int& nextTimestepWrite_start)
{
    mp_NetlistCircuit->writePressuresFlowsAndVolumes(nextTimestepWrite_start);
}

// void ClosedLoopDownstreamSubsection::loadPressuresFlowsAndVolumesOnRestart(const int startingTimeStepIndex)
// {
//     mp_NetlistCircuit->loadPressuresFlowsAndVolumesOnRestart(startingTimeStepIndex);
// }

void ClosedLoopDownstreamSubsection::writePartOfKirchoffEquationIntoClosedLoopSysteMatrix(const boost::shared_ptr<const CircuitData> circuitData, const int multipleIncidentCurrentNode, const int row, const int numberOfHistoryPressures, const int columnOffset)
{
    PetscErrorCode errFlag;
    // Do the equations for the nodes with multiple incident currents (just the part for the upstream boundary condition circuit first...
    // we'll append the currents for the components in the downstream closed loop circuit in a moment...)
    for (size_t component = 0; component < circuitData->numberOfComponents; component++)
    {
      int column = component + circuitData->numberOfPressureNodes + numberOfHistoryPressures + columnOffset;
      bool foundMultipleIncidentCurrentsForEndNode;
      try {
          foundMultipleIncidentCurrentsForEndNode = (circuitData->components.at(component)->endNode->getIndex() == multipleIncidentCurrentNode); 
      } catch (const std::exception& e) {
          std::cout << e.what() << " observed at line " << __LINE__ << " of " << __FILE__ << std::endl;
          throw;
      }
      if (foundMultipleIncidentCurrentsForEndNode)
      {
        errFlag = MatSetValue(m_closedLoopSystemMatrix,row,column,1.0,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
      }

      bool foundMultipleIncidentCurrentsForStartNode = (circuitData->components.at(component)->startNode->getIndex() == multipleIncidentCurrentNode);
      if (foundMultipleIncidentCurrentsForStartNode)
      {
        errFlag = MatSetValue(m_closedLoopSystemMatrix,row,column,-1.0,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
      }
    }
}

int ClosedLoopDownstreamSubsection::getCircuitIndexFromSurfaceIndex(const int upstreamSurfaceIndex) const
{
    int returnValue = -1;
    for (auto upstreamBC = m_upstreamBoundaryConditionCircuits.begin(); upstreamBC != m_upstreamBoundaryConditionCircuits.end(); upstreamBC++)
    {
        if ((*upstreamBC)->surfaceIndexMatches(upstreamSurfaceIndex))
        {
            boost::shared_ptr<NetlistBoundaryCircuitWhenDownstreamCircuitsExist> downcastCircuit = boost::dynamic_pointer_cast<NetlistBoundaryCircuitWhenDownstreamCircuitsExist> (*upstreamBC);
            returnValue = downcastCircuit->getCircuitIndex();
        }
    }
    assert(returnValue != -1);
    return returnValue;
}

void ClosedLoopDownstreamSubsection::enforcePressureEqualityBetweenDuplicatedNodes()
{
    PetscErrorCode errFlag;
    // Vectors to hold the pass-by-reference return from getSharedNodeDownstreamAndUpstreamAndCircuitUpstreamIndices
    std::vector<int> downstreamNodeIndices;
    std::vector<int> upstreamNodeIndices;
    std::vector<int> upstreamSurfaceIndices;
    mp_NetlistCircuit->getSharedNodeDownstreamAndUpstreamAndCircuitUpstreamIndices(downstreamNodeIndices, upstreamNodeIndices, upstreamSurfaceIndices);

    // The upstreamSurfaceIndices (which are solver.inp indices for the surfaces of the 3D model) need to be converted to 
    // the indices of the circuits themselves (i.e. 0th, 1st, 2nd,... circuit; whereas the surfaces may have non-consecutive 
    // arbitrary numbering)
    std::vector<int> upstreamCircuitIndices;
    for (auto upstreamSurfaceIndex = upstreamSurfaceIndices.begin(); upstreamSurfaceIndex != upstreamSurfaceIndices.end(); upstreamSurfaceIndex++)
    {
        int upstreamCircuitIndex = getCircuitIndexFromSurfaceIndex(*upstreamSurfaceIndex);
        upstreamCircuitIndices.push_back(upstreamCircuitIndex);
    }

    for (size_t sharedNodeIndex = 0; sharedNodeIndex < downstreamNodeIndices.size(); sharedNodeIndex++)
    {
        // The entry for the upstream copy of the pressure node:
        {
            assert(m_nextBlankSystemMatrixRow < m_systemSize);
            int currentUpstreamCircuitIndex;
            int upstreamNodeIndex;
            int column;
            try {
                currentUpstreamCircuitIndex = upstreamCircuitIndices.at(sharedNodeIndex);
                upstreamNodeIndex = upstreamNodeIndices.at(sharedNodeIndex);
                column = m_indicesOfFirstColumnOfEachSubcircuitContributionInClosedLoopMatrix.at(currentUpstreamCircuitIndex) + toZeroIndexing(upstreamNodeIndex);
            } catch (const std::exception& e) {
                std::cout << e.what() << " observed at line " << __LINE__ << " of " << __FILE__ << std::endl;
                throw;
            }
            assert(column < m_systemSize);
            errFlag = MatSetValue(m_closedLoopSystemMatrix,m_nextBlankSystemMatrixRow,column,1.0,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
            // std::cout << "just wrote to row " << m_nextBlankSystemMatrixRow << " and column " << column << " in enforcePressure." << std::endl;
        }

        // MAGICAL_DEBUG();
        // std::cout << "Location (see above MAGICAL_DEBUG): System matrix for closed loop " << m_index << ":" << std::endl;
        // errFlag = MatAssemblyBegin(m_closedLoopSystemMatrix,MAT_FINAL_ASSEMBLY); CHKERRABORT(PETSC_COMM_SELF,errFlag);
        // errFlag = MatAssemblyEnd(m_closedLoopSystemMatrix,MAT_FINAL_ASSEMBLY); CHKERRABORT(PETSC_COMM_SELF,errFlag);
        // errFlag = MatView(m_closedLoopSystemMatrix,PETSC_VIEWER_STDOUT_WORLD); CHKERRABORT(PETSC_COMM_SELF,errFlag);

        // The entry for the downstream (i.e. local) copy of the pressure node:
        {
            assert(m_nextBlankSystemMatrixRow < m_systemSize);
            int downstreamNodeIndex;
            try {
                downstreamNodeIndex = downstreamNodeIndices.at(sharedNodeIndex);
            } catch (const std::exception& e) {
                std::cout << e.what() << " observed at line " << __LINE__ << " of " << __FILE__ << std::endl;
                throw;
            }
            int minusOneToIndicateDownstreamCircuitNotBC = -1;
            int column = m_indicesOfFirstColumnOfEachSubcircuitContributionInClosedLoopMatrix.at(minusOneToIndicateDownstreamCircuitNotBC) + toZeroIndexing(downstreamNodeIndex);
            assert(column < m_systemSize);
            errFlag = MatSetValue(m_closedLoopSystemMatrix,m_nextBlankSystemMatrixRow,column,-1.0,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
            // std::cout << "just wrote to row " << m_nextBlankSystemMatrixRow << " and column " << column << " in enforcePressure." << std::endl;
        }
        m_nextBlankSystemMatrixRow++;

        // MAGICAL_DEBUG();
        // std::cout << "Location (see above MAGICAL_DEBUG): System matrix for closed loop " << m_index << ":" << std::endl;
        // errFlag = MatAssemblyBegin(m_closedLoopSystemMatrix,MAT_FINAL_ASSEMBLY); CHKERRABORT(PETSC_COMM_SELF,errFlag);
        // errFlag = MatAssemblyEnd(m_closedLoopSystemMatrix,MAT_FINAL_ASSEMBLY); CHKERRABORT(PETSC_COMM_SELF,errFlag);
        // errFlag = MatView(m_closedLoopSystemMatrix,PETSC_VIEWER_STDOUT_WORLD); CHKERRABORT(PETSC_COMM_SELF,errFlag);
    } 
}

void ClosedLoopDownstreamSubsection::markLinearSystemAsNeedingBuildingAgain()
{
    m_linearSystemAlreadyBuiltAndSolvedOnThisTimestep = false;
}

void ClosedLoopDownstreamSubsection::markLinearSystemAsNeedingUpdatingAgain()
{
    m_linearSystemAlreadyUpdatedOnThisTimestep = false;
}

std::pair<double,double> ClosedLoopDownstreamSubsection::getImplicitCoefficients(const int boundaryConditionIndex) const
{
    PetscErrorCode errFlag;
    // assert the linear system has been solved:
    assert(m_linearSystemAlreadyBuiltAndSolvedOnThisTimestep);

    // The linear system is solved, so we can just extract the necessary values from the resulting solution
    // vector and inverted system matrix:
    int rowToGet[] = {m_indicesOfFirstColumnOfEachSubcircuitContributionInClosedLoopMatrix.at(boundaryConditionIndex)};
    const int numberOfValuesToGet = 1;
    PetscScalar valueFromInverseOfSystemMatrix;
    
    std::pair<double,double> implicitCoefficientsToReturn;
    int columnIndexOf3DInterfaceFlow;
    try {
        columnIndexOf3DInterfaceFlow = m_columnIndicesOf3DInterfacePrescribedFlowsInUpstreamLinearSystems.at(boundaryConditionIndex) + 
                                            m_indicesOfFirstRowOfEachSubcircuitContributionInClosedLoopMatrix.at(boundaryConditionIndex);
    } catch (const std::exception& e) {
        std::cout << e.what() << " observed at line " << __LINE__ << " of " << __FILE__ << std::endl;
        throw;
    }

    errFlag = MatGetValues(m_inverseOfClosedLoopMatrix,numberOfValuesToGet,rowToGet,numberOfValuesToGet,&columnIndexOf3DInterfaceFlow,&valueFromInverseOfSystemMatrix);CHKERRABORT(PETSC_COMM_SELF,errFlag);
    implicitCoefficientsToReturn.first = getSignForPrescribed3DInterfaceFlow(boundaryConditionIndex) * valueFromInverseOfSystemMatrix;

    PetscScalar valueFromRHS;
    errFlag = VecGetValues(m_closedLoopRHS,numberOfValuesToGet,&columnIndexOf3DInterfaceFlow,&valueFromRHS);CHKERRABORT(PETSC_COMM_SELF,errFlag);
    
    PetscScalar valueFromSolutionVector;
    errFlag = VecGetValues(m_solutionVector,numberOfValuesToGet,rowToGet,&valueFromSolutionVector);CHKERRABORT(PETSC_COMM_SELF,errFlag);
    
    implicitCoefficientsToReturn.second = valueFromSolutionVector - valueFromInverseOfSystemMatrix * valueFromRHS;//\todo make dynamic

    // std::cout << "and just set 2 " << implicitCoefficientsToReturn.first << " " <<implicitCoefficientsToReturn.second << std::endl;
    
    return implicitCoefficientsToReturn;
}

double ClosedLoopDownstreamSubsection::getSignForPrescribed3DInterfaceFlow(const int boundaryConditionIndex) const {
    // std::cout << "requested data for BC with index " << boundaryConditionIndex << std::endl;
    return m_signForPrescribed3DInterfaceFlow.at(boundaryConditionIndex);
}

double ClosedLoopDownstreamSubsection::getComputedInterfacePressure(const int boundaryConditionIndex) const
{
    // Note that in some cases, this "computed interface pressure" will actually have been prescribed,
    // and so the value this function returns is not strictly "computed", but rather, just imposed
    // and passed through the linear system without change.
    PetscErrorCode errFlag;
    // assert the linear system has been solved:
    assert(m_linearSystemAlreadyBuiltAndSolvedOnThisTimestep);

    const int numberOfValuesToGet = 1;
    int vectorIndexOf3DInterfacePressure;
    try {
        vectorIndexOf3DInterfacePressure = m_indicesOf3DInterfaceComputedPressuresInUpstreamSolutionVectors.at(boundaryConditionIndex) + 
                                            m_indicesOfFirstColumnOfEachSubcircuitContributionInClosedLoopMatrix.at(boundaryConditionIndex);
    } catch (const std::exception& e) {
        std::cout << e.what() << " observed at line " << __LINE__ << " of " << __FILE__ << std::endl;
        throw;
    }

    PetscScalar pressureFromRHS;
    errFlag = VecGetValues(m_solutionVector,numberOfValuesToGet,&vectorIndexOf3DInterfacePressure,&pressureFromRHS);CHKERRABORT(PETSC_COMM_SELF,errFlag);

    return pressureFromRHS;
}

double ClosedLoopDownstreamSubsection::getComputedInterfaceFlow(const int boundaryConditionIndex) const
{
    // Note that in some cases, this "computed interface flow" will actually have been prescribed,
    // and so the value this function returns is not strictly "computed", but rather, just imposed
    // and passed through the linear system without change.
    PetscErrorCode errFlag;
    // assert the linear system has been solved:
    assert(m_linearSystemAlreadyBuiltAndSolvedOnThisTimestep);

    const int numberOfValuesToGet = 1;
    int vectorIndexOf3DInterfaceFlow;
    try {
        vectorIndexOf3DInterfaceFlow = m_indicesOf3DInterfaceComputedFlowsInUpstreamSolutionVectors.at(boundaryConditionIndex) + 
                                            m_indicesOfFirstColumnOfEachSubcircuitContributionInClosedLoopMatrix.at(boundaryConditionIndex);
    } catch (const std::exception& e) {
        std::cout << e.what() << " observed at line " << __LINE__ << " of " << __FILE__ << std::endl;
        throw;
    }

    PetscScalar flowFromRHS;
    errFlag = VecGetValues(m_solutionVector,numberOfValuesToGet,&vectorIndexOf3DInterfaceFlow,&flowFromRHS);CHKERRABORT(PETSC_COMM_SELF,errFlag);

    return flowFromRHS * getSignForPrescribed3DInterfaceFlow(boundaryConditionIndex);
}

void ClosedLoopDownstreamSubsection::createContiguousIntegerRange(const int startingInteger, const int numberOfIntegers, PetscInt* const arrayToFill)
{
    for (int ii = 0; ii < numberOfIntegers; ii++)
    {
        arrayToFill[ii] = startingInteger + ii;
    }
}

int ClosedLoopDownstreamSubsection::getIndexOfClosedLoop_zeroIndexed() const
{
    return toZeroIndexing(m_index);
}

int ClosedLoopDownstreamSubsection::getIndexOfClosedLoop_oneIndexed() const
{
    return m_index;
}

boost::shared_ptr<NetlistClosedLoopDownstreamCircuit> ClosedLoopDownstreamSubsection::getNetlistCircuit() const
{
    return mp_NetlistCircuit;
}

void ClosedLoopDownstreamSubsection::finalizeLPNAtEndOfTimestep()
{
    mp_NetlistCircuit->finalizeLPNAtEndOfTimestep();
}