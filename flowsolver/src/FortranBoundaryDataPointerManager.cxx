#include "FortranBoundaryDataPointerManager.hxx"
#include <iostream>


extern "C" void giveflowpointertocpp(int& surfaceIndex,double*& flowPointer) {
	FortranBoundaryDataPointerManager::Get()->setBoundaryFlows(surfaceIndex,flowPointer);
}

extern "C" void givepressurepointertocpp(int& surfaceIndex,double*& pressPointer) {
	FortranBoundaryDataPointerManager::Get()->setBoundaryPressures(surfaceIndex,pressPointer);
}

void FortranBoundaryDataPointerManager::setBoundaryFlows(int surfaceIndex, double* flowPointer)
{
	// This was a bug fix, so this is a guard against regression.
	assert(pointerNotDuplicated(flowPointer));

	m_boundaryFlows.insert(std::pair<int,double*>(surfaceIndex,flowPointer));
	m_hasBoundaryFlows = true;
}
void FortranBoundaryDataPointerManager::setBoundaryPressures(int surfaceIndex, double* pressPointer)
{
	// This was a bug fix, so this is a guard against regression.
	assert(pointerNotDuplicated(pressPointer));

	m_boundaryPressures.insert(std::pair<int,double*>(surfaceIndex,pressPointer));
	m_hasBoundaryPressures = true;
}

double* FortranBoundaryDataPointerManager::getBoundaryFlows(int surfaceIndex)
{
	assert(m_hasBoundaryFlows);
	try {
		return m_boundaryFlows.at(surfaceIndex);
	} catch (const std::exception& e) {
	    std::cout << e.what() << " observed at line " << __LINE__ << " of " << __FILE__ << std::endl;
	    std::cout << "This may mean multidomain.dat is missing." << std::endl;
	    throw;
	}
}
double* FortranBoundaryDataPointerManager::getBoundaryPressures(int surfaceIndex)
{
	assert(m_hasBoundaryPressures);
	try {
		return m_boundaryPressures.at(surfaceIndex);
	} catch (const std::exception& e) {
	    std::cout << e.what() << " observed at line " << __LINE__ << " of " << __FILE__ << std::endl;
	    throw;
	}
}

bool FortranBoundaryDataPointerManager::pointerNotDuplicated(double* pointerAboutToBeStored)
{
	// Change this in a moment if the pointerAboutToBeStored is already owned by FortranBoundaryDataPointerManager:
	bool pointerIsNotDuplicated = true;

	// Ensure it's not in m_boundaryFlows
	for (auto storedPointer = m_boundaryFlows.begin(); storedPointer != m_boundaryFlows.end(); storedPointer++)
	{
		if (storedPointer->second == pointerAboutToBeStored)
		{
			pointerIsNotDuplicated = false;
		}
	}

	// Ensure it's not in m_boundaryPressures:
	for (auto storedPointer = m_boundaryPressures.begin(); storedPointer != m_boundaryPressures.end(); storedPointer++)
	{
		if (storedPointer->second == pointerAboutToBeStored)
		{
			pointerIsNotDuplicated = false;
		}
	}
	
	return pointerIsNotDuplicated;
}