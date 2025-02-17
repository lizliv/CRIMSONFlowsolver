from CRIMSONPython import *
from math import pi, cos

# The parameter controller must have exactly this name
class parameterController(abstractParameterController):

	def __init__(self, baseNameOfThisScriptAndOfRelatedFlowOrPressureDatFile, MPIRank):
		abstractParameterController.__init__(self,baseNameOfThisScriptAndOfRelatedFlowOrPressureDatFile, MPIRank)
		self.m_periodicTime = 0.0; #\todo think about this for restarts!
		self.m_heartPeriod = 0.86;
		self.controlSignal = 0.0 #an intiial value
		self.controllerPriority = 1
		self.finishSetup()

	# This method must have exactly this name
	def setFirstTimestepBroadcastValues(self):
		self.clearBroadcastData()
		self.addBroadcastVariable('masterControlSignal', self.controlSignal)
		
	# This method must have exactly this name
	def updateControl(self, delt):
		# print "Master Controller Reporting!"
		# self.printAllRecievedData()
		self.updatePeriodicTime(delt)

		self.controlSignal = cos(self.m_periodicTime)

		self.clearBroadcastData()
		self.addBroadcastVariable('masterControlSignal', self.controlSignal)

	def updatePeriodicTime(self, delt):
	
		self.m_periodicTime = self.m_periodicTime + delt
		# Keep m_periodicTime in the range [0,m_heartPeriod):
		if self.m_periodicTime >= self.m_heartPeriod:
			self.m_periodicTime = self.m_periodicTime - self.m_heartPeriod
