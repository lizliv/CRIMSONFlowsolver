from CRIMSONPython import *
from math import pi, cos

# The parameter controller must have exactly this name
class parameterController(abstractParameterController):

	def __init__(self, baseNameOfThisScriptAndOfRelatedFlowOrPressureDatFile, MPIRank):
		abstractParameterController.__init__(self,baseNameOfThisScriptAndOfRelatedFlowOrPressureDatFile, MPIRank)
		self.controllerPriority = 3
		self.m_periodicTime = 0.0; #\todo think about this for restarts!
		self.m_heartPeriod = 0.86;
		self.finishSetup()

	# This method must have exactly this name
	def setFirstTimestepBroadcastValues(self):
		self.clearBroadcastData()
		self.addBroadcastVariable('three', 3) # just a non-functional example broadcast
		self.addBroadcastVariable('four', 4) # just a non-functional example broadcast

	# This method must have exactly this name
	def updateControl(self, currentParameterValue, delt, dictionaryOfPressuresByNodeIndex, dictionaryOfFlowsByComponentIndex, dictionaryOfVolumesByComponentIndex):
		self.clearBroadcastData()
		self.addBroadcastVariable('three', 3) # just a non-functional example broadcast
		self.addBroadcastVariable('four', 4) # just a non-functional example broadcast
		# print "Node Controller Reporting!"
		# self.printAllRecievedData()

		# Only update the time if this controller is receiving the (otherwise-unused)
		# broadcasts from other controllers. The only purpose of this is to
		# make the test fail (due to the time not being updated properly) if
		# there is a problem with the broadcast reception.
		if self.getRecievedBroadcastValue('elastanceController2','six') == 6:
			self.updatePeriodicTime(delt)

		pressure = self.getRecievedBroadcastValue('masterController','masterControlSignal')
		# print "in nodecontroller:", pressure

		# for key in dictionaryOfPressuresByNodeIndex:
		# 	print "Pressure ", key, " was ", dictionaryOfPressuresByNodeIndex[key]
		# for key in dictionaryOfFlowsByComponentIndex:
		# 	print "Flow ", key, " was ", dictionaryOfFlowsByComponentIndex[key]

		return pressure

	def updatePeriodicTime(self, delt):
	
		self.m_periodicTime = self.m_periodicTime + delt
		# Keep m_periodicTime in the range [0,m_heartPeriod):
		if self.m_periodicTime >= self.m_heartPeriod:
			self.m_periodicTime = self.m_periodicTime - self.m_heartPeriod
