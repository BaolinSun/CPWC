# -*- coding: utf8 -*-
import sys
#sys.path.append("../../Lib")
sys.path.append("D:/Cosmos/v1/RND/SEQ_Data/Lib")
#sys.path.append("/Users/omrindal/Development/Alpinion/SEQ_Data/Lib")
#sys.path.append("/Users/omrindal/Development/Alpinion/Lib_hospital")
from numpy import *
import numpy as np
import TransLib
import SysLib
import FilterLib

transNumElement = 128
numberOfAngles = 11; #The final number is 2*numberOfAngles-1
maxAngle = 16;
depthCm = 5;
numFrames = 10;
center_frequency = 6.0;

# Calculate Angles
angles = np.linspace(0,maxAngle,numberOfAngles);
print(angles);
steerAngleA = [0];
angle_idx = 1;
for i in range(1,numberOfAngles*2-1,2):
	steerAngleA.append(angles[angle_idx]);
	steerAngleA.append(-angles[angle_idx]);
	angle_idx += 1;
spcdCnt = size(steerAngleA)
acqNum = 2*spcdCnt

#(System , Resource, Tw, Tgc, Sequence, Process, Tx, RcvBuffer, Bf, Rx, SeqControl) = SysLib_v8.init_object(objectNum)

System = SysLib.init_object("System")
################################################################
# Specify transducer array
################################################################
System['Transducer'] = TransLib.loadTransducer("L3-8")

################################################################
# Specify system parameters
################################################################
System['Parameters']['transmitNum'] = transNumElement
System['Parameters']['receiveNum'] = 64
System['Parameters']['cineDataType'] = 0 #Channel=0, BeamFormed=1, IQ=2

################################################################
# Specify filter parameters
################################################################ 
Parameter = SysLib.init_object("Parameter", 1)
Parameter[0]['dtgc'] = SysLib.init_object("Parameter.dtgc", 1)
Parameter[0]['dtgc'][0]['dtgcDb'] = ""#[0,0,0,0,0,0,0,0]
SysLib.convertToDtgc(Parameter[0])

Parameter[0]['demod'] = SysLib.init_object("Parameter.demod", 1)
Parameter[0]['demod'][0]['demodMHz'] = ""#[10.5, 10.5, 10.5, 10.0, 10.0, 10.0, 9.5, 9.5]
SysLib.convertToDemod(Parameter[0], System)

Parameter[0]['lpf'] = SysLib.init_object("Parameter.lpf", 1)
Parameter[0]['lpf'][0]['lpfMHz'] = ""#[7.0,7.0,7.0,7.0,7.0,7.0,7.0,7.0]
SysLib.convertToLpf(Parameter[0], System)

Parameter[0]['rxApod'] = SysLib.init_object("Parameter.rxApod", 1)
Parameter[0]['rxApod'][0]['WindowFunc_FuncName'] = ""#"Hamming"
SysLib.convertToRxApod(Parameter[0], 128)

Parameter[0]['speedOfSoundMps'] = 1540

#####################################################
# Specify ROI structure
#####################################################

roiStartElement = 0
roiEndElement = System['Transducer']['elementCnt']-1

Roi = SysLib.init_object("Roi", 1)
Roi[0] = SysLib.computeRoi("B", "2D", System['Transducer'], [roiStartElement, roiEndElement], [0, depthCm], steerAngleA)
Roi[0]['bfProc'] = 'PWI_SPCD'

#####################################################
# Specify High Voltage structure
#####################################################
Hv = SysLib.init_object("Hv", 1)
Hv[0]['highVoltage'] = 40

#####################################################
# Specify Transmit waveform structure
#####################################################
Tw = SysLib.init_object("Tw",1)
Tw[0]['waveformType'] = "Parameter"
Tw[0]['freqMHz'] = center_frequency;
Tw[0]['burstCount'] = 1
Tw[0]['pulseLengthDuty'] = 50
Tw[0]['hvNo'] = 0


#####################################################
# Specify TX structure array
#####################################################
Tx = SysLib.init_object("Tx",spcdCnt)
for i in range(0, spcdCnt):
    Tx[i]['twNo'] = 0
    Tx[i]['steerAngleDeg'] = steerAngleA[i]
    Tx[i]['focalPosCm'] = 0
    Tx[i]['fNum'] = 0
    Tx[i]['apod'] = ones(transNumElement)
    Tx[i]['delay'] = SysLib.computeTXDelay(Tx[i],System,Parameter[0])
    Tx[i]['apertureCtrlF'] = 1
    Tx[i]['aperture'] = zeros(transNumElement)

#####################################################
# Specify Rx buffer
#####################################################
index = 0

RcvBuffer = SysLib.init_object("RcvBuffer", acqNum*numFrames)
for i in range(0,numFrames,1):
    for j in range(0,acqNum,1):
        RcvBuffer[index]['index'] = index
        RcvBuffer[index]['acqNo'] = j
        RcvBuffer[index]['frameNo'] = i
        index = index + 1

#####################################################
# Specify Tgc
#####################################################
Tgc = SysLib.init_object("Tgc", 1)
Tgc[0]['atgc'] = "" # SysLib.convertToAtgc(Tgc[0], array([128,192,228,256,256,256,256,256]))

#####################################################
# Specify Rx controls
#####################################################
index = 0
Rx = SysLib.init_object("Rx",acqNum)
for i in range(0, spcdCnt, 1):
    for j in range(0, acqNum/spcdCnt, 1):
        Rx[index]['acqNum'] = acqNum
        Rx[index]['tgcNo'] = 0
        Rx[index]['numFrames'] = numFrames
        Rx[index]['originElement'] = 32 + (64*j)
        Rx[index]['aperture'] = SysLib.computeRXAperture(Rx[index],System['Transducer'])
        Rx[index]['roiNo'] = 0
        Rx[index]['angleRoiNo'] = i
        index = index + 1        

###################################
# specify Sequence structure array 
###################################
seqControlCnt = 2
SeqControl = SysLib.init_object("SeqControl", seqControlCnt)

SeqControl[0]['waitTimeToNext']	= 0
SeqControl[0]['jump']			= -1
SeqControl[0]['loopCount']		= 0
SeqControl[0]['stop']			= 0
SeqControl[0]['pause']			= 0
SeqControl[0]['triggerOut']	= 0
SeqControl[0]['display']		= 0
SeqControl[1]['waitTimeToNext']	= 0
SeqControl[1]['jump']			= -1
SeqControl[1]['loopCount']		= 0
SeqControl[1]['stop']			= 0
SeqControl[1]['pause']			= 0
SeqControl[1]['triggerOut']	= 0
SeqControl[1]['display']		= 1

DisplayIntv = 10

SeqItem = SysLib.init_object("SeqItem", acqNum*numFrames )

index = 0 # SeqItem number
for i in range(0,numFrames,1):
    for j in range(0,spcdCnt,1):
        for k in range(0, acqNum/spcdCnt, 1):        
            SeqItem[index]['txNo'] = j
            SeqItem[index]['rxNo'] = j*(acqNum/spcdCnt) + k
            SeqItem[index]['bufNo'] = index
            SeqItem[index]['ctrlNo'] = (0, 1)[(i+1)%DisplayIntv == 0 and (j+1)==spcdCnt and (k+1)==acqNum/spcdCnt]
            index = index + 1

    
###################################
# specify Layout structure array 
###################################
Layout = SysLib.init_object("Layout", 1)
Layout[0]['filter'] = FilterLib.createFilterGraph('Unfocused_B')


savedict = dict( 
                SeqItem = SeqItem, 
                RcvBuffer = RcvBuffer, 
                Parameter = Parameter,
                Roi = Roi,
                Hv = Hv,
                Tw = Tw, 
                Tx = Tx,
                System = System,
                Tgc = Tgc,
                Rx = Rx,
                SeqControl = SeqControl,
                Layout = Layout
                )

SysLib.saveMat("L3-8_OMHR_PW_"+str(spcdCnt),savedict)
