# -*- coding: utf-8 -*-
# Modification of original script by Alpinio script by Alpinion 
# @author: Ole Marius Hoel Rindal


import sys
if sys.platform == "win64" or sys.platform == "win32": # If windows we are on the Alpinion scanner
	sys.path.append("D:/Cosmos/v1/RND/SEQ_Data/Lib")
elif sys.platform == "darwin": #If mac we are on my laptop
	sys.path.append("/Users/omrindal/Development/Alpinion/Lib")

from numpy import *
import numpy as np
import TransLib
import SysLib
import FilterLib
import pdb

transNumElement = 128	#Number of transmit elements (is it really 128, though??)
numFrames = 40		#Number of frames to save
acqNum = 128	#Number of scanlines, 384 is maximum
transmit_focus = 3     #in centimeters
center_frequency = 6   #in MHz
depthCm = 6

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
Parameter[0]['dtgc'][0]['dtgcDb'] = [0,0,0,0,0,0,0,0]
SysLib.convertToDtgc(Parameter[0])

Parameter[0]['demod'] = SysLib.init_object("Parameter.demod", 1)
Parameter[0]['demod'][0]['demodMHz'] = [10.5, 10.5, 10.5, 10.0, 10.0, 10.0, 9.5, 9.5]
SysLib.convertToDemod(Parameter[0], System)

Parameter[0]['lpf'] = SysLib.init_object("Parameter.lpf", 1)
Parameter[0]['lpf'][0]['lpfMHz'] = [7.0,7.0,7.0,7.0,7.0,7.0,7.0,7.0]
SysLib.convertToLpf(Parameter[0], System)

Parameter[0]['rxApod'] = SysLib.init_object("Parameter.rxApod", 1)
Parameter[0]['rxApod'][0]['WindowFunc_FuncName'] = "Hamming"
SysLib.convertToRxApod(Parameter[0], 128)

Parameter[0]['speedOfSoundMps'] = 1540

#####################################################
# Specify ROI structure
#####################################################

roiStartElement = 0
roiEndElement = System['Transducer']['elementCnt']-1

Roi = SysLib.init_object("Roi", 1)
Roi[0] = SysLib.computeRoi("B", "2D", System['Transducer'], [roiStartElement, roiEndElement], [0, depthCm])
Roi[0]['maxScanlineNum'] = acqNum;
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
Tw[0]['freqMHz'] = center_frequency
Tw[0]['burstCount'] = 1
Tw[0]['pulseLengthDuty'] = 50
Tw[0]['hvNo'] = 0

#####################################################
# Specify TX structure array
#####################################################
elementCnt = System['Transducer']['elementCnt']
elementPitchCm = System['Transducer']['elementPitchCm']

Tx = SysLib.init_object("Tx",acqNum)
for i in range(0,acqNum,1):
    Tx[i]['twNo'] = 0
    Tx[i]['steerAngleDeg'] = 0
    Tx[i]['focalPosCm'] = transmit_focus
    Tx[i]['fNum'] = 5.0
    Tx[i]['apod'] = ones(transNumElement)
    Tx[i]['originPosCm'] = np.array([((-elementCnt / 2) + float(elementCnt)/acqNum*i) * elementPitchCm, 0.0, 0.0])
    Tx[i]['originElement'] = float(elementCnt) / acqNum * i
    Tx[i]['delay'] = SysLib.computeTXDelay(Tx[i],System,Parameter[0])  
    Tx[i]['aperture'] = SysLib.computeTXAperture(Tx[i],System['Transducer'])


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
Rx = SysLib.init_object("Rx",acqNum)
for i in range(0,acqNum,1):
    Rx[i]['acqNum'] = acqNum
    Rx[i]['tgcNo'] = 0
    Rx[i]['numFrames'] = numFrames
    Rx[i]['roiNo'] = 0

###################################
# specify Sequence structure array 
###################################
seqControlCnt = 2
SeqControl = SysLib.init_object("SeqControl", seqControlCnt)

SeqControl[0]['waitTimeToNext']	= 0
SeqControl[0]['jump']			= -1
SeqControl[0]['loopCount']		     = 0
SeqControl[0]['stop']			= 0
SeqControl[0]['pause']			= 0
SeqControl[0]['triggerOut']		= 0
SeqControl[1]['waitTimeToNext']	= 0
SeqControl[1]['jump']			= -1
SeqControl[1]['loopCount']		= 0
SeqControl[1]['stop']			= 0
SeqControl[1]['pause']			= 0
SeqControl[1]['triggerOut']		= 0


SeqItem = SysLib.init_object("SeqItem", acqNum*numFrames )

index = 0 # SeqItem number
for i in range(0,numFrames,1):
    for j in range(0,acqNum,1):
        SeqItem[index]['txNo'] = j
        SeqItem[index]['rxNo'] = j
        SeqItem[index]['bufNo'] = acqNum*i+j  
        SeqItem[index]['ctrlNo'] = 0   
        index = index + 1

    
###################################
# specify Layout structure array 
###################################
Layout = SysLib.init_object("Layout", 1)
Layout[0]['filter'] = FilterLib.createFilterGraph('Focused_B')


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

#pdb.set_trace(); #set debug point
                
SysLib.saveMat("L3-8_OMHR_focus_at_"+str(transmit_focus)+"cm_"+str(acqNum)+"_scanLines",savedict)
