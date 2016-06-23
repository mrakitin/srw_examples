#############################################################################
# SRWLIB Example#6: Calculating spectral flux of undulator radiation 
# by finite-emittance electron beam collected through a finite aperture
# and power density distribution of this radiation (integrated over all photon energies)
# v 0.02
#############################################################################

from __future__ import print_function  # Python 2.7 compatibility
from srwlib import *
import os
import sys
import copy



# Read data comumns from ASCII file:
def AuxReadInDataColumns(filePath, nCol, strSep):
    return srwl_uti_read_data_cols(filePath, strSep)


def AuxTransmAddSurfHeightProfile(optSlopeErr, heightProfData, dim, ang):
    argHeightProfData = heightProfData[0]
    valHeightProfData = heightProfData[1]
    sinAng = sin(ang)
    npData = len(heightProfData[0])

    # xStep = optSlopeErr.rx/(optSlopeErr.nx - 1)
    # yStep = optSlopeErr.ry/(optSlopeErr.ny - 1)
    # y = optSlopeErr.y - 0.5*optSlopeErr.ry

    auxMesh = optSlopeErr.mesh
    xStep = (auxMesh.xFin - auxMesh.xStart) / (auxMesh.nx - 1)
    yStep = (auxMesh.yFin - auxMesh.yStart) / (auxMesh.ny - 1)

    y = auxMesh.yStart
    hApprox = 0
    ipStart = 0
    # for iy in range(optSlopeErr.ny):
    for iy in range(auxMesh.ny):
        if ('y' in dim):
            hApprox = 0
            y1 = argHeightProfData[ipStart] * sinAng
            for i in range(ipStart + 1, npData):
                y2 = argHeightProfData[i] * sinAng
                if ((y1 <= y) and (y < y2)):
                    hApprox = ((valHeightProfData[i] - valHeightProfData[i - 1]) / (
                        (argHeightProfData[i] - argHeightProfData[i - 1]) * sinAng)) * (y - y1) + valHeightProfData[i - 1]
                    # sys.stdout.write(ipStart, i, iy, y1, y, y2, argHeightProfData[i-1], argHeightProfData[i], valHeightProfData[i-1], valHeightProfData[i], hApprox)
                    ipStart = i - 1
                    break
                y1 = y2

        # x = optSlopeErr.x - 0.5*optSlopeErr.rx
        x = auxMesh.xStart

        # for ix in range(optSlopeErr.nx):
        for ix in range(auxMesh.nx):
            if ('x' in dim):
                if (ix == 0): ipStart = 0
                hApprox = 0
                x1 = argHeightProfData[ipStart] * sinAng
                for i in range(ipStart + 1, npData):
                    x2 = argHeightProfData[i] * sinAng
                    if ((x1 <= x) and (x < x2)):
                        hApprox = ((valHeightProfData[i] - valHeightProfData[i - 1]) / (
                            (argHeightProfData[i] - argHeightProfData[i - 1]) * sinAng)) * (x - x1) + valHeightProfData[
                                      i - 1]
                        ipStart = i - 1
                        break
                    x1 = x2
            # ofst = 2*ix + (2*optSlopeErr.nx)*iy
            ofst = 2 * ix + (2 * auxMesh.nx) * iy

            optSlopeErr.arTr[ofst] = 1.  # Amplitude Transmission
            optSlopeErr.arTr[ofst + 1] = 0.  # Optical Path Difference
            if (hApprox != 0):
                optSlopeErr.arTr[ofst + 1] = -2 * sinAng * hApprox  # Optical Path Difference (to check sign!)
                # sys.stdout.write(ix, iy, optSlopeErr.arTr[ofst + 1])
            x += xStep
        y += yStep


def drift_ebeam(ebeam, dist):
    ebeam.partStatMom1.z += dist
    ebeam.arStatMom2[0] += ebeam.arStatMom2[1] * dist * 2 + ebeam.arStatMom2[2] * dist * dist
    ebeam.arStatMom2[1] += ebeam.arStatMom2[2] * dist
    ebeam.arStatMom2[3] += ebeam.arStatMom2[4] * dist * 2 + ebeam.arStatMom2[5] * dist * dist
    ebeam.arStatMom2[4] += ebeam.arStatMom2[5] * dist
    return ebeam


beamline = 'ES2'  # 'ES1' 'ES2'
BMmode = 'Norm'  # 'Norm'  'LowDiv'
bump = True  # False, True

strBump = '_bump' if bump else ''

print('modeling SMI beamline ' + beamline + ' with bump = ' + repr(bump) + ' BMmode = ' + BMmode)

print('Calculating spectral flux of undulator radiation by finite-emittance electron beam collected through a finite aperture')

# **********************Input Parameters:
strExDataFolderName = 'smi204crlb'  # example data sub-folder name
strIntSourSE_OutFileName = os.path.join(os.getcwd(), strExDataFolderName,
                                        'sour_se' + beamline + strBump + '.dat')  # file name for output UR intensity Single Electron
strIntPropSE_OutFileName = os.path.join(os.getcwd(), strExDataFolderName,
                                        'prop_se' + beamline + strBump + '.dat')  # file name for output intensity propagated Single Electron
strIntPropME_OutFileName = os.path.join(os.getcwd(), strExDataFolderName,
                                        'prop_me' + beamline + strBump + '.dat')  # file name for output intensity propagated Multi Electron
# strTrajOutFileName = "res_trj.dat"
strProfileData = os.path.join(os.getcwd(), strExDataFolderName, 'Si_heat204.dat')
strProfileDataHFM = os.path.join(os.getcwd(), strExDataFolderName, 'HFM_SESO.dat')  # HFM_SESO.dat - 803 points 0.5 m
# strProfileDataHFM        = os.path.join(os.getcwd(),strExDataFolderName,'HFM08rms.dat') #HFM08rms.dat - 464 points 0.5 m
strProfileDataVFM = os.path.join(os.getcwd(), strExDataFolderName, 'VFM_SESO.dat')  # VFM_SESO.dat - 288 points 0.4 m
# strProfileDataVFM        = os.path.join(os.getcwd(),strExDataFolderName,'VFM03rms.dat') #VFM03rms.dat - 464 points 0.5 m
strProfileDataVM = os.path.join(os.getcwd(), strExDataFolderName,
                                'VM03rms.dat')  # VM03rmsL.dat - 127 points, 0.5 m; #VM03rms.dat - 500 points 0.5 m
strProfileDataHKB = os.path.join(os.getcwd(), strExDataFolderName,
                                 'FMsine015.dat')  # FMsine015.dat - 127 points 0.5166 m
strProfileDataVKB = os.path.join(os.getcwd(), strExDataFolderName, 'VM03rms.dat')  # VM03rms.dat - 500 points 0.5 m

# ***********Undulator
numPer = 121.5  # Number of ID Periods (without counting for terminations
undPer = 0.023  # Period Length [m]
By = 0.955  # 0.577 #0.955 #Peak Vertical field [T]
phBy = 0  # Initial Phase of the Vertical field component
sBy = -1  # Symmetry of the Vertical field component vs Longitudinal position
xcID = 0  # Transverse Coordinates of Undulator Center [m]
ycID = 0
zcID = 0.6  # 0 #Longitudinal Coordinate of Undulator Center [m]

sys.stdout.write('   Setup Magnetic Field for Undulator ... ')
sys.stdout.flush()
und = SRWLMagFldU([SRWLMagFldH(1, 'v', By, phBy, sBy, 1)], undPer, numPer)  # Planar Undulator
magFldCnt = SRWLMagFldC([und], array('d', [xcID]), array('d', [ycID]), array('d', [zcID]))  # Container of all Field Elements
sys.stdout.write('done\n')

# ***********Electron Beam
sys.stdout.write('   Setup Electron Beam for Undulator ... ')
sys.stdout.flush()
elecBeam = SRWLPartBeam()
elecBeam.Iavg = 0.5  # Average Current [A]
elecBeam.partStatMom1.x = 0.  # 200e-06 #0. #Initial Transverse Coordinates (initial Longitudinal Coordinate will be defined later on) [m]
elecBeam.partStatMom1.y = 0.  # 30e-06 #0. #-0.00025
elecBeam.partStatMom1.z = -0.9  # 0. #-0.5*undPer*(numPer + 4) #Initial Longitudinal Coordinate (set before the ID)
elecBeam.partStatMom1.xp = 0.  # 10e-06 #0. #Initial Relative Transverse Velocities
elecBeam.partStatMom1.yp = 0.
elecBeam.partStatMom1.gamma = 3. / 0.51099890221e-03  # Relative Energy
# 2nd order statistical moments
elecBeam.arStatMom2[0] = (137.113e-06) ** 2  # <(x-x0)^2>
elecBeam.arStatMom2[1] = -0.0388489e-09  # <(x-x0)*(x'-x'0)>
elecBeam.arStatMom2[2] = (6.57004e-06) ** 2  # <(x'-x'0)^2>

elecBeam.arStatMom2[3] = (5.39499e-06) ** 2  # (15.4091e-06)**2 #<(y-y0)^2>
elecBeam.arStatMom2[4] = -0.00211765e-09  # <(y-y0)*(y'-y'0)>
elecBeam.arStatMom2[5] = (1.53393e-06) ** 2  # <(y'-y'0)^2>
elecBeam.arStatMom2[10] = (0.89e-03) ** 2  # <(E-E0)^2>/E0^2

# elecBeam = drift_ebeam(elecBeam,-0.5)
sys.stdout.write('done\n')

meth = 1  # SR calculation method: 0- "manual", 1- "auto-undulator", 2- "auto-wiggler"
relPrec = 0.008  # 0.001 #relative precision
zStartInteg = 0  # longitudinal position to start integration (effective if < zEndInteg)
zEndInteg = 0  # longitudinal position to finish integration (effective if > zStartInteg)
npTraj = 50000  # Number of points for trajectory calculation
useTermin = 1  # Use "terminating terms" (i.e. asymptotic expansions at zStartInteg and zEndInteg) or not (1 or 0 respectively)
sampFactNxNyForProp = 0.0  # 1.0 #0.2 #0.1 #sampling factor for adjusting nx, ny (effective if > 0) ###TODO fix that
arPrecParSpec = [meth, relPrec, zStartInteg, zEndInteg, npTraj, useTermin, sampFactNxNyForProp]

# ***********UR Wavefront Parameters (mesh)
sys.stdout.write('   Setup Stokes mesh ... ');
sys.stdout.flush()
wfr = SRWLWfr()  # for spectral flux vs photon energy
wfr.allocate(1, 81, 81)  # numbers of points vs photon energy, horizontal and vertical positions
wfr.mesh.zStart = 29.5  # longitudinal position [m] at which UR has to be calculated
wfr.mesh.eStart = 20358.  # initial photon energy [eV]
wfr.mesh.eFin = 20358.  # final photon energy [eV]wfr.mesh.xStart = -0.0006 #initial horizontal position [m]
wfr.mesh.xStart = -0.0006  # initial horizontal position [m]
wfr.mesh.xFin = 0.0006  # final horizontal position [m]
wfr.mesh.yStart = -0.0006  # initial vertical position [m]
wfr.mesh.yFin = 0.0006  # final vertical position [m]
wfr.partBeam = elecBeam
meshInitPartCoh = deepcopy(wfr.mesh)
sys.stdout.write('done\n')

# ========================used for 1D CRLs======================================
deltaV = 8.21692879E-07  # 1.30747685E-06 #Al # 8.21692879E-07 #Be @ 20.4KeV
attenLenV = 28544.7e-06  # 1231.33e-6 # Al # 28544.7e-06 #[m] #Be @20.4KeV
diamCRLV = 1.e-03  # CRL diameter
rMinCRLV = 100e-06  # CRL radius at the tip of parabola [m]
nCRLV = 38  # number of lenses
wallThickCRLV = 50e-06  # CRL wall thickness [m]

deltaH = 8.21692879E-07  # 1.30747685E-06 #Al # 8.21692879E-07 #Be @ 20.4KeV
attenLenH = 28544.7e-06  # 1231.33e-6 # Al # 28544.7e-06 #[m] #Be @20.4KeV
diamCRLH = 1.e-03  # CRL diameter
rMinCRLH = 100e-06  # CRL radius at the tip of parabola [m]
nCRLH = 46  # number of lenses
wallThickCRLH = 50e-06  # CRL wall thickness [m]
# Generating a perfect 2D parabolic CRL:
#   
# CRLV = srwl_opt_setup_CRL(2, deltaV, attenLenV, 1, diamCRLV, diamCRLV, rMinCRLV, nCRLV, wallThickCRLV, 0, 0)
# CRLH = srwl_opt_setup_CRL(1, deltaH, attenLenH, 1, diamCRLH, diamCRLH, rMinCRLH, nCRLH, wallThickCRLH, 0, 0)
# ==========================end of 1D CRLs======================================

delta = 8.21692879E-07  # Be @ 20.4KeV
attenLen = 28544.7e-06  # [m] #20.4KeV
diamCRL = 1.e-03  # CRL diameter
rMinCRL = 50e-06  # CRL radius at the tip of parabola [m]
nCRL = 23  # number of lenses
wallThickCRL = 32.4e-06  # CRL wall thickness [m]

CRL = srwl_opt_setup_CRL(3, delta, attenLen, 1, diamCRL, diamCRL, rMinCRL, nCRL, wallThickCRL, 0, 0)
# Beamline OEs
# APE='aperture',
# MOAT='first mirror of Monocromator error shape',
# VFML='Vertical Focusing Mirror (Spherical) Lens',
# VFMT='Vertical Focusing Mirror (Spherical) error shape'
D_APE_MOA = SRWLOptD(2.44)

# ================introducing Si(111) heat load==============
heightProf = AuxReadInDataColumns(strProfileData, 2, '\t')
# MOAT        = SRWLOptT(100, 2001, 2.0e-02, 3.0e-02*sin(1.223866));
MOAT = SRWLOptT(100, 500, 2.0e-02, 16e-3 * sin(0.09727))
opMOAT = srwl_opt_setup_surf_height_2d(heightProf, 'y', _ang=0.09727, _nx=100, _ny=500, _size_x=2.0e-02, _size_y=16e-3 * sin(0.09727))
AuxTransmAddSurfHeightProfile(MOAT, heightProf, 'y', 0.09727)  # incident angle is 70.122373 deg => 1.223866 rad
opPathDifMOAT = MOAT.get_data(3, 3)
srwl_uti_save_intens_ascii(opPathDifMOAT, MOAT.mesh, os.path.join(os.getcwd(), strExDataFolderName, 'res_er_mono.dat'),
                           0, ['', 'Horizontal Position', 'Vertical Position', 'Opt. Path Diff.'],
                           _arUnits=['', 'm', 'm', 'm'])
# ================finished Si (111) heat load================

D_MOA_HFM = SRWLOptD(2.94244)
D_APE_HFM = SRWLOptD(2.44 + 2.94244)  # used for "no bump" option
if BMmode == 'Norm':
    HFML = SRWLOptL(_Fx=1. / (1. / (29.5 + 2.44 + 2.94244) + 1. / (3.42 + 8.7 + 3.9)))  # to focus at ES1
if BMmode == 'LowDiv':
    HFML = SRWLOptL(_Fx=1. / (1. / (29.5 + 2.44 + 2.94244) + 1. / (3.42 + 8.7 + 3.9 + 8.1 - 0.3)))  # to focus at ES2 with a low divergence

# ================introducing HFM slope error==============
heightProfHFM = AuxReadInDataColumns(strProfileDataHFM, 2, '\t')
HFMT = SRWLOptT(803, 200, 0.5 * sin(3.1415927e-03), 6.0e-03)
# HFMT        = SRWLOptT(464, 200, 0.5*sin(3.1415927e-03), 6.0e-03) #sinusoidal 0.1, 1.8e-08 both 'h' 'v', angle 3.1415927e-03 rad to correct for horizontal.
AuxTransmAddSurfHeightProfile(HFMT, heightProfHFM, 'x', 3.1415927e-03)
opPathDifHFMT = HFMT.get_data(3, 3)
srwl_uti_save_intens_ascii(opPathDifHFMT, HFMT.mesh, os.path.join(os.getcwd(), strExDataFolderName, 'res_er_HFM.dat'),
                           0, ['', 'Horizontal Position', 'Vertical Position', 'Opt. Path Diff.'],
                           _arUnits=['', 'm', 'm', 'm'])
# ================finished HFM slope error=================
# HFMT = SRWLOptT()
D_HFM_VFM = SRWLOptD(3.42)
if BMmode == 'Norm':
    VFML = SRWLOptL(_Fy=1. / (1. / (29.5 + 2.44 + 2.94244 + 3.42 - 0.6) + 1. / (8.7 + 3.9 + 0.3)))  # focus at ES1; if using Bump, VFM must be 3.9+0.3 m (to compensate bump which moves focus 0.2 m upstream)
if BMmode == 'LowDiv':
    VFML = SRWLOptL(_Fy=1. / (1. / (29.5 + 2.44 + 2.94244 + 3.42 - 0.6) + 1. / (8.7 + 3.9 - 5.7 + 8.1)))  # focus at ES2 with a low divergence

# ================introducing VFM slope error==============
heightProfVFM = AuxReadInDataColumns(strProfileDataVFM, 2, '\t')
VFMT = SRWLOptT(200, 288, 6.0e-03, 0.4 * sin(
    3.1415927e-03))  # sinusoidal equal to HFM. the origina spec is 0.1, 6.75e-09 both 'h' 'v', angle 6.1086524e-03 rad to correct for vertical.
AuxTransmAddSurfHeightProfile(VFMT, heightProfVFM, 'y', 3.1415927e-03)
opPathDifVFMT = VFMT.get_data(3, 3)
srwl_uti_save_intens_ascii(opPathDifVFMT, VFMT.mesh, os.path.join(os.getcwd(), strExDataFolderName, 'res_er_VFM.dat'),
                           0, ['', 'Horizontal Position', 'Vertical Position', 'Opt. Path Diff.'],
                           _arUnits=['', 'm', 'm', 'm'])
# ================finished VFM slope error=================
#
# ================introducing VM slope error==============
heightProfVM = AuxReadInDataColumns(strProfileDataVM, 2, '\t')
VMT = SRWLOptT(200, 500, 6.0e-03, 0.5 * sin(
    3.1415927e-03))  # sinusoidal equal to HFM. the origina spec is 0.1, 6.75e-09 both 'h' 'v', angle 6.1086524e-03 rad to correct for vertical.
AuxTransmAddSurfHeightProfile(VMT, heightProfVM, 'y', 3.1415927e-03)
opPathDifVMT = VMT.get_data(3, 3)
srwl_uti_save_intens_ascii(opPathDifVMT, VMT.mesh, os.path.join(os.getcwd(), strExDataFolderName, 'res_er_VM.dat'), 0,
                           ['', 'Horizontal Position', 'Vertical Position', 'Opt. Path Diff.'],
                           _arUnits=['', 'm', 'm', 'm'])
# ================finished VFM slope error=================
# VFMT = SRWLOptT()
# D_VFM_SSA   = SRWLOptD(8.7+3.9) #TODO Where is the SSA?
D_VFM_VM = SRWLOptD(0.7)
D_VM_SSA = SRWLOptD(8.0)
D_VFM_SSA = SRWLOptD(8.7)
if beamline == 'ES1' and BMmode == 'Norm':
    SSA = SRWLOptA('r', 'a', 0.4e-03, 0.4e-03)  # 0.4, 0.4 for NOT low divergence mode;

if beamline == 'ES2' and BMmode == 'Norm':
    SSA = SRWLOptA('r', 'a', 0.9e-03, 0.9e-03)  # 0.4x0.15 for KB

if beamline == 'ES2' and BMmode == 'LowDiv':
    SSA = SRWLOptA('r', 'a', 0.9e-03, 0.9e-03)  # 0.4, 0.4 for low divergence mode;

D_SSA_ES1 = SRWLOptD(3.9)  # TODO Where is the SSA?

# D_SSA_CRLV   = SRWLOptD(9.8903) #D_SSA_CRL   = SRWLOptD(9.9)
# D_CRLV_CRLH  = SRWLOptD(0.4918) #D_CRL_ES2   = SRWLOptD(2.1)
# D_CRLH_ES2   = SRWLOptD(1.6178)

D_SSA_CRL = SRWLOptD(10.33492)
D_CRL_ES2 = SRWLOptD(1.66508)

# ApCRLV = SRWLOptA('c','a',1.0e-3)
# ApCRLH = SRWLOptA('c','a',1.0e-3)

ApCRL = SRWLOptA('c', 'a', 1.0e-3)
# CRL = SRWLOptL(_Fx=1./(1./(6.0)+1./(2.1)), _Fy=1./(1./(6.0)+1./(2.1)))

# angHKB = 3.14e-03 #[rad]
# HKB = SRWLOptMirEl(_p=6.7, _q=1.40, _ang_graz=angHKB, _r_sag=1.e+40, _size_tang=0.3, _nvx=cos(angHKB), _nvy=0, _nvz=-sin(angHKB), _tvx=-sin(angHKB), _tvy=0, _x=0, _y=0, _treat_in_out=1) #HKB Ellipsoidal Mirror
# ================introducing HKB slope error==============
# heightProfHKB  = AuxReadInDataColumns(strProfileDataHKB, 2, '\t')
# HKBT        = SRWLOptT(127, 200, 3.0e-01*sin(3.1415927e-03), 6.0e-03) #sinusoidal 0.1, 1.8e-08 both 'h' 'v', angle 3.1415927e-03 rad to correct for horizontal.
# AuxTransmAddSurfHeightProfile(HKBT, heightProfHKB, 'x', 3.1415927e-03)
# opPathDifHKBT = HKBT.get_data(3, 3)
# srwl_uti_save_intens_ascii(opPathDifHKBT, HKBT.mesh, os.path.join(os.getcwd(), strExDataFolderName, 'res_er_HKB.dat'), 0, ['', 'Horizontal Position', 'Vertical Position', 'Opt. Path Diff.'], _arUnits=['', 'm', 'm', 'm'])
# ================finished HKB slope error=================
# angVKB = 3.14e-03 #grazing angle at VKB center [rad]
# VKB = SRWLOptMirEl(_p=7.1, _q=1., _ang_graz=angVKB, _r_sag=1.e+40, _size_tang=0.35, _nvx=0, _nvy=cos(angVKB), _nvz=-sin(angVKB), _tvx=0, _tvy=-sin(angVKB), _x=0, _y=0, _treat_in_out=1) #VKB Ellipsoidal Mirror
# ================introducing VKB slope error==============
# heightProfVKB = AuxReadInDataColumns(strProfileDataVKB, 2, '\t')
# VKBT        = SRWLOptT(200, 500, 6.0e-03, 0.5*sin(3.1415927e-03))#sinusoidal equal to HFM. the origina spec is 0.1, 6.75e-09 both 'h' 'v', angle 6.1086524e-03 rad to correct for vertical.
# AuxTransmAddSurfHeightProfile(VKBT, heightProfVKB, 'y', 3.1415927e-03)
# opPathDifVKBT = VKBT.get_data(3, 3)
# srwl_uti_save_intens_ascii(opPathDifVKBT, VKBT.mesh, os.path.join(os.getcwd(), strExDataFolderName, 'res_er_VKB.dat'), 0, ['', 'Horizontal Position', 'Vertical Position', 'Opt. Path Diff.'], _arUnits=['', 'm', 'm', 'm'])
# ================finished VKB slope error=================

# D_HKB_VKB = SRWLOptD(0.4) #Distance between centers of Vertically and Horizontally focusing K-B mirrors
# D_VKB_ES2 = SRWLOptD(1.0)
# D_ES1_ES2 = SRWLOptD(8.1)

#             [ 0] [1] [2]  [3] [4] [5]  [6]  [7]  [8]  [9] [10] [11] 
if beamline == 'ES1':
    # ppD_APE_MOA = [ 0,  0, 1.0,  2,  0, 2.0, 6.0, 2.0, 6.0,  0,  0,  0 ]   #settings for single electron emission
    ppD_APE_MOA = [0, 0, 1.0, 2, 0, 8.0, 3.0, 2.0, 3.0, 0, 0, 0]  # settings for single electron emission

    # ppD_APE_MOA = [ 0,  0, 1.0,  2,  0, 8.0, 6.0, 2.0, 10.0,  0,  0,  0 ]  #settings for multi-electron emission
if beamline == 'ES2':
    ppD_APE_MOA = [0, 0, 1.0, 2, 0, 8.0, 3.0, 3.0, 4.0, 0, 0, 0]  # settings for single electron emission
    # ppD_APE_MOA = [ 0,  0, 1.0,  2,  0, 8.0, 2.0, 2.0, 8.0,  0,  0,  0 ]   #settings for multi-electron emission

ppMOAT = [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0]
ppD_MOA_HFM = [0, 0, 1.0, 2, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0]
ppD_APE_HFM = ppD_APE_MOA
ppHFML = [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0]
ppHFMT = [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0]
ppD_HFM_VFM = [0, 0, 1.0, 2, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0]
ppVFML = [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0]
ppVFMT = [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0]
ppVMT = [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0]
ppD_VFM_SSA = [0, 0, 1.0, 2, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0]
ppD_VFM_VM = [0, 0, 1.0, 2, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0]
ppD_VM_SSA = [0, 0, 1.0, 2, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0]
ppSSA = [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0]
ppD_SSA_CRLV = [0, 0, 1.0, 2, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0]
ppD_CRLV_CRLH = [0, 0, 1.0, 2, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0]
ppApCRLH = [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0]
ppCRLH = [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0]
ppApCRLV = [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0]
ppCRLV = [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0]
ppD_CRLH_ES2 = [0, 0, 1.0, 2, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0]
ppD_SSA_CRL = [0, 0, 1.0, 2, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0]
ppApCRL = [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0]
ppCRL = [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0]
ppD_CRL_ES2 = [0, 0, 1.0, 2, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0]
ppD_SSA_ES1 = [0, 0, 1.0, 2, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0]
ppD_ES1_ES2 = [0, 0, 1.0, 2, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0]
# downstream of ES1
ppD_ES1_HKB = [0, 0, 1.0, 2, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0]
ppD_SSA_HKB = [0, 0, 1.0, 2, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0]
ppHKB = [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0]
ppHKBT = [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0]
ppD_HKB_VKB = [0, 0, 1.0, 2, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0]
ppVKB = [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0]
ppVKBT = [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0]
ppD_VKB_ES2 = [0, 0, 1.0, 2, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0]

# [ 0]: Auto-Resize (1) or not (0) Before propagation
# [ 1]: Auto-Resize (1) or not (0) After propagation
# [ 2]: Relative Precision for propagation with Auto-Resizing (1. is nominal)
# [ 3]: Allow (1) or not (0) for semi-analytical treatment of the quadratic (leading) phase terms at the propagation (2) semi-analytical treatment + autoresizing
# [ 4]: Do any Resizing on Fourier side, using FFT, (1) or not (0)
# [ 5]: Horizontal Range modification factor at Resizing (1. means no modification)
# [ 6]: Horizontal Resolution modification factor at Resizing
# [ 7]: Vertical Range modification factor at Resizing
# [ 8]: Vertical Resolution modification factor at Resizing
# [ 9]: Type of wavefront Shift before Resizing (not yet implemented)
# [10]: New Horizontal wavefront Center position after Shift (not yet implemented)
# [11]: New Vertical wavefront Center position after Shift (not yet implemented)

# with transission for error shape
# BML         = SRWLOptC( [   D_APE_MOA,   MOAT,   D_MOA_HFM,   HFML,   HFMT,   D_HFM_VFM,   VFML,   VFMT,   D_VFM_SSA ],
#                        [ ppD_APE_MOA, ppMOAT, ppD_MOA_HFM, ppHFML, ppHFMT, ppD_HFM_VFM, ppVFML, ppVFMT, ppD_VFM_SSA ] )
# without error shape
if beamline == 'ES1' and not bump:
    BML = SRWLOptC([D_APE_HFM, HFML, D_HFM_VFM, VFML, D_VFM_SSA, SSA, D_SSA_ES1],
                   [ppD_APE_HFM, ppHFML, ppD_HFM_VFM, ppVFML, ppD_VFM_SSA, ppSSA, ppD_SSA_ES1])
# if beamline=='ES1' and bump:
#   BML     = SRWLOptC( [   D_APE_MOA,   MOAT,   D_MOA_HFM,   HFML,   D_HFM_VFM,   VFML,   D_VFM_SSA,   SSA,   D_SSA_ES1 ],
#                      [ ppD_APE_MOA, ppMOAT, ppD_MOA_HFM, ppHFML, ppD_HFM_VFM, ppVFML, ppD_VFM_SSA, ppSSA, ppD_SSA_ES1 ])
if beamline == 'ES1' and bump and BMmode == 'Norm':
    # BML     = SRWLOptC( [   D_APE_MOA,   MOAT,   D_MOA_HFM,   HFML,   HFMT,   D_HFM_VFM,   VFML,   D_VFM_SSA,   SSA,   D_SSA_ES1 ],
    #                   [ ppD_APE_MOA, ppMOAT, ppD_MOA_HFM, ppHFML, ppHFMT, ppD_HFM_VFM, ppVFML, ppD_VFM_SSA, ppSSA, ppD_SSA_ES1 ])   
    # BML     = SRWLOptC( [   D_APE_MOA,   MOAT,   D_MOA_HFM,   HFML,   HFMT,   D_HFM_VFM,   VFML,   VFMT,   D_VFM_SSA,   SSA,   D_SSA_ES1 ],
    #                    [ ppD_APE_MOA, ppMOAT, ppD_MOA_HFM, ppHFML, ppHFMT, ppD_HFM_VFM, ppVFML, ppVFMT, ppD_VFM_SSA, ppSSA, ppD_SSA_ES1 ])
    BML = SRWLOptC(
        [D_APE_MOA, MOAT, D_MOA_HFM, HFML, HFMT, D_HFM_VFM, VFML, VFMT, D_VFM_VM, VMT, D_VM_SSA, SSA, D_SSA_ES1],
        [ppD_APE_MOA, ppMOAT, ppD_MOA_HFM, ppHFML, ppHFMT, ppD_HFM_VFM, ppVFML, ppVFMT, ppD_VFM_VM, ppVMT, ppD_VM_SSA,
         ppSSA, ppD_SSA_ES1])

# if beamline=='ES2' and not bump:
#   BML     = SRWLOptC( [   D_APE_HFM,   HFML,   HFMT,   D_HFM_VFM,   VFML,   VFMT,   D_VFM_SSA,   SSA,   D_SSA_HKB,   HKB,   D_HKB_VKB,   VKB,   D_VKB_ES2 ],
#                      [ ppD_APE_HFM, ppHFML, ppHFMT, ppD_HFM_VFM, ppVFML, ppVFMT, ppD_VFM_SSA, ppSSA, ppD_SSA_HKB, ppHKB, ppD_HKB_VKB, ppVKB, ppD_VKB_ES2 ] )

# if beamline =='ES2' and bump and BMmode == 'LowDiv': #focus ES2 without Kb with low divergence
#   BML     = SRWLOptC( [   D_APE_MOA,   MOAT,   D_MOA_HFM,   HFML,   HFMT,  D_HFM_VFM,   VFML,    VFMT,   D_VFM_VM,   VMT,   D_VM_SSA,   SSA,   D_SSA_CRL,   ApCRL,   CRL,   D_CRL_ES2 ],
#                      [ ppD_APE_MOA, ppMOAT, ppD_MOA_HFM, ppHFML, ppHFMT, ppD_HFM_VFM, ppVFML, ppVFMT, ppD_VFM_VM, ppVMT, ppD_VM_SSA, ppSSA, ppD_SSA_CRL, ppApCRL, ppCRL, ppD_CRL_ES2 ] )
if beamline == 'ES2' and bump and BMmode == 'LowDiv':  # focus ES2 without Kb with low divergence
    BML = SRWLOptC(
        [D_APE_MOA, MOAT, D_MOA_HFM, HFML, HFMT, D_HFM_VFM, VFML, VFMT, D_VFM_VM, VMT, D_VM_SSA, SSA, D_SSA_CRL, ApCRL,
         CRL, D_CRL_ES2],
        [ppD_APE_MOA, ppMOAT, ppD_MOA_HFM, ppHFML, ppHFMT, ppD_HFM_VFM, ppVFML, ppVFMT, ppD_VFM_VM, ppVMT, ppD_VM_SSA,
         ppSSA, ppD_SSA_CRL, ppApCRL, ppCRL, ppD_CRL_ES2])

if beamline == 'ES2' and bump and BMmode == 'Norm':
    #    BML     = SRWLOptC( [   D_APE_MOA,   MOAT,   D_MOA_HFM,   HFML,   HFMT,   D_HFM_VFM,   VFML,   VFMT,   D_VFM_VM,   VMT,   D_VM_SSA,   SSA,   D_SSA_CRLV,   ApCRLV,   CRLV,   D_CRLV_CRLH,   ApCRLH,   CRLH,   D_CRLH_ES2 ],
    #                       [ ppD_APE_MOA, ppMOAT, ppD_MOA_HFM, ppHFML, ppHFMT, ppD_HFM_VFM, ppVFML, ppVFMT, ppD_VFM_VM, ppVMT, ppD_VM_SSA, ppSSA, ppD_SSA_CRLV, ppApCRLV, ppCRLV, ppD_CRLV_CRLH, ppApCRLH, ppCRLH, ppD_CRLH_ES2 ] )
    BML = SRWLOptC(
        [D_APE_MOA, MOAT, D_MOA_HFM, HFML, HFMT, D_HFM_VFM, VFML, VFMT, D_VFM_VM, VMT, D_VM_SSA, SSA, D_SSA_CRL, ApCRL,
         CRL, D_CRL_ES2],
        [ppD_APE_MOA, ppMOAT, ppD_MOA_HFM, ppHFML, ppHFMT, ppD_HFM_VFM, ppVFML, ppVFMT, ppD_VFM_VM, ppVMT, ppD_VM_SSA,
         ppSSA, ppD_SSA_CRL, ppApCRL, ppCRL, ppD_CRL_ES2])

# CALCULATION

# **********************Calculation (SRWLIB function calls)

sys.stdout.write('   Performing Single Electron calculation ... ');
sys.stdout.flush()
srwl.CalcElecFieldSR(wfr, 0, magFldCnt, arPrecParSpec)
sys.stdout.write('done\n')

sys.stdout.write('   Saving Single Electron UR Intensity ... ');
sys.stdout.flush()
arI = array('f', [0] * wfr.mesh.nx * wfr.mesh.ny)  # "flat" array to take 2D intensity data
srwl.CalcIntFromElecField(arI, wfr, 6, 0, 3, wfr.mesh.eStart, 0, 0)
srwl.CalcIntFromElecField(arI, wfr, 6, 1, 3, wfr.mesh.eStart, 0, 0)

srwl_uti_save_intens_ascii(arI, wfr.mesh, strIntSourSE_OutFileName)
sys.stdout.write('done\n')

# sys.exit(0)

sys.stdout.write('   Performing Single Electron Radiation Propagation ... ');
sys.stdout.flush()
srwl.PropagElecField(wfr, BML)
sys.stdout.write('done\n')

sys.stdout.write('   Saving Single Electron Propagated Intensity ... ');
sys.stdout.flush()
arI = array('f', [0] * wfr.mesh.nx * wfr.mesh.ny)  # "flat" array to take 2D intensity data
srwl.CalcIntFromElecField(arI, wfr, 6, 0, 3, wfr.mesh.eStart, 0, 0)
srwl_uti_save_intens_ascii(arI, wfr.mesh, strIntPropSE_OutFileName)
sys.stdout.write('done\n')

print('Switching from Coordinate to Angular Representation ... ', end='')
srwl.SetRepresElecField(wfr, 'a');
print('done')

print('Extracting Intensity from the Propagated Electric Field in Angular Representation  ... ', end='')
arIa = array('f', [0] * wfr.mesh.nx * wfr.mesh.ny)  # "flat" 2D array to take intensity data
srwl.CalcIntFromElecField(arIa, wfr, 6, 0, 3, wfr.mesh.eStart, 0, 0)

# Converting intensity in angular representation to [photons/s/.1%bw/mrad^2] vs rad
meshAng = copy.deepcopy(wfr.mesh)
wavelength = 1.239842e-06 / meshAng.eStart  # wavelength in [m] - converiosn / multiplier for argument
invSqWavelength = 1 / (wavelength * wavelength)  # converiosn / multiplier for intensity values
meshAng.xStart *= wavelength
meshAng.xFin *= wavelength
meshAng.yStart *= wavelength
meshAng.yFin *= wavelength
for i in range(len(arIa)):
    arIa[i] *= invSqWavelength

srwl_uti_save_intens_ascii(arIa, meshAng, os.path.join(os.getcwd(), strExDataFolderName, "prop_se_ang.dat"))
print('done')

sys.exit(0)

sys.stdout.write('   Starting simulation of Partially-Coherent Wavefront Propagation (takes a lot of time)... ');
sys.stdout.flush()
nMacroElec = 2000000  # Total number of Macro-Electrons (Wavefronts)
nMacroElecAvgOneProc = 5  # Number of Macro-Electrons (Wavefronts) to average on each node (for MPI calculations)
nMacroElecSavePer = 5  # Saving periodicity (in terms of Macro-Electrons) for the Resulting Intensity
# arPrecParSpec[6] = sampFactNxNyForProp
radStokesProp = srwl_wfr_emit_prop_multi_e(elecBeam, magFldCnt, meshInitPartCoh, 1, 0.01, nMacroElec,
                                           nMacroElecAvgOneProc, nMacroElecSavePer,
                                           strIntPropME_OutFileName, sampFactNxNyForProp,
                                           BML)  # BML, 1) for angluar representation
sys.stdout.write('done\n')
