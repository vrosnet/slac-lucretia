#include "lucretiaManager.hh"
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <string.h>

using namespace std;

lucretiaManager::lucretiaManager(int* blele, int* bunchno, struct Bunch* ThisBunch, double L)
{
  fPrimIndex=NULL;
  fSecondaryPrimaryID = NULL ;
  fSecondariesPerThisPrimary = NULL ;
  Material = NULL ;
  GeomType = NULL ;
  fSecondaryBunch_x = NULL ;
  fPrimaryRegenID = NULL ;
  printf("Do Initialize...\n");
  Initialize(blele, bunchno, ThisBunch, L) ;
  printf("Done.\n");
}

lucretiaManager::~lucretiaManager()
{
}

void lucretiaManager::freeMem()
{
}

void lucretiaManager::Initialize(int* blele, int* bunchno, struct Bunch* ThisBunch, double L)
{
  uint32_T i ;
  fBunch=ThisBunch;
  fBunchNo=bunchno;
  Status=0;
  Lcut=L;
  fEle=blele;
  fRayCount=0;
  fRayGetPtr=0;
  fNumRaysResumed=0;
  fBunchNo=bunchno;
  // Clear any previously assigned memory that isn't required by Matlab session
  if ( fSecondariesPerThisPrimary != NULL )
    free(fSecondariesPerThisPrimary) ;
  if (fPrimIndex!=NULL)
    free(fPrimIndex) ;
  if (Material!=NULL) {
    free(Material) ;
    free(EMInterpMethod) ;
  }
  if (GeomType!=NULL)
    free(GeomType) ;
  // Get required extProcess properties
  /* ===== Get EM field data ===== */
  // B/E field values (scalar or 3D vectors)
  pBx = GetExtProcessData(blele,"Bx") ;
  pBy = GetExtProcessData(blele,"By") ;
  pBz = GetExtProcessData(blele,"Bz") ;
  pEx = GetExtProcessData(blele,"Ex") ;
  pEy = GetExtProcessData(blele,"Ey") ;
  pEz = GetExtProcessData(blele,"Ez") ;
  // check number of dimensions either 1 or 3
  if (mxGetNumberOfElements(pBx)>1 && mxGetNumberOfDimensions(pBx)!=3) {
    Status=1;
    return;
  }
  if (mxGetNumberOfElements(pBy)>1 && mxGetNumberOfDimensions(pBy)!=3) {
    Status=1;
    return;
  }
  if (mxGetNumberOfElements(pBz)>1 && mxGetNumberOfDimensions(pBz)!=3) {
    Status=1;
    return;
  }
  if (mxGetNumberOfElements(pEx)>1 && mxGetNumberOfDimensions(pEx)!=3) {
    Status=1;
    return;
  }
  if (mxGetNumberOfElements(pEy)>1 && mxGetNumberOfDimensions(pEy)!=3) {
    Status=1;
    return;
  }
  if (mxGetNumberOfElements(pEz)>1 && mxGetNumberOfDimensions(pEz)!=3) {
    Status=1;
    return;
  }
  EnableEM = *(char*) mxGetData(GetExtProcessData(blele,"EnableEM")) ;
  EMisUniform = 0;
  if (mxGetNumberOfElements(pBx)==1 && mxGetNumberOfElements(pBy)==1 && mxGetNumberOfElements(pBz)==1 &&
          mxGetNumberOfElements(pEx)==1 && mxGetNumberOfElements(pEy)==1 && mxGetNumberOfElements(pEz)==1) {
    if ( (*mxGetPr(pBx)+*mxGetPr(pBy)+*mxGetPr(pBz))!=0  )
      EMisUniform=1;
    if ( (*mxGetPr(pEx)+*mxGetPr(pEy)+*mxGetPr(pEz))!=0  && EMisUniform==0)
      EMisUniform=2;
    else if ( (*mxGetPr(pEx)+*mxGetPr(pEy)+*mxGetPr(pEz))!=0  && EMisUniform==1)
      EMisUniform = 0;
  }
  // Force GlobalField map calculation, uniform field implementation causes calculation to hang for unknown reason...
  EMisUniform = 0;
  EMStepSize = *mxGetPr(GetExtProcessData(blele,"StepPrec")) ;
  EMDeltaOneStep = *mxGetPr(GetExtProcessData(blele,"DeltaOneStep")) ;
  EMDeltaIntersection = *mxGetPr(GetExtProcessData(blele,"DeltaIntersection")) ;
  EMDeltaChord = *mxGetPr(GetExtProcessData(blele,"DeltaChord")) ;
  EMEpsMin = *mxGetPr(GetExtProcessData(blele,"EpsMin")) ;
  EMEpsMax = *mxGetPr(GetExtProcessData(blele,"EpsMax")) ;
  mxArray* pIM = GetExtProcessData(blele,"Interpolator") ;
  if (pIM == NULL) {
    Status = 1;
    return;
  }
  int buflen = mxGetN(pIM)*sizeof(mxChar)+1;
  EMInterpMethod = (char*) malloc(buflen);
  mxGetString(pIM, EMInterpMethod, buflen) ;
  mxArray* pStep = GetExtProcessData(blele,"StepMethod") ;
  if (pStep == NULL) {
    Status = 1;
    return;
  }
  buflen = mxGetN(pStep)*sizeof(mxChar)+1;
  EMStepperMethod = (char*) malloc(buflen);
  mxGetString(pStep, EMStepperMethod, buflen) ;
  // ====================================
  // - Geometry type (Ellipse or Rectangle)
  mxArray* pGeomType = GetExtProcessData(blele,"GeometryType") ;
  if (pGeomType == NULL) {
    Status = 1;
    return;
  }
  buflen = mxGetN(pGeomType)*sizeof(mxChar)+1;
  GeomType = (char*) malloc(buflen);
  mxGetString(pGeomType, GeomType, buflen) ;
  // - Material
  mxArray* pMaterial = GetExtProcessData(blele,"Material") ;
  if (pMaterial == NULL) {
    Status = 1;
    return;
  }
  buflen = mxGetN(pMaterial)*sizeof(mxChar)+1;
  Material = (char*) malloc(buflen);
  mxGetString(pMaterial, Material, buflen) ;
  // - X and Y Apertures
  mxArray* pAperX = GetExtProcessData(blele,"AperX") ;
  if (pAperX == NULL) {
    Status = 1;
    return;
  }
  AperX = *mxGetPr( pAperX ) ;
  mxArray* pAperY = GetExtProcessData(blele,"AperY") ;
  if (pAperY == NULL) {
    Status = 1;
    return;
  }
  AperY = *mxGetPr( pAperY ) ;
  // - Energy cut for returning tracks
  mxArray* pEcut = GetExtProcessData(blele,"Ecut") ;
  if (pEcut == NULL) {
    Status = 1;
    return;
  }
  Ecut = *mxGetPr( pEcut ) ;
  // - Max Number of primary particles to track
  mxArray* pMax = GetExtProcessData(blele,"MaxPrimaryParticles") ;
  if (pMax == NULL) {
    Status = 1;
    return;
  }
  fMaxPrimaryParticles = *(uint32_T*)mxGetData( pMax ) ;
  fPrimIndex = (int*) malloc(fMaxPrimaryParticles*sizeof(int)) ;
  // - Material thickness
  mxArray* pThickness = GetExtProcessData(blele,"Thickness") ;
  if (pThickness == NULL) {
    Status = 1;
    return;
  }
  Thickness = *mxGetPr( pThickness ) ;
  // - Verbosity control
  mxArray* pVerbose = GetExtProcessData(blele,"Verbose") ;
  Verbose = 0 ;
  if (pVerbose != NULL) {
    Verbose = (int) *mxGetPr( pVerbose ) ;
  }
  // - Primaries ordering
  mxArray* pOrder = mxGetCell(GetExtProcessData(blele,"PrimarySampleOrder"), *bunchno) ;
  if (pOrder == NULL)
    fPrimOrder = NULL ;
  else if ( mxGetM(pOrder)*mxGetN(pOrder) < fMaxPrimaryParticles )
    fPrimOrder = NULL ;
  else {
    fPrimOrder = (uint32_T*) mxGetData(pOrder) ;
  }
  // - vector of re-generated (unstopped) primary particles
  if (fPrimaryRegenID!=NULL)
    free(fPrimaryRegenID) ;
  fPrimaryRegenID = (uint32_T*) malloc(fMaxPrimaryParticles*sizeof(uint32_T)) ;
  for (i=0; i< fMaxPrimaryParticles; i++)
    fPrimaryRegenID[i]=0 ;
  // - Secondaries info
  fSecondariesCounter=0 ;
  //
  mxArray* pMaxSecondaries = GetExtProcessData(blele,"MaxSecondaryParticles") ;
  if ( pMaxSecondaries == NULL)
    fMaxSecondaryParticles = 0 ;
  else
    fMaxSecondaryParticles = *(uint32_T*)mxGetData( pMaxSecondaries ) ;
  //
  mxArray* pMaxSecondariesPerPrimary = GetExtProcessData(blele,"MaxSecondaryParticlesPerPrimary") ;
  if ( pMaxSecondariesPerPrimary == NULL)
    fMaxSecondaryParticlesPerPrimary = 0 ;
  else
    fMaxSecondaryParticlesPerPrimary = *(uint32_T*)mxGetData( pMaxSecondariesPerPrimary ) ;
  //
  fSecondariesPerThisPrimary = (int*) malloc(fMaxPrimaryParticles*sizeof(int)) ;
  // - Create secondary return bunch structure
  if (fMaxSecondaryParticles>0) {
    if (fSecondaryBunch_x!=NULL) {
      free(fSecondaryPrimaryID);
      free(fSecondaryBunch_x);
    }
    fSecondaryBunch_x = (double*) malloc(sizeof(double)*fMaxSecondaryParticles*6) ;
    fSecondaryPrimaryID = (uint32_T*) malloc(sizeof(uint32_T)*fMaxSecondaryParticles) ;
  }
  fTypeCellPtr = mxCreateCellMatrix(1,fMaxSecondaryParticles) ;
}

void lucretiaManager::GetUniformField(double uField[3])
{
  if (EMisUniform==1) {
    uField[0]=*mxGetPr(pBx);
    uField[1]=*mxGetPr(pBy);
    uField[2]=*mxGetPr(pBz);
  }
  else if (EMisUniform==2) {
    uField[0]=*mxGetPr(pEx);
    uField[1]=*mxGetPr(pEy);
    uField[2]=*mxGetPr(pEz);
  }
}

double lucretiaManager::interpField(const int fieldno, const double* point)
{
  const mxArray *F = NULL;
  char method = 0;
  
  // Which field are we looking up?
  switch(fieldno) {
    case 0:
      F = pBx;
      break;
    case 1:
      F = pBy;
      break;
    case 2:
      F = pBz;
      break;
    case 3:
      F = pEx;
      break;
    case 4:
      F = pEy;
      break;
    case 5:
      F = pEz;
      break;
  }
  // If just single field value then return that
  if (mxGetNumberOfElements(F)<2)
    return *mxGetPr(F);
  
  // Interpolation method (1,2,3 == nearest, linear, cubic)
  if (!strcmp(EMInterpMethod,"nearest"))
    method=0;
  else if (!strcmp(EMInterpMethod,"linear"))
    method=1;
  else if (!strcmp(EMInterpMethod,"cubic"))
    method=2;
  
  const mwSize *F_dims = mxGetDimensions(F);
  
  const mwSize M=F_dims[0];
  const mwSize N=F_dims[1];
  const mwSize O=F_dims[2];
  
  const double *pF = mxGetPr(F);
  const double *pX = &point[0];
  const double *pY = &point[1];
  const double *pZ = &point[2];
  double pO ;
  
  // Interpolation is based on x/y/z co-ordinates internal to G4 volume (X/Y/Z given in mm)
  const double x_low = -(AperX+Thickness)*1e3; const double x_high = (AperX+Thickness)*1e3;
  const double y_low = -(AperY+Thickness)*1e3; const double y_high = (AperY+Thickness)*1e3;
  const double z_low = -Lcut*1e3; const double z_high = Lcut*1e3;
  //if (fieldno==0)
    //printf("X: %g Y: %g Z: %g\n",*pX,*pY,*pZ);
  const double s_x = (double(1)-double(N))/(x_low - x_high);
  const double s_y = (double(1)-double(M))/(y_low - y_high);
  const double s_z = (double(1)-double(O))/(z_low - z_high);
  
  const double o_x = double(1)-x_low*s_x;
  const double o_y = double(1)-y_low*s_y;
  const double o_z = double(1)-z_low*s_z;
  
  
  // Do the interpolation
  switch(method) {
    case 0:
      interpolate_nearest(&pO, pF, pX, pY, pZ, 1, M, N, O, 1, s_x, o_x, s_y, o_y, s_z, o_z);
      break;
    case 1:
      interpolate_linear(&pO, pF, pX, pY, pZ, 1, M, N, O, 1, s_x, o_x, s_y, o_y, s_z, o_z);
      break;
    case 2:
      interpolate_bicubic(&pO, pF, pX, pY, pZ, 1, M, N, O, 1, s_x, o_x, s_y, o_y, s_z, o_z);
      break;
    default:
      mexErrMsgTxt("Unimplemented interpolation method.");
  }
  if (fieldno==0)
    //printf("B: %g (BMAX: %g)\n",pO,pF[0]);
  // Return the interpolated field value
  return pO ;
}

/* ======== Interpolation Routines ========= */

int lucretiaManager::access(int M, int N, int O, int x, int y, int z) {
  if (x<0) x=0; else if (x>=N) x=N-1;
  if (y<0) y=0; else if (y>=M) y=M-1;
  if (z<0) z=0; else if (z>=O) z=O-1;
  return y + M*(x + N*z);
}

int lucretiaManager::access_unchecked(int M, int N, int O, int x, int y, int z) {
  return y + M*(x + N*z);
}

void lucretiaManager::indices_linear(
        int &f000_i,
        int &f100_i,
        int &f010_i,
        int &f110_i,
        int &f001_i,
        int &f101_i,
        int &f011_i,
        int &f111_i,
        const int x, const int y, const int z,
        const mwSize &M, const mwSize &N, const mwSize &O) {
  if (x<=1 || y<=1 || z<=1 || x>=N-2 || y>=M-2 || z>=O-2) {
    f000_i = access(M,N,O, x,   y  , z);
    f100_i = access(M,N,O, x+1, y  , z);
    
    f010_i = access(M,N,O, x,   y+1, z);
    f110_i = access(M,N,O, x+1, y+1, z);
    
    f001_i = access(M,N,O, x,   y  , z+1);
    f101_i = access(M,N,O, x+1, y  , z+1);
    
    f011_i = access(M,N,O, x,   y+1, z+1);
    f111_i = access(M,N,O, x+1, y+1, z+1);
  } else {
    f000_i = access_unchecked(M,N,O, x,   y  , z);
    f100_i = access_unchecked(M,N,O, x+1, y  , z);
    
    f010_i = access_unchecked(M,N,O, x,   y+1, z);
    f110_i = access_unchecked(M,N,O, x+1, y+1, z);
    
    f001_i = access_unchecked(M,N,O, x,   y  , z+1);
    f101_i = access_unchecked(M,N,O, x+1, y  , z+1);
    
    f011_i = access_unchecked(M,N,O, x,   y+1, z+1);
    f111_i = access_unchecked(M,N,O, x+1, y+1, z+1);
  }
}

void lucretiaManager::indices_cubic(
        int f_i[64],
        const int x, const int y, const int z,
        const mwSize &M, const mwSize &N, const mwSize &O) {
  if (x<=2 || y<=2 || z<=2 || x>=N-3 || y>=M-3 || z>=O-3) {
    for (int i=0; i<4; ++i)
      for (int j=0; j<4; ++j)
        for (int k=0; k<4; ++k)
          f_i[i+4*(j+4*k)] = access(M,N,O, x+i-1, y+j-1, z+k-1);
  } else {
    for (int i=0; i<4; ++i)
      for (int j=0; j<4; ++j)
        for (int k=0; k<4; ++k)
          f_i[i+4*(j+4*k)] = access_unchecked(M,N,O, x+i-1, y+j-1, z+k-1);
  }
}

void lucretiaManager::interpolate_nearest(double *pO, const double *pF,
        const double *pX, const double *pY, const double *pZ,
        const mwSize ND, const mwSize M, const mwSize N, const mwSize O, const mwSize P,
        const double s_x, const double o_x,
        const double s_y, const double o_y,
        const double s_z, const double o_z) {
  const mwSize LO = M*N*O;
  for (mwSize i=0; i<ND; ++i) {
    const double &x = pX[i];
    const double &y = pY[i];
    const double &z = pZ[i];
    
    const int x_round = int(round(s_x*x+o_x))-1;
    const int y_round = int(round(s_y*y+o_y))-1;
    const int z_round = int(round(s_z*z+o_z))-1;
    
    const int f00_i = access(M,N,O, x_round,y_round,z_round);
    for (mwSize j=0; j<P; ++j) {
      pO[i + j*ND] = pF[f00_i + j*LO];
    }
  }
}

void lucretiaManager::interpolate_linear(double *pO, const double *pF,
        const double *pX, const double *pY, const double *pZ,
        const mwSize ND, const mwSize M, const mwSize N, const mwSize O, const mwSize P,
        const double s_x, const double o_x,
        const double s_y, const double o_y,
        const double s_z, const double o_z) {
  const mwSize LO = M*N*O;
  for (mwSize i=0; i<ND; ++i) {
    const double &x_ = pX[i];
    const double &y_ = pY[i];
    const double &z_ = pZ[i];
    
    const double x = s_x*x_+o_x;
    const double y = s_y*y_+o_y;
    const double z = s_z*z_+o_z;
    
    const double x_floor = floor(x);
    const double y_floor = floor(y);
    const double z_floor = floor(z);
    
    const double dx = x-x_floor;
    const double dy = y-y_floor;
    const double dz = z-z_floor;
    
    const double wx0 = 1.0-dx;
    const double wx1 = dx;
    
    const double wy0 = 1.0-dy;
    const double wy1 = dy;
    
    const double wz0 = 1.0-dz;
    const double wz1 = dz;
    
    int f000_i, f100_i, f010_i, f110_i;
    int f001_i, f101_i, f011_i, f111_i;
    
    // TODO: Use openmp
    indices_linear(
            f000_i, f100_i, f010_i, f110_i,
            f001_i, f101_i, f011_i, f111_i,
            int(x_floor-1), int(y_floor-1), int(z_floor-1), M, N, O);
    
    for (mwSize j=0; j<P; ++j) {
      
      pO[i + j*ND] =
              wz0*(
              wy0*(wx0 * pF[f000_i + j*LO] + wx1 * pF[f100_i + j*LO]) +
              wy1*(wx0 * pF[f010_i + j*LO] + wx1 * pF[f110_i + j*LO])
              )+
              wz1*(
              wy0*(wx0 * pF[f001_i + j*LO] + wx1 * pF[f101_i + j*LO]) +
              wy1*(wx0 * pF[f011_i + j*LO] + wx1 * pF[f111_i + j*LO])
              );
    }
    
  }
}

void lucretiaManager::interpolate_bicubic(double *pO, const double *pF,
        const double *pX, const double *pY, const double *pZ,
        const mwSize ND, const mwSize M, const mwSize N, const mwSize O, const mwSize P,
        const double s_x, const double o_x,
        const double s_y, const double o_y,
        const double s_z, const double o_z)
{
  const mwSize LO = M*N*O;
  for (mwSize i=0; i<ND; ++i) {
    const double &x_ = pX[i];
    const double &y_ = pY[i];
    const double &z_ = pZ[i];
    
    const double x = s_x*x_+o_x;
    const double y = s_y*y_+o_y;
    const double z = s_z*z_+o_z;
    
    const double x_floor = floor(x);
    const double y_floor = floor(y);
    const double z_floor = floor(z);
    
    const double dx = x-x_floor;
    const double dy = y-y_floor;
    const double dz = z-z_floor;
    
    const double dxx = dx*dx;
    const double dxxx = dxx*dx;
    
    const double dyy = dy*dy;
    const double dyyy = dyy*dy;
    
    const double dzz = dz*dz;
    const double dzzz = dzz*dz;
    
    const double wx0 = 0.5 * (    - dx + 2.0*dxx -       dxxx);
    const double wx1 = 0.5 * (2.0      - 5.0*dxx + 3.0 * dxxx);
    const double wx2 = 0.5 * (      dx + 4.0*dxx - 3.0 * dxxx);
    const double wx3 = 0.5 * (         -     dxx +       dxxx);
    
    const double wy0 = 0.5 * (    - dy + 2.0*dyy -       dyyy);
    const double wy1 = 0.5 * (2.0      - 5.0*dyy + 3.0 * dyyy);
    const double wy2 = 0.5 * (      dy + 4.0*dyy - 3.0 * dyyy);
    const double wy3 = 0.5 * (         -     dyy +       dyyy);
    
    const double wz0 = 0.5 * (    - dz + 2.0*dzz -       dzzz);
    const double wz1 = 0.5 * (2.0      - 5.0*dzz + 3.0 * dzzz);
    const double wz2 = 0.5 * (      dz + 4.0*dzz - 3.0 * dzzz);
    const double wz3 = 0.5 * (         -     dzz +       dzzz);
    
    int f_i[64];
    
    indices_cubic(
            f_i,
            int(x_floor-1), int(y_floor-1), int(z_floor-1), M, N, O);
    
    for (mwSize j=0; j<P; ++j) {
      
      pO[i + j*ND] =
              wz0*(
              wy0*(wx0 * pF[f_i[0+4*(0+4*0)] + j*LO] + wx1 * pF[f_i[1+4*(0+4*0)] + j*LO] +  wx2 * pF[f_i[2+4*(0+4*0)] + j*LO] + wx3 * pF[f_i[3+4*(0+4*0)] + j*LO]) +
              wy1*(wx0 * pF[f_i[0+4*(1+4*0)] + j*LO] + wx1 * pF[f_i[1+4*(1+4*0)] + j*LO] +  wx2 * pF[f_i[2+4*(1+4*0)] + j*LO] + wx3 * pF[f_i[3+4*(1+4*0)] + j*LO]) +
              wy2*(wx0 * pF[f_i[0+4*(2+4*0)] + j*LO] + wx1 * pF[f_i[1+4*(2+4*0)] + j*LO] +  wx2 * pF[f_i[2+4*(2+4*0)] + j*LO] + wx3 * pF[f_i[3+4*(2+4*0)] + j*LO]) +
              wy3*(wx0 * pF[f_i[0+4*(3+4*0)] + j*LO] + wx1 * pF[f_i[1+4*(3+4*0)] + j*LO] +  wx2 * pF[f_i[2+4*(3+4*0)] + j*LO] + wx3 * pF[f_i[3+4*(3+4*0)] + j*LO])
              ) +
              wz1*(
              wy0*(wx0 * pF[f_i[0+4*(0+4*1)] + j*LO] + wx1 * pF[f_i[1+4*(0+4*1)] + j*LO] +  wx2 * pF[f_i[2+4*(0+4*1)] + j*LO] + wx3 * pF[f_i[3+4*(0+4*1)] + j*LO]) +
              wy1*(wx0 * pF[f_i[0+4*(1+4*1)] + j*LO] + wx1 * pF[f_i[1+4*(1+4*1)] + j*LO] +  wx2 * pF[f_i[2+4*(1+4*1)] + j*LO] + wx3 * pF[f_i[3+4*(1+4*1)] + j*LO]) +
              wy2*(wx0 * pF[f_i[0+4*(2+4*1)] + j*LO] + wx1 * pF[f_i[1+4*(2+4*1)] + j*LO] +  wx2 * pF[f_i[2+4*(2+4*1)] + j*LO] + wx3 * pF[f_i[3+4*(2+4*1)] + j*LO]) +
              wy3*(wx0 * pF[f_i[0+4*(3+4*1)] + j*LO] + wx1 * pF[f_i[1+4*(3+4*1)] + j*LO] +  wx2 * pF[f_i[2+4*(3+4*1)] + j*LO] + wx3 * pF[f_i[3+4*(3+4*1)] + j*LO])
              ) +
              wz2*(
              wy0*(wx0 * pF[f_i[0+4*(0+4*2)] + j*LO] + wx1 * pF[f_i[1+4*(0+4*2)] + j*LO] +  wx2 * pF[f_i[2+4*(0+4*2)] + j*LO] + wx3 * pF[f_i[3+4*(0+4*2)] + j*LO]) +
              wy1*(wx0 * pF[f_i[0+4*(1+4*2)] + j*LO] + wx1 * pF[f_i[1+4*(1+4*2)] + j*LO] +  wx2 * pF[f_i[2+4*(1+4*2)] + j*LO] + wx3 * pF[f_i[3+4*(1+4*2)] + j*LO]) +
              wy2*(wx0 * pF[f_i[0+4*(2+4*2)] + j*LO] + wx1 * pF[f_i[1+4*(2+4*2)] + j*LO] +  wx2 * pF[f_i[2+4*(2+4*2)] + j*LO] + wx3 * pF[f_i[3+4*(2+4*2)] + j*LO]) +
              wy3*(wx0 * pF[f_i[0+4*(3+4*2)] + j*LO] + wx1 * pF[f_i[1+4*(3+4*2)] + j*LO] +  wx2 * pF[f_i[2+4*(3+4*2)] + j*LO] + wx3 * pF[f_i[3+4*(3+4*2)] + j*LO])
              ) +
              wz3*(
              wy0*(wx0 * pF[f_i[0+4*(0+4*3)] + j*LO] + wx1 * pF[f_i[1+4*(0+4*3)] + j*LO] +  wx2 * pF[f_i[2+4*(0+4*3)] + j*LO] + wx3 * pF[f_i[3+4*(0+4*3)] + j*LO]) +
              wy1*(wx0 * pF[f_i[0+4*(1+4*3)] + j*LO] + wx1 * pF[f_i[1+4*(1+4*3)] + j*LO] +  wx2 * pF[f_i[2+4*(1+4*3)] + j*LO] + wx3 * pF[f_i[3+4*(1+4*3)] + j*LO]) +
              wy2*(wx0 * pF[f_i[0+4*(2+4*3)] + j*LO] + wx1 * pF[f_i[1+4*(2+4*3)] + j*LO] +  wx2 * pF[f_i[2+4*(2+4*3)] + j*LO] + wx3 * pF[f_i[3+4*(2+4*3)] + j*LO]) +
              wy3*(wx0 * pF[f_i[0+4*(3+4*3)] + j*LO] + wx1 * pF[f_i[1+4*(3+4*3)] + j*LO] +  wx2 * pF[f_i[2+4*(3+4*3)] + j*LO] + wx3 * pF[f_i[3+4*(3+4*3)] + j*LO])
              );
    }
    
  }
}

int lucretiaManager::GetNextX()
{
  // Return pointer to next 6D ray co-ordinates that has a stopped flag at this element #
  uint32_T iray, useray, i ;
  //cout << "fRayCount: " << fRayCount << "\n" ;
  if (fRayGetPtr >= (fBunch->nray-1) || fRayCount >= fMaxPrimaryParticles )
    return -1 ;
  for (i=fRayGetPtr; i<fBunch->nray; i++) {
    if (fPrimOrder!=NULL)
      iray=fPrimOrder[i] ;
    else
      iray=i ;
    //cout << "iray: " << iray << "\n" ;
    if ( fBunch->stop[iray] == *fEle+1 ) {
      fRayGetPtr = iray+1 ;
      // return if ray inside geant tracking volume
      useray=1;
      if ( (fBunch->x[iray*6]>(AperX+Thickness)) || (fBunch->x[iray*6+2]>(AperY+Thickness)) ||
              (fBunch->x[iray*6+5]<Ecut) )
        useray=0;
      if (useray==1) {
        fPrimIndex[fRayCount]=iray;
        fSecondariesPerThisPrimary[fRayCount]=0; // reset secondaries count
        fRayCount++;
        /*cout << "GET_ID = " << fRayCount-1 << " iray: " << iray << " X/Y/Z: " <<
         * fBunch->x[6*iray] << " / " << fBunch->x[6*iray+2] << " / " << fBunch->x[6*iray+4] <<
         * " X'/Y' : " << fBunch->x[6*iray+1] << " / " << fBunch->x[6*iray+3] << " E: " <<
         * fBunch->x[6*iray+5] << "\n" ;*/
        return iray ;
      }
    }
  }
  return -1 ;
}

void lucretiaManager::SetLucretiaData()
{
  // Send object data back into Lucretia data structures
  if (fNumRaysResumed>0) {
    mxArray* pPrimaryRegenID = mxCreateNumericMatrix(1, fNumRaysResumed, mxUINT32_CLASS, mxREAL) ;
    uint32_T* dPrimaryRegenID = (uint32_T*) mxMalloc( sizeof(uint32_T)*fNumRaysResumed ) ;
    memcpy(dPrimaryRegenID, fPrimaryRegenID, sizeof(uint32_T)*fNumRaysResumed );
    mxSetData( pPrimaryRegenID, dPrimaryRegenID) ;
    mxSetProperty( GetExtProcessPrimariesData(fEle), *fBunchNo, "regeneratedID", pPrimaryRegenID ) ;
  }
  //
  if (fSecondariesCounter>0) {
    mxArray* pSecondaryBunch_x = mxCreateNumericMatrix(6, fSecondariesCounter, mxDOUBLE_CLASS, mxREAL) ;
    double* dSecondaryBunch_x = (double*) mxMalloc(fSecondariesCounter*sizeof(double)*6) ;
    memcpy( dSecondaryBunch_x, fSecondaryBunch_x, fSecondariesCounter*sizeof(double)*6) ;
    mxSetPr( pSecondaryBunch_x, dSecondaryBunch_x) ;
    mxSetProperty( GetExtProcessSecondariesData(fEle), *fBunchNo, "Pos", pSecondaryBunch_x) ;
    //
    mxArray* pSecondaryPrimaryID = mxCreateNumericMatrix(1, fSecondariesCounter, mxUINT32_CLASS, mxREAL) ;
    uint32_T* dSecondaryPrimaryID = (uint32_T*) mxMalloc(fSecondariesCounter*sizeof(uint32_T)) ;
    memcpy( dSecondaryPrimaryID, fSecondaryPrimaryID, fSecondariesCounter*sizeof(uint32_T)) ;
    mxSetData( pSecondaryPrimaryID, dSecondaryPrimaryID) ;
    mxSetProperty( GetExtProcessSecondariesData(fEle), *fBunchNo, "PrimaryID", pSecondaryPrimaryID) ;
    //
    mxSetProperty( GetExtProcessSecondariesData(fEle), *fBunchNo, "ParticleType", fTypeCellPtr) ;
    //
    mxArray* pNumSec = mxCreateNumericMatrix(1, 1, mxUINT32_CLASS, mxREAL) ;
    uint32_T* dNumSec = (uint32_T*) mxMalloc(sizeof(uint32_T)) ;
    *dNumSec = fSecondariesCounter ;
    mxSetData( pNumSec, dNumSec );
    mxSetProperty( GetExtProcessSecondariesData(fEle), *fBunchNo, "NumStored", pNumSec) ;
  }
}

void lucretiaManager::SetNextSecondary(double x[6], int id, const char* type)
{
  int icoord ;
  int iray=fPrimIndex[id] ;
  if (fMaxSecondaryParticles==0)
    return;
  if (fSecondaryBunch_x != NULL ) {
    for (icoord=0; icoord<6; icoord++)
      fSecondaryBunch_x[fSecondariesCounter*6+icoord] = x[icoord] ;
    mxSetCell(fTypeCellPtr,fSecondariesCounter,mxCreateString(type)) ;
    fSecondaryPrimaryID[fSecondariesCounter] = iray+1 ;
  }
  fSecondariesCounter++ ;
}

void lucretiaManager::SetNextX(double x[6], int id, int doresume)
{
  // Set next stopped particle processed here to GEANT tracked co-ordinates and reset stop (resume tracking) if requested
  int iray, icoord ;
  iray=fPrimIndex[id];
  for (icoord=0; icoord<6; icoord++)
    if (icoord != 4) // Leave z co-ordinate alone
      fBunch->x[6*iray+icoord] = x[icoord] ;
  if (doresume == 1) {
    fBunch->stop[iray] = 0 ;
    fPrimaryRegenID[fNumRaysResumed] = iray+1 ;
    //cout << "iray: " << iray+1 << " stop: " << fBunch->stop[iray] << "\n" ;
    fNumRaysResumed++;
  }
  else { // coordinate should be where particle stopped
    fBunch->x[6*iray+4] = x[4] ;
  }
  /*cout << "SET_ID = " << id << " iray: " << iray << " X/Y/Z: " <<
   * fBunch->x[6*iray] << " / " << fBunch->x[6*iray+2] << " / " << fBunch->x[6*iray+4] <<
   * " X'/Y' : " << fBunch->x[6*iray+1] << " / " << fBunch->x[6*iray+3] << " E: " <<
   * fBunch->x[6*iray+5] << " Resume: " << doresume << "\n" ;*/
  return ;
}
