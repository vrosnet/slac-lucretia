#include "lucretiaManager.hh"
#include <iostream>
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
  Initialize(blele, bunchno, ThisBunch, L) ;
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
  if (Material!=NULL)
    free(Material) ;
  if (GeomType!=NULL)
    free(GeomType) ;
  // Get required extProcess properties
  // - Geometry type (Ellipse or Rectangle)
  mxArray* pGeomType = GetExtProcessData(blele,"GeometryType") ;
  if (pGeomType == NULL) {
    Status = 1;
    return;
  }
  int buflen = mxGetN(pGeomType)*sizeof(mxChar)+1;
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

int lucretiaManager::GetNextX()
{
  // Return pointer to next 6D ray co-ordinates that has a stopped flag at this element #
  uint32_T iray, useray, i ;
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
         fBunch->x[6*iray] << " / " << fBunch->x[6*iray+2] << " / " << fBunch->x[6*iray+4] <<
         " X'/Y' : " << fBunch->x[6*iray+1] << " / " << fBunch->x[6*iray+3] << " E: " <<
         fBunch->x[6*iray+5] << "\n" ;*/
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
    fNumRaysResumed++;
  }
  /*cout << "SET_ID = " << id << " iray: " << iray << " X/Y/Z: " <<
   * fBunch->x[6*iray] << " / " << fBunch->x[6*iray+2] << " / " << fBunch->x[6*iray+4] <<
   * " X'/Y' : " << fBunch->x[6*iray+1] << " / " << fBunch->x[6*iray+3] << " E: " <<
   * fBunch->x[6*iray+5] << " Resume: " << doresume << "\n" ;*/
  return ;
}
