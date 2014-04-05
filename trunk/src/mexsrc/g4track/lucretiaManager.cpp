#include "mex.h"
#include "lucretiaManager.hh"
#include <iostream>

lucretiaManager::lucretiaManager(int* blele, struct Bunch* ThisBunch)
: fBunch(ThisBunch),
        Status(0),
        fEle(blele),
        fRayCount(0),
        fRayGetPtr(0)
{
  // Get required extProcess properties
  // - Geometry type (Ellipse or Rectangle)
  mxArray* pGeomType = GetExtProcessData(blele,"GeometryType") ;
  if (pGeomType == NULL) {
    Status = 1;
    return;
  }
  int buflen = mxGetN(pGeomType)*sizeof(mxChar)+1;
  GeomType = (char*) mxMalloc(buflen);
  mxGetString(pGeomType, GeomType, buflen) ;
  // - Material
  mxArray* pMaterial = GetExtProcessData(blele,"Material") ;
  if (pMaterial == NULL) {
    Status = 1;
    return;
  }
  buflen = mxGetN(pMaterial)*sizeof(mxChar)+1;
  Material = (char*) mxMalloc(buflen);
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
  fMaxPrimaryParticles = *mxGetPr( pMax ) ;
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
}


lucretiaManager::~lucretiaManager()
{
  free(fPrimIndex) ;
  mxFree(Material) ;
  mxFree(GeomType) ;
}

int lucretiaManager::GetNextX()
{
  // Return pointer to next 6D ray co-ordinates that has a stopped flag at this element #
  int iray, useray ;
  if (fRayGetPtr >= (fBunch->nray-1) || fRayCount >= fMaxPrimaryParticles )
    return -1 ;
  for (iray=fRayGetPtr; iray<fBunch->nray; iray++) {
    if ( fBunch->stop[iray] == *fEle+1 ) {
      fRayGetPtr = iray+1 ;
      // return if ray inside geant tracking volume
      useray=1;
      if ( (fBunch->x[iray*6]>(AperX+Thickness)) || (fBunch->x[iray*6+2]>(AperY+Thickness)) ||
           (fBunch->x[iray*6+5]<Ecut) )
          useray=0;
      if (useray==1) {
        fPrimIndex[fRayCount]=iray;
        fRayCount++;
        /*cout << "GET_ID = " << fRayCount-1 << " iray: " << iray << " X/Y/Z: " <<
          fBunch->x[6*iray] << " / " << fBunch->x[6*iray+2] << " / " << fBunch->x[6*iray+4] <<
          " X'/Y' : " << fBunch->x[6*iray+1] << " / " << fBunch->x[6*iray+3] << " E: " <<
          fBunch->x[6*iray+5] << "\n" ;*/
        return 6*iray ;
      }
    }
  }
  return -1 ;
}

void lucretiaManager::SetNextX(double x[6], int id, int doresume)
{
  // Set next stopped particle processed here to GEANT tracked co-ordinates and reset stop (resume tracking) if requested
  int iray, icoord ;
  iray=fPrimIndex[id];
  for (icoord=0; icoord<6; icoord++)
    if (icoord != 4) // Leave z co-ordinate alone
      fBunch->x[6*iray+icoord] = x[icoord] ;
  if (doresume == 1)
    fBunch->stop[iray] = 0 ;
  /*cout << "SET_ID = " << id << " iray: " << iray << " X/Y/Z: " <<
          fBunch->x[6*iray] << " / " << fBunch->x[6*iray+2] << " / " << fBunch->x[6*iray+4] <<
          " X'/Y' : " << fBunch->x[6*iray+1] << " / " << fBunch->x[6*iray+3] << " E: " <<
          fBunch->x[6*iray+5] << " Resume: " << doresume << "\n" ;*/
  return ;
}
