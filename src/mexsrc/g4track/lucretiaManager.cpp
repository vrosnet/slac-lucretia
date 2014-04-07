#include "mex.h"
#include "lucretiaManager.hh"
#include <iostream>

using namespace std;

lucretiaManager::lucretiaManager(int* blele, int* bunchno, struct Bunch* ThisBunch, double L)
{
  fPrimIndex=NULL;
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
  fBunch=ThisBunch;
  fBunchNo=bunchno;
  Status=0;
  Lcut=L;
  fEle=blele;
  fRayCount=0;
  fRayGetPtr=0;
  fNumRaysResumed=0;
  fBunchNo=bunchno;
  // Get required extProcess properties
  // - Geometry type (Ellipse or Rectangle)
  mxArray* pGeomType = GetExtProcessData(blele,"GeometryType") ;
  if (pGeomType == NULL) {
    Status = 1;
    return;
  }
  int buflen = mxGetN(pGeomType)*sizeof(mxChar)+1;
  if (GeomType!=NULL)
    free(GeomType) ;
  GeomType = (char*) malloc(buflen);
  mxGetString(pGeomType, GeomType, buflen) ;
  // - Material
  mxArray* pMaterial = GetExtProcessData(blele,"Material") ;
  if (pMaterial == NULL) {
    Status = 1;
    return;
  }
  buflen = mxGetN(pMaterial)*sizeof(mxChar)+1;
  if (Material!=NULL)
    free(Material) ;
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
  fMaxPrimaryParticles = *mxGetPr( pMax ) ;
  if (fPrimIndex!=NULL)
    free(fPrimIndex) ;
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
  mxArray* pOrder = GetExtProcessData(blele,"PrimaryOrder") ;
  if (pOrder == NULL)
    fPrimOrder = NULL ;
  else if ( mxGetM(pOrder)*mxGetN(pOrder) < fMaxPrimaryParticles )
    fPrimOrder = NULL ;
  else
    fPrimOrder = mxGetPr(pOrder) ;
  // - Secondaries info
  fSecondariesCounter=0 ;
  mxArray* pMaxSecondaries = GetExtProcessData(blele,"MaxSecondaryParticles") ;
  if ( pMaxSecondaries == NULL)
    fMaxSecondaryParticles = 0 ;
  else
    fMaxSecondaryParticles = *mxGetPr( pMaxSecondaries ) ;
  mxArray* pMaxSecondariesPerPrimary = GetExtProcessData(blele,"MaxSecondaryParticlesPerPrimary") ;
  if ( pMaxSecondariesPerPrimary == NULL)
    fMaxSecondaryParticlesPerPrimary = 0 ;
  else
    fMaxSecondaryParticlesPerPrimary = *mxGetPr( pMaxSecondariesPerPrimary ) ;
  fSecondariesPerThisPrimary = (int*) malloc(fMaxPrimaryParticles*sizeof(int)) ;
  mxArray* pSecondaryBeam = GetExtProcessData(blele,"SecondaryBeam") ;
  mxArray* pSecondaryBunch ;
  if ( pSecondaryBeam !=NULL ) {
    pSecondaryBunch = mxGetField( pSecondaryBeam, *bunchno, "Bunch" ) ;
  }
  else
    pSecondaryBunch = NULL ;
  if ( pSecondaryBunch !=NULL ) {
    fSecondaryBunch_x = mxGetPr(mxGetField( pSecondaryBunch, *bunchno, "x" )) ;
    fSecondaryBunch_q = mxGetPr(mxGetField( pSecondaryBunch, *bunchno, "Q" )) ;
    fSecondaryBunch_type = mxGetField( pSecondaryBunch, *bunchno, "type" ) ;
  }
  else {
    fSecondaryBunch_x = NULL ;
    fSecondaryBunch_q = NULL ;
    fSecondaryBunch_type = NULL ;
  }
  fNumSecondariesStored = mxGetPr( GetExtProcessData(blele,"NumSecondariesStored") );
}

int lucretiaManager::GetNextX()
{
  // Return pointer to next 6D ray co-ordinates that has a stopped flag at this element #
  int iray, useray, i ;
  if (fRayGetPtr >= (fBunch->nray-1) || fRayCount >= fMaxPrimaryParticles )
    return -1 ;
  for (i=fRayGetPtr; i<fBunch->nray; i++) {
    if (fPrimOrder!=NULL)
      iray=(int) fPrimOrder[i] ;
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

void lucretiaManager::SetNextSecondary(double x[6], const char* type)
{
  int icoord ;
  double tval ;
  if (fSecondaryBunch_x != NULL ) {
    for (icoord=0; icoord<6; icoord++)
      fSecondaryBunch_x[6*fSecondariesCounter+icoord] = x[icoord] ;\
    if ( fSecondaryBunch_q != NULL )
      fSecondaryBunch_q[fSecondariesCounter] = fBunch->Q[fSecondariesCounter] ;
    if ( fSecondaryBunch_type !=NULL )
      mxSetCell(fSecondaryBunch_type, fSecondariesCounter, mxCreateString(type)) ;
  }
  fSecondariesCounter++ ;
  *fNumSecondariesStored = (double) fSecondariesCounter ;
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
    fNumRaysResumed++;
  }
  /*cout << "SET_ID = " << id << " iray: " << iray << " X/Y/Z: " <<
   * fBunch->x[6*iray] << " / " << fBunch->x[6*iray+2] << " / " << fBunch->x[6*iray+4] <<
   * " X'/Y' : " << fBunch->x[6*iray+1] << " / " << fBunch->x[6*iray+3] << " E: " <<
   * fBunch->x[6*iray+5] << " Resume: " << doresume << "\n" ;*/
  return ;
}
