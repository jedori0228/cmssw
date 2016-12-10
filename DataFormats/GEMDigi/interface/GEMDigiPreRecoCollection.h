#ifndef GEMDigi_GEMDigiPreRecoCollection_h
#define GEMDigi_GEMDigiPreRecoCollection_h
/** \class GEMDigiPreRecoCollection
 *  
 *  \author Marcello Maggi
 */

#include <DataFormats/MuonDetId/interface/GEMDetId.h>
#include <DataFormats/GEMDigi/interface/GEMDigiPreReco.h>
#include <DataFormats/MuonData/interface/MuonDigiCollection.h>

typedef MuonDigiCollection<GEMDetId, GEMDigiPreReco> GEMDigiPreRecoCollection;

#endif
