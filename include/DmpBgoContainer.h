#ifndef DMPBGOCONTAINER_H
#define DMPBGOCONTAINER_H

#include <vector>
#include <memory>

#include "DmpEvtBgoHits.h"

class DmpBgoContainer
{
public:
    DmpBgoContainer(int m_size) : layerBarIndex(m_size, std::vector<short>()),
                                  layerBarNumber(m_size, std::vector<short>()),
                                  idxBarMaxLayer(m_size, -1),
                                  iMaxLayer(m_size, -1),
                                  rmsLayer(m_size, 0),
                                  fracLayer(m_size, 0),
                                  eLayer(m_size, 0),
                                  eCoreLayer(m_size, 0),
                                  eCoreCoord(m_size, 0)
    {
    }
    ~DmpBgoContainer(){};
    void scanBGOHits(
        const std::shared_ptr<DmpEvtBgoHits> bgohits,
        const double bgoTotalE,
        const int DAMPE_bgo_nLayers);
    std::vector<short> GetSingleLayerBarNumber(int nLayer);
    std::vector<std::vector<short>> GetLayerBarNumber();
    std::vector<double> GetRmsLayer();
    std::vector<double> GetFracLayer();
    std::vector<int> GetIdxBarMaxLayer();
    const int GetIdxBarMaxSingleLayer(const int layedIdx);
    std::vector<double> GetELayer();
    const double GetESingleLayer(const int layedIdx);
    std::vector<int> GetiMaxLayer();
    const int GetiMaxSingleLayer(const int layedIdx);
    const double GetSumRMS();

private:
    std::vector<std::vector<short>> layerBarIndex;  // arrange BGO hits by layer
    std::vector<std::vector<short>> layerBarNumber; // arrange BGO bars by layer
    std::vector<int> idxBarMaxLayer;
    std::vector<double> rmsLayer;
    std::vector<double> fracLayer;
    std::vector<double> eLayer;
    std::vector<double> eCoreLayer;
    std::vector<double> eCoreCoord;
    std::vector<int> iMaxLayer;
    double sumRms = 0;
    double Xtr = 0;
};

#endif