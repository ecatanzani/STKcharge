#include "DmpBgoContainer.h"

void DmpBgoContainer::scanBGOHits(
    const std::shared_ptr<DmpEvtBgoHits> bgohits,
    const double bgoTotalE,
    const int DAMPE_bgo_nLayers)
{
    // Get the number of BGO hits
    int nBgoHits = bgohits->GetHittedBarNumber();

    // Scan BGO hits
    for (int ihit = 0; ihit < nBgoHits; ++ihit)
    {
        // Get layer ID
        auto layerID = bgohits->GetLayerID(ihit);
        // Get bar global ID
        auto iBar = ((bgohits->fGlobalBarID)[ihit] >> 6) & 0x1f;
        layerBarIndex[layerID].push_back(ihit);
        layerBarNumber[layerID].push_back(iBar);
    }

    for (int lay = 0; lay < DAMPE_bgo_nLayers; ++lay)
    {
        // Setting default value for maximum bar index and energy for each layer
        int imax = -1;
        double maxE = 0;

        // Find the maximum of the nergy release in a certain layer, together with the bar ID
        for (auto it = layerBarNumber[lay].begin(); it != layerBarNumber[lay].end(); ++it)
        {
            auto index = std::distance(layerBarNumber[lay].begin(), it);
            int ihit = layerBarIndex[lay][index];
            double hitE = (bgohits->fEnergy)[ihit];
            if (hitE > maxE)
            {
                maxE = hitE;
                imax = ihit;
            }
        }
        iMaxLayer[lay] = imax;
        rmsLayer[lay] = 0;
        fracLayer[lay] = 0;
        eLayer[lay] = 0;

        if (maxE > 0)
        {
            // Find the bar index regarding the maximum energy release in a certain layer
            auto iBarMax = ((bgohits->fGlobalBarID)[imax] >> 6) & 0x1f;
            idxBarMaxLayer[lay] = iBarMax;
            // Register the maximum energy release of a layer
            eCoreLayer[lay] = maxE;
            // Find the coordinate (weighted by the nergy release) of the bar with the biggest energy release in a certain layer
            double coordMax = lay % 2 ? bgohits->GetHitX(imax) : bgohits->GetHitY(imax);
            eCoreCoord[lay] = maxE * coordMax;
            // Consider the nearest bar respect to the max one in order to better interpolate the position
            if (iBarMax > 0 && iBarMax < 21)
            {
                for (auto it = layerBarNumber[lay].begin(); it != layerBarNumber[lay].end(); ++it)
                {
                    auto index = std::distance(layerBarNumber[lay].begin(), it);
                    int ihit = layerBarIndex[lay][index];
                    auto iBar = ((bgohits->fGlobalBarID)[ihit] >> 6) & 0x1f;
                    if (iBar - iBarMax == 1 || iBar - iBarMax == -1)
                    {
                        double hitE = (bgohits->fEnergy)[ihit];
                        double thisCoord = lay % 2 ? bgohits->GetHitX(ihit) : bgohits->GetHitY(ihit);
                        eCoreLayer[lay] += hitE;
                        eCoreCoord[lay] += hitE * thisCoord;
                    }
                }
            }
            // Get the CoG coordinate of the max energy bar cluster
            eCoreCoord[lay] /= eCoreLayer[lay];
            // Get the energy RMS of a layer
            for (auto it = layerBarNumber[lay].begin(); it != layerBarNumber[lay].end(); ++it)
            {
                auto index = std::distance(layerBarNumber[lay].begin(), it);
                int ihit = layerBarIndex[lay][index];
                auto hitE = (bgohits->fEnergy)[ihit];
                auto thisCoord = lay % 2 ? bgohits->GetHitX(ihit) : bgohits->GetHitY(ihit);
                eLayer[lay] += hitE;
                rmsLayer[lay] += hitE * (thisCoord - eCoreCoord[lay]) * (thisCoord - eCoreCoord[lay]);
            }
            rmsLayer[lay] = sqrt(rmsLayer[lay] / eLayer[lay]);
            fracLayer[lay] = eLayer[lay] / bgoTotalE;
            if (layerBarNumber[lay].size() <= 1)
                rmsLayer[lay] = 0;
        }
        sumRms += rmsLayer[lay];
    }

    // Build XTR
    Xtr = pow(sumRms, 4) * fracLayer[13] / 8000000.;
}

std::vector<short> DmpBgoContainer::GetSingleLayerBarNumber(int nLayer)
{
    return layerBarNumber[nLayer];
}

std::vector<std::vector<short>> DmpBgoContainer::GetLayerBarNumber()
{
    return layerBarNumber;
}

std::vector<double> DmpBgoContainer::GetRmsLayer()
{
    return rmsLayer;
}

const double DmpBgoContainer::GetSumRMS()
{
    return sumRms;
}

std::vector<double> DmpBgoContainer::GetFracLayer()
{
    return fracLayer;
}

std::vector<int> DmpBgoContainer::GetIdxBarMaxLayer()
{
    return idxBarMaxLayer;
}

const int DmpBgoContainer::GetIdxBarMaxSingleLayer(const int layedIdx)
{
    return idxBarMaxLayer[layedIdx];
}

std::vector<double> DmpBgoContainer::GetELayer()
{
    return eLayer;
}

const double DmpBgoContainer::GetESingleLayer(const int layedIdx)
{
    return eLayer[layedIdx];
}

std::vector<int> DmpBgoContainer::GetiMaxLayer()
{
    return iMaxLayer;
}

const int DmpBgoContainer::GetiMaxSingleLayer(const int layedIdx)
{
    return iMaxLayer[layedIdx];
}