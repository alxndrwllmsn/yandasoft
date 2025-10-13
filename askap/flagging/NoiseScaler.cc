/// @file NoiseScaler.cc
///
/// @copyright (c) 2025 CSIRO
/// Australia Telescope National Facility (ATNF)
/// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
/// PO Box 76, Epping NSW 1710, Australia
/// atnf-enquiries@csiro.au
///
/// This file is part of the ASKAP software distribution.
///
/// The ASKAP software distribution is free software: you can redistribute it
/// and/or modify it under the terms of the GNU General Public License as
/// published by the Free Software Foundation; either version 2 of the License,
/// or (at your option) any later version.
///
/// This program is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
/// GNU General Public License for more details.
///
/// You should have received a copy of the GNU General Public License
/// along with this program; if not, write to the Free Software
/// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
///

// Include own header file first
#include "askap/flagging/NoiseScaler.h"

// Include package level header file
#include "askap/askap_synthesis.h"

using namespace askap;
using namespace casacore;
using namespace askap::synthesis;

NoiseScaler::NoiseScaler(const LOFAR::ParameterSet& parset)
{
    // supply polynomial coefficients highest order first, for freq in GHz
    // default is the ASKAP sensitivity curve (0.8-1.8 GHz) fitted with 12th order polynomial
    itsParams = parset.getDoubleVector("noisePoly",std::vector<double>(
        { 1.45412234e+06, -2.15243951e+07,  1.44698939e+08, -5.84085389e+08,
          1.57652581e+09, -2.99733763e+09,  4.11572297e+09, -4.11249045e+09,
          2.96786192e+09, -1.50871184e+09,  5.12882796e+08, -1.04708447e+08,
          9.71167420e+06})); 
    ASKAPASSERT(itsParams.size()>0);
}

void NoiseScaler::setFrequencies(const casacore::Vector<double>& freqGHz)
{
    size_t nChan = freqGHz.size();
    ASKAPASSERT(nChan > 0);
    itsScale.resize(nChan);
    itsScale = itsParams[0];
    for (size_t i = 1; i < itsParams.size(); i++) {
        itsScale *= freqGHz;
        itsScale += itsParams[i];
    }
    double mean = casacore::mean(itsScale);
    ASKAPASSERT(mean > 0);
    ASKAPASSERT(casacore::allGT(itsScale,0.))
    itsScale /= mean;
}


