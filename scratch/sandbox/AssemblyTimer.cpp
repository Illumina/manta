//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2018 Illumina, Inc.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//

///
/// \author Chris Saunders
///

#include "AssemblyTimer.hh"

#include <iostream>



AssemblyTimer::
~AssemblyTimer()
{
    std::ostream& os (std::cerr);

    os << "[AssemblyTimingInfo]\n";
    os << "total\t" << totalTime.getUserSeconds() << "\n";
    os << "read\t" << readTime.getUserSeconds() << "\n";
    os << "coreAssm\t" << coreAssmTime.getUserSeconds() << "\n";
    os << "graph\t" << graphBuildTime.getUserSeconds() << "\n";
    os << "findContig\t" << findContigTime.getUserSeconds() << "\n";
    os << "align\t" << alignTime.getUserSeconds() << "\n";
}
