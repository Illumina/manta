// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013-2014 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//

///
/// \author Refactored by Richard Shaw from PUMA by Mauricio Varea
///

/*****************************************************************************/

#pragma once

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

#include <boost/throw_exception.hpp>

#include "common/Exceptions.hh"

using namespace illumina::common;

/*****************************************************************************/

/**
 * \brief Genomic coordinates
 */
template< typename Chr, typename Pos >
class VcfLocus
{
public:
    VcfLocus()
        : chromosome_( UndefinedChr() )
        , position_( UndefinedPos() )
    {}
    virtual ~VcfLocus() {}

    Chr getChromosome() const
    {
        return chromosome_;
    }
    void setChromosome(Chr chromosome)
    {
        chromosome_ = chromosome;
    }
    Chr UndefinedChr() const
    {
        return std::numeric_limits<Chr>::max();
    }
    bool chrDefined() const
    {
        return UndefinedChr() != chromosome_;
    }

    Pos getPosition() const
    {
        return position_;
    }
    void setPosition( Pos position)
    {
        position_ = position;
    }
    Pos UndefinedPos() const
    {
        return std::numeric_limits<Pos>::max();
    }
    bool posDefined() const
    {
        return UndefinedPos() != position_;
    }

    // TODO: vinculate position with chromosome, one we bring in the limits of each chromosome
    void incPosition()
    {
        ++position_;
    }
    void decPosition()
    {
        --position_;
    }

    void swap(VcfLocus& rhs)
    {
        std::swap( chromosome_, rhs.chromosome_ );
        std::swap( position_, rhs.position_ );
    }

    VcfLocus& operator=(const VcfLocus& rhs)
    {
        chromosome_ = rhs.chromosome_;
        position_ = rhs.position_;
        return *this;
    }

    bool operator<( const VcfLocus& rhs ) const
    {
        if (!chrDefined())
        {
            BOOST_THROW_EXCEPTION( PreConditionException(
                                       "Undefined CHROMOSOME reference: Could not compare Genomic coordinates." ));
        }
        if (!posDefined())
        {
            BOOST_THROW_EXCEPTION( PreConditionException(
                                       "Undefined POSITION: Could not compare Genomic coordinates." ));
        }
        return this->chromosome_ < rhs.chromosome_  ||
               ( this->chromosome_ == rhs.chromosome_
                 && this->position_ < rhs.position_ );
    }

protected:
    Chr chromosome_;
    Pos position_;
};

/*****************************************************************************/
