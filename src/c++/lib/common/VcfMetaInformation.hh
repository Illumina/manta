// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013-2015 Illumina, Inc.
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

#include <string>
#include <iostream>

#include "String.hh"

/*****************************************************************************/

namespace VCF_VALUE_TYPE
{
enum index_t
{
    NONE=0,
    INTEGER,
    FLOAT,
    FLAG,
    CHARACTER,
    STRING,
    SIZE
};
}

/*****************************************************************************/

class VcfMetaInformation
{
public:
    VcfMetaInformation();

    const char* getKey() const
    {
        return id_.c_str();
    }

    const std::string& getId() const
    {
        return id_.string();
    }

    const std::string& getNumber() const
    {
        return number_;
    }

    const VCF_VALUE_TYPE::index_t& getType() const
    {
        return type_;
    }

    const std::string& getDescription() const
    {
        return description_;
    }

    void setId(const std::string& id)
    {
        id_ = id;
    }

    void setNumber(const std::string& num)
    {
        validateNumber(num);
        number_ = num;
    }

    void setType(VCF_VALUE_TYPE::index_t type)
    {
        type_ = type;
    }

    void setDescription(const std::string& desc)
    {
        validateDescription(desc);
        description_ = desc;
    }

private:
    void validateFlag() const;
    void validateNumber(const std::string& num) const;
    void validateDescription(const std::string& desc) const;

    FastString id_;
    std::string number_;
    /// when type_ == NONE, it means that only id_ and description_ are valid
    VCF_VALUE_TYPE::index_t type_;
    std::string description_;

    friend std::istream& operator>>(std::istream& is,
                                    VcfMetaInformation& info);
};

/*****************************************************************************/

std::istream& operator>>(std::istream& is, VcfMetaInformation& info);
std::ostream& operator<<(std::ostream& os, const VcfMetaInformation& info);

/*****************************************************************************/
