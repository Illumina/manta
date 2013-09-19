// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Manta
// Copyright (c) 2013 Illumina, Inc.
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

#include "String.hh"
#include "common/Exceptions.hh"

/*****************************************************************************/

class ChromosomeMetadata
{
public:
    const char *getKey() const              { return id_.c_str(); }
    const std::string &getId() const        { return id_.string(); }
    size_t getLength() const                { return length_; }
    const std::string &getAssembly() const  { return assembly_; }
    const std::string &getMd5() const       { return md5_; }
    const std::string &getSpecies() const   { return species_; }
    const std::string &getUrl() const       { return url_; }

    void setId(const std::string &id)                   { id_ = id; }
    void setLength(size_t length)                       { length_ = length; }
    void setAssembly(const std::string &assembly)       { assembly_ = assembly; }
    void setMd5(const std::string &md5)                 { md5_ = md5;}
    void setSpecies(const std::string &species)         { species_ = species;}
    void setUrl(const std::string &url)                 { url_ = url;}

private:
    FastString id_;
    size_t length_;
    std::string assembly_;
    std::string md5_;
    std::string species_;
    std::string url_;
};

/*****************************************************************************/
