//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2017 Illumina, Inc.
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

#pragma once


/// convenience base class for polymorphic objects
///
/// This class helps get around some of the boilerplate imposed by
/// c++ spec implicit copy ctor rules for virtual classes. Use
/// this as a base class for any standard virtual object with liberal
/// default copy/move semantics.
///
/// Per suggestion from: http://stackoverflow.com/questions/19997646/no-implicit-copy-constructor-in-polymorphic-class
///
struct PolymorphicObject
{
    PolymorphicObject() = default;
    virtual ~PolymorphicObject() = default;

    explicit
    PolymorphicObject(const PolymorphicObject&) = default;
    PolymorphicObject& operator =(const PolymorphicObject&) = default;

#if ((!defined(_MSC_VER)) || (_MSC_VER > 1800))
    // support moving
    explicit
    PolymorphicObject(PolymorphicObject&&) = default;
    PolymorphicObject& operator=(PolymorphicObject&&) = default;
#endif
};
