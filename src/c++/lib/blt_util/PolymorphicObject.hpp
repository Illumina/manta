//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2025 Illumina, Inc.
//
// This program is licensed under the terms of the Polyform strict license
//
// ***As far as the law allows, the software comes as is, without
// any warranty or condition, and the licensor will not be liable
// to you for any damages arising out of these terms or the use
// or nature of the software, under any kind of legal claim.***
//
// You should have received a copy of the PolyForm Strict License 1.0.0
// along with this program.  If not, see <https://polyformproject.org/licenses/strict/1.0.0>.
//
//

#pragma once

/// \brief Convenience base class for polymorphic objects
///
/// This class helps get around some of the boilerplate imposed by
/// c++11 spec implicit copy ctor rules for virtual classes. Use
/// this as a public base class for any standard virtual object with
/// default copy/move semantics.
///
/// Per suggestion on thread:
/// http://stackoverflow.com/questions/19997646/no-implicit-copy-constructor-in-polymorphic-class
///
struct PolymorphicObject {
  PolymorphicObject()          = default;
  virtual ~PolymorphicObject() = default;

  PolymorphicObject(const PolymorphicObject&) = default;
  PolymorphicObject& operator=(const PolymorphicObject&) = default;

#if ((!defined(_MSC_VER)) || (_MSC_VER > 1800))
  // support moving
  PolymorphicObject(PolymorphicObject&&) = default;
  PolymorphicObject& operator=(PolymorphicObject&&) = default;
#endif
};
