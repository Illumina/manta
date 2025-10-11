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

/// \file
/// \author Chris Saunders
///

///
/// \brief extremely minimal observer pattern
///
/// this class is designed to assist setting up an observer pattern which has zero memory overhead, but
/// has extremely limited abilities, see blt_util/observer.hh for a more general observer pattern support
///

#pragma once

template <typename T>
struct flyweight_notifier;

template <typename T>
struct flyweight_observer {
  friend struct flyweight_notifier<T>;

  flyweight_observer& operator=(const flyweight_observer&) = default;

  virtual ~flyweight_observer() {}

private:
  virtual void recieve_flyweight_notification(const T&) = 0;
};

template <typename T>
struct flyweight_notifier {
  typedef flyweight_observer<T> flyweight_observer_t;

protected:
  void notify_flyweight_observer(flyweight_observer_t* val, const T& msg) const
  {
    val->recieve_flyweight_notification(msg);
  }
};
