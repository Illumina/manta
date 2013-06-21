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
// <https://github.com/downloads/sequencing/licenses/>.
//

///
/// \author Chris Saunders
///

///
/// \brief simple observer/notifier pattern
///
/// see unit test for demonstration, note this is not meant to be used across treads
///

#pragma once

#include "boost/foreach.hpp"

#include <set>



template <typename T>
struct notifier;


template <typename T>
struct observer
{
  friend struct notifier<T>;

  typedef observer self_t;

  observer() {}

  observer(const self_t&)  {} // do not copy notifier set

  virtual ~observer()
  {
      BOOST_FOREACH(typename nots_t::value_type val, _nots)
      {
          val->unregister_observer(this);
      }
  }

protected:
  void
  observe_notifier(const notifier<T>& n)
  {
      n.register_observer(this);
      _nots.insert(&n);
  }

private:
  self_t& operator=(const self_t&);

  virtual void
  recieve_notification(const notifier<T>&,
                       const T&) = 0;

  void
  unregister_notifier(const notifier<T>* n)
  {
    const typename nots_t::iterator i(_nots.find(n));
    if(i != _nots.end()) _nots.erase(i);
  }

  ////////// data:
  typedef typename std::set<const notifier<T>*> nots_t;
  mutable nots_t _nots;
};


template <typename T>
struct notifier {
  friend struct observer<T>;

  typedef notifier self_t;

  notifier() {}

  notifier(const self_t&)  {} // do not copy observer set

  virtual ~notifier()
  {
      BOOST_FOREACH(typename obss_t::value_type val, _obss)
      {
          val->unregister_notifier(this);
      }
  }

protected:
  void
  notify_observers(const T& msg) const
  {
      BOOST_FOREACH(typename obss_t::value_type val, _obss)
      {
          val->recieve_notification(*this, msg);
      }
  }

private:
  self_t& operator=(const self_t&);

  void
  register_observer(observer<T>* n) const { _obss.insert(n); }

  void
  unregister_observer(observer<T>* n) const
  {
    const typename obss_t::iterator i(_obss.find(n));
    if(i != _obss.end()) _obss.erase(i);
  }

  ////////// data:
  typedef typename std::set<observer<T>*> obss_t;
  mutable obss_t _obss;
};

