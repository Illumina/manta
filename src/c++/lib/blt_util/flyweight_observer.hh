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
struct flyweight_observer
{
    friend struct flyweight_notifier<T>;

    virtual ~flyweight_observer() {}

private:
    virtual void
    recieve_flyweight_notification(const T&) = 0;
};


template <typename T>
struct flyweight_notifier
{
    typedef flyweight_observer<T> flyweight_observer_t;

protected:
    void
    notify_flyweight_observer(
        flyweight_observer_t* val,
        const T& msg) const
    {
        val->recieve_flyweight_notification(msg);
    }
};

