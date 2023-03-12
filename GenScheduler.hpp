#pragma once
#include "Scheduler.hpp"
class GenScheduler : private Scheduler
{
public:
    using Scheduler::fft_test;
    using Scheduler::sumset_test;
    using Scheduler::subsetsum_test;
    using Scheduler::maxminskewed_test;
};