#pragma once
#include "GenScheduler.hpp"
#include "ConvScheduler.hpp"
#include "SolvScheduler.hpp"
using namespace std;
void unittests()
{
    GenScheduler schedule = GenScheduler();
    schedule.fft_test();
    schedule.sumset_test();
    schedule.subsetsum_test();
    schedule.maxminskewed_test();
}
void convscheduler()
{
    ConvScheduler schedule = ConvScheduler();
    schedule.insert( 2, 2 );
    schedule.insert( 2, 2 );
    schedule.insert( 4, 1 );
    schedule.insert( 6, 4 );
    cout << schedule.run() << endl;
    schedule.insert(5, 1);
    cout << schedule.run() << endl;
    schedule.insert(3, 1);
    cout << schedule.run() << endl;
    schedule.insert(1, 1);
    cout << schedule.run() << endl;
    schedule.remove(5, 1);
    cout << schedule.run() << endl;
}
void solve()
{
    SolvScheduler schedule = SolvScheduler();
    schedule.insert(2, 2);
    schedule.insert(2, 2);
    schedule.insert(4, 1);
    schedule.insert(6, 4);
    cout << schedule.run() << endl;
    schedule.insert(5, 1);
    cout << schedule.run() << endl;
    schedule.insert(3, 1);
    cout << schedule.run() << endl;
    schedule.insert(1, 1);
    cout << schedule.run() << endl;
    schedule.remove(5, 1);
    cout << schedule.run() << endl;
}
void compare()
{
    ll depth = 100LL;
    ConvScheduler oldScheduler = ConvScheduler();
    SolvScheduler newScheduler = SolvScheduler();
    ofstream fout("data/schedulerdata.txt");
    srand((unsigned int)(time(0)));
    for (ll n = 1LL; n <= 1000LL; n++)
    {
        ll sum1 = 0LL, sum2 = 0LL;
        for (ll j = 0LL; j < depth; j++)
        {
            oldScheduler.clear();
            newScheduler.clear();
            ll P = n;
            while (P)
            {
                ll nextlength = (rand() % P) + 1LL;
                ll nexttime = (rand() % P) + 1LL;
                P -= nextlength;
                oldScheduler.insert(nexttime, nextlength);
                newScheduler.insert(nexttime, nextlength);
            }
            auto start = chrono::high_resolution_clock::now();
            ll ans1 = oldScheduler.run();
            auto stop = chrono::high_resolution_clock::now();
            auto duration1 = chrono::duration_cast<chrono::microseconds>(stop - start);
            start = chrono::high_resolution_clock::now();
            ll ans2 = newScheduler.run();
            stop = chrono::high_resolution_clock::now();
            auto duration2 = chrono::duration_cast<chrono::microseconds>(stop - start);
            assert(ans1 == ans2);
            sum1 += duration1.count();
            sum2 += duration2.count();
        }
        sum1 /= depth;
        sum2 /= depth;
        fout << n << " " << sum1 << " " << sum2 << endl;
        cout << n << endl;
    }
    fout.close();
}
int main()
{
    compare();
}