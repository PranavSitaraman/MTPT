#pragma once
#include "Scheduler.hpp"
class ConvScheduler: private Scheduler
{
private:
    multiset<pair<ll, ll>> J;
public:
    using Scheduler::fft_test;
    using Scheduler::sumset_test;
    using Scheduler::subsetsum_test;
    using Scheduler::maxminskewed_test;
    ConvScheduler()
    {
        J = multiset<pair<ll, ll>>();
    }
    void insert(ll due, ll process)
    {
        J.insert({ due, process });
    }
    void remove(ll due, ll process)
    {
        J.erase({ due, process });
    }
    void clear()
    {
        J.clear();
    }
    ll run()
    {
        vector<pair<ll, multiset<ll>>> X;
        ll P = 0LL;
        for (auto i = J.begin(); i != J.end(); i++)
        {
            auto j = i;
            if (i == J.begin() || (*i).first != (*(--j)).first) X.push_back({ (*i).first, multiset<ll>() });
            X[X.size() - 1LL].second.insert((*i).second);
            P += (*i).second;
        }
        ll dsharp = X.size();
        vector<multiset<ll>> S(dsharp);
        vector<vector<ll>> M(dsharp, vector<ll>(P + 1LL, 0LL));
        for (ll i = 0LL; i < dsharp; i++)
        {
            S[i] = subsetsum(X[i].second);
            for (ll j = 0LL; j <= P; j++)
            {
                if (j == 0LL) M[i][j] = X[i].first;
                else if (S[i].find(j) != S[i].end() && j <= X[i].first) M[i][j] = X[i].first - j;
                else M[i][j] = NINF;
            }
        }
        while (dsharp > 1LL)
        {
            for (ll i = 1LL; i <= dsharp / 2LL; i++)
            {
                M[2LL * i - 2LL][0LL] = M[2LL * i - 1LL][0LL];
                vector<ll> b = M[2LL * i - 2LL], a = M[2LL * i - 1LL];
                for (ll j = 0LL; j <= P; j++) a[j] += j;
                M[i - 1LL] = maxminskewedconvolution(a, b);
                M[i - 1LL].resize(P + 1);
                for (ll j = 0LL; j <= P; j++) if ((M[i - 1LL][j] = M[i - 1LL][j] - j) < 0LL) M[i - 1LL][j] = NINF;
            }
            for (ll i = dsharp / 2LL + 1LL; 2LL * i - 1LL <= dsharp; i++) M[i - 1LL] = M[2LL * i - 2LL];
            dsharp = ll(ceil(dsharp / 2.0L));
        }
        for (ll i = P; i >= 0LL; i--) if (M[0LL][i] >= 0LL) return P - i;
        return -1LL;
    }
    ~ConvScheduler()
    {
        J.clear();
    }
};