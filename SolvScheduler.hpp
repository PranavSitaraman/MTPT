#pragma once
#include "Scheduler.hpp"
class SolvScheduler : private Scheduler
{
private:
    multiset<pair<ll, ll>> J;
public:
    using Scheduler::fft_test;
    using Scheduler::sumset_test;
    using Scheduler::subsetsum_test;
    using Scheduler::maxminskewed_test;
    SolvScheduler()
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
        vector<pair<pair<ll, ll>, pair<ll, multiset<ll>>>> X;
        multiset<ll> T = { 0 };
        ll P = 0LL, dsharp, count = 1LL, sum = 0LL, cur = 0LL, k = 0LL, Pi, size;
        bool next = false;
        for (auto i = J.begin(); i != J.end(); i++)
        {
            auto j = i;
            if (i == J.begin() || (*i).first != (*(--j)).first) X.push_back({ {(*i).first, 0LL}, {-1LL, multiset<ll>()} });
            X[X.size() - 1LL].second.second.insert((*i).second);
            X[X.size() - 1LL].first.second += (*i).second;
            P += (*i).second;
        }
        dsharp = X.size();
        for (ll i = 0LL; i < dsharp; i++)
        {
            if (X[i].first.second > (ll)(pow(P, 0.6L)))
            {
                X[i].second.first = 0LL;
                sum = 0LL;
                if (next) count++;
                continue;
            }
            next = true;
            if ((X[i].second.first = count) && (sum = sum + X[i].first.second) > (ll)(pow(P, 0.6L)))
            {
                sum = X[i].first.second;
                X[i].second.first = ++count;
            }
        }
        vector<multiset<ll>> S(dsharp);
        for (ll i = 0LL; i < dsharp; i++)
        {
            if (X[i].second.first == 0LL)
            {
                T = sumset(T, (S[i] = subsetsum(X[i].second.second)));
                for (auto j = T.begin(); j != T.end(); ++j) if (*j > X[i].first.first) j = --(T.erase(j));
                continue;
            }
            if (X[i].second.first != cur) cur = X[(k = i)].second.first;
            if (i == dsharp - 1LL || X[i].second.first != X[i + 1LL].second.first)
            {
                Pi = 0LL;
                for (ll j = k; j <= i; j++)
                {
                    S[j] = subsetsum(X[j].second.second);
                    Pi += X[j].first.second;
                }
                vector<vector<ll>> M(i - k + 1LL, vector<ll>(2LL * Pi + 1LL, 0LL));
                for (ll j = k; j <= i; j++)
                {
                    for (ll l = 0LL; l <= 2LL * Pi; l++)
                    {
                        if (l == 0LL) M[j - k][l] = X[j].first.first;
                        else if (S[j].find(l) != S[j].end() && l <= X[j].first.first) M[j - k][l] = X[j].first.first - l;
                        else M[j - k][l] = NINF;
                    }
                }
                size = i - k + 1LL;
                while (size > 1LL)
                {
                    for (ll j = 1LL; j <= size / 2LL; j++)
                    {
                        M[2LL * j - 2LL][0LL] = M[2LL * j - 1LL][0LL];
                        vector<ll> b = M[2LL * j - 2LL], a = M[2LL * j - 1LL];
                        for (ll k = 0LL; k <= 2LL * Pi; k++) a[k] += k;
                        M[j - 1LL] = maxminskewedconvolution(a, b);
                        M[j - 1LL].resize(2LL * Pi + 1);
                        for (ll k = 0LL; k <= 2LL * Pi; k++) if ((M[j - 1LL][k] = M[j - 1LL][k] - k) < 0LL) M[j - 1LL][k] = NINF;
                    }
                    for (ll j = size / 2LL + 1LL; 2LL * j - 1LL <= size; j++) M[j - 1LL] = M[2LL * j - 2LL];
                    size = ll(ceil(size / 2.0L));
                }
                multiset<ll> Si, prefix = T, suffix = T;
                for (ll j = 0LL; j <= Pi; j++) if (M[0LL][j] >= 0LL) Si.insert(j);
                prefix.erase(prefix.upper_bound(max(X[k].first.first - Pi, 0LL)), prefix.end());
                suffix.erase(suffix.begin(), suffix.lower_bound(max(X[k].first.first - Pi, 0LL)));
                if (X[k].first.first >= Pi) for (auto j : sumset(prefix, Si)) if (T.find(j) == T.end()) T.insert(j);
                vector<ll> Mprime(2LL * Pi + 1LL, NINF);
                for (auto j : suffix) Mprime[j - max(X[k].first.first - Pi, 0LL)] = 0LL;
                for (ll j = 0LL; j <= 2LL * Pi; j++) M[0LL][j] -= max(0LL, X[k].first.first - Pi);
                vector<ll> b = Mprime, a = M[0LL];
                for (ll j = 0LL; j <= 2LL * Pi; j++) a[j] += j;
                vector<ll> C = maxminskewedconvolution(a, b);
                C.resize(2LL * Pi + 1LL);
                for (ll j = 0LL; j <= 2LL * Pi; j++) if ((C[j] = C[j] - j) < 0LL) C[j] = NINF;
                for (ll j = 0LL; j <= 2LL * Pi; j++) if (C[j] == 0LL) if (T.find(j + max(X[k].first.first - Pi, 0LL)) == T.end()) T.insert(j + max(X[k].first.first - Pi, 0LL));
                for (auto j = T.begin(); j != T.end(); ++j) if (*j > X[i].first.first) j = --(T.erase(j));
            }
        }
        return (T.size() == 0LL) ? -1LL : P - (*(--T.end()));
    }
    ~SolvScheduler()
    {
        J.clear();
    }
};