#pragma once
#include "Scheduling.hpp"
class Scheduler
{
public:
    ll depth = 25LL;
    virtual void fft(vector<cd>& a, bool invert)
    {
        ll n = a.size();
        for (ll i = 1LL, j = 0LL; i < n; i++)
        {
            ll bit = n >> 1LL;
            for (; j & bit; bit >>= 1LL) j ^= bit;
            j ^= bit;
            if (i < j) swap(a[i], a[j]);
        }
        for (ll len = 2LL; len <= n; len <<= 1LL)
        {
            double ang = 2LL * M_PI / len * (invert ? -1LL : 1LL);
            cd wlen(cos(ang), sin(ang));
            for (ll i = 0LL; i < n; i += len)
            {
                cd w(1LL);
                for (ll j = 0LL; j < len / 2LL; j++)
                {
                    cd u = a[i + j], v = a[i + j + len / 2LL] * w;
                    a[i + j] = u + v;
                    a[i + j + len / 2LL] = u - v;
                    w *= wlen;
                }
            }
        }
        if (invert) for (cd& x : a) x /= n;
    }
    virtual vector<ll> multiply(vector<ll>const& a, vector<ll> const& b)
    {
        vector<cd> fa(a.begin(), a.end()), fb(b.begin(), b.end());
        ll n = 1LL;
        while (n < ll(a.size() + b.size())) n <<= 1LL;
        fa.resize(n);
        fb.resize(n);
        fft(fa, false);
        fft(fb, false);
        for (ll i = 0LL; i < n; i++) fa[i] *= fb[i];
        fft(fa, true);
        vector<ll> result(n);
        for (ll i = 0LL; i < n; i++) result[i] = (ll)round(fa[i].real());
        return result;
    }
    virtual multiset<ll> sumset(multiset<ll> const& a, multiset<ll> const& b)
    {
        ll n = max(*(--a.end()), *(--b.end()));
        vector<ll> a1(n + 1LL), b1(n + 1LL);
        for (ll i = 0LL; i <= n; i++)
        {
            multiset<ll>::const_iterator search;
            a1[i] = ((search = a.find(i)) != a.end());
            b1[i] = ((search = b.find(i)) != b.end());
        }
        vector<ll> product = multiply(a1, b1);
        multiset<ll> sums;
        for (ll i = 0LL; i < (ll)product.size(); i++) if (product[i]) sums.insert(i);
        return sums;
    }
    virtual multiset<ll> subsetsum(multiset<ll>& X)
    {
        if (X.size() == 1LL)
        {
            X.insert(0LL);
            return X;
        }
        multiset<ll> a, b;
        multiset<ll>::iterator i = X.begin();
        multiset<ll>::iterator it = X.begin();
        for (ll j = 1LL; j < (ll)X.size() / 2LL + 1LL; j++) ++i;
        for (; it != i; ++it) a.insert(*(it));
        for (; it != X.end(); ++it) b.insert(*(it));
        return sumset(subsetsum(a), subsetsum(b));
    }
    virtual bool ckgreaterthanv(vector<ll> const& a, vector<ll> const& b, vector<ll> const& d, vector<pair<ll, ll>> const& a1, vector<pair<ll, ll>> const& b1, vector<vector<vector<ll>>> const& matrix, ll k, ll v)
    {
        assert(a.size() == b.size());
        ll n = a.size();
        ll i1, j1, i2, j2, i3, j3;
        for (i1 = 0LL; i1 < n; i1++) if (a1[i1].first > v) break;
        for (j1 = 0LL; j1 < n; j1++) if (b1[j1].first > v - d[k]) break;
        if (i1 == n || j1 == n) return false;
        ll p = (ll)floor(pow(n, 2.0L / 3.0L));
        i2 = min((ll)(ceil((i1 + 1.0L) / p) * p), n) - 1LL;
        j2 = min((ll)(ceil((j1 + 1.0L) / p) * p), n) - 1LL;
        i3 = i1 / p;
        j3 = j1 / p;
        if (matrix[i3][j3][k] > 0LL) return true;
        for (ll i = i1; i < i2; i++) if (k - a1[i].second >= 0LL && k - a1[i].second < n && (b[k - a1[i].second]) / (4LL * n) >= (b1[j1].first) / (4LL * n)) return true;
        for (ll j = j1; j < j2; j++) if (k - b1[j].second >= 0LL && k - b1[j].second < n && (a[k - b1[j].second]) / (4LL * n) >= (a1[i1].first) / (4LL * n)) return true;
        return false;
    }
    virtual vector<ll> maxminskewedconvolution(vector<ll>& a, vector<ll>& b)
    {
        assert(a.size() == b.size());
        ll n = a.size();
        vector<pair<ll, ll>> a1(n), b1(n);
        vector<ll> a2(n), b2(n);
        vector<ll> d(2LL * n);
        for (ll i = 0LL; i < n; i++) a[i] = a[i] * 4LL * n + i;
        for (ll i = n; i < 2LL * n; i++) b[i - n] = b[i - n] * 4LL * n + i;
        for (ll i = 0LL; i < 2LL * n; i++) d[i] = 4LL * n * i;
        for (ll i = 0LL; i < n; i++) a1[i] = { a[i], i };
        for (ll i = 0LL; i < n; i++) b1[i] = { b[i], i };
        ll p = (ll)floor(pow(n, 2.0L / 3.0L));
        sort(a1.begin(), a1.end());
        sort(b1.begin(), b1.end());
        vector<vector<vector<ll>>> matrix;
        ll i = 0LL;
        for (ll u = 0LL; u < n; u += p)
        {
            if (u == 0LL) continue;
            matrix.push_back({});
            for (ll w = 0; w < n; w += p)
            {
                if (w == 0LL) continue;
                for (ll i = 0LL; i < n; i++) a2[i] = (a[i] / (4LL * n) >= a1[u - 1LL].first / (4LL * n));
                for (ll i = 0LL; i < n; i++) b2[i] = (b[i] / (4LL * n) >= b1[w - 1LL].first / (4LL * n));
                matrix[i].push_back(multiply(a2, b2));
            }
            i++;
        }
        i = 0LL;
        for (ll u = 0LL; u < n; u += p)
        {
            ll w = n;
            if (u == 0LL) continue;
            for (ll i = 0LL; i < n; i++) a2[i] = (a[i] / (4LL * n) >= a1[u - 1LL].first / (4LL * n));
            for (ll i = 0LL; i < n; i++) b2[i] = (b[i] / (4LL * n) >= b1[w - 1LL].first / (4LL * n));
            matrix[i].push_back(multiply(a2, b2));
            i++;
        }
        matrix.push_back({});
        for (ll w = 0; w < n; w += p)
        {
            ll u = n;
            if (w == 0LL) continue;
            for (ll i = 0LL; i < n; i++) a2[i] = (a[i] / (4LL * n) >= a1[u - 1LL].first / (4LL * n));
            for (ll i = 0LL; i < n; i++) b2[i] = (b[i] / (4LL * n) >= b1[w - 1LL].first / (4LL * n));
            matrix[i].push_back(multiply(a2, b2));
        }
        ll u = n, w = n;
        for (ll i = 0LL; i < n; i++) a2[i] = (a[i] / (4LL * n) >= a1[u - 1LL].first / (4LL * n));
        for (ll i = 0LL; i < n; i++) b2[i] = (b[i] / (4LL * n) >= b1[w - 1LL].first / (4LL * n));
        matrix[i].push_back(multiply(a2, b2));
        vector<ll> c(2LL * n);
        for (ll k = 0LL; k < 2LL * n; k++)
        {
            vector<ll> options;
            ll l = 0LL, r = a1.size() - 1LL;
            while (l < r)
            {
                ll m = l + (r - l) / 2LL;
                if (ckgreaterthanv(a, b, d, a1, b1, matrix, k, a1[m].first)) l = m + 1LL;
                else r = m;
            }
            ll candidateone = a1[l].first;
            l = 0LL;
            r = b1.size() - 1LL;
            while (l < r)
            {
                ll m = l + (r - l) / 2LL;
                if (ckgreaterthanv(a, b, d, a1, b1, matrix, k, b1[m].first + d[k])) l = m + 1LL;
                else r = m;
            }
            ll candidatetwo = b1[r].first + d[k];
            c[k] = min(candidateone, candidatetwo) / (4LL * n);
        }
        return c;
    }
    virtual void fft_test()
    {
        ofstream fout("data/fftunittest.txt");
        srand((unsigned int)(time(0)));
        for (ll n = 1LL; n <= depth; n++)
        {
            vector<ll> a(n), b(n);
            ll sum = 0LL;
            for (ll j = 0LL; j < depth; j++)
            {
                for (ll i = 0LL; i < n; i++)
                {
                    a[i] = rand() % depth;
                    b[i] = rand() % depth;
                }
                auto start = chrono::high_resolution_clock::now();
                vector<ll> ans = multiply(a, b);
                auto stop = chrono::high_resolution_clock::now();
                auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
                sum += duration.count();
            }
            sum /= depth;
            fout << "(" << n << ", " << sum << ")" << endl;
            cout << n << endl;
        }
        fout.close();
    }
    virtual void sumset_test()
    {
        ofstream fout("data/sumsetunittest.txt");
        srand((unsigned int)(time(0)));
        for (ll n = 1LL; n <= depth; n++)
        {
            multiset<ll> a, b;
            ll sum = 0LL;
            for (ll j = 0LL; j < depth; j++)
            {
                for (ll i = 0LL; i < n; i++)
                {
                    a.insert(rand() % n);
                    b.insert(rand() % n);
                }
                auto start = chrono::high_resolution_clock::now();
                multiset<ll> sums = sumset(a, b);
                auto stop = chrono::high_resolution_clock::now();
                auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
                sum += duration.count();
            }
            sum /= depth;
            fout << "(" << n << ", " << sum << ")" << endl;
            cout << n << endl;
        }
        fout.close();
    }
    virtual void subsetsum_test()
    {
        ofstream fout("data/subsetsumunittest.txt");
        srand((unsigned int)(time(0)));
        for (ll n = 1LL; n <= depth; n++)
        {
            multiset<ll> a;
            ll sum = 0LL;
            for (ll j = 0LL; j < depth; j++)
            {
                ll P = n;
                while (P)
                {
                    ll next = (rand() % P) + 1LL;
                    P -= next;
                    a.insert(next);
                }
                auto start = chrono::high_resolution_clock::now();
                multiset<ll> sums = subsetsum(a);
                auto stop = chrono::high_resolution_clock::now();
                auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
                sum += duration.count();
            }
            sum /= depth;
            fout << "(" << n << ", " << sum << ")" << endl;
            cout << n << endl;
        }
        fout.close();
    }
    virtual void maxminskewed_test()
    {
        ofstream fout("data/maxminskewedconvolutionunittest.txt");
        srand((unsigned int)(time(0)));
        for (ll n = 1LL; n <= depth; n++)
        {
            vector<ll> a, b;
            ll sum = 0LL;
            for (ll j = 0LL; j < depth; j++)
            {
                for (ll i = 0LL; i < n; i++)
                {
                    a.push_back(rand() % n);
                    b.push_back(rand() % n);
                }
                auto start = chrono::high_resolution_clock::now();
                vector<ll> sums = maxminskewedconvolution(a, b);
                auto stop = chrono::high_resolution_clock::now();
                auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
                sum += duration.count();
            }
            sum /= depth;
            fout << "(" << n << ", " << sum << ")" << endl;
            cout << n << endl;
        }
        fout.close();
    }
};